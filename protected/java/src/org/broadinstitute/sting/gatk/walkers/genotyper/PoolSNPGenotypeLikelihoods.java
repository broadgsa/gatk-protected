package org.broadinstitute.sting.gatk.walkers.genotyper;


import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;

import static java.lang.Math.log10;
import static java.lang.Math.pow;


/**
 * Stable, error checking version of the pool genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
*/
public class PoolSNPGenotypeLikelihoods extends PoolGenotypeLikelihoods/* implements Cloneable*/ {

    final List<Allele> myAlleles;
    final int[] alleleIndices;
    final boolean useBAQedPileup;
    final byte refByte;
    int mbq;
    //final double[] PofDGivenBase;

    protected static final double[][][] qualLikelihoodCache;
    /**
     * Create a new GenotypeLikelhoods object with given priors and PCR error rate for each pool genotype
     * @param alleles           Alleles associated with this likelihood object
     * @param logLikelihoods     Likelihoods (can be null if no likelihoods known)
     * @param ploidy            Ploidy of sample (# of chromosomes)
     * @param perLaneErrorModels error model objects for each lane
     * @param useBQAedPileup    Use BAQed pileup
     * @param ignoreLaneInformation  If true, lane info is ignored
     */
    public PoolSNPGenotypeLikelihoods(final List<Allele> alleles, final double[] logLikelihoods, final int ploidy,
                                      final HashMap<String, ErrorModel> perLaneErrorModels, final boolean useBQAedPileup,final boolean ignoreLaneInformation) {
        super(alleles, logLikelihoods, ploidy, perLaneErrorModels, ignoreLaneInformation);
        this.useBAQedPileup = useBQAedPileup;

        myAlleles = new ArrayList<Allele>(alleles);

        refByte = alleles.get(0).getBases()[0];  // by construction, first allele in list is always ref!

        if (myAlleles.size() < BaseUtils.BASES.length) {
            // likelihood only defined for subset of possible alleles. Fill then with other alleles to have all possible ones,
            for (byte b : BaseUtils.BASES) {
                // if base is not included in myAlleles, add new allele
                boolean isRef = (b==refByte);
                if (!myAlleles.contains(Allele.create(b,isRef)))
                    myAlleles.add(Allele.create(b,isRef));

            }

        }


        // compute permutation vector to figure out mapping from indices to bases
        int idx = 0;
        alleleIndices = new int[myAlleles.size()];
        for (byte b : BaseUtils.BASES) {
            boolean isRef = (b==refByte);
            alleleIndices[idx++] = myAlleles.indexOf(Allele.create(b,isRef));
        }

    }

    // -------------------------------------------------------------------------------------
    //
    // add() routines.  These are the workhorse routines for calculating the overall genotype
    // likelihoods given observed bases and reads.  Includes high-level operators all the
    // way down to single base and qual functions.
    //
    // -------------------------------------------------------------------------------------

    public int add(ReadBackedPileup pileup, UnifiedArgumentCollection UAC) {
        mbq = UAC.MIN_BASE_QUALTY_SCORE; // record for later use
        return add(pileup, true, true, mbq);
    }

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @param minBaseQual               the minimum base quality at which to consider a base valid
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        int n = 0;

        if ( useBAQedPileup )
            pileup = createBAQedPileup( pileup );

        if (!hasReferenceSampleData) {
            return add(pileup, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual, null);
        }

        for (String laneID : perLaneErrorModels.keySet() ) {
            // get pileup for this lane
            ReadBackedPileup perLanePileup;
            if (ignoreLaneInformation)
                perLanePileup = pileup;
            else
                perLanePileup = pileup.getPileupForLane(laneID);

            if (perLanePileup == null || perLanePileup.isEmpty())
                continue;

            ErrorModel errorModel = perLaneErrorModels.get(laneID);
            n += add(perLanePileup, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual, errorModel);
            if (ignoreLaneInformation)
                break;

        }

        return n;
    }

    /**
     * Calculates the pool's probability for all possible allele counts for all bases. Calculation is based on the error model
     * generated by the reference sample on the same lane. The probability is given by :
     *
     * Pr(ac=jA,jC,jG,jT| pool, errorModel) = sum_over_all_Qs ( Pr(ac=jA,jC,jG,jT) * Pr(errorModel_q) *
     * [jA*(1-eq)/2n + eq/3*(jc+jg+jt)/2N)^nA *   jC*(1-eq)/2n + eq/3*(ja+jg+jt)/2N)^nC *
     * jG*(1-eq)/2n + eq/3*(jc+ja+jt)/2N)^nG * jT*(1-eq)/2n + eq/3*(jc+jg+ja)/2N)^nT
     *
     *  log Pr(ac=jA,jC,jG,jT| pool, errorModel) = logsum( Pr(ac=jA,jC,jG,jT) * Pr(errorModel_q) *
     * [jA*(1-eq)/2n + eq/3*(jc+jg+jt)/2N)^nA *   jC*(1-eq)/2n + eq/3*(ja+jg+jt)/2N)^nC *
     * jG*(1-eq)/2n + eq/3*(jc+ja+jt)/2N)^nG * jT*(1-eq)/2n + eq/3*(jc+jg+ja)/2N)^nT)
     * = logsum(logPr(ac=jA,jC,jG,jT) + log(Pr(error_Model(q)
     * )) + nA*log(jA/2N(1-eq)+eq/3*(2N-jA)/2N) + nC*log(jC/2N(1-eq)+eq/3*(2N-jC)/2N)
     * + log(jG/2N(1-eq)+eq/3*(2N-jG)/2N) + log(jT/2N(1-eq)+eq/3*(2N-jT)/2N)
     *
     * Let Q(j,k) = log(j/2N*(1-e[k]) + (2N-j)/2N*e[k]/3)
     *
     * Then logPr(ac=jA,jC,jG,jT|D,errorModel) = logPR(ac=Ja,jC,jG,jT) + logsum_k( logPr (errorModel[k],
     * nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
     *
     * If pileup data comes from several error models (because lanes can have different error models),
     * Pr(Ac=j|D,E1,E2) = sum(Pr(AC1=j1|D,E1,E2) * Pr(AC2=j-j2|D,E1,E2))
     * = sum(Pr(AC1=j1|D,E1)*Pr(AC2=j-j1|D,E2)) from j=0..2N
     *
     * So, for each lane, build error model and combine lanes.
     * To store model, can do
     * for jA=0:2N
     *  for jC = 0:2N-jA
     *   for jG = 0:2N-jA-jC
     *    for jT = 0:2N-jA-jC-jG
     *      Q(jA,jC,jG,jT)
     *      for k = minSiteQual:maxSiteQual
     *        likelihood(jA,jC,jG,jT) = logsum(logPr (errorModel[k],nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
     *
     *
     *
     * where: nA,nC,nG,nT = counts of bases observed in pileup.
     *
     *
     * @param pileup                            Base pileup
     * @param ignoreBadBases                    Whether to ignore bad bases
     * @param capBaseQualsAtMappingQual         Cap base at mapping qual
     * @param minBaseQual                       Minimum base quality to consider
     * @param errorModel                        Site error model
     * @return                                  Number of bases added
     */
    private int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual, ErrorModel errorModel) {
         // Number of [A C G T]'s in pileup, in that order
        List<Integer> numSeenBases = new ArrayList<Integer>(BaseUtils.BASES.length);
        for (byte b: BaseUtils.BASES)
            numSeenBases.add(0);

        if (hasReferenceSampleData) {
            // count number of elements in pileup
            for (PileupElement elt : pileup) {
                byte obsBase = elt.getBase();
                byte qual = qualToUse(elt, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
                if ( qual == 0 )
                    continue;
    
                int idx = 0;
    
                for (byte base:BaseUtils.BASES) {
                    int cnt = numSeenBases.get(idx); 
                    numSeenBases.set(idx++,cnt + (base == obsBase?1:0));
    
                }
    
            }
            if (VERBOSE)
                System.out.format("numSeenBases: %d %d %d %d\n",numSeenBases.get(0),numSeenBases.get(1),numSeenBases.get(2),numSeenBases.get(3));
        }
        computeLikelihoods(errorModel, myAlleles, numSeenBases, pileup);
        return pileup.getNumberOfElements();
    }

    /**
     * Compute likelihood of current conformation
     *
     * @param ACset       Count to compute
     * @param errorModel    Site-specific error model object
     * @param alleleList    List of alleles
     * @param numObservations Number of observations for each allele in alleleList
      */
    public void getLikelihoodOfConformation(final AlleleFrequencyCalculationModel.ExactACset ACset,
                                            final ErrorModel errorModel,
                                            final List<Allele> alleleList,
                                            final List<Integer> numObservations,
                                            final ReadBackedPileup pileup) {
        final int[] currentCnt = Arrays.copyOf(ACset.ACcounts.counts, BaseUtils.BASES.length);
        final int[] ac = new int[BaseUtils.BASES.length];
        
        for (int k=0; k < BaseUtils.BASES.length; k++ )
            ac[k] = currentCnt[alleleIndices[k]];

        double p1 = 0.0;
        
        if (!hasReferenceSampleData) {
            // no error model: loop throught pileup to compute likalihoods just on base qualities
            for (final PileupElement elt : pileup) {
                final byte obsBase = elt.getBase();
                final byte qual = qualToUse(elt, true, true, mbq);
                if ( qual == 0 )
                    continue;
                final double acc[] = new double[ACset.ACcounts.counts.length];
                for (int k=0; k < acc.length; k++ )
                    acc[k] = qualLikelihoodCache[BaseUtils.simpleBaseToBaseIndex(alleleList.get(k).getBases()[0])][BaseUtils.simpleBaseToBaseIndex(obsBase)][qual] +MathUtils.log10Cache[ACset.ACcounts.counts[k]]
                            - LOG10_PLOIDY;
                p1 += MathUtils.log10sumLog10(acc);
            }
        }
        else {
            final int minQ = errorModel.getMinSignificantQualityScore();
            final int maxQ = errorModel.getMaxSignificantQualityScore();
            final double[] acVec = new double[maxQ - minQ + 1];
    
            final int nA = numObservations.get(0);
            final int nC = numObservations.get(1);
            final int nG = numObservations.get(2);
            final int nT = numObservations.get(3);
    

            for (int k=minQ; k<=maxQ; k++)
                acVec[k-minQ] = nA*logMismatchProbabilityArray[ac[0]][k] +
                        nC*logMismatchProbabilityArray[ac[1]][k] +
                        nG*logMismatchProbabilityArray[ac[2]][k] +
                        nT*logMismatchProbabilityArray[ac[3]][k];
    
            p1 = MathUtils.logDotProduct(errorModel.getErrorModelVector().getProbabilityVector(minQ,maxQ), acVec);
        }
        ACset.log10Likelihoods[0] = p1;
        /*        System.out.println(Arrays.toString(ACset.ACcounts.getCounts())+" "+String.valueOf(p1));
        System.out.println(Arrays.toString(errorModel.getErrorModelVector().getProbabilityVector(minQ,maxQ)));
      */
    }

    public ReadBackedPileup createBAQedPileup( final ReadBackedPileup pileup ) {
        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
        for( final PileupElement PE : pileup ) {
            final PileupElement newPE = new BAQedPileupElement( PE );
            BAQedElements.add( newPE );
        }
        return new ReadBackedPileupImpl( pileup.getLocation(), BAQedElements );
    }

    public class BAQedPileupElement extends PileupElement {
        public BAQedPileupElement( final PileupElement PE ) {
            super(PE.getRead(), PE.getOffset(), PE.isDeletion(), PE.isBeforeDeletedBase(), PE.isAfterDeletedBase(), PE.isBeforeInsertion(), PE.isAfterInsertion(), PE.isNextToSoftClip());
        }

        @Override
        public byte getQual( final int offset ) { return BAQ.calcBAQFromTag(getRead(), offset, true); }
    }


    /**
     * Helper function that returns the phred-scaled base quality score we should use for calculating
     * likelihoods for a pileup element.  May return 0 to indicate that the observation is bad, and may
     * cap the quality score by the mapping quality of the read itself.
     *
     * @param p                            Pileup element
     * @param ignoreBadBases               Flag to ignore bad bases
     * @param capBaseQualsAtMappingQual    Whether to cap base Q at mapping quality
     * @param minBaseQual                  Min qual to use
     * @return                             New phred-scaled base quality
     */
    private static byte qualToUse(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        if ( ignoreBadBases && !BaseUtils.isRegularBase( p.getBase() ) )
            return 0;

        byte qual = p.getQual();

        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBAM(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)qual, p.getMappingQual());
        if ( (int)qual < minBaseQual )
            qual = (byte)0;

        return qual;
    }

    static {
        qualLikelihoodCache = new double[BaseUtils.BASES.length][BaseUtils.BASES.length][1+SAMUtils.MAX_PHRED_SCORE];
        for (byte j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++) {
            for (byte b1:BaseUtils.BASES) {
                for (byte b2:BaseUtils.BASES) {
                    qualLikelihoodCache[BaseUtils.simpleBaseToBaseIndex(b1)][BaseUtils.simpleBaseToBaseIndex(b2)][j] = log10PofObservingBaseGivenChromosome(b1,b2,j);   
                }
            }                
        }

    }

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param qual         base quality
     * @return log10 likelihood
     */
    private static double log10PofObservingBaseGivenChromosome(byte observedBase, byte chromBase, byte qual) {
        final double log10_3 = log10(3.0);
        double logP;

        if ( observedBase == chromBase ) {
            // the base is consistent with the chromosome -- it's 1 - e
            //logP = oneMinusData[qual];
            double e = pow(10, (qual / -10.0));
            logP = log10(1.0 - e);
        } else {
            // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
            logP = qual / -10.0 + (-log10_3);
        }

        //System.out.printf("%c %c %d => %f%n", observedBase, chromBase, qual, logP);
        return logP;
    }
    
}
