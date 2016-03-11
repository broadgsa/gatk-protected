/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.genotyper;


import htsjdk.samtools.SAMUtils;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.ExactACset;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import htsjdk.variant.variantcontext.Allele;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import static java.lang.Math.log10;
import static java.lang.Math.pow;


/**
 * Stable, error checking version of the pool genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
*/
public class GeneralPloidySNPGenotypeLikelihoods extends GeneralPloidyGenotypeLikelihoods/* implements Cloneable*/ {

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
    public GeneralPloidySNPGenotypeLikelihoods(final List<Allele> alleles, final double[] logLikelihoods, final int ploidy,
                                               final HashMap<String, ErrorModel> perLaneErrorModels, final boolean useBQAedPileup, final boolean ignoreLaneInformation) {
        super(alleles, logLikelihoods, ploidy, perLaneErrorModels, ignoreLaneInformation);
        this.useBAQedPileup = useBQAedPileup;

        myAlleles = new ArrayList<Allele>(alleles);

        Allele refAllele = alleles.get(0);
        //sanity check: by construction, first allele should ALWAYS be the reference alleles
        if (!refAllele.isReference())
            throw new ReviewedGATKException("BUG: First allele in list passed to GeneralPloidySNPGenotypeLikelihoods should be reference!");

        refByte = refAllele.getBases()[0];  // by construction, first allele in list is always ref!

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
     * @return                                  Number of bases added - only good bases actually added to GLs are counted.
     */
    private int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual, ErrorModel errorModel) {
         // Number of [A C G T]'s in pileup, in that order
        List<Integer> numSeenBases = new ArrayList<Integer>(BaseUtils.BASES.length);
        for (byte b: BaseUtils.BASES)
            numSeenBases.add(0);

        int nGoodBases = 0;
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
            nGoodBases++;
        }

        if (VERBOSE)
            System.out.format("numSeenBases: %d %d %d %d\n",numSeenBases.get(0),numSeenBases.get(1),numSeenBases.get(2),numSeenBases.get(3));

        computeLikelihoods(errorModel, myAlleles, numSeenBases, pileup);
        return nGoodBases;
    }

    /**
     * Compute likelihood of current conformation
     *
     * @param ACset       Count to compute
     * @param errorModel    Site-specific error model object
     * @param alleleList    List of alleles
     * @param numObservations Number of observations for each allele in alleleList
      */
    public void getLikelihoodOfConformation(final ExactACset ACset,
                                            final ErrorModel errorModel,
                                            final List<Allele> alleleList,
                                            final List<Integer> numObservations,
                                            final ReadBackedPileup pileup) {
        final int[] currentCnt = Arrays.copyOf(ACset.getACcounts().getCounts(), BaseUtils.BASES.length);
        final int[] ac = new int[BaseUtils.BASES.length];
        
        for (int k=0; k < BaseUtils.BASES.length; k++ )
            ac[k] = currentCnt[alleleIndices[k]];

        double p1 = 0.0;
        
        if (!hasReferenceSampleData) {
            // no error model: loop through pileup to compute likelihoods just on base qualities
            // In this case, vector numObservations is not used directly for GL computation
            for (final PileupElement elt : pileup) {
                final byte obsBase = elt.getBase();
                final byte qual = qualToUse(elt, true, true, mbq);
                if ( qual == 0 )
                    continue;
                final double acc[] = new double[ACset.getACcounts().getCounts().length];
                for (int k=0; k < acc.length; k++ )
                    acc[k] = qualLikelihoodCache[BaseUtils.simpleBaseToBaseIndex(alleleList.get(k).getBases()[0])][BaseUtils.simpleBaseToBaseIndex(obsBase)][qual] +
                            MathUtils.Log10Cache.get(ACset.getACcounts().getCounts()[k]) - LOG10_PLOIDY;
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
        ACset.getLog10Likelihoods()[0] = p1;
        /*        System.out.println(Arrays.toString(ACset.ACcounts.getCounts())+" "+String.valueOf(p1));
        System.out.println(Arrays.toString(errorModel.getErrorModelVector().getProbabilityVector(minQ,maxQ)));
      */
    }

    public ReadBackedPileup createBAQedPileup( final ReadBackedPileup pileup ) {
        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
        for( final PileupElement PE : pileup ) {
            final PileupElement newPE = new SNPGenotypeLikelihoodsCalculationModel.BAQedPileupElement( PE );
            BAQedElements.add( newPE );
        }
        return new ReadBackedPileupImpl( pileup.getLocation(), BAQedElements );
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
