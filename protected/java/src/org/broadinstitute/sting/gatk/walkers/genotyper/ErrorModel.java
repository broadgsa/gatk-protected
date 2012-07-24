package org.broadinstitute.sting.gatk.walkers.genotyper;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 7/21/11
 * Time: 2:21 PM
 *
 * This is a site based implementation of an Error Model. The error model is a probability
 * distribution for the site given the phred scaled quality.
 */
public class ErrorModel  {
    private byte maxQualityScore;
    private byte minQualityScore;
    private byte phredScaledPrior;
    private double log10minPower;
    private int refDepth;
    private boolean hasData = false;
    private ProbabilityVector probabilityVector;
    private static final boolean compressRange = false;
    
    private static final double log10MinusE = Math.log10(Math.exp(1.0));

    /**
     * Calculates the probability of the data (reference sample reads) given the phred scaled site quality score.
     * 
     * @param minQualityScore            Minimum site quality score to evaluate
     * @param maxQualityScore            Maximum site quality score to evaluate
     * @param phredScaledPrior           Prior for site quality
     * @param refSamplePileup            Reference sample pileup
     * @param refSampleVC                VC with True alleles in reference sample pileup
     * @param minPower                   Minimum power
     */
    public ErrorModel (byte minQualityScore, byte maxQualityScore, byte phredScaledPrior,
                       ReadBackedPileup refSamplePileup, VariantContext refSampleVC, double minPower) {
        this.maxQualityScore = maxQualityScore;
        this.minQualityScore = minQualityScore;
        this.phredScaledPrior = phredScaledPrior;
        log10minPower = Math.log10(minPower);


        double[] model = new double[maxQualityScore+1];
        Arrays.fill(model,Double.NEGATIVE_INFINITY);

        boolean hasCalledAlleles = false;
        if (refSampleVC != null) {

            for (Allele allele : refSampleVC.getAlleles()) {
                if (allele.isCalled()) {
                    hasCalledAlleles = true;
                    break;
                }
            }


        }
        if (refSamplePileup == null || refSampleVC == null  || !hasCalledAlleles) {
            double p = MathUtils.phredScaleToLog10Probability((byte)(maxQualityScore-minQualityScore));
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                // maximum uncertainty if there's no ref data at site
                model[q] = p;
            }
            this.refDepth = 0;
        }
        else {
            hasData = true;
            int matches = 0;
            int coverage = refSamplePileup.getNumberOfElements();

            Allele refAllele = refSampleVC.getReference();

            for (PileupElement refPileupElement : refSamplePileup) {
                boolean isMatch = false;
                for (Allele allele : refSampleVC.getAlleles())
                    isMatch |= pileupElementMatches(refPileupElement, allele, refAllele);

                matches += (isMatch?1:0);
          //      System.out.format("MATCH:%b\n",isMatch);
            }

            int mismatches = coverage - matches;
            //System.out.format("Cov:%d match:%d mismatch:%d\n",coverage, matches, mismatches);
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                model[q] = log10PoissonProbabilitySiteGivenQual(q,coverage,  mismatches);
            }
            this.refDepth = coverage;
        }
        
        // compress probability vector
        this.probabilityVector = new ProbabilityVector(model, compressRange);
    }


    /**
     * Simple constructor that just takes a given log-probability vector as error model.
     * Only intended for unit testing, not general usage.
     * @param pvector       Given vector of log-probabilities
     *
     */
    public ErrorModel(double[] pvector) {
        this.maxQualityScore = (byte)(pvector.length-1);
        this.minQualityScore = 0;
        this.probabilityVector = new ProbabilityVector(pvector, compressRange);
        this.hasData = true;

    }

    public static boolean pileupElementMatches(PileupElement pileupElement, Allele allele, Allele refAllele) {
 /*       System.out.format("PE: base:%s isNextToDel:%b isNextToIns:%b eventBases:%s eventLength:%d Allele:%s RefAllele:%s\n",
                pileupElement.getBase(), pileupElement.isBeforeDeletionStart(),
                pileupElement.isBeforeInsertion(),pileupElement.getEventBases(),pileupElement.getEventLength(), allele.toString(), refAllele.toString());
   */

        // if test allele is ref, any base mismatch, or any insertion/deletion at start of pileup count as mismatch
        if (allele.isReference()) {
            // for a ref allele, any base mismatch or new indel is a mismatch.
            if(allele.getBases().length>0 && allele.getBases().length == refAllele.getBases().length ) // SNP/MNP case
                return (/*!pileupElement.isBeforeInsertion() && !pileupElement.isBeforeDeletionStart() &&*/ pileupElement.getBase() == allele.getBases()[0]);
            else
                // either null allele to compare, or ref/alt lengths are different (indel by definition).
                // if we have an indel that we are comparing against a REF allele, any indel presence (of any length/content) is a mismatch
                return (!pileupElement.isBeforeInsertion() && !pileupElement.isBeforeDeletionStart());
        }

        if (refAllele.getBases().length == allele.getBases().length)
            // alleles have the same length (eg snp or mnp)
            return pileupElement.getBase() == allele.getBases()[0];

        // for non-ref alleles,
        byte[] alleleBases = allele.getBases();
        int eventLength = alleleBases.length - refAllele.getBases().length;
        if (eventLength < 0 && pileupElement.isBeforeDeletionStart() && pileupElement.getEventLength() == -eventLength)
            return true;

        if (eventLength > 0 && pileupElement.isBeforeInsertion() &&
                Arrays.equals(pileupElement.getEventBases().getBytes(),alleleBases))
            return true;

        return false;
    }


    /**
     * What's the log-likelihood that a site's quality is equal to q? If we see N observations and n mismatches,
     * and assuming each match is independent of each other and that the match probability is just dependent of
     * the site quality, so p = 10.^-q/10.
     * Since we'll normally have relatively high Q sites and deep coverage in reference samples (ie p small, N high),
     * to avoid underflows we'll use the Poisson approximation with lambda = N*p.
     * Hence, the log-likelihood of q i.e. Pr(Nmismatches = n | SiteQ = q) ~ Poisson(n | lambda = p*N) with p as above.
     * @param q                     Desired q to get likelihood from
     * @param coverage              Total coverage
     * @param mismatches            Number of mismatches
     * @return                      Likelihood of observations as a function of q
     */
    @Requires({
            "q >= minQualityScore",
            "q <= maxQualityScore",
            "coverage >= 0",
            "mismatches >= 0",
            "mismatches <= coverage"
    })
    private double log10PoissonProbabilitySiteGivenQual(byte q, int coverage, int mismatches) {
        // same as   log10ProbabilitySiteGivenQual but with Poisson approximation to avoid numerical underflows
        double lambda = MathUtils.phredScaleToProbability(q) * (double )coverage;
        // log10(e^-lambda*lambda^k/k!) = -lambda + k*log10(lambda) - log10factorial(k)
        return Math.log10(lambda)*mismatches - lambda*log10MinusE- MathUtils.log10Factorial(mismatches);
    }

    @Requires({"qual-minQualityScore <= maxQualityScore"})
    public double getSiteLogErrorProbabilityGivenQual (int qual) {
        return probabilityVector.getLogProbabilityForIndex(qual);
    }

    public byte getMaxQualityScore() {
        return maxQualityScore;
    }

    public byte getMinQualityScore() {
        return minQualityScore;
    }

    public int getMinSignificantQualityScore() {
        return new ProbabilityVector(probabilityVector,true).getMinVal();
    }

    public int getMaxSignificantQualityScore() {
        return new ProbabilityVector(probabilityVector,true).getMaxVal();
    }

    public int getReferenceDepth() {
        return refDepth;
    }
    public boolean hasData() {
        return hasData;
    }

    public ProbabilityVector getErrorModelVector() {
        return probabilityVector;
    }

    public String toString() {
        String result = "(";
        boolean skipComma = true;
        for (double v : probabilityVector.getProbabilityVector()) {
            if (skipComma) {
                skipComma = false;
            }
            else {
                result += ",";
            }
            result += String.format("%.4f", v);
        }
        return result + ")";
    }
    
    public static int getTotalReferenceDepth(HashMap<String, ErrorModel>  perLaneErrorModels) {
        int n=0;
        for (ErrorModel e : perLaneErrorModels.values()) {
            n += e.getReferenceDepth();
        }
        return n;
    }
    
   /* 
@Requires({"maxAlleleCount >= 0"})
//todo -- memoize this function
    public boolean hasPowerForMaxAC (int maxAlleleCount) {
        int siteQ = (int) Math.ceil(MathUtils.probabilityToPhredScale((double) 1/maxAlleleCount));
        double log10CumSum = getCumulativeSum(siteQ);
        return log10CumSum < log10minPower;
    }  */
}
