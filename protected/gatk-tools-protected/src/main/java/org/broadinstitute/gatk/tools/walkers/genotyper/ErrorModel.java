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

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.tools.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

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
    private static final boolean DEBUG = false;
    /**
     * Calculates the probability of the data (reference sample reads) given the phred scaled site quality score.
     * 
     * @param UAC                           Argument Collection
     * @param refSamplePileup            Reference sample pileup
     * @param refSampleVC                VC with True alleles in reference sample pileup
     */
    public ErrorModel (final UnifiedArgumentCollection UAC,
                       final ReadBackedPileup refSamplePileup,
                       VariantContext refSampleVC, final ReferenceContext refContext) {
        this.maxQualityScore = UAC.maxQualityScore;
        this.minQualityScore = UAC.minQualityScore;
        this.phredScaledPrior = UAC.phredScaledPrior;
        log10minPower = Math.log10(UAC.minPower);

        PairHMMIndelErrorModel pairModel = null;
        LinkedHashMap<Allele, Haplotype> haplotypeMap = null;
        double[][] perReadLikelihoods = null;

        double[] model = new double[maxQualityScore+1];
        Arrays.fill(model,Double.NEGATIVE_INFINITY);

        boolean hasCalledAlleles = false;

        final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        if (refSampleVC != null) {

            for (Allele allele : refSampleVC.getAlleles()) {
                if (allele.isCalled()) {
                    hasCalledAlleles = true;
                    break;
                }
            }
            haplotypeMap = new LinkedHashMap<Allele, Haplotype>();
            if (refSampleVC.isIndel()) {
                pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                        UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.pairHMM);
                IndelGenotypeLikelihoodsCalculationModel.getHaplotypeMapFromAlleles(refSampleVC.getAlleles(), refContext, refContext.getLocus(), haplotypeMap); // will update haplotypeMap adding elements
            }
        }

        double p = QualityUtils.qualToErrorProbLog10((byte)(maxQualityScore-minQualityScore));
        if (refSamplePileup == null || refSampleVC == null  || !hasCalledAlleles) {
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                // maximum uncertainty if there's no ref data at site
                model[q] = p;
            }
            this.refDepth = 0;
        }
        else {
            hasData = true;
            int matches = 0;
            int coverage = 0;

            Allele refAllele = refSampleVC.getReference();

            if ( refSampleVC.isIndel()) {
                //perReadLikelihoods = new double[readCounts.length][refSampleVC.getAlleles().size()];
                final int eventLength = IndelGenotypeLikelihoodsCalculationModel.getEventLength(refSampleVC.getAlleles());
                if (!haplotypeMap.isEmpty())
                    perReadLikelihoods = pairModel.computeGeneralReadHaplotypeLikelihoods(refSamplePileup,haplotypeMap,refContext, eventLength, perReadAlleleLikelihoodMap);
            }
            int idx = 0;
            for (PileupElement refPileupElement : refSamplePileup) {
                if (DEBUG)
                    System.out.println(refPileupElement.toString());
                boolean isMatch = false;
                for (Allele allele : refSampleVC.getAlleles()) {
                    boolean m = pileupElementMatches(refPileupElement, allele, refAllele, refContext.getBase());
                    if (DEBUG) System.out.println(m);
                    isMatch |= m;
                }
                if (refSampleVC.isIndel() && !haplotypeMap.isEmpty()) {
                    // ignore match/mismatch if reads, as determined by their likelihood, are not informative
                    double[] perAlleleLikelihoods = perReadLikelihoods[idx++];
                    if (!isInformativeElement(perAlleleLikelihoods))
                        matches++;
                    else
                        matches += (isMatch?1:0);

                }   else {
                    matches += (isMatch?1:0);
                }
                coverage++;
            }

            int mismatches = coverage - matches;
            //System.out.format("Cov:%d match:%d mismatch:%d\n",coverage, matches, mismatches);
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                if (coverage==0)
                    model[q] = p;
                else
                    model[q] = log10PoissonProbabilitySiteGivenQual(q,coverage,  mismatches);
            }
            this.refDepth = coverage;
        }
        
        // compress probability vector
        this.probabilityVector = new ProbabilityVector(model, compressRange);
    }


    @Requires("likelihoods.length>0")
    private boolean isInformativeElement(double[] likelihoods) {
        // if likelihoods are the same, they're not informative
        final double thresh = 0.1;
        int maxIdx = MathUtils.maxElementIndex(likelihoods);
        int minIdx = MathUtils.minElementIndex(likelihoods);
        if (likelihoods[maxIdx]-likelihoods[minIdx]< thresh)
            return false;
        else
            return true;
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

    public static boolean pileupElementMatches(PileupElement pileupElement, Allele allele, Allele refAllele, byte refBase) {
        if (DEBUG)
            System.out.format("PE: base:%s isNextToDel:%b isNextToIns:%b eventBases:%s eventLength:%d Allele:%s RefAllele:%s\n",
                pileupElement.getBase(), pileupElement.isBeforeDeletionStart(),
                pileupElement.isBeforeInsertion(),pileupElement.getBasesOfImmediatelyFollowingInsertion(),pileupElement.getLengthOfImmediatelyFollowingIndel(), allele.toString(), refAllele.toString());

        //pileupElement.
        // if test allele is ref, any base mismatch, or any insertion/deletion at start of pileup count as mismatch
        if (allele.isReference()) {
            // for a ref allele, any base mismatch or new indel is a mismatch.
            if(allele.getBases().length>0)
                // todo - can't check vs. allele because allele is not padded so it doesn't include the reference base at this location
                // could clean up/simplify this when unpadding is removed
                return (pileupElement.getBase() == refBase && !pileupElement.isBeforeInsertion() && !pileupElement.isBeforeDeletionStart());
            else
                // either null allele to compare, or ref/alt lengths are different (indel by definition).
                // if we have an indel that we are comparing against a REF allele, any indel presence (of any length/content) is a mismatch
                return (!pileupElement.isBeforeInsertion() && !pileupElement.isBeforeDeletionStart());
        }

        // for non-ref alleles to compare:
        if (refAllele.getBases().length == allele.getBases().length)
            // alleles have the same length (eg snp or mnp)
            return pileupElement.getBase() == allele.getBases()[0];

        // for non-ref alleles,
        byte[] alleleBases = allele.getBases();
        int eventLength = alleleBases.length - refAllele.getBases().length;
        if (eventLength < 0 && pileupElement.isBeforeDeletionStart() && pileupElement.getLengthOfImmediatelyFollowingIndel() == -eventLength)
            return true;

                if (eventLength > 0 && pileupElement.isBeforeInsertion() &&
                Arrays.equals(pileupElement.getBasesOfImmediatelyFollowingInsertion().getBytes(),Arrays.copyOfRange(alleleBases,1,alleleBases.length))) // allele contains ref byte, but pileupElement's event bases doesn't
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
        double lambda = QualityUtils.qualToErrorProb(q) * (double )coverage;
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
        StringBuilder result = new StringBuilder("(");
        boolean skipComma = true;
        for (double v : probabilityVector.getProbabilityVector()) {
            if (skipComma) {
                skipComma = false;
            }
            else {
                result.append(",");
            }
            result.append(String.format("%.4f", v));
        }
        result.append(")");
        return result.toString();
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
