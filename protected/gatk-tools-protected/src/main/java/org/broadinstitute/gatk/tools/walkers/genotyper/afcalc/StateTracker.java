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

package org.broadinstitute.gatk.tools.walkers.genotyper.afcalc;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.MathUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * TODO this class (+AFCalculator) is a bit messy... it seems that it combines "debugging" (unnecessarily adding CPU cost in production)
 * TODO but it also contains important part of the AF calculation state... why mix both!!!? It seems that the second
 * TODO part could be just blend into AFCalculator ... one one hand you want to reduce classes code size ... but these
 * TODO two classes code seems to be quite intertwine and makes it difficult to understand what is going on.
 * in the production setting without much need
 *
 * Keeps track of the state information during the exact model AF calculation.
 *
 * Tracks things like the MLE and MAP AC values, their corresponding likelihood and posterior
 * values, the likelihood of the AF == 0 state, and the number of evaluations needed
 * by the calculation to compute the P(AF == 0)
 */
final class StateTracker {
    protected static final double VALUE_NOT_CALCULATED = Double.NEGATIVE_INFINITY;
    protected final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    /**
     * These variables are intended to contain the MLE and MAP (and their corresponding allele counts)
     * of the site over all alternate alleles
     */
    protected double log10MLE;
    protected double log10MAP;

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     */
    private int[] alleleCountsOfMLE;
    private int[] alleleCountsOfMAP;

    /**
     * A vector of log10 likelihood values seen, for future summation.  When the size of the
     * vector is exceeed -- because we've pushed more posteriors than there's space to hold
     * -- we simply sum up the existing values, make that the first value, and continue.
     */
    private final double[] log10LikelihoodsForAFGt0 = new double[LIKELIHOODS_CACHE_SIZE];
    private static final int LIKELIHOODS_CACHE_SIZE = 5000;
    private int log10LikelihoodsForAFGt0CacheIndex = 0;

    /**
     * The actual sum of the likelihoods.  Null if the sum hasn't been computed yet
     */
    protected Double log10LikelihoodsForAFGt0Sum = null;

    /**
     * Contains the likelihood for the site's being monomorphic (i.e. AF=0 for all alternate alleles)
     */
    private double log10LikelihoodOfAFzero = 0.0;

    /**
     * The number of evaluates we've gone through in the AFCalc
     */
    private int nEvaluations = 0;

    /**
     * The list of alleles actually used in computing the AF
     */
    private List<Allele> allelesUsedInGenotyping = null;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     *
     * @param maxAltAlleles an integer >= 1
     */
    public StateTracker(final int maxAltAlleles) {
        if ( maxAltAlleles < 0 ) throw new IllegalArgumentException("maxAltAlleles must be >= 0, saw " + maxAltAlleles);

        alleleCountsOfMLE = new int[maxAltAlleles];
        alleleCountsOfMAP = new int[maxAltAlleles];

        reset();
    }

    /**
     * Is the likelihood of configuration K too low to consider, related to the
     * maximum likelihood seen already?
     *
     * @param log10LofK the log10 likelihood of the configuration we're considering analyzing
     * @return true if the configuration cannot meaningfully contribute to our likelihood sum
     */
    private boolean tooLowLikelihood(final double log10LofK) {
        return log10LofK < log10MLE - MAX_LOG10_ERROR_TO_STOP_EARLY;
    }

    /**
     * @return true iff all ACs in this object are less than or equal to their corresponding ACs in the provided set
     */
    private boolean isLowerAC(final ExactACcounts otherACs, final boolean otherACsContainsReference) {
        final int[] otherACcounts = otherACs.getCounts();

        final int firstAltAlleleIndex = otherACsContainsReference ? 1 : 0;

        for ( int i = firstAltAlleleIndex; i < otherACcounts.length; i++ ) {
            if ( alleleCountsOfMLE[i - firstAltAlleleIndex] > otherACcounts[i] )
                return false;
        }

        return true;
    }

    /**
     * Should we stop exploring paths from ACs, given it's log10LofK
     *
     * @param log10LofK the log10LofK of these ACs
     * @param ACs the ACs of this state
     * @param exactACcountsContainReference whether the {@code ACs} contains the reference allele count (index == 0) beside all other alternative alleles.
     * @return return true if there's no reason to continue with subpaths of AC, or false otherwise
     */
    protected boolean abort(final double log10LofK, final ExactACcounts ACs, final boolean enforceLowerACs, final boolean exactACcountsContainReference) {
        return tooLowLikelihood(log10LofK) && (!enforceLowerACs || isLowerAC(ACs,exactACcountsContainReference));
    }

    @Ensures("result != null")
    protected int[] getAlleleCountsOfMAP() {
        return alleleCountsOfMAP;
    }

    @Ensures("result >= 0")
    protected int getNEvaluations() {
        return nEvaluations;
    }

    /**
     * @return the likelihoods summed across all AC values for AC > 0
     */
    private double getLog10LikelihoodOfAFNotZero() {
        if ( log10LikelihoodsForAFGt0Sum == null ) {
            if ( log10LikelihoodsForAFGt0CacheIndex == 0 ) // there's nothing to sum up, so make the sum equal to the smallest thing we have
                log10LikelihoodsForAFGt0Sum = MathUtils.LOG10_P_OF_ZERO;
            else
                log10LikelihoodsForAFGt0Sum = MathUtils.log10sumLog10(log10LikelihoodsForAFGt0, 0, log10LikelihoodsForAFGt0CacheIndex);
        }
        return log10LikelihoodsForAFGt0Sum;
    }

    /**
     * @return the log10 likelihood of AF == 0
     */
    private double getLog10LikelihoodOfAFzero() {
        return log10LikelihoodOfAFzero;
    }

    /**
     * Convert this state to an corresponding AFCalcResult.
     *
     * Assumes that the values in this state have been filled in with meaningful values during the calculation.
     * For example, that the allelesUsedInGenotyping has been set, that the alleleCountsOfMLE contains meaningful
     * values, etc.
     *
     * @param log10PriorsByAC the priors by AC
     *
     * @return an AFCalcResult summarizing the final results of this calculation
     */
    @Requires("allelesUsedInGenotyping != null")
    protected AFCalculationResult toAFCalculationResult(final double[] log10PriorsByAC) {
        final int [] subACOfMLE = Arrays.copyOf(alleleCountsOfMLE, allelesUsedInGenotyping.size() - 1);
        //TODO bad calculation of normalized log10 ACeq0 and ACgt0 likelihoods, priors and consequently posteriors calculated in AFCalculationResult constructor.
        final double[] log10Likelihoods = MathUtils.normalizeFromLog10(new double[]{getLog10LikelihoodOfAFzero(), getLog10LikelihoodOfAFNotZero()}, true);
        final double[] log10Priors = MathUtils.normalizeFromLog10(new double[]{log10PriorsByAC[0], MathUtils.log10sumLog10(log10PriorsByAC, 1)}, true);

        final Map<Allele, Double> log10pRefByAllele = new HashMap<Allele, Double>(allelesUsedInGenotyping.size());
        for ( int i = 0; i < subACOfMLE.length; i++ ) {
            final Allele allele = allelesUsedInGenotyping.get(i+1);
            final double log10PRef = alleleCountsOfMAP[i] > 0 ? -10000 : 0; // TODO -- a total hack but in effect what the old behavior was
            log10pRefByAllele.put(allele, log10PRef);
        }

        return new AFCalculationResult(subACOfMLE, nEvaluations, allelesUsedInGenotyping, log10Likelihoods, log10Priors, log10pRefByAllele);
    }

    // --------------------------------------------------------------------------------
    //
    // Protected mutational methods only for use within the calculation models themselves
    //
    // --------------------------------------------------------------------------------

    /**
     * Reset the data in this results object, so that it can be used in a subsequent AF calculation
     *
     * Resetting of the data is done by the calculation model itself, so shouldn't be done by callers any longer
     *
     * @param ensureAltAlleleCapacity indicate the minimum number of alt-alleles that should be supported by the
     *                                tracker.
     */
    protected void reset(final int ensureAltAlleleCapacity) {
        log10MLE = log10MAP = log10LikelihoodOfAFzero = VALUE_NOT_CALCULATED;
        log10LikelihoodsForAFGt0CacheIndex = 0;
        log10LikelihoodsForAFGt0Sum = null;
        allelesUsedInGenotyping = null;
        nEvaluations = 0;
        if (alleleCountsOfMAP.length < ensureAltAlleleCapacity) {
            final int newCapacity = Math.max(ensureAltAlleleCapacity,alleleCountsOfMAP.length << 1);
            alleleCountsOfMAP = new int[newCapacity];
            alleleCountsOfMLE = new int[newCapacity];
        } else {
            Arrays.fill(alleleCountsOfMLE, 0);
            Arrays.fill(alleleCountsOfMAP, 0);
        }
        Arrays.fill(log10LikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
    }

    /**
     * Reset the data in this results object, so that it can be used in a subsequent AF calculation
     *
     * Resetting of the data is done by the calculation model itself, so shouldn't be done by callers any longer
     */
    protected void reset() {
        log10MLE = log10MAP = log10LikelihoodOfAFzero = VALUE_NOT_CALCULATED;
        log10LikelihoodsForAFGt0CacheIndex = 0;
        log10LikelihoodsForAFGt0Sum = null;
        allelesUsedInGenotyping = null;
        nEvaluations = 0;
        Arrays.fill(alleleCountsOfMLE, 0);
        Arrays.fill(alleleCountsOfMAP, 0);
        Arrays.fill(log10LikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
    }


    /**
     * Tell this result we used one more evaluation cycle
     */
    protected void incNEvaluations() {
        nEvaluations++;
    }

    /**
     * Update the maximum log10 likelihoods seen, if log10LofKs is higher, and the corresponding ACs of this state
     *
     * @param log10LofK the likelihood of our current configuration state, cannot be the 0 state
     * @param alleleCountsForK the allele counts for this state
     */
    @Requires({"alleleCountsForK != null", "MathUtils.sum(alleleCountsForK) >= 0"})
    @Ensures("log10MLE == Math.max(log10LofK, log10MLE)")
    protected void updateMLEifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        addToLikelihoodsCache(log10LofK);

        if ( log10LofK > log10MLE ) {
            log10MLE = log10LofK;
            System.arraycopy(alleleCountsForK, 0, alleleCountsOfMLE, 0, alleleCountsForK.length);
        }
    }

    /**
     * Update the maximum log10 posterior seen, if log10PofKs is higher, and the corresponding ACs of this state
     *
     * @param log10PofK the posterior of our current configuration state
     * @param alleleCountsForK the allele counts for this state
     */
    @Requires({"alleleCountsForK != null", "MathUtils.sum(alleleCountsForK) >= 0"})
    @Ensures("log10MAP == Math.max(log10PofK, log10MAP)")
    protected void updateMAPifNeeded(final double log10PofK, final int[] alleleCountsForK) {
        if ( log10PofK > log10MAP ) {
            log10MAP = log10PofK;
            System.arraycopy(alleleCountsForK, 0, alleleCountsOfMAP, 0, alleleCountsForK.length);
        }
    }

    private void addToLikelihoodsCache(final double log10LofK) {
        // add to the cache
        log10LikelihoodsForAFGt0[log10LikelihoodsForAFGt0CacheIndex++] = log10LofK;

        // if we've filled up the cache, then condense by summing up all of the values and placing the sum back into the first cell
        if ( log10LikelihoodsForAFGt0CacheIndex == LIKELIHOODS_CACHE_SIZE) {
            final double temporarySum = MathUtils.log10sumLog10(log10LikelihoodsForAFGt0, 0, log10LikelihoodsForAFGt0CacheIndex);
            Arrays.fill(log10LikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
            log10LikelihoodsForAFGt0[0] = temporarySum;
            log10LikelihoodsForAFGt0CacheIndex = 1;
        }
    }

    protected void setLog10LikelihoodOfAFzero(final double log10LikelihoodOfAFzero) {
        this.log10LikelihoodOfAFzero = log10LikelihoodOfAFzero;
        if ( log10LikelihoodOfAFzero > log10MLE ) {
            log10MLE = log10LikelihoodOfAFzero;
            Arrays.fill(alleleCountsOfMLE, 0);
        }
    }

    @Requires({"MathUtils.goodLog10Probability(log10PosteriorOfAFzero)"})
    protected void setLog10PosteriorOfAFzero(final double log10PosteriorOfAFzero) {
        if ( log10PosteriorOfAFzero > log10MAP ) {
            log10MAP = log10PosteriorOfAFzero;
            Arrays.fill(alleleCountsOfMAP, 0);
        }
    }

    /**
     * Set the list of alleles used in genotyping
     *
     * @param allelesUsedInGenotyping the list of alleles, where the first allele is reference
     */
    @Requires({"allelesUsedInGenotyping != null", "allelesUsedInGenotyping.size() > 1"})
    protected void setAllelesUsedInGenotyping(List<Allele> allelesUsedInGenotyping) {
        if ( allelesUsedInGenotyping == null || allelesUsedInGenotyping.isEmpty() )
            throw new IllegalArgumentException("allelesUsedInGenotyping cannot be null or empty");
        if ( allelesUsedInGenotyping.get(0).isNonReference() )
            throw new IllegalArgumentException("The first element of allelesUsedInGenotyping must be the reference allele");

        this.allelesUsedInGenotyping = allelesUsedInGenotyping;
    }

    public void ensureMaximumAlleleCapacity(final int capacity) {
        if (this.alleleCountsOfMAP.length < capacity)
            reset(capacity);
    }
}
