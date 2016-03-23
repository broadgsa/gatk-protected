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
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import htsjdk.variant.variantcontext.Allele;

import java.util.*;

/**
 * Describes the results of the AFCalc
 *
 * Only the bare essentials are represented here, as all AFCalc models must return meaningful results for
 * all of these fields.
 *
 * Note that all of the values -- i.e. priors -- are checked now that they are meaningful, which means
 * that users of this code can rely on the values coming out of these functions.
 */
public class AFCalculationResult {
    private final static int AF0 = 0;  //index of the AC=0 entry in the log10 array (of priors, likelihoods, or posteriors listed below)
    private final static int AF1p = 1;  //index of the AC>0 entry in the log10 array
    private final static int LOG_10_ARRAY_SIZES = 2;

    private final double[] log10LikelihoodsOfAC;
    private final double[] log10PriorsOfAC;
    private final double[] log10PosteriorsOfAC;

    private final Map<Allele, Double> log10pRefByAllele;

    /**
     * The AC values for all ALT alleles at the MLE
     */
    private final int[] alleleCountsOfMLE;

    int nEvaluations = 0;

    /**
     * The list of alleles actually used in computing the AF
     */
    private List<Allele> allelesUsedInGenotyping = null;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     */
    public AFCalculationResult(final int[] alleleCountsOfMLE,
                               final int nEvaluations,
                               final List<Allele> allelesUsedInGenotyping,
                               final double[] log10LikelihoodsOfAC,
                               final double[] log10PriorsOfAC,
                               final Map<Allele, Double> log10pRefByAllele) {
        if ( allelesUsedInGenotyping == null || allelesUsedInGenotyping.size() < 1 ) throw new IllegalArgumentException("allelesUsedInGenotyping must be non-null list of at least 1 value " + allelesUsedInGenotyping);
        if ( alleleCountsOfMLE == null ) throw new IllegalArgumentException("alleleCountsOfMLE cannot be null");
        if ( alleleCountsOfMLE.length != allelesUsedInGenotyping.size() - 1) throw new IllegalArgumentException("alleleCountsOfMLE.length " + alleleCountsOfMLE.length + " != number of alternate alleles used in genotyping " + (allelesUsedInGenotyping.size() - 1));
        if ( nEvaluations < 0 ) throw new IllegalArgumentException("nEvaluations must be >= 0 but saw " + nEvaluations);
        if ( log10LikelihoodsOfAC.length != 2 ) throw new IllegalArgumentException("log10LikelihoodsOfAC must have length equal 2");
        if ( log10PriorsOfAC.length != 2 ) throw new IllegalArgumentException("log10PriorsOfAC must have length equal 2");
        if ( log10pRefByAllele == null ) throw new IllegalArgumentException("log10pRefByAllele cannot be null");
        if ( log10pRefByAllele.size() != allelesUsedInGenotyping.size() - 1 ) throw new IllegalArgumentException("log10pRefByAllele has the wrong number of elements: log10pRefByAllele " + log10pRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        if ( ! allelesUsedInGenotyping.containsAll(log10pRefByAllele.keySet()) ) throw new IllegalArgumentException("log10pRefByAllele doesn't contain all of the alleles used in genotyping: log10pRefByAllele " + log10pRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        if ( ! MathUtils.goodLog10ProbVector(log10LikelihoodsOfAC, LOG_10_ARRAY_SIZES, false) ) throw new IllegalArgumentException("log10LikelihoodsOfAC are bad " + Utils.join(",", log10LikelihoodsOfAC));
        if ( ! MathUtils.goodLog10ProbVector(log10PriorsOfAC, LOG_10_ARRAY_SIZES, true) ) throw new IllegalArgumentException("log10priors are bad " + Utils.join(",", log10PriorsOfAC));

        this.alleleCountsOfMLE = alleleCountsOfMLE;
        this.nEvaluations = nEvaluations;
        this.allelesUsedInGenotyping = allelesUsedInGenotyping;

        this.log10LikelihoodsOfAC = Arrays.copyOf(log10LikelihoodsOfAC, LOG_10_ARRAY_SIZES);
        this.log10PriorsOfAC = Arrays.copyOf(log10PriorsOfAC, LOG_10_ARRAY_SIZES);
        this.log10PosteriorsOfAC = computePosteriors(log10LikelihoodsOfAC, log10PriorsOfAC);
        this.log10pRefByAllele = new HashMap<Allele, Double>(log10pRefByAllele);
    }

    /**
     * Return a new AFCalcResult with a new prior probability
     *
     * @param log10PriorsOfAC
     * @return
     */
    public AFCalculationResult withNewPriors(final double[] log10PriorsOfAC) {
        return new AFCalculationResult(alleleCountsOfMLE, nEvaluations, allelesUsedInGenotyping, log10LikelihoodsOfAC, log10PriorsOfAC, log10pRefByAllele);
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     *
     * @return a vector with allele counts, not all of which may be meaningful
     */
    @Ensures("result != null")
    public int[] getAlleleCountsOfMLE() {
        return alleleCountsOfMLE;
    }

    /**
     * Returns the AC of allele a la #getAlleleCountsOfMLE
     *
     * @param allele the allele whose AC we want to know.  Error if its not in allelesUsedInGenotyping
     * @throws IllegalStateException if allele isn't in allelesUsedInGenotyping
     * @return the AC of allele
     */
    public int getAlleleCountAtMLE(final Allele allele) {
        return getAlleleCountsOfMLE()[altAlleleIndex(allele)];
    }

    /**
     * Returns the number of cycles used to evaluate the pNonRef for this AF calculation
     *
     * @return the number of evaluations required to produce the answer for this AF calculation
     */
    public int getnEvaluations() {
        return nEvaluations;
    }

    /**
     * Get the list of alleles actually used in genotyping.
     *
     * Due to computational / implementation constraints this may be smaller than
     * the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping, the first of which is the reference allele
     */
    @Ensures({"result != null", "! result.isEmpty()"})
    public List<Allele> getAllelesUsedInGenotyping() {
        return allelesUsedInGenotyping;
    }

    /**
     * Get the log10 normalized -- across all ACs -- posterior probability of AC == 0 for all alleles
     *
     * @return
     */
    @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PosteriorOfAFEq0() {
        return log10PosteriorsOfAC[AF0];
    }

    /**
     * Get the log10 normalized -- across all ACs -- posterior probability of AC > 0 for any alleles
     *
     * @return
     */
    @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PosteriorOfAFGT0() {
        return log10PosteriorsOfAC[AF1p];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- likelihood of AC == 0 for all alleles
     *
     * @return
     */
    @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10LikelihoodOfAFEq0() {
        return log10LikelihoodsOfAC[AF0];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- likelihood of AC > 0 for any alleles
     *
     * @return
     */
    @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10LikelihoodOfAFGT0() {
        return log10LikelihoodsOfAC[AF1p];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- prior probability of AC == 0 for all alleles
     *
     * @return
     */
    @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PriorOfAFEq0() {
        return log10PriorsOfAC[AF0];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- prior probability of AC > 0
     *
     * @return
     */
    @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PriorOfAFGT0() {
        return log10PriorsOfAC[AF1p];
    }

    @Override
    public String toString() {
        final List<String> byAllele = new LinkedList<String>();
        for ( final Allele a : getAllelesUsedInGenotyping() )
            if ( a.isNonReference() ) byAllele.add(String.format("%s => MLE %d / posterior %.2f", a, getAlleleCountAtMLE(a), getLog10PosteriorOfAFEq0ForAllele(a)));
        return String.format("AFCalc%n\t\tlog10PosteriorOfAFGT0=%.2f%n\t\t%s", getLog10LikelihoodOfAFGT0(), Utils.join("\n\t\t", byAllele));
    }

    /**
     * Are we sufficiently confidence in being non-ref that the site is considered polymorphic?
     *
     * We are non-ref if the probability of being non-ref > the emit confidence (often an argument).
     * Suppose posterior AF > 0 is log10: -5 => 10^-5
     * And that log10minPNonRef is -3.
     * We are considered polymorphic since 10^-5 < 10^-3 => -5 < -3
     *
     * Note that log10minPNonRef is really the minimum confidence, scaled as an error rate, so
     * if you want to be 99% confidence, then log10PNonRef should be log10(0.01) = -2.
     *
     * @param log10minPNonRef the log10 scaled min pr of being non-ref to be considered polymorphic
     *
     * @return true if there's enough confidence (relative to log10minPNonRef) to reject AF == 0
     */
    @Requires("MathUtils.goodLog10Probability(log10minPNonRef)")
    public boolean isPolymorphic(final Allele allele, final double log10minPNonRef) {
        return getLog10PosteriorOfAFEq0ForAllele(allele) < log10minPNonRef;
    }

    /**
     * Same as #isPolymorphic but takes a phred-scaled quality score as input
     */
    public boolean isPolymorphicPhredScaledQual(final Allele allele, final double minPNonRefPhredScaledQual) {
        if ( minPNonRefPhredScaledQual < 0 ) throw new IllegalArgumentException("phredScaledQual " + minPNonRefPhredScaledQual + " < 0 ");
        final double log10Threshold = minPNonRefPhredScaledQual / -10;
        return isPolymorphic(allele, log10Threshold);
    }

    /**
     * Are any of the alleles polymorphic w.r.t. #isPolymorphic?
     *
     * @param log10minPNonRef the confidence threshold, in log10 space
     * @return true if any are poly, false otherwise
     */
    public boolean anyPolymorphic(final double log10minPNonRef) {
        for ( final Allele a : getAllelesUsedInGenotyping() )
            if ( a.isNonReference() && isPolymorphic(a, log10minPNonRef) )
                return true;
        return false;
    }

    /**
     * Returns the log10 probability that allele is not segregating
     *
     * Note that this function is p not segregating so that we can store
     * internally the log10 value of AF == 0, which grows very quickly
     * negative and yet has sufficient resolution for high confidence tests.
     * For example, if log10pRef == -100, not an unreasonably high number,
     * if we tried to store log10pNonRef we'd be looking at 1 - 10^-100, which
     * quickly underflows to 1.  So the logic here is backward from what
     * you really want (the p of segregating) but we do that for numerical
     * reasons
     *
     * Unlike the sites-level annotation, this calculation is specific to allele, and can be
     * used to separately determine how much evidence there is that allele is independently
     * segregating as opposed to the site being polymorphic with any allele.  In the bi-allelic
     * case these are obviously the same but for multiple alt alleles there can be lots of
     * evidence for one allele but not so much for any other allele
     *
     * @param allele the allele we're interested in, must be in getAllelesUsedInGenotyping
     * @return the log10 probability that allele is not segregating at this site
     */
    @Ensures("MathUtils.goodLog10Probability(result)")
    public double getLog10PosteriorOfAFEq0ForAllele(final Allele allele) {
        final Double log10pNonRef = log10pRefByAllele.get(allele);
        if ( log10pNonRef == null ) throw new IllegalArgumentException("Unknown allele " + allele);
        return log10pNonRef;
    }

    /**
     * Returns the log10 normalized posteriors given the log10 likelihoods and priors
     *
     * @param log10LikelihoodsOfAC
     * @param log10PriorsOfAC
     *
     * @return freshly allocated log10 normalized posteriors vector
     */
    @Requires("log10LikelihoodsOfAC.length == log10PriorsOfAC.length")
    @Ensures("MathUtils.goodLog10ProbVector(result, LOG_10_ARRAY_SIZES, true)")
    private static double[] computePosteriors(final double[] log10LikelihoodsOfAC, final double[] log10PriorsOfAC) {
        final double[] log10UnnormalizedPosteriors = new double[log10LikelihoodsOfAC.length];
        for ( int i = 0; i < log10LikelihoodsOfAC.length; i++ )
            log10UnnormalizedPosteriors[i] = log10LikelihoodsOfAC[i] + log10PriorsOfAC[i];
        return MathUtils.normalizeFromLog10(log10UnnormalizedPosteriors, true, false);
    }

    /**
     * Computes the offset into linear vectors indexed by alt allele for allele
     *
     * Things like our MLE allele count vector are indexed by alt allele index, with
     * the first alt allele being 0, the second 1, etc.  This function computes the index
     * associated with allele.
     *
     * @param allele the allele whose alt index we'd like to know
     * @throws IllegalArgumentException if allele isn't in allelesUsedInGenotyping
     * @return an index value greater than 0 suitable for indexing into the MLE and other alt allele indexed arrays
     */
    @Requires("allele != null")
    @Ensures({"result >= 0", "result < allelesUsedInGenotyping.size() - 1"})
    private int altAlleleIndex(final Allele allele) {
        if ( allele.isReference() ) throw new IllegalArgumentException("Cannot get the alt allele index for reference allele " + allele);
        final int index = allelesUsedInGenotyping.indexOf(allele);
        if ( index == -1 )
            throw new IllegalArgumentException("could not find allele " + allele + " in " + allelesUsedInGenotyping);
        else
            return index - 1;
    }
}