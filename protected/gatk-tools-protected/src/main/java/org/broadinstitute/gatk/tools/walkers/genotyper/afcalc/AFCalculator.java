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
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingLikelihoods;
import org.broadinstitute.gatk.utils.SimpleTimer;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.SimpleTimer;

import java.io.File;
import java.util.List;


/**
 * Generic interface for calculating the probability of alleles segregating given priors and genotype likelihoods
 */
public abstract class AFCalculator implements Cloneable {
    private static final Logger defaultLogger = Logger.getLogger(AFCalculator.class);
    public static final int MAX_NUM_PL_VALUES_DEFAULT = 100;

    protected Logger logger = defaultLogger;
    protected int maxNumPLValues = MAX_NUM_PL_VALUES_DEFAULT; // if PL vectors longer than this # of elements, don't log them
    protected static int maxNumPLValuesObserved = 0;
    protected static long numTimesMaxNumPLValuesExceeded = 0;

    private SimpleTimer callTimer = new SimpleTimer();
    private StateTracker stateTracker;
    private ExactCallLogger exactCallLogger = null;

    /**
     * Create a new AFCalc object capable of calculating the prob. that alleles are
     * segregating among many samples.
     *
     * <p>
     *    Restrictions in ploidy and number of alternative alleles that a instance can handle will be determined
     *    by its implementation class {@link AFCalculatorImplementation}
     * </p>
     */
    protected AFCalculator() {
    }

    /**
     * Enable exact call logging to file
     *
     * @param exactCallsLog the destination file
     */
    public void enableProcessLog(final File exactCallsLog) {
        exactCallLogger = new ExactCallLogger(exactCallsLog);
    }

    /**
     * Use this logger instead of the default logger
     *
     * @param logger
     */
    public void setLogger(final Logger logger) {
        this.logger = logger;
    }

    /**
     * Set the maximum number of PL values to log.  If the number of PL values exceeds this, no PL values will be logged.
     * @param maxNumPLValues    maximum number of PL values to log
     */
    public AFCalculator setMaxNumPLValues(final int maxNumPLValues) {
        this.maxNumPLValues = maxNumPLValues;
        return this;
    }

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @param log10AlleleFrequencyPriors a prior vector nSamples x 2 in length indicating the Pr(AF = i)
     * @return result (for programming convenience)
     */
    public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles, final double[] log10AlleleFrequencyPriors) {

        if ( vc == null ) throw new IllegalArgumentException("VariantContext cannot be null");
        if ( vc.getNAlleles() == 1 ) throw new IllegalArgumentException("VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);
        if ( log10AlleleFrequencyPriors == null ) throw new IllegalArgumentException("priors vector cannot be null");

        // reset the result, so we can store our new result there
        final StateTracker stateTracker = getStateTracker(true,maximumAlternativeAlleles);

        //TODO All implementations of the reduce-scope seems to employ a bad criterion to
        //TODO decide what alleles to keep. This must be changed eventually.
        //TODO issue {@see https://github.com/broadinstitute/gsa-unstable/issues/1376}
        final VariantContext vcWorking = reduceScope(vc,defaultPloidy, maximumAlternativeAlleles);

        callTimer.start();
        final AFCalculationResult result = computeLog10PNonRef(vcWorking, defaultPloidy, log10AlleleFrequencyPriors, stateTracker);
        final long nanoTime = callTimer.getElapsedTimeNano();

        if ( exactCallLogger != null )
            exactCallLogger.printCallInfo(vcWorking, log10AlleleFrequencyPriors, nanoTime, result);

        return result;
    }

    /**
     * Convert the final state of the state tracker into our result as an AFCalculationResult
     *
     * Assumes that stateTracker has been updated accordingly
     *
     * @param vcWorking the VariantContext we actually used as input to the calc model (after reduction)
     * @param log10AlleleFrequencyPriors the priors by AC vector
     * @return a AFCalculationResult describing the result of this calculation
     */
    @Requires("stateTracker.getnEvaluations() >= 0")
    @Ensures("result != null")
    protected AFCalculationResult getResultFromFinalState(final VariantContext vcWorking, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        stateTracker.setAllelesUsedInGenotyping(vcWorking.getAlleles());
        return stateTracker.toAFCalculationResult(log10AlleleFrequencyPriors);
    }

    // ---------------------------------------------------------------------------
    //
    // Abstract methods that should be implemented by concrete implementations
    // to actually calculate the AF
    //
    // ---------------------------------------------------------------------------

    /**
     * Look at VC and perhaps return a new one of reduced complexity, if that's necessary
     *
     * Used before the call to computeLog10PNonRef to simply the calculation job at hand,
     * if vc exceeds bounds.  For example, if VC has 100 alt alleles this function
     * may decide to only genotype the best 2 of them.
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    @Requires({"vc != null", "vc.getNAlleles() > 1"})
    @Ensures("result != null")
    protected abstract VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles);

    /**
     * Actually carry out the log10PNonRef calculation on vc, storing results in results
     *
     * @param vc                                variant context with alleles and genotype likelihoods,
     *                                          must have at least one alt allele
     * @param log10AlleleFrequencyPriors        priors
     * @return a AFCalcResult object describing the results of this calculation
     */
    @Requires({"vc != null", "log10AlleleFrequencyPriors != null", "vc.getNAlleles() > 1"})
    protected abstract AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                        final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker);

    /**
     * Subset VC to the just allelesToUse, updating genotype likelihoods
     *
     * Must be overridden by concrete subclasses
     *
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     default ploidy to assume in case {@code vc} does not indicate it for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes
     * @return GenotypesContext object
     */
    public abstract GenotypesContext subsetAlleles(final VariantContext vc,
                                                   final int defaultPloidy,
                                                   final List<Allele> allelesToUse,
                                                   final boolean assignGenotypes);

    // ---------------------------------------------------------------------------
    //
    // accessors
    //
    // ---------------------------------------------------------------------------


    /**
     * Retrieves the state tracker.
     *
     * <p>
     *     The tracker will be reset if so requested or if it needs to be resized due to an increase in the
     *     maximum number of alleles is must be able to handle.
     * </p>
     *
     * @param reset make sure the tracker is reset.
     * @param maximumAlternativeAlleleCount the maximum alternative allele count it must be able to handle. Has no effect if
     *                                     the current tracker is able to handle that number.
     *
     * @return never {@code null}
     */
    protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) {
        if (stateTracker == null)
            stateTracker = new StateTracker(maximumAlternativeAlleleCount);
        else if (reset)
            stateTracker.reset(maximumAlternativeAlleleCount);
        else
            stateTracker.ensureMaximumAlleleCapacity(maximumAlternativeAlleleCount);
        return stateTracker;
    }

    /**
     * Used by testing code.
     *
     * Please don't use this method in production.
     *
     * @deprecated
     */
    @Deprecated
    protected int getAltAlleleCountOfMAP(final int allele) {
        return getStateTracker(false,allele + 1).getAlleleCountsOfMAP()[allele];
    }

    /**
     * Strips PLs from the specified GenotypeBuilder if their number exceeds the maximum allowed.  Corresponding counters are updated.
     * @param gb                the GenotypeBuilder to modify
     * @param vc                the VariantContext
     * @param sampleName        the sample name
     * @param newLikelihoods    the PL array
     */
    protected void removePLsIfMaxNumPLValuesExceeded(final GenotypeBuilder gb, final VariantContext vc, final String sampleName, final double[] newLikelihoods) {
        final int numPLValuesFound = newLikelihoods.length;
        if (numPLValuesFound > maxNumPLValues) {
            logMaxNumPLValuesWarning(vc, sampleName, numPLValuesFound);
            numTimesMaxNumPLValuesExceeded++;
            gb.noPL();
            if (numPLValuesFound > maxNumPLValuesObserved) {
                maxNumPLValuesObserved = numPLValuesFound;
            }
        }
    }

    private void logMaxNumPLValuesWarning(final VariantContext vc, final String sampleName, final int numPLValuesFound) {
        final String message = String.format("Maximum allowed number of PLs (%d) exceeded for sample %s at %s:%d-%d with %d possible genotypes. " +
                        "No PLs will be output for these genotypes (which may cause incorrect results in subsequent analyses) " +
                        "unless the --max_num_PL_values argument is increased accordingly",
                maxNumPLValues, sampleName, vc.getContig(), vc.getStart(), vc.getEnd(), numPLValuesFound);

        if ( numTimesMaxNumPLValuesExceeded == 0 ) {
            logger.warn(message + ". Unless the DEBUG logging level is used, this warning message is output just once per run and further warnings are suppressed.");
        } else {
            logger.debug(message);
        }
    }

    /**
     * Logs the number of times the maximum allowed number of PLs was exceeded and the largest number of PLs observed. The corresponding counters are reset.
     */
    public void printFinalMaxNumPLValuesWarning() {
        if ( numTimesMaxNumPLValuesExceeded > 0 ) {
            final String message = String.format("Maximum allowed number of PLs (%d) was exceeded %d time(s); the largest number of PLs found was %d. " +
                            "No PLs will be output for these genotypes (which may cause incorrect results in subsequent analyses) " +
                            "unless the --max_num_PL_values argument is increased accordingly",
                    maxNumPLValues, numTimesMaxNumPLValuesExceeded, maxNumPLValuesObserved);
            logger.warn(message);
        }
        maxNumPLValuesObserved = 0;
        numTimesMaxNumPLValuesExceeded = 0;
    }
}