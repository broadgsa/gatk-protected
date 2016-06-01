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

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.GenotypeCalculationArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.genotyper.StandardCallerArgumentCollection;

/**
 * A single fixed instance AF calculator provider.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class FixedAFCalculatorProvider extends AFCalculatorProvider {

    private final AFCalculator singleton;

    private final boolean verifyRequests;

    private final int maximumAltAlleleCount;

    private final int maximumNumPLValues;

    private final int ploidy;

    /**
     * Constructs a fixed AF Calculator provider.
     * @param configuration the called configuration. This is the source of the fixed ploidy and maximum number of
     *                      supported alleles.
     * @param logger the logger to be used by the AF calculator instance.
     * @param verifyRequests whether this provider will verify that each request for the AF calculator meets the
     *                       initial parameter values (ploidy, sample-count and maximum number of alleles.
     *
     * @throws NullPointerException if {@code configuration} is {@code null}, or it contains invalid values for
     *    sample ploidy and maximum number of alternative alleles, or {@code sampleCount} is less than 0.
     */
    public FixedAFCalculatorProvider(final StandardCallerArgumentCollection configuration,
                                     final Logger logger, final boolean verifyRequests) {
        this(configuration.requestedAlleleFrequencyCalculationModel,configuration.genotypeArgs,logger,verifyRequests);

    }


    /**
     * Constructs a fixed AF Calculator provider.
     *
     * @param configuration the called configuration. This is the source of the fixed ploidy and maximum number of
     *                      supported alleles.
     * @param logger the logger to be used by the AF calculator instance.
     * @param verifyRequests whether this provider will verify that each request for the AF calculator meets the
     *                       initial parameter values (ploidy, sample-count and maximum number of alleles.
     *
     * @throws IllegalArgumentException if {@code configuration} is {@code null}, or it contains invalid values for
     *    sample ploidy and maximum number of alternative alleles, or {@code sampleCount} is less than 0.
     */
    public FixedAFCalculatorProvider(final GenotypeCalculationArgumentCollection configuration,
                                     final Logger logger, final boolean verifyRequests) {
        this(null,configuration,logger,verifyRequests);
    }

    /**
     * Constructs a fixed AF Calculator provider.
     *
     * @param preferred preferred implementation.
     * @param configuration the called configuration. This is the source of the fixed ploidy and maximum number of
     *                      supported alleles.
     * @param logger the logger to be used by the AF calculator instance.
     * @param verifyRequests whether this provider will verify that each request for the AF calculator meets the
     *                       initial parameter values (ploidy, sample-count and maximum number of alleles.
     *
     * @throws IllegalArgumentException if {@code configuration} is {@code null}, or it contains invalid values for
     *    sample ploidy and maximum number of alternative alleles, or {@code sampleCount} is less than 0.
     */
    public FixedAFCalculatorProvider(final AFCalculatorImplementation preferred, final GenotypeCalculationArgumentCollection configuration,
                                      final Logger logger, final boolean verifyRequests) {

        if (configuration == null)
            throw new IllegalArgumentException("null genotype-arguments configuration");
        if (configuration.samplePloidy < 1)
            throw new IllegalArgumentException("invalid sample ploidy " + configuration.samplePloidy);
        if (configuration.MAX_ALTERNATE_ALLELES < 0)
            throw new IllegalArgumentException("invalid maximum number of alleles " + (configuration.MAX_ALTERNATE_ALLELES + 1));
        if (configuration.MAX_NUM_PL_VALUES < 0)
            throw new IllegalArgumentException("invalid maximum number of PL values " + configuration.MAX_NUM_PL_VALUES);

        ploidy = configuration.samplePloidy;
        maximumAltAlleleCount = configuration.MAX_ALTERNATE_ALLELES;
        maximumNumPLValues = configuration.MAX_NUM_PL_VALUES;
        singleton = AFCalculatorImplementation.bestValue(ploidy,maximumAltAlleleCount,preferred).newInstance().setMaxNumPLValues(maximumNumPLValues);
        singleton.setLogger(logger);
        this.verifyRequests = verifyRequests;
    }

    @Override
    public AFCalculator getInstance(final VariantContext vc, final int defaultPloidy, final int maximumAlleleCount) {
        if (verifyRequests)
            // supers implementation will call eventually one of the other methods, so no need to verify anything here.
            return super.getInstance(vc,defaultPloidy,maximumAlleleCount);
        return singleton;
    }

    @Override
    public AFCalculator getInstance(final int ploidy, final int maxAltAlleleCount) {
        if (verifyRequests) {
            if (this.ploidy != AFCalculatorImplementation.UNBOUND_PLOIDY && ploidy != this.ploidy)
                throw new IllegalStateException("non-supported ploidy");
            if (maximumAltAlleleCount != AFCalculatorImplementation.UNBOUND_ALTERNATIVE_ALLELE_COUNT && maxAltAlleleCount > maximumAltAlleleCount)
                throw new IllegalStateException("non-supported alleleCount");
        }
        return singleton;
    }

    /**
     * Creates a fixed AF calculator provider that is thread safe given the engine configuration.
     *
     * @param toolkit the engine.
     * @param config the caller configuration.
     * @param logger reference to the logger for the AF calculator to dump messages to.
     *
     * @throws IllegalArgumentException if any of the input argument is {@code null} or contain invalid configuration
     *   like zero-samples, zero or negative ploidy or negative-zero maximum number of alleles.
     *
     * @return never {@code null}
     */
    public static AFCalculatorProvider createThreadSafeProvider(final GenomeAnalysisEngine toolkit,
                                                                final StandardCallerArgumentCollection config,
                                                                final Logger logger) {
        final GATKArgumentCollection arguments = toolkit.getArguments();
        final boolean isMultithread = arguments.numberOfCPUThreadsPerDataThread > 1 || arguments.numberOfDataThreads > 1;
        return !isMultithread
                ? new FixedAFCalculatorProvider(config,logger,false) :
                new ConcurrentAFCalculatorProvider() {
                    @Override
                    protected AFCalculatorProvider createProvider() {
                        return new FixedAFCalculatorProvider(config,logger,false);
                    }
                };
    }

    /**
     * Creates a fixed AF calculator provider that is thread safe given the engine configuration.
     *
     * @param toolkit the engine.
     * @param config the caller configuration.
     * @param logger reference to the logger for the AF calculator to dump messages to.
     * @param verifyRequests whether each request should be check for a compatible ploidy and max-alt-allele count.
     *       A non-compliant request would result in an exception when requesting the AF calculator instance.
     *
     * @throws IllegalArgumentException if any of the input argument is {@code null} or contain invalid configuration
     *   like zero-samples, zero or negative ploidy or negative-zero maximum number of alleles.
     *
     * @return never {@code null}.
     */
    @SuppressWarnings("unused")
    public static AFCalculatorProvider createThreadSafeProvider(final GenomeAnalysisEngine toolkit,
                                                                  final StandardCallerArgumentCollection config,
                                                                  final Logger logger, final boolean verifyRequests) {
        final GATKArgumentCollection arguments = toolkit.getArguments();
        final boolean isMultithread = arguments.numberOfCPUThreadsPerDataThread > 1 || arguments.numberOfDataThreads > 1;
        return !isMultithread
                ? new FixedAFCalculatorProvider(config,logger,verifyRequests) :
                new ConcurrentAFCalculatorProvider() {
                    @Override
                    protected AFCalculatorProvider createProvider() {
                        return new FixedAFCalculatorProvider(config,logger,verifyRequests);
                    }
                };
    }


    /**
     * Creates a fixed AF calculator provider that is thread safe given the engine configuration.
     *
     * @param toolkit the engine.
     * @param config the caller configuration.
     * @param logger reference to the logger for the AF calculator to dump messages to.
     *
     * @throws IllegalArgumentException if any of the input argument is {@code null} or contain invalid configuration
     *   like zero-samples, zero or negative ploidy or negative-zero maximum number of alleles.
     *
     * @return never {@code null}
     */
    public static AFCalculatorProvider createThreadSafeProvider(final GenomeAnalysisEngine toolkit,
                                                                final GenotypeCalculationArgumentCollection config,
                                                                        final Logger logger, final boolean verifyRequests) {
        return toolkit.getArguments().numberOfCPUThreadsPerDataThread <= 1
                ? new FixedAFCalculatorProvider(config,logger,false) :
                new ConcurrentAFCalculatorProvider() {
                    @Override
                    protected AFCalculatorProvider createProvider() {
                        return new FixedAFCalculatorProvider(config,logger,verifyRequests);
                    }
                };
    }
}
