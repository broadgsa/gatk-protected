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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.baq.BAQ;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.gga.GenotypingGivenAllelesUtils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.*;

/**
 * {@link UnifiedGenotyper}'s genotyping strategy implementation.
 */
public class UnifiedGenotypingEngine extends GenotypingEngine<UnifiedArgumentCollection> {

    private static final int SNP_MODEL = 0;
    private static final int INDEL_MODEL = 1;

    // the model used for calculating genotypes
    private ThreadLocal<Map<String, GenotypeLikelihoodsCalculationModel>> glcm;
    private final List<GenotypeLikelihoodsCalculationModel.Model> modelsToUse = new ArrayList<>(2);


    // the various loggers and writers
    private PrintStream verboseWriter;

    private final boolean BAQEnabledOnCMDLine;

    private boolean doAlleleSpecificCalcs = false;

    // ---------------------------------------------------------------------------------------------------------
    //
    // Public interface functions
    //
    // ---------------------------------------------------------------------------------------------------------


    /**
     * Creates a new unified genotyping given the UG configuration parameters and the GA engine.
     * @param configuration the UG configuration.
     * @param toolkit the GA engine.
     *
     * @throws NullPointerException if either {@code configuration} or {@code toolkit} is {@code null}.
     */
    public UnifiedGenotypingEngine(final UnifiedArgumentCollection configuration, final AFCalculatorProvider afCalculatorProvider,
                                   final GenomeAnalysisEngine toolkit) {
        this(configuration,toolkit.getSampleList(),toolkit.getGenomeLocParser(),afCalculatorProvider,toolkit.getArguments().BAQMode);
    }

    /**
     * Creates a new unified genotyping given the UG configuration parameters, the targeted set of samples and
     * a genome location parser.
     *
     * @param configuration the UG configuration.
     * @param samples {@inheritDoc}
     * @param baqCalculationMode the BAQ calculation mode.
     *
     * @throws NullPointerException if any of {@code configuration}, {@code samples} or {@code genomeLocParser} is {@code null}.
     *
     * @throws IllegalArgumentException if {@code baqCalculationMode} is {@code null}.
     */
    public UnifiedGenotypingEngine(final UnifiedArgumentCollection configuration,
                                   final SampleList samples, final GenomeLocParser genomeLocParser, final AFCalculatorProvider afCalculatorProvider,
                                   final BAQ.CalculationMode baqCalculationMode) {
        this(configuration, samples, genomeLocParser, afCalculatorProvider, baqCalculationMode, false);
    }


    /**
     * Creates a new unified genotyping given the UG configuration parameters, the targeted set of samples and
     * a genome location parser.
     *
     * @param configuration the UG configuration.
     * @param samples {@inheritDoc}
     * @param baqCalculationMode the BAQ calculation mode.
     *
     * @throws NullPointerException if any of {@code configuration}, {@code samples} or {@code genomeLocParser} is {@code null}.
     *
     * @throws IllegalArgumentException if {@code baqCalculationMode} is {@code null}.
     */
    public UnifiedGenotypingEngine(final UnifiedArgumentCollection configuration,
                                   final SampleList samples, final GenomeLocParser genomeLocParser, final AFCalculatorProvider afCalculatorProvider,
                                   final BAQ.CalculationMode baqCalculationMode,
                                   final boolean doAlleleSpecificCalcs) {

        super(configuration,samples,genomeLocParser,afCalculatorProvider, doAlleleSpecificCalcs);

        if (baqCalculationMode == null)
            throw new IllegalArgumentException("the BAQ calculation mode cannot be null");

        this.BAQEnabledOnCMDLine = baqCalculationMode != BAQ.CalculationMode.OFF;
        this.doAlleleSpecificCalcs = doAlleleSpecificCalcs;
        determineGLModelsToUse();
        initializeGenotypeLikelihoodsCalculationModels();
    }



    /**
     * Changes the verbose output writer for this engine.
     *
     * @param writer the new writer; it can be {@code null}.
     */
    public void setVerboseWriter(final PrintStream writer) {
        verboseWriter = writer;
    }

    /**
     * Initialize {@link #glcm}.
     */
    private void initializeGenotypeLikelihoodsCalculationModels() {
        glcm = new ThreadLocal<Map<String, GenotypeLikelihoodsCalculationModel>>() {

            @Override
            protected Map<String,GenotypeLikelihoodsCalculationModel> initialValue() {
                return getGenotypeLikelihoodsCalculationObject(logger,UnifiedGenotypingEngine.this.configuration);
            }
        };
    }

    /**
     * Compute full calls at a given locus. Entry point for engine calls from the UnifiedGenotyper.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public List<VariantCallContext> calculateLikelihoodsAndGenotypes(final RefMetaDataTracker tracker,
                                                                     final ReferenceContext refContext,
                                                                     final AlignmentContext rawContext) {
        final List<VariantCallContext> results = new ArrayList<>(2);

        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, rawContext);

        final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = new HashMap<>();

        final VariantCallContext defaultResult = configuration.outputMode == OutputMode.EMIT_ALL_SITES
                && configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
                ? emptyCallContext(tracker, refContext, rawContext)
                : null;

        if ( models.isEmpty() )
            results.add(defaultResult);
        else {
            for ( final GenotypeLikelihoodsCalculationModel.Model model : models ) {
                perReadAlleleLikelihoodMap.clear();
                final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(refContext, rawContext, model);
                if ( stratifiedContexts == null )
                    results.add(defaultResult);
                else {
                    final VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model, perReadAlleleLikelihoodMap);
                    if ( vc != null )
                        results.add(calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, true, perReadAlleleLikelihoodMap));
// todo - uncomment if we want to also emit a null ref call (with no QUAL) if there's no evidence for REF and if EMIT_ALL_SITES is set
//                    else if (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES)
//                        results.add(generateEmptyContext(tracker, refContext, null, rawContext));

                }
            }        
        }

        return results;
    }



    private static Map<String, GenotypeLikelihoodsCalculationModel> getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {

        final Map<String, GenotypeLikelihoodsCalculationModel> glcm = new HashMap<>();
        final List<Class<? extends GenotypeLikelihoodsCalculationModel>> glmClasses = new PluginManager<GenotypeLikelihoodsCalculationModel>(GenotypeLikelihoodsCalculationModel.class).getPlugins();

        for (final Class<? extends GenotypeLikelihoodsCalculationModel> glmClass : glmClasses) {
            final String key = glmClass.getSimpleName().replaceAll("GenotypeLikelihoodsCalculationModel","").toUpperCase();
            try {
                final Object args[] = new Object[]{UAC,logger};
                final Constructor c = glmClass.getDeclaredConstructor(UnifiedArgumentCollection.class, Logger.class);
                glcm.put(key, (GenotypeLikelihoodsCalculationModel)c.newInstance(args));
            }
            catch (Exception e) {
                throw new UserException("The likelihoods model provided for the -glm argument (" + UAC.GLmodel + ") is not a valid option: " + e.getMessage());
            }
        }

        return glcm;
    }

    /**
     * Compute GLs at a given locus. Entry point for engine calls from UGCalcLikelihoods.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param perReadAlleleLikelihoodMap    Map to store per-sample, per-read, per-allele likelihoods (only used for indels)
     * @return the VariantContext object
     */
    public VariantContext calculateLikelihoods(final RefMetaDataTracker tracker,
                                               final ReferenceContext refContext,
                                               final AlignmentContext rawContext,
                                               final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, rawContext);
        if ( models.isEmpty() ) {
            return null;
        }

        for ( final GenotypeLikelihoodsCalculationModel.Model model : models ) {
            final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(refContext, rawContext, model);
            // return the first valid one we encounter
            if ( stratifiedContexts != null )
                return calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model, perReadAlleleLikelihoodMap);

        }

        return null;
    }

    /**
     * Compute genotypes at a given locus. Entry point for engine calls from UGCallVariants.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param vc         the GL-annotated variant context
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker,
                                                 final ReferenceContext refContext,
                                                 final AlignmentContext rawContext,
                                                 final VariantContext vc) {
        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, rawContext);
        if ( models.isEmpty() ) {
            return null;
        }

        // return the first one
        final GenotypeLikelihoodsCalculationModel.Model model = models.get(0);
        final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(refContext, rawContext, model);
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, null, false);
    }

    /**
     * Compute genotypes at a given locus.
     *
     * @param vc         the GL-annotated variant context
     * @return the VariantCallContext object (can be null)
     */
    public VariantCallContext calculateGenotypes(VariantContext vc) {
        final VariantContext.Type type = vc.getType();
        final GenotypeLikelihoodsCalculationModel.Model model;
        /**
         * Query the VariantContext for the appropriate model.  If type == MIXED, one would want to use model = BOTH.
         * However GenotypingEngine.getAlleleFrequencyPriors throws an exception if you give it anything but a SNP or INDEL model.
         */
        if ( type == VariantContext.Type.INDEL ) {
            model = GenotypeLikelihoodsCalculationModel.Model.INDEL;
        } else {
            model = GenotypeLikelihoodsCalculationModel.Model.SNP;
        }

        return calculateGenotypes(null, null, null, null, vc, model, null, doAlleleSpecificCalcs);
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Private implementation helpers
    //
    // ---------------------------------------------------------------------------------------------------------

    // private method called by both UnifiedGenotyper and UGCalcLikelihoods entry points into the engine
    private VariantContext calculateLikelihoods(final RefMetaDataTracker tracker,
                                                final ReferenceContext refContext,
                                                final Map<String, AlignmentContext> stratifiedContexts,
                                                final AlignmentContextUtils.ReadOrientation type,
                                                final List<Allele> alternateAllelesToUse,
                                                final boolean useBAQedPileup,
                                                final GenotypeLikelihoodsCalculationModel.Model model,
                                                final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        return glcm.get().get(model.name()).getLikelihoods(tracker, refContext, stratifiedContexts, type, alternateAllelesToUse, useBAQedPileup && BAQEnabledOnCMDLine,
                genomeLocParser != null || refContext == null ? genomeLocParser : refContext.getGenomeLocParser(), perReadAlleleLikelihoodMap);
    }


    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        return calculateGenotypes(null, null, null, null, vc, model, null, false);
    }

    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                 final AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final boolean inheritAttributesFromInputVC,
                                                 final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, inheritAttributesFromInputVC, perReadAlleleLikelihoodMap, false);
    }

    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker,
                                                 final ReferenceContext refContext,
                                                 final AlignmentContext rawContext,
                                                 final Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc,
                                                 final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                                 final boolean useAlleleSpecificCalcs) {
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, false, perReadAlleleLikelihoodMap, useAlleleSpecificCalcs);
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
                || configuration.annotateAllSitesWithPLs;
    }

    @Override
    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                 final AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final boolean inheritAttributesFromInputVC,
                                                 final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                                 final boolean useAlleleSpecificCalcs) {
        boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;
        final VariantCallContext result = super.calculateGenotypes(tracker,refContext,rawContext,stratifiedContexts,vc,model,inheritAttributesFromInputVC,perReadAlleleLikelihoodMap, useAlleleSpecificCalcs);
        if ( verboseWriter != null && !limitedContext )
            printVerboseData(refContext.getLocus().toString(), vc, model);
        return result;
    }

    @Override
    protected String callSourceString() {
        return "UG_call";
    }

    @Override
    protected boolean forceSiteEmission() {
        return configuration.outputMode == OutputMode.EMIT_ALL_SITES;
    }


    @Override
    protected Map<String,Object> composeCallAttributes(final boolean inheritAttributesFromInputVC, final VariantContext vc,
                                                       final AlignmentContext rawContext, final Map<String, AlignmentContext> stratifiedContexts, final RefMetaDataTracker tracker, final ReferenceContext refContext, final List<Integer> alleleCountsofMLE, final boolean bestGuessIsRef,
                                                       final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes,
                                                       final GenotypeLikelihoodsCalculationModel.Model model, final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                                       final boolean useAlleleSpecificCalcs) {
        final Map<String,Object> result = super.composeCallAttributes(inheritAttributesFromInputVC, vc,rawContext,stratifiedContexts,tracker,refContext,alleleCountsofMLE,bestGuessIsRef,
                                    AFresult,allAllelesToUse,genotypes,model,perReadAlleleLikelihoodMap, useAlleleSpecificCalcs);

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        if ( configuration.COMPUTE_SLOD && !limitedContext && !bestGuessIsRef ) {
            final double strandScore = calculateSLOD(stratifiedContexts, tracker, refContext, AFresult, allAllelesToUse, model, perReadAlleleLikelihoodMap);
            if ( !Double.isNaN(strandScore) )
                result.put("SB", strandScore);
        }
        return result;
    }

    private double calculateSLOD(final Map<String, AlignmentContext> stratifiedContexts,
                                 final RefMetaDataTracker tracker,
                                 final ReferenceContext refContext, final AFCalculationResult AFresult,
                                 final List<Allele> allAllelesToUse,
                                 final GenotypeLikelihoodsCalculationModel.Model model,
                                 final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        // the overall lod
        //double overallLog10PofNull = AFresult.log10AlleleFrequencyPosteriors[0];
        final double overallLog10PofF = AFresult.getLog10LikelihoodOfAFGT0();
        //if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

        // the forward lod
        final AFCalculationResult forwardAFresult = getDirectionalAfCalcResult(AlignmentContextUtils.ReadOrientation.FORWARD,stratifiedContexts, tracker, refContext, allAllelesToUse, model, perReadAlleleLikelihoodMap);
        final double forwardLog10PofNull = forwardAFresult.getLog10LikelihoodOfAFEq0();
        final double forwardLog10PofF = forwardAFresult.getLog10LikelihoodOfAFGT0();
        //if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

        // the reverse lod
        final AFCalculationResult reverseAFresult = getDirectionalAfCalcResult(AlignmentContextUtils.ReadOrientation.REVERSE,stratifiedContexts, tracker, refContext, allAllelesToUse, model, perReadAlleleLikelihoodMap);
        final double reverseLog10PofNull = reverseAFresult.getLog10LikelihoodOfAFEq0();
        final double reverseLog10PofF = reverseAFresult.getLog10LikelihoodOfAFGT0();
        //if ( DEBUG_SLOD ) System.out.println("reverseLog10PofNull=" + reverseLog10PofNull + ", reverseLog10PofF=" + reverseLog10PofF);

        final double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
        final double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
        //if ( DEBUG_SLOD ) System.out.println("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

        // strand score is max bias between forward and reverse strands
        double strandScore = Math.max(forwardLod, reverseLod);
        // rescale by a factor of 10
        strandScore *= 10.0;
        //logger.debug(String.format("SLOD=%f", strandScore));
        return strandScore;
    }

    private AFCalculationResult getDirectionalAfCalcResult(final AlignmentContextUtils.ReadOrientation orientation,
                                                    final Map<String, AlignmentContext> stratifiedContexts,
                                                    final RefMetaDataTracker tracker,
                                                    final ReferenceContext refContext, List<Allele> allAllelesToUse,
                                                    final GenotypeLikelihoodsCalculationModel.Model model,
                                                    final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts,orientation,
                allAllelesToUse, false, model, perReadAlleleLikelihoodMap);
        final int defaultPloidy = configuration.genotypeArgs.samplePloidy;
        final int maxAltAlleles = configuration.genotypeArgs.MAX_ALTERNATE_ALLELES;
        final int maxNumPLValues = configuration.genotypeArgs.MAX_NUM_PL_VALUES;
        return afCalculatorProvider.getInstance(vc, defaultPloidy, maxAltAlleles).setMaxNumPLValues(maxNumPLValues)
                .getLog10PNonRef(vc, defaultPloidy, maxAltAlleles, getAlleleFrequencyPriors(vc,defaultPloidy,model));
    }

    private Map<String, AlignmentContext> getFilteredAndStratifiedContexts(final ReferenceContext refContext,
                                                                           final AlignmentContext rawContext,
                                                                           final GenotypeLikelihoodsCalculationModel.Model model) {

        if ( !BaseUtils.isRegularBase(refContext.getBase()) )
            return null;


        switch (model) {
            case INDEL:
            case GENERALPLOIDYINDEL:

                final ReadBackedPileup pileup = rawContext.getBasePileup().getMappingFilteredPileup(configuration.MIN_BASE_QUALTY_SCORE);
                // don't call when there is no coverage
                if ( pileup.getNumberOfElements() == 0 && configuration.outputMode != OutputMode.EMIT_ALL_SITES  )
                    return null;

                // stratify the AlignmentContext and cut by sample
                return AlignmentContextUtils.splitContextBySampleName(pileup);
            case SNP:
            case GENERALPLOIDYSNP:

                if ( !(configuration.outputMode == OutputMode.EMIT_ALL_SITES && configuration.genotypingOutputMode != GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) ) {
                    int numDeletions = 0;
                    for ( final PileupElement p : rawContext.getBasePileup() ) {
                        if ( p.isDeletion() )
                            numDeletions++;
                    }
                    if ( ((double) numDeletions) / ((double) rawContext.getBasePileup().depthOfCoverage()) > configuration.MAX_DELETION_FRACTION ) {
                        return null;
                    }
                }
                // stratify the AlignmentContext and cut by sample
                return AlignmentContextUtils.splitContextBySampleName(rawContext.getBasePileup());
            default :
                throw new IllegalStateException("unexpected model: " + model);
        }
    }

    protected void printVerboseData(final String pos, final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        Allele refAllele = null, altAllele = null;
        for ( Allele allele : vc.getAlleles() ) {
            if ( allele.isReference() )
                refAllele = allele;
            else
                altAllele = allele;
        }

        for (int i = 0; i <= numberOfGenomes; i++) {
            StringBuilder AFline = new StringBuilder("AFINFO\t");
            AFline.append(pos);
            AFline.append('\t');
            AFline.append(refAllele);
            AFline.append('\t');
            if ( altAllele != null )
                AFline.append(altAllele);
            else
                AFline.append("N/A");
            AFline.append('\t');
            AFline.append(i).append('/').append(numberOfGenomes).append('\t');
            AFline.append(String.format("%.2f\t", ((float) i) / numberOfGenomes));
            AFline.append(String.format("%.8f\t", getAlleleFrequencyPriors(vc,configuration.genotypeArgs.samplePloidy,model)[i]));
            verboseWriter.println(AFline.toString());
        }

        verboseWriter.println("Qscore = " + vc.getLog10PError());
        verboseWriter.println();
    }

    /**
     * Determine the models to be use for genotype likelihood calculation.
     *
     * <p>
     * Whether to use the general ones or the concrete diploid ones need to depend on what the user has said
     * in its parameter glm. Iff he selected GENERALPLOIDYINDEL or GENERALPLOIDYSNP is the general set, otherwise
     * </p>
     *
     * the standard set (SNP and INDEL).
     * <p>
     *     Even if the user did not select to use both models, GGA force the inclusion of both: snp and indels.
     * </p>
     * <p>
     *     Also, we must fail
     * </p>
     *
     * The models are added to the field {@link #modelsToUse}.
     */
    private void determineGLModelsToUse() {

        modelsToUse.clear();

        boolean useGeneral = false;
        boolean useSNP = false;
        boolean useINDEL = false;

        switch (configuration.GLmodel) {
            case BOTH:
                useSNP = useINDEL = true; break;
            case SNP:
                useSNP = true; break;
            case INDEL:
                useINDEL = true; break;
            case GENERALPLOIDYINDEL:
                useINDEL = useGeneral = true; break;
            case GENERALPLOIDYSNP:
                useSNP = useGeneral = true; break;
            default: //Paranoia
                throw new IllegalStateException("unexpected genotype likelihood model " + configuration.GLmodel);
        }

        // Force the use of both models in GGA:
        if (configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
            useSNP = useINDEL = true;

        // The non-general models only support Diploid so need to go to general if not the default_ploidy == 2.
        if (configuration.genotypeArgs.samplePloidy != GATKVariantContextUtils.DEFAULT_PLOIDY)
            useGeneral = true;

        // If annotateAllSitesWithPLs requested , SNP model must be used.
        if (!useSNP && configuration.annotateAllSitesWithPLs)
            throw new UserException.BadArgumentValue("glm","Invalid genotype likelihood model specification: " +
                    "only diploid SNP model can be used in conjunction with option allSitePLs");

        // Finally add the relevant model
        if (useSNP)
            modelsToUse.add(useGeneral ? GenotypeLikelihoodsCalculationModel.Model.GENERALPLOIDYSNP :
                        GenotypeLikelihoodsCalculationModel.Model.SNP);
        if (useINDEL)
            modelsToUse.add(useGeneral ? GenotypeLikelihoodsCalculationModel.Model.GENERALPLOIDYINDEL :
                    GenotypeLikelihoodsCalculationModel.Model.INDEL);
    }

    // decide whether we are currently processing SNPs, indels, neither, or both
    private List<GenotypeLikelihoodsCalculationModel.Model> getGLModelsToUse(final RefMetaDataTracker tracker,
                                                                             final AlignmentContext rawContext) {
        if ( configuration.genotypingOutputMode != GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES )
            return modelsToUse;

        if ( modelsToUse.size() != 2 )
            throw new IllegalStateException("GGA mode assumes that we have initialized both the SNP and indel models but found " + modelsToUse);

        // if we're genotyping given alleles then we need to choose the model corresponding to the variant type requested
        final VariantContext vcInput = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(tracker, rawContext.getLocation(), false, logger, configuration.alleles);

        if ( vcInput == null ) {
            return Collections.emptyList(); // no work to be done
        } else if ( vcInput.isSNP() )  {
            return Collections.singletonList(modelsToUse.get(SNP_MODEL));
        } else if ( vcInput.isIndel() || vcInput.isMixed() ) {
            return Collections.singletonList(modelsToUse.get(INDEL_MODEL));
        } else {
            return Collections.emptyList(); // No support for other types yet
        }
    }


}
