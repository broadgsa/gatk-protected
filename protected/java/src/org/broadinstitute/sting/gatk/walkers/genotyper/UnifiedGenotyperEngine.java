/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.genotyper;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalc;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcFactory;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcResult;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.variant.variantcontext.*;

import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.*;

public class UnifiedGenotyperEngine {
    public static final String LOW_QUAL_FILTER_NAME = "LowQual";
    private static final String GPSTRING = "GENERALPLOIDY";

    public static final String NUMBER_OF_DISCOVERED_ALLELES_KEY = "NDA";

    public static final double HUMAN_SNP_HETEROZYGOSITY = 1e-3;
    public static final double HUMAN_INDEL_HETEROZYGOSITY = 1e-4;

    public enum OUTPUT_MODE {
        /** produces calls only at variant sites */
        EMIT_VARIANTS_ONLY,
        /** produces calls at variant sites and confident reference sites */
        EMIT_ALL_CONFIDENT_SITES,
        /** produces calls at any callable site regardless of confidence; this argument is intended only for point
         * mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by
         * no means produce a comprehensive set of indels in DISCOVERY mode */
        EMIT_ALL_SITES
    }

    // the unified argument collection
    private final UnifiedArgumentCollection UAC;
    public UnifiedArgumentCollection getUAC() { return UAC; }

    // the annotation engine
    private final VariantAnnotatorEngine annotationEngine;

    // the model used for calculating genotypes
    private ThreadLocal<Map<String, GenotypeLikelihoodsCalculationModel>> glcm = new ThreadLocal<Map<String, GenotypeLikelihoodsCalculationModel>>();
    private final List<GenotypeLikelihoodsCalculationModel.Model> modelsToUse = new ArrayList<GenotypeLikelihoodsCalculationModel.Model>(2);

    // the model used for calculating p(non-ref)
    private ThreadLocal<AFCalc> afcm = new ThreadLocal<AFCalc>();

    // because the allele frequency priors are constant for a given i, we cache the results to avoid having to recompute everything
    private final double[] log10AlleleFrequencyPriorsSNPs;
    private final double[] log10AlleleFrequencyPriorsIndels;

    // samples in input
    private final Set<String> samples;

    // the various loggers and writers
    private final Logger logger;
    private final PrintStream verboseWriter;

    // number of chromosomes (ploidy * samples) in input
    private final int ploidy;
    private final int N;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);

    private final GenomeLocParser genomeLocParser;
    private final boolean BAQEnabledOnCMDLine;

    // ---------------------------------------------------------------------------------------------------------
    //
    // Public interface functions
    //
    // ---------------------------------------------------------------------------------------------------------
    @Requires({"toolkit != null", "UAC != null"})
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        this(toolkit, UAC, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, SampleUtils.getSAMFileSamples(toolkit.getSAMFileHeader()), GATKVariantContextUtils.DEFAULT_PLOIDY);
    }

    protected UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, Set<String> samples, UnifiedArgumentCollection UAC) {
        this(toolkit, UAC, Logger.getLogger(UnifiedGenotyperEngine.class), null, null, samples, GATKVariantContextUtils.DEFAULT_PLOIDY);
    }

    @Requires({"toolkit != null", "UAC != null", "logger != null", "samples != null && samples.size() > 0","ploidy>0"})
    public UnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, Set<String> samples, int ploidy) {
        this.BAQEnabledOnCMDLine = toolkit.getArguments().BAQMode != BAQ.CalculationMode.OFF;
        genomeLocParser = toolkit.getGenomeLocParser();
        this.samples = new TreeSet<String>(samples);
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.UAC = UAC;

        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        this.ploidy = ploidy;
        this.N = samples.size() * ploidy;
        log10AlleleFrequencyPriorsSNPs = new double[N+1];
        log10AlleleFrequencyPriorsIndels = new double[N+1];
        computeAlleleFrequencyPriors(N, log10AlleleFrequencyPriorsSNPs, UAC.heterozygosity);
        computeAlleleFrequencyPriors(N, log10AlleleFrequencyPriorsIndels, UAC.INDEL_HETEROZYGOSITY);

        filter.add(LOW_QUAL_FILTER_NAME);

        determineGLModelsToUse();
    }

    /**
     * @see #calculateLikelihoodsAndGenotypes(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext, java.util.Set)
     *
     * same as the full call but with allSamples == null
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public List<VariantCallContext> calculateLikelihoodsAndGenotypes(final RefMetaDataTracker tracker,
                                                                     final ReferenceContext refContext,
                                                                     final AlignmentContext rawContext) {
        return calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext, null);
    }


    /**
     * Compute full calls at a given locus. Entry point for engine calls from the UnifiedGenotyper.
     *
     * If allSamples != null, then the output variantCallContext is guarenteed to contain a genotype
     * for every sample in allSamples.  If it's null there's no such guarentee.  Providing this
     * argument is critical when the resulting calls will be written to a VCF file.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param allSamples set of all sample names that we might call (i.e., those in the VCF header)
     * @return the VariantCallContext object
     */
    public List<VariantCallContext> calculateLikelihoodsAndGenotypes(final RefMetaDataTracker tracker,
                                                                     final ReferenceContext refContext,
                                                                     final AlignmentContext rawContext,
                                                                     final Set<String> allSamples) {
        final List<VariantCallContext> results = new ArrayList<VariantCallContext>(2);

        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, refContext, rawContext);

        final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = new HashMap<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap>();

        if ( models.isEmpty() ) {
            results.add(UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(tracker, refContext, null, rawContext) : null);
        }
        else {
            for ( final GenotypeLikelihoodsCalculationModel.Model model : models ) {
                perReadAlleleLikelihoodMap.clear();
                final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
                if ( stratifiedContexts == null ) {
                    results.add(UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(tracker, refContext, null, rawContext) : null);
                }
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
                                               final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, refContext, rawContext);
        if ( models.isEmpty() ) {
            return null;
        }

        for ( final GenotypeLikelihoodsCalculationModel.Model model : models ) {
            final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
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
        final List<GenotypeLikelihoodsCalculationModel.Model> models = getGLModelsToUse(tracker, refContext, rawContext);
        if ( models.isEmpty() ) {
            return null;
        }

        // return the first one
        final GenotypeLikelihoodsCalculationModel.Model model = models.get(0);
        final Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, null);
    }

    /**
     * Compute genotypes at a given locus.
     *
     * @param vc         the GL-annotated variant context
     * @return the VariantCallContext object
     */
    public VariantCallContext calculateGenotypes(VariantContext vc) {
        return calculateGenotypes(null, null, null, null, vc, GenotypeLikelihoodsCalculationModel.Model.valueOf("SNP"), null);
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
                                                final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        // initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(getGenotypeLikelihoodsCalculationObject(logger, UAC));
        }

        return glcm.get().get(model.name()).getLikelihoods(tracker, refContext, stratifiedContexts, type, alternateAllelesToUse, useBAQedPileup && BAQEnabledOnCMDLine, genomeLocParser, perReadAlleleLikelihoodMap);
    }

    private VariantCallContext generateEmptyContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, AlignmentContext rawContext) {
        VariantContext vc;
        if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            VariantContext vcInput = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, rawContext.getLocation(), false, logger, UAC.alleles);
            if ( vcInput == null )
                return null;
            vc = new VariantContextBuilder("UG_call", ref.getLocus().getContig(), vcInput.getStart(), vcInput.getEnd(), vcInput.getAlleles()).make();
        } else {
            // deal with bad/non-standard reference bases
            if ( !Allele.acceptableAlleleBases(new byte[]{ref.getBase()}) )
                return null;

            Set<Allele> alleles = new HashSet<Allele>();
            alleles.add(Allele.create(ref.getBase(), true));
            vc = new VariantContextBuilder("UG_call", ref.getLocus().getContig(), ref.getLocus().getStart(), ref.getLocus().getStart(), alleles).make();
        }
        
        if ( annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadBackedPileup pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

            vc = annotationEngine.annotateContext(tracker, ref, stratifiedContexts, vc);
        }

        return new VariantCallContext(vc, false);
    }

    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model, final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return calculateGenotypes(null, null, null, null, vc, model, perReadAlleleLikelihoodMap);
    }

    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        return calculateGenotypes(null, null, null, null, vc, model, null);
    }

    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker,
                                                 final ReferenceContext refContext,
                                                 final AlignmentContext rawContext,
                                                 final Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc,
                                                 final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, vc, model, false,perReadAlleleLikelihoodMap);
    }

    /**
     * Main entry function to calculate genotypes of a given VC with corresponding GL's
     * @param tracker                            Tracker
     * @param refContext                         Reference context
     * @param rawContext                         Raw context
     * @param stratifiedContexts                 Stratified alignment contexts
     * @param vc                                 Input VC
     * @param model                              GL calculation model
     * @param inheritAttributesFromInputVC       Output VC will contain attributes inherited from input vc
     * @return                                   VC with assigned genotypes
     */
    public VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                 final AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final boolean inheritAttributesFromInputVC,
                                                 final Map<String, org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        // initialize the data for this thread if that hasn't been done yet
        if ( afcm.get() == null ) {
            afcm.set(AFCalcFactory.createAFCalc(UAC, N, logger));
        }

        // estimate our confidence in a reference call and return
        if ( vc.getNSamples() == 0 ) {
            if ( limitedContext )
                return null;
            return (UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES ?
                    estimateReferenceConfidence(vc, stratifiedContexts, getTheta(model), false, 1.0) :
                    generateEmptyContext(tracker, refContext, stratifiedContexts, rawContext));
        }

        AFCalcResult AFresult = afcm.get().getLog10PNonRef(vc, getAlleleFrequencyPriors(model));

        // is the most likely frequency conformation AC=0 for all alternate alleles?
        boolean bestGuessIsRef = true;

        // determine which alternate alleles have AF>0
        final List<Allele> myAlleles = new ArrayList<Allele>(vc.getAlleles().size());
        final List<Integer> alleleCountsofMLE = new ArrayList<Integer>(vc.getAlleles().size());
        myAlleles.add(vc.getReference());
        for ( int i = 0; i < AFresult.getAllelesUsedInGenotyping().size(); i++ ) {
            final Allele alternateAllele = AFresult.getAllelesUsedInGenotyping().get(i);
            if ( alternateAllele.isReference() )
                continue;

            // Compute if the site is considered polymorphic with sufficient confidence relative to our
            // phred-scaled emission QUAL
            final boolean isNonRef = AFresult.isPolymorphicPhredScaledQual(alternateAllele, UAC.STANDARD_CONFIDENCE_FOR_EMITTING);

            // if the most likely AC is not 0, then this is a good alternate allele to use
            if ( isNonRef ) {
                myAlleles.add(alternateAllele);
                alleleCountsofMLE.add(AFresult.getAlleleCountAtMLE(alternateAllele));
                bestGuessIsRef = false;
            }
            // if in GENOTYPE_GIVEN_ALLELES mode, we still want to allow the use of a poor allele
            else if ( UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
                myAlleles.add(alternateAllele);
                alleleCountsofMLE.add(AFresult.getAlleleCountAtMLE(alternateAllele));
            }
        }

        final double PoFGT0 = Math.pow(10, AFresult.getLog10PosteriorOfAFGT0());

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        final double phredScaledConfidence =
                Math.abs(! bestGuessIsRef || UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
                        ? -10 * AFresult.getLog10PosteriorOfAFEq0()
                        : -10 * AFresult.getLog10PosteriorOfAFGT0());

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES && !passesEmitThreshold(phredScaledConfidence, bestGuessIsRef) ) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            return limitedContext ? null : estimateReferenceConfidence(vc, stratifiedContexts, getTheta(model), true, PoFGT0);
        }

        // start constructing the resulting VC
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(vc);
        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), myAlleles);
        builder.log10PError(phredScaledConfidence/-10.0);
        if ( ! passesCallThreshold(phredScaledConfidence) )
            builder.filters(filter);

        // create the genotypes
        final GenotypesContext genotypes = afcm.get().subsetAlleles(vc, myAlleles, true,ploidy);
        builder.genotypes(genotypes);

        // print out stats if we have a writer
        if ( verboseWriter != null && !limitedContext )
            printVerboseData(refContext.getLocus().toString(), vc, PoFGT0, phredScaledConfidence, model);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        final HashMap<String, Object> attributes = new HashMap<String, Object>();

        // inherit attributed from input vc if requested
        if (inheritAttributesFromInputVC)
            attributes.putAll(vc.getAttributes());
        // if the site was downsampled, record that fact
        if ( !limitedContext && rawContext.hasPileupBeenDownsampled() )
            attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

        if ( UAC.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED )
            attributes.put(NUMBER_OF_DISCOVERED_ALLELES_KEY, vc.getAlternateAlleles().size());

        // add the MLE AC and AF annotations
        if ( alleleCountsofMLE.size() > 0 ) {
            attributes.put(VCFConstants.MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
            final int AN = builder.make().getCalledChrCount();
            final ArrayList<Double> MLEfrequencies = new ArrayList<Double>(alleleCountsofMLE.size());
            // the MLEAC is allowed to be larger than the AN (e.g. in the case of all PLs being 0, the GT is ./. but the exact model may arbitrarily choose an AC>1)
            for ( int AC : alleleCountsofMLE )
                MLEfrequencies.add(Math.min(1.0, (double)AC / (double)AN));
            attributes.put(VCFConstants.MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
        }

        if ( UAC.COMPUTE_SLOD && !limitedContext && !bestGuessIsRef ) {
            //final boolean DEBUG_SLOD = false;

            // the overall lod
            //double overallLog10PofNull = AFresult.log10AlleleFrequencyPosteriors[0];
            double overallLog10PofF = AFresult.getLog10LikelihoodOfAFGT0();
            //if ( DEBUG_SLOD ) System.out.println("overallLog10PofF=" + overallLog10PofF);

            List<Allele> allAllelesToUse = builder.make().getAlleles();
            
            // the forward lod
            VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.FORWARD, allAllelesToUse, false, model, perReadAlleleLikelihoodMap);
            AFresult = afcm.get().getLog10PNonRef(vcForward, getAlleleFrequencyPriors(model));
            //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(AFresult.log10AlleleFrequencyPosteriors, true);
            double forwardLog10PofNull = AFresult.getLog10LikelihoodOfAFEq0();
            double forwardLog10PofF = AFresult.getLog10LikelihoodOfAFGT0();
            //if ( DEBUG_SLOD ) System.out.println("forwardLog10PofNull=" + forwardLog10PofNull + ", forwardLog10PofF=" + forwardLog10PofF);

            // the reverse lod
            VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.REVERSE, allAllelesToUse, false, model, perReadAlleleLikelihoodMap);
            AFresult = afcm.get().getLog10PNonRef(vcReverse, getAlleleFrequencyPriors(model));
            //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(AFresult.log10AlleleFrequencyPosteriors, true);
            double reverseLog10PofNull = AFresult.getLog10LikelihoodOfAFEq0();
            double reverseLog10PofF = AFresult.getLog10LikelihoodOfAFGT0();
            //if ( DEBUG_SLOD ) System.out.println("reverseLog10PofNull=" + reverseLog10PofNull + ", reverseLog10PofF=" + reverseLog10PofF);

            double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
            double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
            //if ( DEBUG_SLOD ) System.out.println("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

            // strand score is max bias between forward and reverse strands
            double strandScore = Math.max(forwardLod, reverseLod);
            // rescale by a factor of 10
            strandScore *= 10.0;
            //logger.debug(String.format("SLOD=%f", strandScore));

            if ( !Double.isNaN(strandScore) )
                attributes.put("SB", strandScore);
        }

        // finish constructing the resulting VC
        builder.attributes(attributes);
        VariantContext vcCall = builder.make();

        // if we are subsetting alleles (either because there were too many or because some were not polymorphic)
        // then we may need to trim the alleles (because the original VariantContext may have had to pad at the end).
        if ( myAlleles.size() != vc.getAlleles().size() && !limitedContext ) // limitedContext callers need to handle allele trimming on their own to keep their perReadAlleleLikelihoodMap alleles in sync
            vcCall = GATKVariantContextUtils.reverseTrimAlleles(vcCall);

        if ( annotationEngine != null && !limitedContext ) { // limitedContext callers need to handle annotations on their own by calling their own annotationEngine
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadBackedPileup pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

            vcCall = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall, perReadAlleleLikelihoodMap);
        }

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PoFGT0));
    }

    private Map<String, AlignmentContext> getFilteredAndStratifiedContexts(UnifiedArgumentCollection UAC, ReferenceContext refContext, AlignmentContext rawContext, final GenotypeLikelihoodsCalculationModel.Model model) {

        if ( !BaseUtils.isRegularBase(refContext.getBase()) )
            return null;

        Map<String, AlignmentContext> stratifiedContexts = null;

        if ( model.name().contains("INDEL") ) {

            final ReadBackedPileup pileup = rawContext.getBasePileup().getMappingFilteredPileup(UAC.MIN_BASE_QUALTY_SCORE);
            // don't call when there is no coverage
            if ( pileup.getNumberOfElements() == 0 && UAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES  )
                return null;

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

        } else if ( model.name().contains("SNP") ) {

            // stratify the AlignmentContext and cut by sample
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(rawContext.getBasePileup());

            if ( !(UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES) ) {
                int numDeletions = 0;
                for ( final PileupElement p : rawContext.getBasePileup() ) {
                    if ( p.isDeletion() )
                        numDeletions++;
                }
                if ( ((double) numDeletions) / ((double) rawContext.getBasePileup().getNumberOfElements()) > UAC.MAX_DELETION_FRACTION ) {
                    return null;
                }
            }
        }

        return stratifiedContexts;
    }

    private final static double[] binomialProbabilityDepthCache = new double[10000];
    private final static double REF_BINOMIAL_PROB_LOG10_0_5 = Math.log10(0.5);

    static {
        for ( int i = 1; i < binomialProbabilityDepthCache.length; i++ ) {
            binomialProbabilityDepthCache[i] = MathUtils.log10BinomialProbability(i, 0, REF_BINOMIAL_PROB_LOG10_0_5);
        }
    }

    private final double getRefBinomialProbLog10(final int depth) {
        if ( depth < binomialProbabilityDepthCache.length )
            return binomialProbabilityDepthCache[depth];
        else
            return MathUtils.log10BinomialProbability(depth, 0, REF_BINOMIAL_PROB_LOG10_0_5);
    }

    private VariantCallContext estimateReferenceConfidence(VariantContext vc, Map<String, AlignmentContext> contexts, double theta, boolean ignoreCoveredSamples, double initialPofRef) {
        if ( contexts == null )
            return null;

        double log10POfRef = Math.log10(initialPofRef);

        // for each sample that we haven't examined yet
        for ( String sample : samples ) {
            final AlignmentContext context = contexts.get(sample);
            if ( ignoreCoveredSamples && context != null )
                continue;
            final int depth = context == null ? 0 : context.getBasePileup().depthOfCoverage();
            log10POfRef += estimateLog10ReferenceConfidenceForOneSample(depth, theta);
        }

        return new VariantCallContext(vc, QualityUtils.phredScaleLog10CorrectRate(log10POfRef) >= UAC.STANDARD_CONFIDENCE_FOR_CALLING, false);
    }

    /**
     * Compute the log10 probability of a sample with sequencing depth and no alt allele is actually truly homozygous reference
     *
     * Assumes the sample is diploid
     *
     * @param depth the depth of the sample
     * @param theta the heterozygosity of this species (between 0 and 1)
     * @return a valid log10 probability of the sample being hom-ref
     */
    @Requires({"depth >= 0", "theta >= 0.0 && theta <= 1.0"})
    @Ensures("MathUtils.goodLog10Probability(result)")
    protected double estimateLog10ReferenceConfidenceForOneSample(final int depth, final double theta) {
        final double log10PofNonRef = Math.log10(theta / 2.0) + getRefBinomialProbLog10(depth);
        return MathUtils.log10OneMinusX(Math.pow(10.0, log10PofNonRef));
    }

    protected void printVerboseData(String pos, VariantContext vc, double PofF, double phredScaledConfidence, final GenotypeLikelihoodsCalculationModel.Model model) {
        Allele refAllele = null, altAllele = null;
        for ( Allele allele : vc.getAlleles() ) {
            if ( allele.isReference() )
                refAllele = allele;
            else
                altAllele = allele;
        }

        for (int i = 0; i <= N; i++) {
            StringBuilder AFline = new StringBuilder("AFINFO\t");
            AFline.append(pos);
            AFline.append("\t");
            AFline.append(refAllele);
            AFline.append("\t");
            if ( altAllele != null )
                AFline.append(altAllele);
            else
                AFline.append("N/A");
            AFline.append("\t");
            AFline.append(i + "/" + N + "\t");
            AFline.append(String.format("%.2f\t", ((float)i)/N));
            AFline.append(String.format("%.8f\t", getAlleleFrequencyPriors(model)[i]));
            verboseWriter.println(AFline.toString());
        }

        verboseWriter.println("P(f>0) = " + PofF);
        verboseWriter.println("Qscore = " + phredScaledConfidence);
        verboseWriter.println();
    }

    protected boolean passesEmitThreshold(double conf, boolean bestGuessIsRef) {
        return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) && conf >= Math.min(UAC.STANDARD_CONFIDENCE_FOR_CALLING, UAC.STANDARD_CONFIDENCE_FOR_EMITTING);
    }

    protected boolean passesCallThreshold(double conf) {
        return conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    protected boolean confidentlyCalled(double conf, double PofF) {
        return conf >= UAC.STANDARD_CONFIDENCE_FOR_CALLING ||
                (UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES && QualityUtils.phredScaleErrorRate(PofF) >= UAC.STANDARD_CONFIDENCE_FOR_CALLING);
    }

    private void determineGLModelsToUse() {

        String modelPrefix = "";
        if ( !UAC.GLmodel.name().contains(GPSTRING) && UAC.samplePloidy != GATKVariantContextUtils.DEFAULT_PLOIDY )
            modelPrefix = GPSTRING;

        if ( UAC.GLmodel.name().toUpperCase().contains("BOTH") ) {
            modelPrefix += UAC.GLmodel.name().toUpperCase().replaceAll("BOTH","");
            modelsToUse.add(GenotypeLikelihoodsCalculationModel.Model.valueOf(modelPrefix+"SNP"));
            modelsToUse.add(GenotypeLikelihoodsCalculationModel.Model.valueOf(modelPrefix+"INDEL"));
        }
        else {
            modelsToUse.add(GenotypeLikelihoodsCalculationModel.Model.valueOf(modelPrefix+UAC.GLmodel.name().toUpperCase()));
        }
    }

    // decide whether we are currently processing SNPs, indels, neither, or both
    private List<GenotypeLikelihoodsCalculationModel.Model> getGLModelsToUse(final RefMetaDataTracker tracker,
                                                                             final ReferenceContext refContext,
                                                                             final AlignmentContext rawContext) {

        if ( UAC.GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES )
            return modelsToUse;

        // if we're genotyping given alleles then we need to choose the model corresponding to the variant type requested
        final List<GenotypeLikelihoodsCalculationModel.Model> GGAmodel = new ArrayList<GenotypeLikelihoodsCalculationModel.Model>(1);
        final VariantContext vcInput = getVCFromAllelesRod(tracker, refContext, rawContext.getLocation(), false, logger, UAC.alleles);
        if ( vcInput == null )
            return GGAmodel; // no work to be done

        if ( vcInput.isSNP() )  {
            // use the SNP model unless the user chose INDEL mode only
            if ( modelsToUse.size() == 2 || modelsToUse.get(0).name().endsWith("SNP") )
                GGAmodel.add(modelsToUse.get(0));
        }
        else if ( vcInput.isIndel() || vcInput.isMixed() ) {
            // use the INDEL model unless the user chose SNP mode only
            if ( modelsToUse.size() == 2 )
                GGAmodel.add(modelsToUse.get(1));
            else if ( modelsToUse.get(0).name().endsWith("INDEL") )
                GGAmodel.add(modelsToUse.get(0));
        }
        // No support for other types yet

        return GGAmodel;
    }

    public static void computeAlleleFrequencyPriors(final int N, final double[] priors, final double theta) {

        double sum = 0.0;

        // for each i
        for (int i = 1; i <= N; i++) {
            final double value = theta / (double)i;
            priors[i] = Math.log10(value);
            sum += value;
        }

        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
    }

    protected double[] getAlleleFrequencyPriors( final GenotypeLikelihoodsCalculationModel.Model model ) {
        if (model.name().toUpperCase().contains("SNP"))
            return log10AlleleFrequencyPriorsSNPs;
        else if (model.name().toUpperCase().contains("INDEL"))
            return log10AlleleFrequencyPriorsIndels;
        else
            throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);

    }

    protected double getTheta( final GenotypeLikelihoodsCalculationModel.Model model ) {
        if( model.name().contains("SNP") )
            return HUMAN_SNP_HETEROZYGOSITY;
        if( model.name().contains("INDEL") )
            return HUMAN_INDEL_HETEROZYGOSITY;
        else throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
    }

    private static Map<String, GenotypeLikelihoodsCalculationModel> getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {

        final Map<String, GenotypeLikelihoodsCalculationModel> glcm = new HashMap<String, GenotypeLikelihoodsCalculationModel>();
        final List<Class<? extends GenotypeLikelihoodsCalculationModel>> glmClasses = new PluginManager<GenotypeLikelihoodsCalculationModel>(GenotypeLikelihoodsCalculationModel.class).getPlugins();

        for (int i = 0; i < glmClasses.size(); i++) {
            final Class<? extends GenotypeLikelihoodsCalculationModel> glmClass = glmClasses.get(i);
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

    public static VariantContext getVCFromAllelesRod(RefMetaDataTracker tracker, ReferenceContext ref, GenomeLoc loc, boolean requireSNP, Logger logger, final RodBinding<VariantContext> allelesBinding) {
        if ( tracker == null || ref == null || logger == null )
            return null;
        VariantContext vc = null;

        // search for usable record
        for ( final VariantContext vc_input : tracker.getValues(allelesBinding, loc) ) {
            if ( vc_input != null && ! vc_input.isFiltered() && (! requireSNP || vc_input.isSNP() )) {
                if ( vc == null ) {
                    vc = vc_input;
                } else {
                    logger.warn("Multiple valid VCF records detected in the alleles input file at site " + ref.getLocus() + ", only considering the first record");
                }
            }
        }

        return vc;
    }
}
