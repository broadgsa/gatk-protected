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

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.utils.genotyper.IndexedSampleList;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.gatk.utils.baq.BAQ;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Call SNPs and indels on a per-locus basis
 *
 * <p>
 * This tool uses a Bayesian genotype likelihood model to estimate simultaneously the most likely genotypes and
 * allele frequency in a population of N samples, emitting a genotype for each sample. The system can either emit
 * just the variant sites or complete genotypes (which includes homozygous reference calls) satisfying some
 * phred-scaled confidence value.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * The read data from which to make variant calls.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A raw, unfiltered, highly sensitive callset in VCF format.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <h4>Multi-sample SNP calling</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T UnifiedGenotyper \
 *   -R reference.fasta \
 *   -I sample1.bam [-I sample2.bam ...] \
 *   --dbsnp dbSNP.vcf \
 *   -o snps.raw.vcf \
 *   -stand_call_conf [50.0] \
 *   -stand_emit_conf 10.0 \
 *   [-L targets.interval_list]
 * </pre>
 *
 * <h4>Generate calls at all sites</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T UnifiedGenotyper \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -o raw_variants.vcf \
 *   --output_mode EMIT_ALL_SITES
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The caller can be very aggressive in calling variants in order to be very sensitive, so the raw output will
 * contain many false positives. We use extensive post-calling filters to eliminate most of these FPs. See the documentation on filtering (especially by Variant Quality Score Recalibration) for more details.</li>
 * <li><b>This tool has been deprecated in favor of HaplotypeCaller, a much more sophisticated variant caller that
 * produces much better calls, especially on indels, and includes features that allow it to scale to much larger
 * cohort sizes.</b></li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle almost any ploidy (except very high ploidies in large pooled experiments); the ploidy
 * can be specified using the -ploidy argument for non-diploid organisms.</p>
 *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = ReadTransformer.ApplicationTime.ON_INPUT)
@ReadFilters( {BadMateFilter.class, MappingQualityUnavailableFilter.class} )
@Reference(window=@Window(start=-200,stop=200))
@PartitionBy(value = PartitionType.LOCUS, includeUnmapped = false)
@By(DataSource.REFERENCE)
// TODO -- When LocusIteratorByState gets cleaned up, we should enable multiple @By sources:
// TODO -- @By( {DataSource.READS, DataSource.REFERENCE_ORDERED_DATA} )
@Downsample(by=DownsampleType.BY_SAMPLE, toCoverage=250)
public class UnifiedGenotyper extends LocusWalker<List<VariantCallContext>, UnifiedGenotyper.UGStatistics> implements TreeReducible<UnifiedGenotyper.UGStatistics>, AnnotatorCompatible, NanoSchedulable {

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    /**
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).
     * Records that are filtered in the comp track will be ignored.
     * Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Input(fullName="comp", shortName = "comp", doc="Comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    //@Gather(className = "org.broadinstitute.gatk.queue.extensions.gatk.CatVariantsGatherer")
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter writer = null;

    @Advanced
    @Argument(fullName = "onlyEmitSamples", shortName = "onlyEmitSamples", doc = "If provided, only these samples will be emitted into the VCF, regardless of which samples are present in the BAM file", required = false)
    protected Set<String> onlyEmitSamples = Collections.emptySet();

    @Hidden
    @Argument(fullName = "debug_file", shortName = "debug_file", doc = "File to print all of the annotated and detailed debugging output", required = false)
    protected PrintStream verboseWriter = null;

    @Hidden
    @Argument(fullName = "metrics_file", shortName = "metrics", doc = "File to print any relevant callability metrics output", required = false)
    protected PrintStream metricsWriter = null;

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<>();

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<>();

    /**
     * If specified, all available annotations in the group will be applied. See the VariantAnnotator -list argument to view available groups.
     * Keep in mind that RODRequiringAnnotations are not intended to be used as a group, because they require specific ROD inputs.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls.  The single value 'none' removes the default group", required=false)
    protected String[] annotationClassesToUse = { "Standard", "StandardUG" };

    // the calculation arguments
    private UnifiedGenotypingEngine genotypingEngine = null;

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * Inner class for collecting output statistics from the UG
     */
    public static class UGStatistics {
        /** The total number of passes examined -- i.e., the number of map calls */
        long nBasesVisited = 0;

        /** The number of bases that were potentially callable -- i.e., those not at excessive coverage or masked with N */
        long nBasesCallable = 0;

        /** The number of bases called confidently (according to user threshold), either ref or other */
        long nBasesCalledConfidently = 0;

        /** The number of bases for which calls were emitted */
        long nCallsMade = 0;

        /** The total number of extended events encountered */
        long nExtendedEvents = 0;

        double percentCallableOfAll()    { return (100.0 * nBasesCallable) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfAll()      { return (100.0 * nBasesCalledConfidently) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / (nBasesCallable); }
    }

    /**
     * Initialize the samples, output, and genotype calculation model
     *
     **/
    public void initialize() {
        super.initialize();

        final GenomeAnalysisEngine toolkit = getToolkit();
        final Set<String> sampleNameSet;
        if ( UAC.TREAT_ALL_READS_AS_SINGLE_POOL ) {
            sampleNameSet = Collections.singleton(GenotypeLikelihoodsCalculationModel.DUMMY_SAMPLE_NAME);
        } else {
            // get all of the unique sample names
            sampleNameSet = ReadUtils.getSAMFileSamples(toolkit.getSAMFileHeader());
            if ( UAC.referenceSampleName != null )
                sampleNameSet.remove(UAC.referenceSampleName);
        }
        final SampleList samples = new IndexedSampleList(sampleNameSet);

        if ( UAC.CONTAMINATION_FRACTION_FILE != null )
            UAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(UAC.CONTAMINATION_FRACTION_FILE, UAC.CONTAMINATION_FRACTION, sampleNameSet, logger));

        // check for a bad max alleles value
        if ( UAC.genotypeArgs.MAX_ALTERNATE_ALLELES > GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
            throw new UserException.BadArgumentValue("max_alternate_alleles", "the maximum possible value is " + GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);

        // warn the user for misusing EMIT_ALL_SITES
        if ( UAC.outputMode == OutputMode.EMIT_ALL_SITES &&
                UAC.genotypingOutputMode == GenotypingOutputMode.DISCOVERY &&
                UAC.GLmodel != GenotypeLikelihoodsCalculationModel.Model.SNP )
            logger.warn("WARNING: note that the EMIT_ALL_SITES option is intended only for point mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce a comprehensive set of indels in DISCOVERY mode");
        
         // initialize the verbose writer
        if ( verboseWriter != null )
            verboseWriter.println("AFINFO\tLOC\tREF\tALT\tMAF\tF\tAFprior\tMLE\tMAP");

        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

        final AFCalculatorProvider afCalcAFCalculatorProvider = FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(),UAC,logger);
        genotypingEngine = new UnifiedGenotypingEngine(UAC, samples, toolkit.getGenomeLocParser(), afCalcAFCalculatorProvider, toolkit.getArguments().BAQMode);
        genotypingEngine.setVerboseWriter(verboseWriter);
        genotypingEngine.setAnnotationEngine(annotationEngine);

        // initialize the header
        Set<VCFHeaderLine> headerInfo = getHeaderInfo(UAC, annotationEngine, dbsnp);
        headerInfo.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // invoke initialize() method on each of the annotation classes, allowing them to add their own header lines
        // and perform any necessary initialization/validation steps
        annotationEngine.invokeAnnotationInitializationMethods(headerInfo);

        final Set<String> samplesForHeader;
        if ( ! onlyEmitSamples.isEmpty() ) {
            // make sure that onlyEmitSamples is a subset of samples
            if ( ! sampleNameSet.containsAll(onlyEmitSamples) )
                throw new UserException.BadArgumentValue("onlyEmitSamples", "must be a strict subset of the samples in the BAM files but is wasn't");
            samplesForHeader = onlyEmitSamples;
        } else {
            samplesForHeader = sampleNameSet;
        }
        writer.writeHeader(new VCFHeader(headerInfo, samplesForHeader));
    }



    public static Set<VCFHeaderLine> getHeaderInfo(final UnifiedArgumentCollection UAC,
                                                   final VariantAnnotatorEngine annotationEngine,
                                                   final DbsnpArgumentCollection dbsnp) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        // all annotation fields from VariantAnnotatorEngine
        if ( annotationEngine != null )
            headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        if ( UAC.COMPUTE_SLOD )
            VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true, VCFConstants.STRAND_BIAS_KEY);

        // add the pool values for each genotype
        if (UAC.genotypeArgs.samplePloidy != GATKVariantContextUtils.DEFAULT_PLOIDY) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MLE_PER_SAMPLE_ALLELE_COUNT_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.MLE_PER_SAMPLE_ALLELE_FRACTION_KEY));
        }
        if (UAC.referenceSampleName != null) {
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.REFSAMPLE_DEPTH_KEY));
        }

        if (UAC.annotateAllSitesWithPLs) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PL_FOR_ALL_SNP_ALLELES_KEY));
        }
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.DOWNSAMPLED_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));

        // also, check to see whether comp rods were included
        if ( dbsnp != null && dbsnp.dbsnp.isBound() )
            VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true, VCFConstants.DBSNP_KEY);

        // FORMAT fields
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        return headerInfo;
    }

    /**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public List<VariantCallContext> map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        return genotypingEngine.calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
    }

    public UGStatistics reduceInit() { return new UGStatistics(); }

    public UGStatistics treeReduce(UGStatistics lhs, UGStatistics rhs) {
        lhs.nBasesCallable += rhs.nBasesCallable;
        lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
        lhs.nBasesVisited += rhs.nBasesVisited;
        lhs.nCallsMade += rhs.nCallsMade;
        return lhs;
    }

    public UGStatistics reduce(List<VariantCallContext> calls, UGStatistics sum) {
        // we get a point for reaching reduce
        sum.nBasesVisited++;

        boolean wasCallable = false;
        boolean wasConfidentlyCalled = false;

        for ( VariantCallContext call : calls ) {
            if ( call == null )
                continue;

            // A call was attempted -- the base was callable
            wasCallable = true;

            // was the base confidently callable?
            wasConfidentlyCalled = call.confidentlyCalled;

            if ( call.shouldEmit ) {
                try {
                    // we are actually making a call
                    sum.nCallsMade++;
                    writer.add(subsetToEmitSamples(call));
                } catch (IllegalArgumentException e) {
                    throw new IllegalArgumentException(e.getMessage());
                }
            }
        }

        if ( wasCallable )
            sum.nBasesCallable++;

        if ( wasConfidentlyCalled )
            sum.nBasesCalledConfidently++;

        return sum;
    }

    /**
     * Subset the VariantContext down to just the emitting samples, if onlyEmitSamples has been provided
     * @param fullVC the VariantContext containing calls for all samples in the BAM files
     * @return a VariantContext that has been appropriately reduced to a subset of samples, if required
     */
    private VariantContext subsetToEmitSamples(final VariantContext fullVC) {
        if ( onlyEmitSamples.isEmpty() ) {
            return fullVC;
        } else {
            return GATKVariantContextUtils.trimAlleles(fullVC.subContextFromSamples(onlyEmitSamples, false), false, true);
        }
    }

    public void onTraversalDone(UGStatistics sum) {
        if ( metricsWriter != null ) {
            metricsWriter.println(String.format("Visited bases                                %d", sum.nBasesVisited));
            metricsWriter.println(String.format("Callable bases                               %d", sum.nBasesCallable));
            metricsWriter.println(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
            metricsWriter.println(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
            metricsWriter.println(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
            metricsWriter.println(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
            metricsWriter.println(String.format("Actual calls made                            %d", sum.nCallsMade));
        }
    }
}
