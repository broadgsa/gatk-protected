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

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.io.PrintStream;
import java.util.*;

/**
 * A variant caller which unifies the approaches of several disparate callers -- Works for single-sample and multi-sample data.
 *
 * <p>
 * The GATK Unified Genotyper is a multiple-sample, technology-aware SNP and indel caller. It uses a Bayesian genotype
 * likelihood model to estimate simultaneously the most likely genotypes and allele frequency in a population of N samples,
 * emitting an accurate posterior probability of there being a segregating variant allele at each locus as well as for the
 * genotype of each sample. The system can either emit just the variant sites or complete genotypes (which includes
 * homozygous reference calls) satisfying some phred-scaled confidence value. The genotyper can make accurate calls on
 * both single sample data and multi-sample data.
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * The read data from which to make variant calls.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A raw, unfiltered, highly sensitive callset in VCF format.
 * </p>
 *
 * <h2>Example generic command for multi-sample SNP calling</h2>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -T UnifiedGenotyper \
 *   -I sample1.bam [-I sample2.bam ...] \
 *   --dbsnp dbSNP.vcf \
 *   -o snps.raw.vcf \
 *   -stand_call_conf [50.0] \
 *   -stand_emit_conf 10.0 \
 *   -dcov [50 for 4x, 200 for >30x WGS or Whole exome] \
 *   [-L targets.interval_list]
 * </pre>
 *
 * <p>
 * The above command will call all of the samples in your provided BAM files [-I arguments] together and produce a VCF file
 * with sites and genotypes for all samples. The easiest way to get the dbSNP file is from the GATK resource bundle (see Guide FAQs for details). Several
 * arguments have parameters that should be chosen based on the average coverage per sample in your data. See the detailed
 * argument descriptions below.
 * </p>
 *
 * <h2>Example command for generating calls at all sites</h2>
 * <pre>
 * java -jar /path/to/GenomeAnalysisTK.jar \
 *   -l INFO \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -T UnifiedGenotyper \
 *   -I /DCC/ftp/pilot_data/data/NA12878/alignment/NA12878.SLX.maq.SRP000031.2009_08.bam \
 *   -o my.vcf \
 *   --output_mode EMIT_ALL_SITES
 * </pre>
 *
 * <h2>Caveats</h2>
 * <ul>
 * <li>The system is under active and continuous development. All outputs, the underlying likelihood model, arguments, and
 * file formats are likely to change.</li>
 * <li>The system can be very aggressive in calling variants. In the 1000 genomes project for pilot 2 (deep coverage of ~35x)
 * we expect the raw Qscore > 50 variants to contain at least ~10% FP calls. We use extensive post-calling filters to eliminate
 * most of these FPs. Variant Quality Score Recalibration is a tool to perform this filtering.</li>
 * <li>The generalized ploidy model can be used to handle non-diploid or pooled samples (see the -ploidy argument in the table below).</li>
 * </ul>
 *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@BAQMode(QualityMode = BAQ.QualityMode.ADD_TAG, ApplicationTime = ReadTransformer.ApplicationTime.ON_INPUT)
@ReadFilters( {BadMateFilter.class, MappingQualityUnavailableFilter.class} )
@Reference(window=@Window(start=-200,stop=200))
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
    @Input(fullName="comp", shortName = "comp", doc="comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    //@Gather(className = "org.broadinstitute.sting.queue.extensions.gatk.CatVariantsGatherer")
    @Output(doc="File to which variants should be written",required=true)
    protected VariantContextWriter writer = null;

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
    protected List<String> annotationsToUse = new ArrayList<String>();

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<String>();

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private Set<String> samples;

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

        if ( UAC.TREAT_ALL_READS_AS_SINGLE_POOL ) {
            samples.add(GenotypeLikelihoodsCalculationModel.DUMMY_SAMPLE_NAME);
        } else {
            // get all of the unique sample names
            samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
            if ( UAC.referenceSampleName != null )
                samples.remove(UAC.referenceSampleName);
        }
        if ( UAC.CONTAMINATION_FRACTION_FILE != null )
            UAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(UAC.CONTAMINATION_FRACTION_FILE, UAC.CONTAMINATION_FRACTION, samples, logger));

        // check for a bad max alleles value
        if ( UAC.MAX_ALTERNATE_ALLELES > GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
            throw new UserException.BadArgumentValue("max_alternate_alleles", "the maximum possible value is " + GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);

        // warn the user for misusing EMIT_ALL_SITES
        if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_ALL_SITES &&
                UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY &&
                UAC.GLmodel != GenotypeLikelihoodsCalculationModel.Model.SNP )
            logger.warn("WARNING: note that the EMIT_ALL_SITES option is intended only for point mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce a comprehensive set of indels in DISCOVERY mode");
        
         // initialize the verbose writer
        if ( verboseWriter != null )
            verboseWriter.println("AFINFO\tLOC\tREF\tALT\tMAF\tF\tAFprior\tMLE\tMAP");

        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, verboseWriter, annotationEngine, samples, UAC.samplePloidy);

        // initialize the header
        Set<VCFHeaderLine> headerInfo = getHeaderInfo(UAC, annotationEngine, dbsnp);

        // invoke initialize() method on each of the annotation classes, allowing them to add their own header lines
        // and perform any necessary initialization/validation steps
        annotationEngine.invokeAnnotationInitializationMethods(headerInfo);

        writer.writeHeader(new VCFHeader(headerInfo, samples));


    }

    public static Set<VCFHeaderLine> getHeaderInfo(final UnifiedArgumentCollection UAC,
                                                   final VariantAnnotatorEngine annotationEngine,
                                                   final DbsnpArgumentCollection dbsnp) {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
        if ( annotationEngine != null )
            headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        if ( UAC.COMPUTE_SLOD )
            VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true, VCFConstants.STRAND_BIAS_KEY);

        if ( UAC.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED )
            headerInfo.add(new VCFInfoHeaderLine(UnifiedGenotyperEngine.NUMBER_OF_DISCOVERED_ALLELES_KEY, 1, VCFHeaderLineType.Integer, "Number of alternate alleles discovered (but not necessarily genotyped) at this site"));

        // add the pool values for each genotype
        if (UAC.samplePloidy != GATKVariantContextUtils.DEFAULT_PLOIDY) {
            headerInfo.add(new VCFFormatHeaderLine(VCFConstants.MLE_PER_SAMPLE_ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Maximum likelihood expectation (MLE) for the alternate allele count, in the same order as listed, for each individual sample"));
            headerInfo.add(new VCFFormatHeaderLine(VCFConstants.MLE_PER_SAMPLE_ALLELE_FRACTION_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Maximum likelihood expectation (MLE) for the alternate allele fraction, in the same order as listed, for each individual sample"));
        }
        if (UAC.referenceSampleName != null) {
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.REFSAMPLE_DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Total reference sample depth"));
        }

        VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true,
                VCFConstants.DOWNSAMPLED_KEY,
                VCFConstants.MLE_ALLELE_COUNT_KEY,
                VCFConstants.MLE_ALLELE_FREQUENCY_KEY);

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
        headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));

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
        return UG_engine.calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext, samples);
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
                    writer.add(call);
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
