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

package org.broadinstitute.gatk.tools.walkers.cancer.m2;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.io.DirectOutputTracker;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.*;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.gatk.utils.activeregion.ActivityProfileState;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.fragments.FragmentUtils;
import org.broadinstitute.gatk.utils.genotyper.*;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.DroppedReadsTracker;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pairhmm.PairHMM;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.*;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;

import java.io.FileNotFoundException;
import java.util.*;

import static java.lang.Math.pow;

/**
 * Call somatic SNPs and indels via local re-assembly of haplotypes
 *
 * <p>MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (<a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>) with the assembly-based machinery of HaplotypeCaller. The basic operation of MuTect2 proceeds similarly to that of the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>.</p>
 *
 * <h3>Differences from HaplotypeCaller</h3>
 * <p>While the HaplotypeCaller relies on a ploidy assumption (diploid by default) to inform its genotype likelihood and variant quality calculations, MuTect2 allows for a varying allelic fraction for each variant, as is often seen in tumors with purity less than 100%, multiple subclones, and/or copy number variation (either local or aneuploidy). MuTect2 also differs from the HaplotypeCaller in that it does apply some hard filters to variants before producing output.</p>
 * <p>Note that the GVCF generation capabilities of HaplotypeCaller are NOT available in MuTect2, even though some of the relevant arguments are listed below. There are currently no plans to make GVCF calling available in MuTect2.</p>
 *
 * <h3>Usage examples</h3>
 * <p>These are example commands that show how to run MuTect2 for typical use cases. Square brackets ("[ ]") indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Tumor/Normal variant calling</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -T MuTect2 \
 *     -R reference.fasta \
 *     -I:tumor tumor.bam \
 *     -I:normal normal.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     [-L targets.interval_list] \
 *     -o output.vcf
 * </pre>
 *
 * <h4>Normal-only calling for panel of normals creation</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -T MuTect2 \
 *     -R reference.fasta \
 *     -I:tumor normal1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     --artifact_detection_mode \
 *     [-L targets.interval_list] \
 *     -o output.normal1.vcf
 * </pre>
 * <br />
 * For full PON creation, call each of your normals separately in artifact detection mode as shown above. Then use CombineVariants to output only sites where a variant was seen in at least two samples:
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -T CombineVariants \
 *     -R reference.fasta \
 *     -V output.normal1.vcf -V output.normal2.vcf [-V output.normal2.vcf ...] \
 *     -minN 2 \
 *     --setKey "null" \
 *     --filteredAreUncalled \
 *     --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
 *     [-L targets.interval_list] \
 *     -o MuTect2_PON.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>As noted in several places in the documentation, MuTect2 is currently released under BETA status; it is NOT recommended for production work and is NOT available for commercial/for-profit licensing.</li>
 * <li>MuTect2 currently only supports the calling of a single tumor-normal pair at a time.</li>
 * <li>Tumor-only variant calling is possible but it is NOT supported and we will not answer any questions about it until it becomes a supported feature.</li>
 * <li>Some of the arguments listed below are not functional; they are exclusive to HaplotypeCaller and are listed here due to technical entanglements in the code. This will be resolved in the upcoming GATK 4 release.</li>
 * </ul>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ActiveRegionTraversalParameters(extension=100, maxRegion=300)
public class MuTect2 extends ActiveRegionWalker<List<VariantContext>, Integer> implements AnnotatorCompatible, NanoSchedulable {
    public static final String BAM_TAG_TUMOR = "tumor";
    public static final String BAM_TAG_NORMAL = "normal";

    protected Set<SAMReaderID> tumorSAMReaderIDs = new HashSet<>();
    protected Set<SAMReaderID> normalSAMReaderIDs = new HashSet<>();
    protected String tumorSampleName;
    protected String normalSampleName;

    protected SampleList samplesList;
    protected boolean printTCGAsampleHeader = false;

    protected CachingIndexedFastaSequenceFile referenceReader;
    protected LocalAssemblyEngine assemblyEngine = null;
    protected ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;
    protected SomaticGenotypingEngine genotypingEngine = null;
    private HaplotypeBAMWriter haplotypeBAMWriter;

    private byte MIN_TAIL_QUALITY;
    private double log10GlobalReadMismappingRate;
    private static final int REFERENCE_PADDING = 500;
    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @ArgumentCollection
    protected ReadThreadingAssemblerArgumentCollection RTAC = new ReadThreadingAssemblerArgumentCollection();

    @ArgumentCollection
    protected LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();

    @Argument(fullName = "debug_read_name", required = false, doc="trace this read name through the calling process")
    protected String DEBUG_READ_NAME = null;

    @Hidden
    @Advanced
    @Argument(fullName = "MQ_filtering_level", shortName = "MQthreshold", required = false, doc="Set an alternate MQ threshold for debugging")
    final public int MQthreshold = 20;

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     *  as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).
     *  Records that are filtered in the comp track will be ignored.
     *  Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Input(fullName="comp", shortName = "comp", doc="comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"DepthPerAlleleBySample", "BaseQualitySumPerAlleleBySample", "TandemRepeatAnnotator", "OxoGReadCounts"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Advanced
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<>(Arrays.asList(new String[]{"SpanningDeletions"}));

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { };

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    /**
     * Active region trimmer reference.
     */
    @ArgumentCollection
    protected ActiveRegionTrimmer trimmer = new ActiveRegionTrimmer();

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use read from this read group when making calls (but use all reads to build the assembly)", required = false)
    protected String keepRG = null;

    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public byte MIN_BASE_QUALTY_SCORE = 10;

    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.HOSTILE;

    @Hidden
    @Argument(fullName="errorCorrectReads", shortName="errorCorrectReads", doc = "Use an exploratory algorithm to error correct the kmers used during assembly.  May cause fundamental problems with the assembly graph itself", required=false)
    protected boolean errorCorrectReads = false;

    @Hidden
    @Argument(fullName="captureAssemblyFailureBAM", shortName="captureAssemblyFailureBAM", doc="If specified, we will write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", required = false)
    protected boolean captureAssemblyFailureBAM = false;

    @Advanced
    @Argument(fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", doc="If specified, we will not analyze soft clipped bases in the reads", required = false)
    protected boolean dontUseSoftClippedBases = false;

    @Hidden
    @Argument(fullName="justDetermineActiveRegions", shortName="justDetermineActiveRegions", doc = "If specified, the HC won't actually do any assembly or calling, it'll just run the upfront active region determination code.  Useful for benchmarking and scalability testing", required=false)
    protected boolean justDetermineActiveRegions = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName="doNotRunPhysicalPhasing", shortName="doNotRunPhysicalPhasing", doc="Disable physical phasing", required = false)
    protected boolean doNotRunPhysicalPhasing = false;

    /**
     * When downsampling, level the coverage of the reads in each sample to no more than maxReadsInRegionPerSample reads,
     * not reducing coverage at any read start to less than minReadsPerAlignmentStart
     */
    @Argument(fullName = "maxReadsInRegionPerSample", shortName = "maxReadsInRegionPerSample", doc="Maximum reads in an active region", required = false)
    protected int maxReadsInRegionPerSample = 1000;

    @Argument(fullName = "minReadsPerAlignmentStart", shortName = "minReadsPerAlignStart", doc="Minimum number of reads sharing the same alignment start for each genomic location in an active region", required = false)
    protected int minReadsPerAlignmentStart = 5;

    public RodBinding<VariantContext> getDbsnpRodBinding() { return MTAC.dbsnp.dbsnp; }

    @Override
    public void initialize() {
        super.initialize();

        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader())));

        // MUTECT: check that we have at least one tumor bam
        for(final SAMReaderID id : getToolkit().getReadsDataSource().getReaderIDs()) {
            if (id.getTags().getPositionalTags().size() == 0) {
                throw new RuntimeException("BAMs must be tagged as either 'tumor' or 'normal'");
            }

            // only supports single-sample BAMs (ie first read group is representative)
            final String bamSampleName = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();

            for(final String tag : id.getTags().getPositionalTags()) {
                if (BAM_TAG_TUMOR.equalsIgnoreCase(tag)) {
                    tumorSAMReaderIDs.add(id);
                    if (tumorSampleName == null) {
                        tumorSampleName = bamSampleName;
                    } else {
                        if (!tumorSampleName.equals(bamSampleName)) {
                            throw new UserException.BadInput("Found more than one tumor sample name in read data");
                        }
                    }
                } else if (BAM_TAG_NORMAL.equalsIgnoreCase(tag)) {
                    normalSAMReaderIDs.add(id);
                    if (normalSampleName == null) {
                        normalSampleName = bamSampleName;
                    } else {
                        if (!normalSampleName.equals(bamSampleName)) {
                            throw new UserException.BadInput("Found more than one normal sample name in read data");
                        }
                    }
                } else {
                    throw new RuntimeException("Unknown BAM tag '" + tag + "' must be either 'tumor' or 'normal'");
                }
            }
        }

        //If the samples specified are exactly one normal and one tumor, use the TCGA VCF sample header format
        if (samplesList.sampleCount() == 2 && normalSampleName != null && tumorSampleName != null && ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader()).size() == 2)
            printTCGAsampleHeader = true;

        final VariantAnnotatorEngine annotationEngine = initializeVCFOutput();

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        // create and setup the assembler
        assemblyEngine = new ReadThreadingAssembler(RTAC.maxNumHaplotypesInPopulation, RTAC.kmerSizes, RTAC.dontIncreaseKmerSizesForCycles, RTAC.allowNonUniqueKmersInRef, RTAC.numPruningSamples);

        assemblyEngine.setErrorCorrectKmers(RTAC.errorCorrectKmers);
        assemblyEngine.setPruneFactor(RTAC.MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(MTAC.DEBUG);
        assemblyEngine.setDebugGraphTransformations(RTAC.debugGraphTransformations);
        assemblyEngine.setAllowCyclesInKmerGraphToGeneratePaths(RTAC.allowCyclesInKmerGraphToGeneratePaths);
        assemblyEngine.setRecoverDanglingBranches(!RTAC.doNotRecoverDanglingBranches);
        assemblyEngine.setMinBaseQualityToUseInAssembly(MIN_BASE_QUALTY_SCORE);

        MIN_TAIL_QUALITY = (byte)(MIN_BASE_QUALTY_SCORE - 1);

        if ( RTAC.graphWriter != null ) assemblyEngine.setGraphWriter(RTAC.graphWriter);

        // setup the likelihood calculation engine
        if ( LEAC.phredScaledGlobalReadMismappingRate < 0 ) LEAC.phredScaledGlobalReadMismappingRate = -1;

        // configure the global mismapping rate
        if ( LEAC.phredScaledGlobalReadMismappingRate < 0 ) {
            log10GlobalReadMismappingRate = - Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(LEAC.phredScaledGlobalReadMismappingRate);
            logger.info("Using global mismapping rate of " + LEAC.phredScaledGlobalReadMismappingRate + " => " + log10GlobalReadMismappingRate + " in log10 likelihood units");
        }

        //static member function - set number of threads
        PairHMM.setNumberOfThreads(getToolkit().getTotalNumberOfThreads());
        // create our likelihood calculation engine
        likelihoodCalculationEngine = createLikelihoodCalculationEngine();

        final MergeVariantsAcrossHaplotypes variantMerger = new MergeVariantsAcrossHaplotypes();

        final GenomeAnalysisEngine toolkit = getToolkit();

        genotypingEngine = new SomaticGenotypingEngine( MTAC, samplesList, toolkit.getGenomeLocParser(), !doNotRunPhysicalPhasing, MTAC,
                tumorSampleName, normalSampleName, DEBUG_READ_NAME);
        genotypingEngine.setCrossHaplotypeEventMerger(variantMerger);
        genotypingEngine.setAnnotationEngine(annotationEngine);


        if ( MTAC.bamWriter != null ) {
            // we currently do not support multi-threaded BAM writing, so exception out
            if ( getToolkit().getTotalNumberOfThreads() > 1 )
                throw new UserException.BadArgumentValue("bamout", "Currently cannot emit a BAM file from the HaplotypeCaller in multi-threaded mode.");
            haplotypeBAMWriter = HaplotypeBAMWriter.create(MTAC.bamWriterType, MTAC.bamWriter, getToolkit().getSAMFileHeader());
        }

        // why isn't this a constructor (instead of initialize)?  Since the method is package-friendly
        trimmer.initialize(getToolkit().getGenomeLocParser(), MTAC.DEBUG,
                MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES, false);

        // KCIBUL: what's the right way to set this sensible default for somatic mutation calling from here?
        trimmer.snpPadding = 50;

        samplesList = toolkit.getReadSampleList();
        final Set<String> sampleSet = SampleListUtils.asSet(samplesList);

        if( MTAC.CONTAMINATION_FRACTION_FILE != null )
            MTAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(MTAC.CONTAMINATION_FRACTION_FILE, MTAC.CONTAMINATION_FRACTION, sampleSet, logger));

    }

    private VariantAnnotatorEngine initializeVCFOutput() {
        // initialize the output VCF header
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            annotationsToUse.add("ClusteredReadPosition");
        }
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        headerInfo.addAll(getM2HeaderLines());
        headerInfo.addAll(getSampleHeaderLines());

        // if printTCGAsampleHeader, we already checked for exactly 1 tumor and 1 normal in printTCGAsampleHeader assignment in initialize()
        final List<String> outputSampleNames = printTCGAsampleHeader ? Arrays.asList("TUMOR", "NORMAL") : SampleListUtils.asList(samplesList);

        vcfWriter.writeHeader(new VCFHeader(headerInfo, outputSampleNames));

        return annotationEngine;
    }

    private Set<VCFHeaderLine> getM2HeaderLines(){
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NORMAL_LOD_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TUMOR_LOD_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.PANEL_OF_NORMALS_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.HAPLOTYPE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY));

        if (MTAC.ENABLE_STRAND_ARTIFACT_FILTER){
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TLOD_FWD_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TLOD_REV_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY));
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY));
        }

        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ALLELE_FRACTION_KEY));

        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.PON_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.ALT_ALLELE_IN_NORMAL_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.TUMOR_LOD_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.TRIALLELIC_SITE_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME));
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.CLUSTERED_READ_POSITION_FILTER_NAME));


        if ( ! doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }
        return headerInfo;
    }

    private Set<VCFHeaderLine> getSampleHeaderLines(){
        final Set<VCFHeaderLine> sampleLines = new HashSet<>();
        if (printTCGAsampleHeader) {
            //NOTE: This will only list the first bam file for each tumor/normal sample if there is more than one
            final Map<String, String> normalSampleHeaderAttributes = new HashMap<>();
            normalSampleHeaderAttributes.put("ID", "NORMAL");
            normalSampleHeaderAttributes.put("SampleName", normalSampleName);
            if (normalSAMReaderIDs.iterator().hasNext() && !getToolkit().getArguments().disableCommandLineInVCF)
                normalSampleHeaderAttributes.put("File", normalSAMReaderIDs.iterator().next().getSamFilePath());
            final VCFSimpleHeaderLine normalSampleHeader = new VCFSimpleHeaderLine("SAMPLE", normalSampleHeaderAttributes);
            final Map<String, String> tumorSampleHeaderAttributes = new HashMap<>();
            tumorSampleHeaderAttributes.put("ID", "TUMOR");
            tumorSampleHeaderAttributes.put("SampleName", tumorSampleName);
            if (tumorSAMReaderIDs.iterator().hasNext() && !getToolkit().getArguments().disableCommandLineInVCF)
                tumorSampleHeaderAttributes.put("File", tumorSAMReaderIDs.iterator().next().getSamFilePath());
            final VCFSimpleHeaderLine tumorSampleHeader = new VCFSimpleHeaderLine("SAMPLE", tumorSampleHeaderAttributes);

            sampleLines.add(normalSampleHeader);
            sampleLines.add(tumorSampleHeader);
        }
        return sampleLines;
    }

    @Override
    public ActivityProfileState isActive(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if( context == null || context.getBasePileup().isEmpty() )
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getLocus(), 0.0);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final AlignmentContext tumorContext = splitContexts.get(tumorSampleName);
        final AlignmentContext normalContext = splitContexts.get(normalSampleName);

        // if there are no tumor reads... there is no activity!
        if (tumorContext == null) {
            return new ActivityProfileState(ref.getLocus(), 0);
        }

        // KCIBUL -- this method was inlined and modified from ReferenceConfidenceModel
        final ReadBackedPileup tumorPileup = tumorContext.getBasePileup().getMappingFilteredPileup(MQthreshold);
        final double[] tumorGLs = calcGenotypeLikelihoodsOfRefVsAny(tumorPileup, ref.getBase(), MIN_BASE_QUALTY_SCORE);
        final double tumorLod = tumorGLs[1] - tumorGLs[0];

        // NOTE: do I want to convert to a probability (or just keep this as a LOD score)

        // also at this point, we should throw out noisy sites (hence the nonRefInNormalCheck) but this is non-optimal
        double prob = 0;
        if (tumorLod > MTAC.INITIAL_TUMOR_LOD_THRESHOLD) {

            // TODO: should we even do this performance optimization?
            // in any case, we have to handle the case where there is no normal (and thus no normal context) which is
            // different than having a normal but having no reads (where we should not enter the active region)
            if (normalSampleName != null && normalContext != null) {
                final int nonRefInNormal = getCountOfNonRefEvents(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE);

                final double[] normalGLs = calcGenotypeLikelihoodsOfRefVsAny(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE, 0.5f);
                final double normalLod = normalGLs[0] - normalGLs[1];

                // TODO: parameterize these
                if (normalLod > 1.0 && nonRefInNormal < 4) {
                    prob = 1;
                    logger.debug("At " + ref.getLocus().toString() + " tlod: " + tumorLod + " nlod: " + normalLod + " with normal non-ref of " + nonRefInNormal);
                }
            } else {
                prob = 1;
                logger.debug("At " + ref.getLocus().toString() + " tlod: " + tumorLod + " and no-normal calling");
            }
        }

        return new ActivityProfileState( ref.getLocus(), prob, ActivityProfileState.Type.NONE, null);
    }

    private final static List<VariantContext> NO_CALLS = Collections.emptyList();
    @Override
    public List<VariantContext> map( final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker ) {
        if ( justDetermineActiveRegions ) {
            return NO_CALLS;
        } else if( !originalActiveRegion.isActive() ) {
            return referenceModelForNoVariation(originalActiveRegion, true);
        } else if( originalActiveRegion.size() == 0 ) {
            return referenceModelForNoVariation(originalActiveRegion, true);
        }
        logReadInfo(DEBUG_READ_NAME, originalActiveRegion.getReads(), "Present in original active region");

        // create the assembly using just high quality reads (Q20 or higher).  We want to use lower
        // quality reads in the PairHMM (and especially in the normal) later, so we can't use a ReadFilter
        final ActiveRegion assemblyActiveRegion = new ActiveRegion(originalActiveRegion.getLocation(), originalActiveRegion.getSupportingStates(),originalActiveRegion.isActive(), getToolkit().getGenomeLocParser(), originalActiveRegion.getExtension());
        for (final GATKSAMRecord rec : originalActiveRegion.getReads()) {
            if (rec.getMappingQuality() >= MQthreshold ) {
                assemblyActiveRegion.add(rec);
            }
        }
        logReadInfo(DEBUG_READ_NAME, assemblyActiveRegion.getReads(), "Present in assembly active region");

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(assemblyActiveRegion, Collections.EMPTY_LIST);
        final TreeSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        final ActiveRegionTrimmer.Result trimmingResult = trimmer.trim(originalActiveRegion,allVariationEvents);
        if (!trimmingResult.isVariationPresent()) {
            return referenceModelForNoVariation(originalActiveRegion, false);
        }
        logReadInfo(DEBUG_READ_NAME, trimmingResult.getCallableRegion().getReads(), "Present in trimming result");

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        // it is conceivable that the region is active yet has no events upon assembling only the well-mapped reads
        if( ! assemblyResult.isVariationPresent() ) {
            return referenceModelForNoVariation(originalActiveRegion, false);
        }

        final ActiveRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping");

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.TRIMMMED, originalActiveRegion.getReads(), regionForGenotyping.getReads());
        }

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKSAMRecord> filteredReads = filterNonPassingReads(regionForGenotyping);

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReads(DroppedReadsTracker.Reason.FILTERED, filteredReads);
        }

        final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample(filteredReads);
        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping after filtering reads");

        // evaluate each sample's reads against all haplotypes

        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKSAMRecord>> reads = splitReadsBySample( regionForGenotyping.getReads() );
        for (final List<GATKSAMRecord> rec : reads.values()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after splitting assemblyResult by sample");
        }

        //TODO: this should be written as
        //TODO final Map<String, Integer> ARreads_origNormalMQ = regionForGenotyping.getReads().stream()
        //TODO        .collect(Collectors.toMap(GATKSAMRecord::getReadName, GATKSAMRecord::getMappingQuality));
        //TODO: but for some reason sometimes streamifying Mutect2 code causes a MalformedWalkerArgumentsException
        //TODO: this will probably not occur after the port to GATK4
        final HashMap<String, Integer> ARreads_origNormalMQ = new HashMap<>();
        for (final GATKSAMRecord read : regionForGenotyping.getReads()) {
            ARreads_origNormalMQ.put(read.getReadName(), read.getMappingQuality());
        }

        // modify MAPQ scores in normal to be high so that we don't do any base quality score capping
        for(final GATKSAMRecord rec : regionForGenotyping.getReads()) {
            if (isReadFromNormal(rec)) {
                rec.setMappingQuality(60);
            }
        }

        logger.debug("Computing read likelihoods with " + regionForGenotyping.getReads().size() + " reads against " + haplotypes.size() + " haplotypes across region " + assemblyResult.getRegionForGenotyping().toString());
        final ReadLikelihoods<Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);

        final Map<GATKSAMRecord,GATKSAMRecord> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(
                    DroppedReadsTracker.Reason.REALIGNMENT_FAILURE,
                    regionForGenotyping.getReads(),
                    readRealignments.values());
        }

        readLikelihoods.changeReads(readRealignments);
        for (final GATKSAMRecord rec : readRealignments.keySet()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after computing read likelihoods");
        }

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods,
                ARreads_origNormalMQ,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getLocation(),
                metaDataTracker);

        if ( MTAC.bamWriter != null ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (MTAC.disableOptimizations)
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    assemblyResult.getPaddedReferenceLoc(),
                    haplotypes,
                    calledHaplotypeSet,
                    readLikelihoods);

            if ( MTAC.emitDroppedReads ) {
                haplotypeBAMWriter.writeDroppedReads();
            }
        }

        if( MTAC.DEBUG ) { logger.info("----------------------------------------------------------------------------------"); }
        return annotateVCs(calledHaplotypes, metaDataTracker);
    }

    private Set<String> calculateFilters(final RefMetaDataTracker metaDataTracker, final VariantContext vc, final Map<String, Object> eventDistanceAttributes) {
        final Set<String> filters = new HashSet<>();

        final Integer eventCount = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY);
        final Integer maxEventDistance = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY);

        final Collection<VariantContext> panelOfNormalsVC = metaDataTracker.getValues(MTAC.normalPanelRod,
                getToolkit().getGenomeLocParser().createGenomeLoc(vc.getChr(), vc.getStart()));
        final VariantContext ponVc = panelOfNormalsVC.isEmpty()?null:panelOfNormalsVC.iterator().next();

        if (ponVc != null) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }

        // FIXME: how do we sum qscores here?
        // FIXME: parameterize thresholds
        // && sum of alt likelihood scores > 20

        // TODO: make the change to have only a single normal sample (but multiple tumors is ok...)
        int normalAltCounts = 0;
        double normalF = 0;
        int normalAltQualityScoreSum = 0;
        if (hasNormal()) {
            final Genotype normalGenotype = vc.getGenotype(normalSampleName);

            // NOTE: how do we get the non-ref depth here?
            normalAltCounts = normalGenotype.getAD()[1];
            normalF = (Double) normalGenotype.getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);

            final Object qss = normalGenotype.getExtendedAttribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY);
            if (qss != null) {
                normalAltQualityScoreSum = (Integer) ((Object[]) qss)[1];
            } else {
                logger.error("Null qss at " + vc.getStart());
            }
        }

        if ( (normalAltCounts > MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || normalF > MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION ) && normalAltQualityScoreSum > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
            filters.add(GATKVCFConstants.ALT_ALLELE_IN_NORMAL_FILTER_NAME);
        } else if ( eventCount > 1 && normalAltCounts >= 1) {
            filters.add(GATKVCFConstants.MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME);
        } else if (eventCount >= 3) {
            filters.add(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME);
        }

        // STR contractions, that is the deletion of one repeat unit of a short repeat (>1bp repeat unit)
        // such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we
        // hard filter them out by default
        if (vc.isIndel()) {
            final ArrayList rpa = (ArrayList) vc.getAttribute(GATKVCFConstants.REPEATS_PER_ALLELE_KEY);
            final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa != null && rpa.size() > 1 && ru.length() > 1) {
                final int refCount = (Integer) rpa.get(0);
                final int altCount = (Integer) rpa.get(1);

                if (refCount - altCount == 1) {
                    filters.add(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME);
                }
            }
        }

        // NOTE: what if there is a 3bp indel followed by a snp... we are comparing starts
        // so it would be thrown but it's really an adjacent event
        if ( eventCount >= 2 && maxEventDistance >= 3) {
            filters.add(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME);
        }

        // clustered read position filter
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER){
            final Double tumorFwdPosMedian = (Double) vc.getAttribute(GATKVCFConstants.MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMedian = (Double) vc.getAttribute(GATKVCFConstants.MEDIAN_RIGHT_OFFSET_KEY);
            final Double tumorFwdPosMAD = (Double) vc.getAttribute(GATKVCFConstants.MAD_MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMAD = (Double) vc.getAttribute(GATKVCFConstants.MAD_MEDIAN_RIGHT_OFFSET_KEY);
            //If the variant is near the read end (median threshold) and the positions are very similar (MAD threshold) then filter
            if ( (tumorFwdPosMedian != null && tumorFwdPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorFwdPosMAD != null && tumorFwdPosMAD <= MTAC.PIR_MAD_THRESHOLD) ||
                    (tumorRevPosMedian != null && tumorRevPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorRevPosMAD != null && tumorRevPosMAD <= MTAC.PIR_MAD_THRESHOLD))
                filters.add(GATKVCFConstants.CLUSTERED_READ_POSITION_FILTER_NAME);
        }
        // TODO: Move strand bias filter here

        return filters;
    }

    private final static byte REF_MODEL_DELETION_QUAL = (byte) 30;
    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @return genotype likelihoods of [AA,AB]
     */
    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual, final double f) {
        final double[] genotypeLikelihoods = new double[2];
        final int AA = 0;
        final int AB=1;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {

                // TODO: why not use base qualities here?
                //double pobs = QualityUtils.qualToErrorProbLog10(qual);
                final double pobs = 1.0d - pow(10, (30 / -10.0));
                if( isNonRef(refBase, p)) {
                    genotypeLikelihoods[AB] += Math.log10(f*pobs + (1-f)*pobs/3.0d);
                    genotypeLikelihoods[AA] += Math.log10((1-pobs)/3);
                } else {
                    genotypeLikelihoods[AB] += Math.log10(f*(1-pobs)/3.0d + (1-f)*pobs);
                    genotypeLikelihoods[AA] += Math.log10(pobs);
                }
            }
        }

        return genotypeLikelihoods;
    }

    private boolean hasNormal() {
        return (normalSampleName != null);
    }

    //TODO: streamify in GATK4
    protected int getCountOfNonRefEvents(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual) {
        int i=0;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    i++;
                }
            }
        }
        return i;
    }

    protected double[] calcGenotypeLikelihoodsOfRefVsAny(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual) {
        final double f = calculateF(pileup, refBase, minBaseQual);
        return calcGenotypeLikelihoodsOfRefVsAny(pileup, refBase, minBaseQual, f);
    }

    private double calculateF(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual) {
        int refCount = 0, altCount = 0;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());

            // only consider deletions AND sites of sufficient quality
            if( p.isDeletion() || qual > minBaseQual ) {
                if( isNonRef(refBase, p)) {
                    altCount++;
                } else {
                    refCount++;
                }
            }
        }
        return (double) altCount / (refCount + altCount);
    }

    private boolean isNonRef(final byte refBase, final PileupElement p) {
        return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
    }

    int MIN_READ_LENGTH = 30; // private in superclass

    protected Set<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            //TODO: Takuto points out that this is questionable.  Let's think hard abut it.
            // KCIBUL: only perform read quality filtering on tumor reads...
            if (isReadFromNormal(rec)) {
                if( rec.getReadLength() < MIN_READ_LENGTH ) {
                    readsToRemove.add(rec);
                }
            } else if( rec.getReadLength() < MIN_READ_LENGTH || rec.getMappingQuality() < MQthreshold || BadMateFilter.hasBadMate(rec) ||
                    (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.maxNumberOfThreads = LEAC.maxNumThreadsHMM;
        args.useDoublePrecision = LEAC.useDoublePrecisionHMM;
        return new PairHMMLikelihoodCalculationEngine( (byte)LEAC.gcpHMM, LEAC.pairHMM, args, log10GlobalReadMismappingRate, LEAC.noFpga, pcrErrorModel );
    }

    /**
     * FROM HC
     *
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    protected List<VariantContext> referenceModelForNoVariation(final ActiveRegion region, final boolean needsToBeFinalized) {
            return NO_CALLS;
    }

    protected Map<String, List<GATKSAMRecord>> splitReadsBySample( final Collection<GATKSAMRecord> reads ) {
        return HaplotypeCaller.splitReadsBySample(samplesList, reads);
    }

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary and extended reads in the active region
    @Override
    public EnumSet<ActiveRegionReadState> desiredReadStates() {
            return EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY, ActiveRegionReadState.EXTENDED);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(final List<VariantContext> callsInRegion, final Integer numCalledRegions) {
        for( final VariantContext call : callsInRegion ) {
            vcfWriter.add( call );
        }
        return (callsInRegion.isEmpty() ? 0 : 1) + numCalledRegions;
    }

    @Override
    public void onTraversalDone(final Integer result) {
//        if ( SCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) ((GVCFWriter)vcfWriter).close(false); // GROSS -- engine forces us to close our own VCF writer since we wrapped it
//        referenceConfidenceModel.close();
        //TODO remove the need to call close here for debugging, the likelihood output stream should be managed
        //TODO (open & close) at the walker, not the engine.
        likelihoodCalculationEngine.close();
        logger.info("Ran local assembly on " + result + " active regions");
    }

    // The following are not used but are required by the AnnotatorCompatible interface
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * High-level function that runs the assembler on the active region reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param activeRegion the region we should assemble
     * @param giveAlleles additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    protected AssemblyResultSet assembleReads(final ActiveRegion activeRegion, final List<VariantContext> giveAlleles) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeActiveRegion(activeRegion); // handle overlapping fragments, clip adapter and low qual tails

        final byte[] fullReferenceWithPadding = activeRegion.getActiveRegionReference(referenceReader, REFERENCE_PADDING);
        final GenomeLoc paddedReferenceLoc = getPaddedLoc(activeRegion);
        final Haplotype referenceHaplotype = createReferenceHaplotype(activeRegion, paddedReferenceLoc);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        final ReadErrorCorrector readErrorCorrector = errorCorrectReads ? new ReadErrorCorrector(RTAC.kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION,
                RTAC.minObservationsForKmerToBeSolid, MTAC.DEBUG, fullReferenceWithPadding) : null;
        try {
            final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly( activeRegion, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, giveAlleles,readErrorCorrector );
            assemblyResultSet.debugDump(logger);
            return assemblyResultSet;
        } catch ( final Exception e ) {
            // Capture any exception that might be thrown, and write out the assembly failure BAM if requested
            if ( captureAssemblyFailureBAM ) {
                final SAMFileWriter writer = SAMFileWriterStub.createSAMFileWriter("assemblyFailure.bam", getToolkit());
                new DirectOutputTracker().addOutput((SAMFileWriterStub) writer);
                for ( final GATKSAMRecord read : activeRegion.getReads() ) {
                    writer.addAlignment(read);
                }
                writer.close();
            }
            throw e;
        }
    }

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if (activeRegion.isFinalized()) return;

        if( MTAC.DEBUG ) { logger.info("Assembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:    (with overlap region = " + activeRegion.getExtendedLoc() + ")"); }

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKSAMRecord> readsToUse = new ArrayList<>(activeRegion.getReads().size());
        for( final GATKSAMRecord myRead : activeRegion.getReads() ) {
            GATKSAMRecord clippedRead;
            if (errorCorrectReads)
                clippedRead = ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION );
            else  // default case: clip low qual ends of reads
                clippedRead= ReadClipper.hardClipLowQualEnds( myRead, MIN_TAIL_QUALITY );

            if ( dontUseSoftClippedBases || ! ReadUtils.hasWellDefinedFragmentSize(clippedRead) ) {
                // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
                clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
            } else {
                // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                // TODO -- reference haplotype start must be removed
                clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
            }

            clippedRead = ( clippedRead.getReadUnmappedFlag() ? clippedRead : ReadClipper.hardClipAdaptorSequence( clippedRead ) );
            if( !clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0 ) {
                clippedRead = ReadClipper.hardClipToRegion(clippedRead, activeRegion.getExtendedLoc().getStart(), activeRegion.getExtendedLoc().getStop());
                if( activeRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                    readsToUse.add(clippedRead);
                }
            }
        }

        final List<GATKSAMRecord> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.DOWNSAMPLED, activeRegion.getReads(), downsampledReads);
        }

        activeRegion.clearReads();
        activeRegion.addAll(downsampledReads);
        activeRegion.setFinalized(true);
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getExtendedLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getExtendedLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getExtendedLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getExtendedLoc().getContig(), padLeft, padRight);
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param activeRegion the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the GenomeLoc which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    private Haplotype createReferenceHaplotype(final ActiveRegion activeRegion, final GenomeLoc paddedReferenceLoc) {
        return ReferenceConfidenceModel.createReferenceHaplotype(activeRegion, activeRegion.getActiveRegionReference(referenceReader), paddedReferenceLoc);
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private void cleanOverlappingReadPairs(final List<GATKSAMRecord> reads, final Set<String> normalSampleNames) {
        final Map<String, List<GATKSAMRecord>> data = splitReadsBySample(reads);
        for ( final String sampleName : data.keySet() ) {
            final boolean isTumor = !normalSampleNames.contains(sampleName);
            final List<GATKSAMRecord> perSampleReadList = data.get(sampleName);

            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create(perSampleReadList);
            for ( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() )
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
        }
    }

    public static void logReadInfo(final String readName, final Collection<GATKSAMRecord> records, final String message) {
        if (readName != null) {
            for (final GATKSAMRecord rec : records) {
                logReadInfo(readName, rec, message);
            }
        }
    }

    public static void logReadInfo(final String readName, final GATKSAMRecord rec, final String message) {
        if (readName != null && rec != null && readName.equals(rec.getReadName())) {
            logger.info("Found " + rec.toString() + " - " + message);
        }
    }

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    // TODO: migrate from HC -> HCUtils Class and share it!
    private Map<GATKSAMRecord,GATKSAMRecord> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final GenomeLoc paddedReferenceLoc) {

        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKSAMRecord,GATKSAMRecord> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKSAMRecord originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKSAMRecord realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead,realignedRead);
        }
        return result;
    }

    private boolean isReadFromNormal(final GATKSAMRecord rec) {
        return normalSampleName != null && normalSampleName.equals(rec.getReadGroup().getSample());
    }

    public String getTumorSampleName(){
        return tumorSampleName;
    }

    final List<VariantContext> annotateVCs(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes, final RefMetaDataTracker metaDataTracker) {
        final int eventCount = calledHaplotypes.getCalls().size();
        final Map<String, Object> eventDistanceAttributes = new HashMap<>();    //TODO: should be Map<String, Integer> -- see TODO below
        eventDistanceAttributes.put(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount);
        if (eventCount > 1) {
            final int lastPosition = calledHaplotypes.getCalls().get(0).getStart();
            final int[] eventDistances = new int[calledHaplotypes.getCalls().size() - 1];
            for (int n = 0; n < eventDistances.length; n++) {
                eventDistances[n] = Math.abs(calledHaplotypes.getCalls().get(n + 1).getStart() - lastPosition);
            }
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY, MathUtils.arrayMin(eventDistances));
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY, MathUtils.arrayMax(eventDistances));
        } else {    //TODO: putting null is a hack -- we should remove this and update the integration test md5s
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY, null);
            eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY, null);
        }

        final List<VariantContext> annotatedCalls = new ArrayList<>();
        // can we do this with the Annotation classes instead?
        for (final VariantContext originalVC : calledHaplotypes.getCalls()) {
            final VariantContextBuilder vcb = new VariantContextBuilder(originalVC);

            final Map<String, Object> attributes = new HashMap<>(originalVC.getAttributes());
            attributes.putAll(eventDistanceAttributes);
            vcb.attributes(attributes);

            final Set<String> filters = new HashSet<>(originalVC.getFilters());

            final double tumorLod = originalVC.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, -1);
            if (tumorLod < MTAC.TUMOR_LOD_THRESHOLD) {
                filters.add(GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
            }

            // if we are in artifact detection mode, apply the thresholds for the LOD scores
            if (!MTAC.ARTIFACT_DETECTION_MODE) {
                filters.addAll(calculateFilters(metaDataTracker, originalVC, eventDistanceAttributes));
            }

            vcb.filters(filters.isEmpty() ? VariantContext.PASSES_FILTERS : filters);

            if (printTCGAsampleHeader) {
                final Genotype tumorGenotype = new GenotypeBuilder(originalVC.getGenotype(tumorSampleName)).name("TUMOR").make();
                final Genotype normalGenotype = new GenotypeBuilder(originalVC.getGenotype(normalSampleName)).name("NORMAL").make();
                vcb.genotypes(Arrays.asList(tumorGenotype, normalGenotype));
            }
            annotatedCalls.add(vcb.make());
        }
        return annotatedCalls;
    }
}


