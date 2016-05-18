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

import java.io.FileNotFoundException;
import java.util.*;

import static java.lang.Math.pow;

/**
 * Call somatic SNPs and indels via local re-assembly of haplotypes
 *
 * <p>MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (<a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>) with the assembly-based machinery of HaplotypeCaller.</p>
 *
 * <p>The basic operation of MuTect2 proceeds similarly to that of the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>   </p>
 *
 * <h3>Differences from HaplotypeCaller</h3>
 * <p>While the HaplotypeCaller relies on a ploidy assumption (diploid by default) to inform its genotype likelihood and
 * variant quality calculations, MuTect2 allows for a varying allelic fraction for each variant, as is often seen in tumors with purity less
 * than 100%, multiple subclones, and/or copy number variation (either local or aneuploidy). MuTect2 also differs from the HaplotypeCaller in that it does apply some hard filters
 * to variants before producing output.</p>
 *
 * <h3>Usage examples</h3>
 * <p>These are example commands that show how to run MuTect2 for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the
 * Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Tumor/Normal variant calling</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar \
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
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T MuTect2
 *     -R reference.fasta
 *     -I:tumor normal1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     [--cosmic COSMIC.vcf] \
 *     --artifact_detection_mode \
 *     [-L targets.interval_list] \
 *     -o output.normal1.vcf
 * </pre>
 * <br />
 * For full PON creation, call each of your normals separately in artifact detection mode. Then use CombineVariants to
 * output only sites where a variant was seen in at least two samples:
 * <pre>
 * java -jar GenomeAnalysisTK.jar
 *     -T CombineVariants
 *     -R reference.fasta
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
 * <li>MuTect2 currently only supports the calling of a single tumor-normal pair at a time</li>
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

    // fasta reference reader to supplement the edges of the reference sequence
    protected CachingIndexedFastaSequenceFile referenceReader;

    // the assembly engine
    protected LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    protected ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    protected HaplotypeCallerGenotypingEngine genotypingEngine = null;


    private byte MIN_TAIL_QUALITY;
    private double log10GlobalReadMismappingRate;



    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @ArgumentCollection
    protected ReadThreadingAssemblerArgumentCollection RTAC = new ReadThreadingAssemblerArgumentCollection();

    @ArgumentCollection
    protected LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();


    @Argument(fullName = "debug_read_name", required = false, doc="trace this read name through the calling process")
    public String DEBUG_READ_NAME = null;

    @Hidden
    @Advanced
    @Argument(fullName = "MQ_filtering_level", shortName = "MQthreshold", required = false, doc="Set an alternate MQ threshold for debugging")
    final public int MQthreshold = 20;


    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    /**
     * MuTect2 has the ability to use COSMIC data in conjunction with dbSNP to adjust the threshold for evidence of a variant
     * in the normal.  If a variant is present in dbSNP, but not in COSMIC, then more evidence is required from the normal
     * sample to prove the variant is not present in germline.
     */
    @Input(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", required=false)
    public List<RodBinding<VariantContext>> cosmicRod = Collections.emptyList();

    /**
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Input(fullName="normal_panel", shortName = "PON", doc="VCF file of sites observed in normal", required=false)
    public List<RodBinding<VariantContext>> normalPanelRod = Collections.emptyList();

    private HaplotypeBAMWriter haplotypeBAMWriter;

    @Override
    public void initialize() {
        super.initialize();

        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader())));

        // MUTECT: check that we have at least one tumor bam
        for(SAMReaderID id : getToolkit().getReadsDataSource().getReaderIDs()) {
            if (id.getTags().getPositionalTags().size() == 0) {
                throw new RuntimeException("BAMs must be tagged as either 'tumor' or 'normal'");
            }

            // only supports single-sample BAMs (ie first read group is representative)
            String bamSampleName = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();

            for(String tag : id.getTags().getPositionalTags()) {
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
        final GenomeLocParser genomeLocParser = toolkit.getGenomeLocParser();

        genotypingEngine = new SomaticGenotypingEngine( MTAC, samplesList, genomeLocParser, FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(), MTAC, logger), !doNotRunPhysicalPhasing, MTAC);

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
        Set<String> sampleSet = SampleListUtils.asSet(samplesList);

        if( MTAC.CONTAMINATION_FRACTION_FILE != null )
            MTAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(MTAC.CONTAMINATION_FRACTION_FILE, MTAC.CONTAMINATION_FRACTION, sampleSet, logger));

    }

    private VariantAnnotatorEngine initializeVCFOutput() {
        // initialize the output VCF header
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

        Set<VCFHeaderLine> headerInfo = new HashSet<>();

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

        List<String> outputSampleNames = getOutputSampleNames();

        vcfWriter.writeHeader(new VCFHeader(headerInfo, outputSampleNames));

        return annotationEngine;
    }

    private Set<VCFHeaderLine> getM2HeaderLines(){
        Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NORMAL_LOD_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.TUMOR_LOD_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.PANEL_OF_NORMALS_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.HAPLOTYPE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY));

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

        if ( ! doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }
        return headerInfo;
    }

    private Set<VCFHeaderLine> getSampleHeaderLines(){
        Set<VCFHeaderLine> sampleLines = new HashSet<>();
        if (printTCGAsampleHeader) {
            //NOTE: This will only list the first bam file for each tumor/normal sample if there is more than one
            Map<String, String> normalSampleHeaderAttributes = new HashMap<>();
            normalSampleHeaderAttributes.put("ID", "NORMAL");
            normalSampleHeaderAttributes.put("SampleName", normalSampleName);
            if (normalSAMReaderIDs.iterator().hasNext() && !getToolkit().getArguments().disableCommandLineInVCF)
                normalSampleHeaderAttributes.put("File", normalSAMReaderIDs.iterator().next().getSamFilePath());
            VCFSimpleHeaderLine normalSampleHeader = new VCFSimpleHeaderLine("SAMPLE", normalSampleHeaderAttributes);
            Map<String, String> tumorSampleHeaderAttributes = new HashMap<>();
            tumorSampleHeaderAttributes.put("ID", "TUMOR");
            tumorSampleHeaderAttributes.put("SampleName", tumorSampleName);
            if (tumorSAMReaderIDs.iterator().hasNext() && !getToolkit().getArguments().disableCommandLineInVCF)
                tumorSampleHeaderAttributes.put("File", tumorSAMReaderIDs.iterator().next().getSamFilePath());
            VCFSimpleHeaderLine tumorSampleHeader = new VCFSimpleHeaderLine("SAMPLE", tumorSampleHeaderAttributes);

            sampleLines.add(normalSampleHeader);
            sampleLines.add(tumorSampleHeader);
        }
        return sampleLines;
    }

    private List<String> getOutputSampleNames(){
        if (printTCGAsampleHeader) {
         //Already checked for exactly 1 tumor and 1 normal in printTCGAsampleHeader assignment in initialize()
            List<String> sampleNamePlaceholders = new ArrayList<>(2);
            sampleNamePlaceholders.add("TUMOR");
            sampleNamePlaceholders.add("NORMAL");
            return sampleNamePlaceholders;
        }
        else {
            return SampleListUtils.asList(samplesList);
        }
    }

    @Override
    public ActivityProfileState isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if( context == null || context.getBasePileup().isEmpty() )
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getLocus(), 0.0);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        AlignmentContext tumorContext = splitContexts.get(tumorSampleName);
        AlignmentContext normalContext = splitContexts.get(normalSampleName);

        // if there are no tumor reads... there is no activity!
        if (tumorContext == null) {
            return new ActivityProfileState(ref.getLocus(), 0);
        }

        // KCIBUL -- this method was inlined and modified from ReferenceConfidenceModel
        ReadBackedPileup tumorPileup = tumorContext.getBasePileup().getMappingFilteredPileup(MQthreshold);
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
                int nonRefInNormal = getCountOfNonRefEvents(normalContext.getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE);

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
        if ( justDetermineActiveRegions )
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;

        if( !originalActiveRegion.isActive() )
            // Not active so nothing to do!
            return referenceModelForNoVariation(originalActiveRegion, true);

        // No reads here so nothing to do!
        if( originalActiveRegion.size() == 0 ) { return referenceModelForNoVariation(originalActiveRegion, true); }

        logReadInfo(DEBUG_READ_NAME, originalActiveRegion.getReads(), "Present in original active region");

        // create the assembly using just high quality reads (Q20 or higher).  We want to use lower
        // quality reads in the PairHMM (and especially in the normal) later, so we can't use a ReadFilter
        ActiveRegion assemblyActiveRegion = new ActiveRegion(originalActiveRegion.getLocation(), originalActiveRegion.getSupportingStates(),originalActiveRegion.isActive(), getToolkit().getGenomeLocParser(), originalActiveRegion.getExtension());
        for (GATKSAMRecord rec : originalActiveRegion.getReads()) {
            if (rec.getMappingQuality() >= MQthreshold ) {
                assemblyActiveRegion.add(rec);
            }
        }

        logReadInfo(DEBUG_READ_NAME, assemblyActiveRegion.getReads(), "Present in assembly active region");

        // run the local assembler, getting back a collection of information on how we should proceed
        final List<VariantContext> givenAlleles = new ArrayList<>();
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(assemblyActiveRegion, givenAlleles);


        final TreeSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        // TODO - line bellow might be unecessary : it might be that assemblyResult will always have those alleles anyway
        // TODO - so check and remove if that is the case:
        allVariationEvents.addAll(givenAlleles);

        final ActiveRegionTrimmer.Result trimmingResult = trimmer.trim(originalActiveRegion,allVariationEvents);


        // Stop the trimming madness!!!
        if (!trimmingResult.isVariationPresent())
            return referenceModelForNoVariation(originalActiveRegion,false);

        logReadInfo(DEBUG_READ_NAME, trimmingResult.getCallableRegion().getReads(), "Present in trimming result");

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

//        final AssemblyResultSet assemblyResult = untrimmedAssemblyResult;

        // after talking to Ryan -- they grab the reads out of the assembly (and trim then) to pass into the PairHMM
        // because at one point they were trying error correcting of the reads based on the haplotypes.. but that is not
        // working out, so it's safe for us just to take the reads
//
        final ActiveRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        logReadInfo(DEBUG_READ_NAME, regionForGenotyping.getReads(), "Present in region for genotyping");

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.TRIMMMED, originalActiveRegion.getReads(), regionForGenotyping.getReads());
        }

//
//        final ActiveRegion regionForGenotyping = trimmingResult.getCallableRegion();

//        final ActiveRegion regionForGenotyping = originalActiveRegion;

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

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( ! assemblyResult.isVariationPresent() )
            return referenceModelForNoVariation(originalActiveRegion, false);

        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( regionForGenotyping.size() == 0 ) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(originalActiveRegion, false);
        }

        // evaluate each sample's reads against all haplotypes

        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKSAMRecord>> reads = splitReadsBySample( regionForGenotyping.getReads() );
        for (List<GATKSAMRecord> rec : reads.values()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after splitting assemblyResult by sample");
        }

        final HashMap<String, Integer> ARreads_origNormalMQ = new HashMap<>();
        for (GATKSAMRecord read : regionForGenotyping.getReads()) {
            ARreads_origNormalMQ.put(read.getReadName(), read.getMappingQuality());
        }

        // modify MAPQ scores in normal to be high so that we don't do any base quality score capping
        for(GATKSAMRecord rec : regionForGenotyping.getReads()) {
            if (isReadFromNormal(rec)) {
                rec.setMappingQuality(60);
            }
        }

        logger.debug("Computing read likelihoods with " + regionForGenotyping.getReads().size() + " reads against " + haplotypes.size() + " haplotypes across region " + assemblyResult.getRegionForGenotyping().toString());


        // Calculate the likelihoods: CPU intensive part.
        final ReadLikelihoods<Haplotype> readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);

        // Realign reads to their best haplotype.
        // KCIBUL: this is new stuff -- review it!
        final Map<GATKSAMRecord,GATKSAMRecord> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(
                    DroppedReadsTracker.Reason.REALIGNMENT_FAILURE,
                    regionForGenotyping.getReads(),
                    readRealignments.values());
        }

        readLikelihoods.changeReads(readRealignments);

        for (GATKSAMRecord rec : readRealignments.keySet()) {
            logReadInfo(DEBUG_READ_NAME, rec, "Present after computing read likelihoods");
        }

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = ((SomaticGenotypingEngine)genotypingEngine).callMutations(
                haplotypes,
                readLikelihoods,
                ARreads_origNormalMQ,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getLocation(),
                getToolkit().getGenomeLocParser(),
                metaDataTracker,
                givenAlleles, false ,
                tumorSampleName,
                normalSampleName,
                dbsnp.dbsnp,
                cosmicRod,
                DEBUG_READ_NAME
                );

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


        List<VariantContext> annotatedCalls = new ArrayList<>();
        int eventCount = calledHaplotypes.getCalls().size();
        Integer minEventDistance = null;
        Integer maxEventDistance = null;
        Integer lastPosition = null;
        for (VariantContext vc : calledHaplotypes.getCalls()) {
            if (lastPosition == null) {
                lastPosition = vc.getStart();
            } else {
                int dist = Math.abs(vc.getStart() - lastPosition);
                if (maxEventDistance == null || dist > maxEventDistance) {
                    maxEventDistance = dist;
                }
                if (minEventDistance == null || dist < minEventDistance) {
                    minEventDistance = dist;
                }
            }
        }
        Map<String, Object> eventDistanceAttributes = new HashMap<>();
        eventDistanceAttributes.put(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount);
        eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MIN_KEY, minEventDistance);
        eventDistanceAttributes.put(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY, maxEventDistance);


        // can we do this with the Annotation classes instead?
        for (VariantContext originalVC : calledHaplotypes.getCalls()) {
            VariantContextBuilder vcb = new VariantContextBuilder(originalVC);

            Map<String, Object> attributes = new HashMap<>(originalVC.getAttributes());
            attributes.putAll(eventDistanceAttributes);
            vcb.attributes(attributes);

            Set<String> filters = new HashSet<>(originalVC.getFilters());

            double tumorLod = originalVC.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, -1);
            if (tumorLod < MTAC.TUMOR_LOD_THRESHOLD) {
                filters.add(GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
            }

            // if we are in artifact detection mode,  apply the thresholds for the LOD scores
            if (!MTAC.ARTIFACT_DETECTION_MODE) {
                 filters.addAll(calculateFilters(metaDataTracker, originalVC, eventDistanceAttributes));
            }

            if (filters.size() > 0) {
                vcb.filters(filters);
            } else {
                vcb.passFilters();
            }

            if (printTCGAsampleHeader) {
                GenotypesContext genotypesWithBamSampleNames = originalVC.getGenotypes();
                List<Genotype> renamedGenotypes = new ArrayList<>();
                GenotypeBuilder GTbuilder = new GenotypeBuilder(genotypesWithBamSampleNames.get(tumorSampleName));
                GTbuilder.name("TUMOR");
                renamedGenotypes.add(GTbuilder.make());
                GTbuilder = new GenotypeBuilder(genotypesWithBamSampleNames.get(normalSampleName));
                GTbuilder.name("NORMAL");
                renamedGenotypes.add(GTbuilder.make());
                vcb.genotypes(renamedGenotypes);
            }

            annotatedCalls.add(vcb.make());
        }






        return annotatedCalls;
    }

    private Set<String> calculateFilters(RefMetaDataTracker metaDataTracker, VariantContext vc, Map<String, Object> eventDistanceAttributes) {
        Set<String> filters = new HashSet<>();

        Integer eventCount = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY);
        Integer maxEventDistance = (Integer) eventDistanceAttributes.get(GATKVCFConstants.EVENT_DISTANCE_MAX_KEY);

        Collection<VariantContext> panelOfNormalsVC = metaDataTracker.getValues(normalPanelRod,
                getToolkit().getGenomeLocParser().createGenomeLoc(vc.getChr(), vc.getStart()));
        VariantContext ponVc = panelOfNormalsVC.isEmpty()?null:panelOfNormalsVC.iterator().next();

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
            Genotype normalGenotype = vc.getGenotype(normalSampleName);

            // NOTE: how do we get the non-ref depth here?
            normalAltCounts = normalGenotype.getAD()[1];
            normalF = (Double) normalGenotype.getExtendedAttribute(GATKVCFConstants.ALLELE_FRACTION_KEY);

            Object qss = normalGenotype.getExtendedAttribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY);
            if (qss != null) {
                normalAltQualityScoreSum = (Integer) ((Object[]) qss)[1];
            } else {
                logger.error("Null qss at " + vc.getStart());
            }
        }

        if ( (normalAltCounts > MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || normalF > MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION ) && normalAltQualityScoreSum > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
            filters.add(GATKVCFConstants.ALT_ALLELE_IN_NORMAL_FILTER_NAME);
        } else {

            // NOTE: does normal alt counts presume the normal had all these events in CIS?
            if ( eventCount > 1 && normalAltCounts >= 1) {
                filters.add(GATKVCFConstants.MULTI_EVENT_ALT_ALLELE_IN_NORMAL_FILTER_NAME);
            } else if (eventCount >= 3) {
                filters.add(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME);
            }

        }

        // STR contractions, that is the deletion of one repeat unit of a short repeat (>1bp repeat unit)
        // such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we
        // hard filter them out by default
        if (vc.isIndel()) {
            ArrayList rpa = (ArrayList) vc.getAttribute(GATKVCFConstants.REPEATS_PER_ALLELE_KEY);
            String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa != null && rpa.size() > 1 && ru.length() > 1) {
                int refCount = (Integer) rpa.get(0);
                int altCount = (Integer) rpa.get(1);

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
        int AA = 0, AB=1;
        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {

                // TODO: why not use base qualities here?
                //double pobs = QualityUtils.qualToErrorProbLog10(qual);
                double pobs = 1.0d - pow(10, (30 / -10.0));
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
        double f = calculateF(pileup, refBase, minBaseQual);
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
        double f = (double) altCount / ((double) refCount + (double) altCount);
        return f;
    }

    private boolean isNonRef(byte refBase, PileupElement p) {
        return p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
    }

    int MIN_READ_LENGTH = 30; // private in superclass

    protected Set<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {

            // KCIBUL: only perform read quality filtering on tumor reads...
            if (isReadFromNormal(rec)) {

                if( rec.getReadLength() < MIN_READ_LENGTH ) {
                    readsToRemove.add(rec);
                }

            } else {


                if( rec.getReadLength() < MIN_READ_LENGTH ||
                    rec.getMappingQuality() < MQthreshold ||
                    BadMateFilter.hasBadMate(rec) ||

                    (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                    readsToRemove.add(rec);
                }
            }
        }
        activeRegion.removeAll(readsToRemove);
        return readsToRemove;
    }

    private static GATKSAMRecord findReadByName(Collection<GATKSAMRecord> reads, String name) {
        for(GATKSAMRecord read : reads) {
            if (name.equals(read.getReadName())) return read;
        }
        return null;
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        return new PairHMMLikelihoodCalculationEngine( (byte)LEAC.gcpHMM, LEAC.pairHMM, LEAC.pairHMMSub, LEAC.alwaysLoadVectorLoglessPairHMMLib, log10GlobalReadMismappingRate, LEAC.noFpga, pcrErrorModel );
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
//        if ( includeUnmappedReads )
//            throw new UserException.BadArgumentValue("includeUnmappedReads", "is not yet functional");
//        else
            return EnumSet.of(
                    ActiveRegionReadState.PRIMARY,
                    ActiveRegionReadState.NONPRIMARY,
                    ActiveRegionReadState.EXTENDED);
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
    public Integer reduce(List<VariantContext> callsInRegion, Integer numCalledRegions) {
        for( final VariantContext call : callsInRegion ) {
            vcfWriter.add( call );
        }
        return (callsInRegion.isEmpty() ? 0 : 1) + numCalledRegions;
    }

    @Override
    public void onTraversalDone(Integer result) {
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
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP overlap is only used to require more evidence of absence in the normal if the variant in question has been seen before in germline.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

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
//    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"ClippingRankSumTest", "DepthPerSampleHC"}));
    //protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"DepthPerAlleleBySample", "BaseQualitySumPerAlleleBySample", "TandemRepeatAnnotator",
    //    "RMSMappingQuality","MappingQualityRankSumTest","FisherStrand","StrandOddsRatio","ReadPosRankSumTest","QualByDepth", "Coverage"}));
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
    //protected String[] annotationGroupsToUse = { StandardAnnotation.class.getSimpleName() };
    protected String[] annotationClassesToUse = { };

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
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




    /**
     * The minimum confidence needed for a given base for it to be used in variant calling.
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public byte MIN_BASE_QUALTY_SCORE = 10;



// PAIR-HMM-Related Goodness

//    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE;
//    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.AGGRESSIVE;
    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.HOSTILE;

    // Parameters to control read error correction
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




    // reference base padding size
    private static final int REFERENCE_PADDING = 500;

    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;
    private final static int maxReadsInRegionPerSample = 1000; // TODO -- should be an argument
    private final static int minReadsPerAlignmentStart = 5; // TODO -- should be an argument



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
        ReadErrorCorrector readErrorCorrector = null;
        if (errorCorrectReads)
            readErrorCorrector = new ReadErrorCorrector(RTAC.kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION, RTAC.minObservationsForKmerToBeSolid, MTAC.DEBUG, fullReferenceWithPadding);

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
                    //logger.info("Keeping read " + clippedRead + " start " + clippedRead.getAlignmentStart() + " end " + clippedRead.getAlignmentEnd());
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.

        final List<GATKSAMRecord> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        if ( MTAC.bamWriter != null && MTAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.DOWNSAMPLED, activeRegion.getReads(), downsampledReads);
        }

        // handle overlapping read pairs from the same fragment
        // KC: commented out as we handle overlapping read pairs in a different way...
        //cleanOverlappingReadPairs(downsampledReads, normalSampleNames);

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
    private void cleanOverlappingReadPairs(final List<GATKSAMRecord> reads, Set<String> normalSampleNames) {
        Map<String, List<GATKSAMRecord>> data = splitReadsBySample(reads);
        for ( String sampleName : data.keySet() ) {
            final boolean isTumor = !normalSampleNames.contains(sampleName);
            final List<GATKSAMRecord> perSampleReadList = data.get(sampleName);

            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create(perSampleReadList);
            for ( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() )

                // in MuTect -- right now we compare the
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);


        }
    }

    public static void logReadInfo(String readName, Collection<GATKSAMRecord> records, String message) {
        if (readName != null) {
            for (GATKSAMRecord rec : records) {
                logReadInfo(readName, rec, message);
            }

        }
    }

    public static void logReadInfo(String readName, GATKSAMRecord rec, String message) {
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

    private boolean isReadFromNormal(GATKSAMRecord rec) {
        return normalSampleName != null && normalSampleName.equals(rec.getReadGroup().getSample());

    }
    // KCIBUL: new stuff -- read up on this!!
    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName="doNotRunPhysicalPhasing", shortName="doNotRunPhysicalPhasing", doc="Disable physical phasing", required = false)
    protected boolean doNotRunPhysicalPhasing = false;

}


