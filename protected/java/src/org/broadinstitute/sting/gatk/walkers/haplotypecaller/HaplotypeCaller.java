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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardCallerArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingUtils;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcFactory;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileState;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.haplotype.*;
import org.broadinstitute.sting.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.pairhmm.PairHMM;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.*;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Call SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. Haplotypes are evaluated using an affine gap penalty Pair HMM.
 *
 * <h3>Input</h3>
 * <p>
 * Input bam file(s) from which to make calls
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * VCF file with raw, unrecalibrated SNP and indel calls.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T HaplotypeCaller
 *     -R reference/human_g1k_v37.fasta
 *     -I sample1.bam [-I sample2.bam ...] \
 *     --dbsnp dbSNP.vcf \
 *     -stand_call_conf [50.0] \
 *     -stand_emit_conf 10.0 \
 *     [-L targets.interval_list]
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The system is under active and continuous development. All outputs, the underlying likelihood model, and command line arguments are likely to change often.</li>
 * </ul>
 *
 * @author rpoplin
 * @since 8/22/11
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ActiveRegionTraversalParameters(extension=100, maxRegion=300)
@ReadFilters({HCMappingQualityFilter.class})
@Downsample(by= DownsampleType.BY_SAMPLE, toCoverage=250)
public class HaplotypeCaller extends ActiveRegionWalker<List<VariantContext>, Integer> implements AnnotatorCompatible, NanoSchedulable {
    // -----------------------------------------------------------------------------------------------
    // general haplotype caller arguments
    // -----------------------------------------------------------------------------------------------

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false, defaultToStdout = false)
    protected PrintStream graphWriter = null;

    /**
     * The assembled haplotypes will be written as BAM to this file if requested.  Really for debugging purposes only.
     * Note that the output here does not include uninformative reads so that not every input read is emitted to the bam.
     *
     * Turning on this mode may result in serious performance cost for the HC.  It's really only appropriate to
     * use in specific areas where you want to better understand why the HC is making specific calls.
     *
     * The reads are written out containing a HC tag (integer) that encodes which haplotype each read best matches
     * according to the haplotype caller's likelihood calculation.  The use of this tag is primarily intended
     * to allow good coloring of reads in IGV.  Simply go to Color Alignments By > Tag and enter HC to more
     * easily see which reads go with these haplotype.
     *
     * Note that the haplotypes (called or all, depending on mode) are emitted as single reads covering the entire
     * active region, coming from read HC and a special read group.
     *
     * Note that only reads that are actually informative about the haplotypes are emitted.  By informative we mean
     * that there's a meaningful difference in the likelihood of the read coming from one haplotype compared to
     * its next best haplotype.
     *
     * The best way to visualize the output of this mode is with IGV.  Tell IGV to color the alignments by tag,
     * and give it the HC tag, so you can see which reads support each haplotype.  Finally, you can tell IGV
     * to group by sample, which will separate the potential haplotypes from the reads.  All of this can be seen
     * in the following screenshot: https://www.dropbox.com/s/xvy7sbxpf13x5bp/haplotypecaller%20bamout%20for%20docs.png
     *
     */
    @Advanced
    @Output(fullName="bamOutput", shortName="bamout", doc="File to which assembled haplotypes should be written", required = false, defaultToStdout = false)
    protected StingSAMFileWriter bamWriter = null;
    private HaplotypeBAMWriter haplotypeBAMWriter;

    /**
     * The type of BAM output we want to see.
     */
    @Advanced
    @Argument(fullName="bamWriterType", shortName="bamWriterType", doc="How should haplotypes be written to the BAM?", required = false)
    public HaplotypeBAMWriter.Type bamWriterType = HaplotypeBAMWriter.Type.CALLED_HAPLOTYPES;

    /**
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
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

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>(Arrays.asList(new String[]{"ClippingRankSumTest"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Advanced
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<String>(Arrays.asList(new String[]{"SpanningDeletions", "TandemRepeatAnnotator"}));

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    @ArgumentCollection
    private StandardCallerArgumentCollection SCAC = new StandardCallerArgumentCollection();

    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the debruijn assembler
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName="useDebruijnAssembler", shortName="useDebruijnAssembler", doc="If specified, we will use the old DeBruijn assembler.  Depreciated as of 2.6", required = false)
    protected boolean useDebruijnAssembler = false;

    @Advanced
    @Argument(fullName="minKmerForDebruijnAssembler", shortName="minKmerForDebruijnAssembler", doc="Minimum kmer length to use in the debruijn assembly graph", required = false)
    protected int minKmerForDebruijnAssembler = 11;

    @Advanced
    @Argument(fullName="onlyUseKmerSizeForDebruijnAssembler", shortName="onlyUseKmerSizeForDebruijnAssembler", doc="If specified, we will only build kmer graphs with this kmer size in the debruijn", required = false)
    protected int onlyUseKmerSizeForDebruijnAssembler = -1;

    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the read threading assembler
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName="kmerSize", shortName="kmerSize", doc="Kmer size to use in the read threading assembler", required = false)
    protected List<Integer> kmerSizes = Arrays.asList(10, 25);

    /**
     * Assembly graph can be quite complex, and could imply a very large number of possible haplotypes.  Each haplotype
     * considered requires N PairHMM evaluations if there are N reads across all samples.  In order to control the
     * run of the haplotype caller we only take maxPathsPerSample * nSample paths from the graph, in order of their
     * weights, no matter how many paths are possible to generate from the graph.  Putting this number too low
     * will result in dropping true variation because paths that include the real variant are not even considered.
     */
    @Advanced
    @Argument(fullName="maxPathsPerSample", shortName="maxPathsPerSample", doc="Max number of paths to consider for the read threading assembler per sample.", required = false)
    protected int maxPathsPerSample = 10;

    /**
     * The minimum number of paths to advance forward for genotyping, regardless of the
     * number of samples
     */
    private final static int MIN_PATHS_PER_GRAPH = 128;

    @Hidden
    @Argument(fullName="dontRecoverDanglingTails", shortName="dontRecoverDanglingTails", doc="Should we disable dangling tail recovery in the read threading assembler?", required = false)
    protected boolean dontRecoverDanglingTails = false;

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName="minPruning", shortName="minPruning", doc = "The minimum allowed pruning factor in assembly graph. Paths with <= X supporting kmers are pruned from the graph", required = false)
    protected int MIN_PRUNE_FACTOR = 2;

    @Advanced
    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="Flat gap continuation penalty for use in the Pair HMM", required = false)
    protected int gcpHMM = 10;

    /**
     * If this flag is provided, the haplotype caller will include unmapped reads in the assembly and calling
     * when these reads occur in the region being analyzed.  Typically, for paired end analyses, one pair of the
     * read can map, but if its pair is too divergent then it may be unmapped and placed next to its mate, taking
     * the mates contig and alignment start.  If this flag is provided the haplotype caller will see such reads,
     * and may make use of them in assembly and calling, where possible.
     */
    @Hidden
    @Argument(fullName="includeUmappedReads", shortName="unmapped", doc="If provided, unmapped reads with chromosomal coordinates (i.e., those placed to their maps) will be included in the assembly and calling", required = false)
    protected boolean includeUnmappedReads = false;

    @Advanced
    @Argument(fullName="useAllelesTrigger", shortName="allelesTrigger", doc = "If specified, use additional trigger on variants found in an external alleles file", required=false)
    protected boolean USE_ALLELES_TRIGGER = false;

    @Advanced
    @Argument(fullName="useFilteredReadsForAnnotations", shortName="useFilteredReadsForAnnotations", doc = "If specified, use the contamination-filtered read maps for the purposes of annotating variants", required=false)
    protected boolean USE_FILTERED_READ_MAP_FOR_ANNOTATIONS = false;

    /**
     * The phredScaledGlobalReadMismappingRate reflects the average global mismapping rate of all reads, regardless of their
     * mapping quality.  This term effects the probability that a read originated from the reference haplotype, regardless of
     * its edit distance from the reference, in that the read could have originated from the reference haplotype but
     * from another location in the genome.  Suppose a read has many mismatches from the reference, say like 5, but
     * has a very high mapping quality of 60.  Without this parameter, the read would contribute 5 * Q30 evidence
     * in favor of its 5 mismatch haplotype compared to reference, potentially enough to make a call off that single
     * read for all of these events.  With this parameter set to Q30, though, the maximum evidence against the reference
     * that this (and any) read could contribute against reference is Q30.
     *
     * Set this term to any negative number to turn off the global mapping rate
     */
    @Advanced
    @Argument(fullName="phredScaledGlobalReadMismappingRate", shortName="globalMAPQ", doc="The global assumed mismapping rate for reads", required = false)
    protected int phredScaledGlobalReadMismappingRate = 60;

    @Advanced
    @Argument(fullName="maxNumHaplotypesInPopulation", shortName="maxNumHaplotypesInPopulation", doc="Maximum number of haplotypes to consider for your population. This number will probably need to be increased when calling organisms with high heterozygosity.", required = false)
    protected int maxNumHaplotypesInPopulation = 25;

    @Advanced
    @Argument(fullName="mergeVariantsViaLD", shortName="mergeVariantsViaLD", doc="If specified, we will merge variants together into block substitutions that are in strong local LD", required = false)
    protected boolean mergeVariantsViaLD = false;

    // -----------------------------------------------------------------------------------------------
    // arguments for debugging / developing the haplotype caller
    // -----------------------------------------------------------------------------------------------
    /**
     * The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Hidden
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for genotype likelihood calculations", required = false)
    public PairHMM.HMM_IMPLEMENTATION pairHMM = PairHMM.HMM_IMPLEMENTATION.LOGLESS_CACHING;

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use read from this read group when making calls (but use all reads to build the assembly)", required = false)
    protected String keepRG = null;

    @Hidden
    @Argument(fullName="justDetermineActiveRegions", shortName="justDetermineActiveRegions", doc = "If specified, the HC won't actually do any assembly or calling, it'll just run the upfront active region determination code.  Useful for benchmarking and scalability testing", required=false)
    protected boolean justDetermineActiveRegions = false;

    @Hidden
    @Argument(fullName="dontGenotype", shortName="dontGenotype", doc = "If specified, the HC will do any assembly but won't do calling.  Useful for benchmarking and scalability testing", required=false)
    protected boolean dontGenotype = false;

    @Hidden
    @Argument(fullName="errorCorrectKmers", shortName="errorCorrectKmers", doc = "Use an exploratory algorithm to error correct the kmers used during assembly.  May cause fundamental problems with the assembly graph itself", required=false)
    protected boolean errorCorrectKmers = false;

    @Advanced
    @Argument(fullName="debug", shortName="debug", doc="If specified, print out very verbose debug information about each triggering active region", required = false)
    protected boolean DEBUG;

    @Hidden
    @Argument(fullName="debugGraphTransformations", shortName="debugGraphTransformations", doc="If specified, we will write DOT formatted graph files out of the assembler for only this graph size", required = false)
    protected boolean debugGraphTransformations = false;

    @Hidden // TODO -- not currently useful
    @Argument(fullName="useLowQualityBasesForAssembly", shortName="useLowQualityBasesForAssembly", doc="If specified, we will include low quality bases when doing the assembly", required = false)
    protected boolean useLowQualityBasesForAssembly = false;

    @Hidden
    @Argument(fullName="dontTrimActiveRegions", shortName="dontTrimActiveRegions", doc="If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping", required = false)
    protected boolean dontTrimActiveRegions = false;

    @Hidden
    @Argument(fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", doc="If specified, we will not analyze soft clipped bases in the reads", required = false)
    protected boolean dontUseSoftClippedBases = false;

    @Hidden
    @Argument(fullName="allowCyclesInKmerGraphToGeneratePaths", shortName="allowCyclesInKmerGraphToGeneratePaths", doc="If specified, we will allow cycles in the kmer graphs to generate paths with multiple copies of the path sequenece rather than just the shortest paths", required = false)
    protected boolean allowCyclesInKmerGraphToGeneratePaths = false;

    // -----------------------------------------------------------------------------------------------
    // done with Haplotype caller parameters
    // -----------------------------------------------------------------------------------------------

    // the UG engines
    private UnifiedGenotyperEngine UG_engine = null;
    private UnifiedGenotyperEngine UG_engine_simple_genotyper = null;

    // the assembly engine
    private LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    private LikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    private GenotypingEngine genotypingEngine = null;

    private VariantAnnotatorEngine annotationEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    private CachingIndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 500;

    // include at least this many bases around an event for calling it
    private final static int PADDING_AROUND_SNPS_FOR_CALLING = 20;
    private final static int PADDING_AROUND_OTHERS_FOR_CALLING = 150;

    // the maximum extent into the full active region extension that we're willing to go in genotyping our events
    private final static int MAX_GENOTYPING_ACTIVE_REGION_EXTENSION = 25;

    private ActiveRegionTrimmer trimmer = null;

    private final static int maxReadsInRegionPerSample = 1000; // TODO -- should be an argument
    private final static int minReadsPerAlignmentStart = 5; // TODO -- should be an argument

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 20;

    // the minimum length of a read we'd consider using for genotyping
    private final static int MIN_READ_LENGTH = 10;

    private List<String> samplesList = new ArrayList<String>();
    private final static double LOG_ONE_HALF = -Math.log10(2.0);
    private final static double LOG_ONE_THIRD = -Math.log10(3.0);
    private final List<VariantContext> allelesToGenotype = new ArrayList<VariantContext>();

    private final static Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private final static Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        super.initialize();

        if ( SCAC.AFmodel == AFCalcFactory.Calculation.EXACT_GENERAL_PLOIDY )
            throw new UserException.BadArgumentValue("pnrm", "HaplotypeCaller doesn't currently support " + SCAC.AFmodel);

        // get all of the unique sample names
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        samplesList.addAll( samples );
        final int nSamples = samples.size();
        // initialize the UnifiedGenotyper Engine which is used to call into the exact model
        final UnifiedArgumentCollection UAC = new UnifiedArgumentCollection( SCAC ); // this adapter is used so that the full set of unused UG arguments aren't exposed to the HC user
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples, GATKVariantContextUtils.DEFAULT_PLOIDY);

        // create a UAC but with the exactCallsLog = null, so we only output the log for the HC caller itself, if requested
        UnifiedArgumentCollection simpleUAC = new UnifiedArgumentCollection(UAC);
        simpleUAC.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;
        simpleUAC.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY;
        simpleUAC.STANDARD_CONFIDENCE_FOR_CALLING = Math.min( 4.0, UAC.STANDARD_CONFIDENCE_FOR_CALLING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.STANDARD_CONFIDENCE_FOR_EMITTING = Math.min( 4.0, UAC.STANDARD_CONFIDENCE_FOR_EMITTING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.CONTAMINATION_FRACTION = 0.0;
        simpleUAC.CONTAMINATION_FRACTION_FILE=null;
        simpleUAC.exactCallsLog = null;
        UG_engine_simple_genotyper = new UnifiedGenotyperEngine(getToolkit(), simpleUAC, logger, null, null, samples, GATKVariantContextUtils.DEFAULT_PLOIDY);

        // Currently, per-sample contamination level is only implemented for UG
        if( UAC.CONTAMINATION_FRACTION_FILE !=null) {
            throw new UserException("Per-Sample contamination level not supported in Haplotype Caller at this point");
        }

        // when we do implement per-sample contamination for HC, this will probably be needed.
        // UAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(UAC.CONTAMINATION_FRACTION_FILE, samples, logger));

        // initialize the output VCF header
        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());
        // all callers need to add these standard annotation header lines
        VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true,
                VCFConstants.DOWNSAMPLED_KEY,
                VCFConstants.MLE_ALLELE_COUNT_KEY,
                VCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));

        vcfWriter.writeHeader(new VCFHeader(headerInfo, samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        // create and setup the assembler
        final int maxAllowedPathsForReadThreadingAssembler = Math.max(maxPathsPerSample * nSamples, MIN_PATHS_PER_GRAPH);
        assemblyEngine = useDebruijnAssembler
                ? new DeBruijnAssembler(minKmerForDebruijnAssembler, onlyUseKmerSizeForDebruijnAssembler)
                : new ReadThreadingAssembler(maxAllowedPathsForReadThreadingAssembler, kmerSizes);

        assemblyEngine.setErrorCorrectKmers(errorCorrectKmers);
        assemblyEngine.setPruneFactor(MIN_PRUNE_FACTOR);
        assemblyEngine.setDebug(DEBUG);
        assemblyEngine.setDebugGraphTransformations(debugGraphTransformations);
        assemblyEngine.setAllowCyclesInKmerGraphToGeneratePaths(allowCyclesInKmerGraphToGeneratePaths);
        assemblyEngine.setRecoverDanglingTails(!dontRecoverDanglingTails);

        if ( graphWriter != null ) assemblyEngine.setGraphWriter(graphWriter);
        if ( useLowQualityBasesForAssembly ) assemblyEngine.setMinBaseQualityToUseInAssembly((byte)1);

        // setup the likelihood calculation engine
        if ( phredScaledGlobalReadMismappingRate < 0 ) phredScaledGlobalReadMismappingRate = -1;

        // configure the global mismapping rate
        final double log10GlobalReadMismappingRate;
        if ( phredScaledGlobalReadMismappingRate < 0 ) {
            log10GlobalReadMismappingRate = - Double.MAX_VALUE;
        } else {
            log10GlobalReadMismappingRate = QualityUtils.qualToErrorProbLog10(phredScaledGlobalReadMismappingRate);
            logger.info("Using global mismapping rate of " + phredScaledGlobalReadMismappingRate + " => " + log10GlobalReadMismappingRate + " in log10 likelihood units");
        }

        // create our likelihood calculation engine
        likelihoodCalculationEngine = new LikelihoodCalculationEngine( (byte)gcpHMM, DEBUG, pairHMM, log10GlobalReadMismappingRate );

        final MergeVariantsAcrossHaplotypes variantMerger = mergeVariantsViaLD ? new LDMerger(DEBUG, 10, 1) : new MergeVariantsAcrossHaplotypes();

        genotypingEngine = new GenotypingEngine( DEBUG, annotationEngine, USE_FILTERED_READ_MAP_FOR_ANNOTATIONS, variantMerger );

        if ( bamWriter != null )
            haplotypeBAMWriter = HaplotypeBAMWriter.create(bamWriterType, bamWriter, getToolkit().getSAMFileHeader());

        trimmer = new ActiveRegionTrimmer(DEBUG, PADDING_AROUND_SNPS_FOR_CALLING, PADDING_AROUND_OTHERS_FOR_CALLING,
                MAX_GENOTYPING_ACTIVE_REGION_EXTENSION, getToolkit().getGenomeLocParser());
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // isActive
    //
    //---------------------------------------------------------------------------------------------------------------

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary and extended reads in the active region
    @Override
    public EnumSet<ActiveRegionReadState> desiredReadStates() {
        if ( includeUnmappedReads ) {
            throw new UserException.BadArgumentValue("includeUmappedReads", "is not yet functional");
//            return EnumSet.of(
//                    ActiveRegionReadState.PRIMARY,
//                    ActiveRegionReadState.NONPRIMARY,
//                    ActiveRegionReadState.EXTENDED,
//                    ActiveRegionReadState.UNMAPPED
//            );
        } else
            return EnumSet.of(
                    ActiveRegionReadState.PRIMARY,
                    ActiveRegionReadState.NONPRIMARY,
                    ActiveRegionReadState.EXTENDED
            );
    }

    @Override
    @Ensures({"result.isActiveProb >= 0.0", "result.isActiveProb <= 1.0"})
    public ActivityProfileState isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {

        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vcFromAllelesRod = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, ref.getLocus(), false, logger, UG_engine.getUAC().alleles);
            if( vcFromAllelesRod != null ) {
                allelesToGenotype.add(vcFromAllelesRod); // save for later for processing during the ActiveRegion's map call. Should be folded into a RefMetaDataTracker object
                return new ActivityProfileState(ref.getLocus(), 1.0);
            }
        }

        if( USE_ALLELES_TRIGGER ) {
            return new ActivityProfileState( ref.getLocus(), tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null || context.getBasePileup().isEmpty() )
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getLocus(), 0.0);

        final List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
        noCall.add(Allele.NO_CALL);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        for( final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet() ) {
            final double[] genotypeLikelihoods = new double[3]; // ref versus non-ref (any event)
            Arrays.fill(genotypeLikelihoods, 0.0);

            for( final PileupElement p : sample.getValue().getBasePileup() ) {
                final byte qual = p.getQual();
                if( p.isDeletion() || qual > (byte) 18) {
                    int AA = 0; final int AB = 1; int BB = 2;
                    if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ) {
                        AA = 2;
                        BB = 0;
                        if( p.isNextToSoftClip() ) {
                            averageHQSoftClips.add(AlignmentUtils.calcNumHighQualitySoftClips(p.getRead(), (byte) 28));
                        }
                    }
                    genotypeLikelihoods[AA] += p.getRepresentativeCount() * QualityUtils.qualToProbLog10(qual);
                    genotypeLikelihoods[AB] += p.getRepresentativeCount() * MathUtils.approximateLog10SumLog10( QualityUtils.qualToProbLog10(qual) + LOG_ONE_HALF, QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD + LOG_ONE_HALF );
                    genotypeLikelihoods[BB] += p.getRepresentativeCount() * QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD;
                }
            }
            genotypes.add( new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make() );
        }

        final List<Allele> alleles = new ArrayList<Allele>();
        alleles.add( FAKE_REF_ALLELE );
        alleles.add( FAKE_ALT_ALLELE );
        final VariantCallContext vcOut = UG_engine_simple_genotyper.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.INDEL);
        final double isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb( vcOut.getPhredScaledQual() );

        return new ActivityProfileState( ref.getLocus(), isActiveProb, averageHQSoftClips.mean() > 6.0 ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean() );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    private final static List<VariantContext> NO_CALLS = Collections.emptyList();
    @Override
    public List<VariantContext> map( final ActiveRegion originalActiveRegion, final RefMetaDataTracker metaDataTracker ) {
        if ( justDetermineActiveRegions )
            // we're benchmarking ART and/or the active region determination code in the HC, just leave without doing any work
            return NO_CALLS;

        if( !originalActiveRegion.isActive() ) { return NO_CALLS; } // Not active so nothing to do!

        final List<VariantContext> activeAllelesToGenotype = new ArrayList<>();
        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for( final VariantContext vc : allelesToGenotype ) {
                if( originalActiveRegion.getLocation().overlapsP( getToolkit().getGenomeLocParser().createGenomeLoc(vc) ) ) {
                    activeAllelesToGenotype.add(vc); // do something with these VCs during GGA mode
                }
            }
            allelesToGenotype.removeAll( activeAllelesToGenotype );
            // No alleles found in this region so nothing to do!
            if ( activeAllelesToGenotype.isEmpty() ) { return NO_CALLS; }
        } else {
            if( originalActiveRegion.size() == 0 ) { return NO_CALLS; } // No reads here so nothing to do!
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResult assemblyResult = assembleReads(originalActiveRegion, activeAllelesToGenotype);

        // abort early if something is out of the acceptable range
        if( ! assemblyResult.isVariationPresent() ) { return NO_CALLS; } // only the reference haplotype remains so nothing else to do!
        if (dontGenotype) return NO_CALLS; // user requested we not proceed

        // filter out reads from genotyping which fail mapping quality based criteria
        final List<GATKSAMRecord> filteredReads = filterNonPassingReads( assemblyResult.regionForGenotyping );
        final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );

        if( assemblyResult.regionForGenotyping.size() == 0 ) { return NO_CALLS; } // no reads remain after filtering so nothing else to do!

        // evaluate each sample's reads against all haplotypes
        //logger.info("Computing read likelihoods with " + assemblyResult.regionForGenotyping.size() + " reads");
        final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = likelihoodCalculationEngine.computeReadLikelihoods( assemblyResult.haplotypes, splitReadsBySample( assemblyResult.regionForGenotyping.getReads() ) );

        // subset down to only the best haplotypes to be genotyped in all samples ( in GGA mode use all discovered haplotypes )
        final List<Haplotype> bestHaplotypes = selectBestHaplotypesForGenotyping(assemblyResult.haplotypes, stratifiedReadMap);

        final GenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.assignGenotypeLikelihoods( UG_engine,
                bestHaplotypes,
                stratifiedReadMap,
                perSampleFilteredReadList,
                assemblyResult.fullReferenceWithPadding,
                assemblyResult.paddedReferenceLoc,
                assemblyResult.regionForGenotyping.getLocation(),
                getToolkit().getGenomeLocParser(),
                activeAllelesToGenotype );

        // TODO -- must disable if we are doing NCT, or set the output type of ! presorted
        if ( bamWriter != null ) {
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(assemblyResult.haplotypes, assemblyResult.paddedReferenceLoc,
                    bestHaplotypes,
                    calledHaplotypes.getCalledHaplotypes(),
                    stratifiedReadMap);
        }

        if( DEBUG ) { logger.info("----------------------------------------------------------------------------------"); }

        return calledHaplotypes.getCalls();
    }

    private final static class AssemblyResult {
        final List<Haplotype> haplotypes;
        final ActiveRegion regionForGenotyping;
        final byte[] fullReferenceWithPadding;
        final GenomeLoc paddedReferenceLoc;
        final boolean variationPresent;

        private AssemblyResult(List<Haplotype> haplotypes, ActiveRegion regionForGenotyping, byte[] fullReferenceWithPadding, GenomeLoc paddedReferenceLoc, boolean variationPresent) {
            this.haplotypes = haplotypes;
            this.regionForGenotyping = regionForGenotyping;
            this.fullReferenceWithPadding = fullReferenceWithPadding;
            this.paddedReferenceLoc = paddedReferenceLoc;
            this.variationPresent = variationPresent;
        }

        public boolean isVariationPresent() {
            return variationPresent && haplotypes.size() > 1;
        }
    }

    /**
     * High-level function that runs the assembler on the active region reads,
     * returning a data structure with the resulting information needed
     * for further HC steps
     *
     * @param activeRegion the region we should assemble
     * @param activeAllelesToGenotype additional alleles we might need to genotype (can be empty)
     * @return the AssemblyResult describing how to proceed with genotyping
     */
    protected AssemblyResult assembleReads(final ActiveRegion activeRegion, final List<VariantContext> activeAllelesToGenotype) {
        // Create the reference haplotype which is the bases from the reference that make up the active region
        finalizeActiveRegion(activeRegion); // merge overlapping fragments, clip adapter and low qual tails

        final Haplotype referenceHaplotype = new Haplotype(activeRegion.getActiveRegionReference(referenceReader), true);
        final byte[] fullReferenceWithPadding = activeRegion.getActiveRegionReference(referenceReader, REFERENCE_PADDING);
        final GenomeLoc paddedReferenceLoc = getPaddedLoc(activeRegion);

        final List<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( activeRegion, referenceHaplotype, fullReferenceWithPadding, paddedReferenceLoc, activeAllelesToGenotype );

        if ( ! dontTrimActiveRegions ) {
            return trimActiveRegion(activeRegion, haplotypes, fullReferenceWithPadding, paddedReferenceLoc);
        } else {
            // we don't want to trim active regions, so go ahead and use the old one
            return new AssemblyResult(haplotypes, activeRegion, fullReferenceWithPadding, paddedReferenceLoc, true);
        }
    }

    /**
     * Trim down the active region to just enough to properly genotype the events among the haplotypes
     *
     * @param originalActiveRegion our full active region
     * @param haplotypes the list of haplotypes we've created from assembly
     * @param fullReferenceWithPadding the reference bases over the full padded location
     * @param paddedReferenceLoc the span of the reference bases
     * @return an AssemblyResult containing the trimmed active region with all of the reads we should use
     *         trimmed down as well, and a revised set of haplotypes.  If trimming down the active region results
     *         in only the reference haplotype over the non-extended active region, returns null.
     */
    private AssemblyResult trimActiveRegion(final ActiveRegion originalActiveRegion,
                                            final List<Haplotype> haplotypes,
                                            final byte[] fullReferenceWithPadding,
                                            final GenomeLoc paddedReferenceLoc) {
        if ( DEBUG ) logger.info("Trimming active region " + originalActiveRegion + " with " + haplotypes.size() + " haplotypes");

        EventMap.buildEventMapsForHaplotypes(haplotypes, fullReferenceWithPadding, paddedReferenceLoc, DEBUG);
        final TreeSet<VariantContext> allVariantsWithinFullActiveRegion = EventMap.getAllVariantContexts(haplotypes);
        final ActiveRegion trimmedActiveRegion = trimmer.trimRegion(originalActiveRegion, allVariantsWithinFullActiveRegion);

        if ( trimmedActiveRegion == null ) {
            // there were no variants found within the active region itself, so just return null
            if ( DEBUG ) logger.info("No variation found within the active region, skipping the region :-)");
            return new AssemblyResult(haplotypes, originalActiveRegion, fullReferenceWithPadding, paddedReferenceLoc, false);
        }

        // trim down the haplotypes
        final Set<Haplotype> haplotypeSet = new HashSet<Haplotype>(haplotypes.size());
        for ( final Haplotype h : haplotypes ) {
            final Haplotype trimmed = h.trim(trimmedActiveRegion.getExtendedLoc());
            if ( trimmed != null ) {
                haplotypeSet.add(trimmed);
            } else if ( DEBUG ) {
                logger.info("Throwing out haplotype " + h + " with cigar " + h.getCigar() + " because it starts with or ends with an insertion or deletion when trimmed to " + trimmedActiveRegion.getExtendedLoc());
            }
        }

        // create the final list of trimmed haplotypes
        final List<Haplotype> trimmedHaplotypes = new ArrayList<Haplotype>(haplotypeSet);

        // sort haplotypes to take full advantage of haplotype start offset optimizations in PairHMM
        Collections.sort( trimmedHaplotypes, new HaplotypeBaseComparator() );

        if ( DEBUG ) logger.info("Trimmed region to " + trimmedActiveRegion.getLocation() + " size " + trimmedActiveRegion.getLocation().size() + " reduced number of haplotypes from " + haplotypes.size() + " to only " + trimmedHaplotypes.size());
        if ( DEBUG ) {
            for ( final Haplotype remaining: trimmedHaplotypes ) {
                logger.info("  Remains: " + remaining + " cigar " + remaining.getCigar());
            }
        }


        // trim down the reads and add them to the trimmed active region
        final List<GATKSAMRecord> trimmedReads = new ArrayList<GATKSAMRecord>(originalActiveRegion.getReads().size());
        for( final GATKSAMRecord read : originalActiveRegion.getReads() ) {
            final GATKSAMRecord clippedRead = ReadClipper.hardClipToRegion( read, trimmedActiveRegion.getExtendedLoc().getStart(), trimmedActiveRegion.getExtendedLoc().getStop() );
            if( trimmedActiveRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                trimmedReads.add(clippedRead);
            }
        }
        trimmedActiveRegion.clearReads();
        trimmedActiveRegion.addAll(ReadUtils.sortReadsByCoordinate(trimmedReads));

        return new AssemblyResult(trimmedHaplotypes, trimmedActiveRegion, fullReferenceWithPadding, paddedReferenceLoc, true);
    }

    /**
     * Select the best N haplotypes according to their likelihoods, if appropriate
     *
     * @param haplotypes a list of haplotypes to consider
     * @param stratifiedReadMap a map from samples -> read likelihoods
     * @return the list of haplotypes to genotype
     */
    protected List<Haplotype> selectBestHaplotypesForGenotyping(final List<Haplotype> haplotypes, final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap) {
        if ( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            return haplotypes;
        } else {
            return likelihoodCalculationEngine.selectBestHaplotypesFromEachSample(haplotypes, stratifiedReadMap, maxNumHaplotypesInPopulation);
        }
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
            // TODO -- uncomment this line once ART-based walkers have a proper RefMetaDataTracker.
            // annotationEngine.annotateDBs(metaDataTracker, getToolkit().getGenomeLocParser().createGenomeLoc(call),  call);
            vcfWriter.add( call );
        }
        return (callsInRegion.isEmpty() ? 0 : 1) + numCalledRegions;
    }

    @Override
    public void onTraversalDone(Integer result) {
        logger.info("Ran local assembly on " + result + " active regions");
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // private helper functions
    //
    //---------------------------------------------------------------------------------------------------------------

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if( DEBUG ) { logger.info("Assembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:    (with overlap region = " + activeRegion.getExtendedLoc() + ")"); }
        final List<GATKSAMRecord> finalizedReadList = new ArrayList<>();
        final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create( activeRegion.getReads() );
        activeRegion.clearReads();

        // Join overlapping paired reads to create a single longer read
        finalizedReadList.addAll( fragmentCollection.getSingletonReads() );
        for( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
            finalizedReadList.addAll( FragmentUtils.mergeOverlappingPairedFragments(overlappingPair) );
        }

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKSAMRecord> readsToUse = new ArrayList<>(finalizedReadList.size());
        for( final GATKSAMRecord myRead : finalizedReadList ) {
            final GATKSAMRecord postAdapterRead = ( myRead.getReadUnmappedFlag() ? myRead : ReadClipper.hardClipAdaptorSequence( myRead ) );
            if( postAdapterRead != null && !postAdapterRead.isEmpty() && postAdapterRead.getCigar().getReadLength() > 0 ) {
                GATKSAMRecord clippedRead = useLowQualityBasesForAssembly ? postAdapterRead : ReadClipper.hardClipLowQualEnds( postAdapterRead, MIN_TAIL_QUALITY );

                if ( dontUseSoftClippedBases ) {
                    // uncomment to remove hard clips from consideration at all
                    clippedRead = ReadClipper.hardClipSoftClippedBases(clippedRead);
                } else {
                    // revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
                    // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
                    // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
                    // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
                    // TODO -- reference haplotype start must be removed
                    clippedRead = ReadClipper.revertSoftClippedBases(clippedRead);
                }

                clippedRead = ReadClipper.hardClipToRegion( clippedRead, activeRegion.getExtendedLoc().getStart(), activeRegion.getExtendedLoc().getStop() );
                if( activeRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                    //logger.info("Keeping read " + clippedRead + " start " + clippedRead.getAlignmentStart() + " end " + clippedRead.getAlignmentEnd());
                    readsToUse.add(clippedRead);
                }
            }
        }

        activeRegion.addAll(DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart));
    }

    private List<GATKSAMRecord> filterNonPassingReads( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion ) {
        final List<GATKSAMRecord> readsToRemove = new ArrayList<>();
//        logger.info("Filtering non-passing regions: n incoming " + activeRegion.getReads().size());
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getReadLength() < MIN_READ_LENGTH || rec.getMappingQuality() < 20 || BadMateFilter.hasBadMate(rec) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
//                logger.info("\tremoving read " + rec + " len " + rec.getReadLength());
            }
        }
        activeRegion.removeAll( readsToRemove );
//        logger.info("Filtered non-passing regions: n remaining " + activeRegion.getReads().size());
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getExtendedLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getExtendedLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getExtendedLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getExtendedLoc().getContig(), padLeft, padRight);
    }

    private Map<String, List<GATKSAMRecord>> splitReadsBySample( final List<GATKSAMRecord> reads ) {
        final Map<String, List<GATKSAMRecord>> returnMap = new HashMap<String, List<GATKSAMRecord>>();
        for( final String sample : samplesList) {
            List<GATKSAMRecord> readList = returnMap.get( sample );
            if( readList == null ) {
                readList = new ArrayList<>();
                returnMap.put(sample, readList);
            }
        }
        for( final GATKSAMRecord read : reads ) {
            returnMap.get(read.getReadGroup().getSample()).add(read);
        }

        return returnMap;
    }


}