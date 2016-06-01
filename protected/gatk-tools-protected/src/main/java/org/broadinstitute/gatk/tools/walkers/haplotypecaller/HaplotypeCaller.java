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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.filters.BadMateFilter;
import org.broadinstitute.gatk.engine.io.DirectOutputTracker;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterStub;
import org.broadinstitute.gatk.engine.io.stubs.VariantContextWriterStub;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardHCAnnotation;
import org.broadinstitute.gatk.tools.walkers.genotyper.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
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
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.fragments.FragmentUtils;
import org.broadinstitute.gatk.utils.genotyper.*;
import org.broadinstitute.gatk.utils.gga.GenotypingGivenAllelesUtils;
import org.broadinstitute.gatk.utils.gvcf.GVCFWriter;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.DroppedReadsTracker;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pairhmm.PairHMM;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Call germline SNPs and indels via local re-assembly of haplotypes
 *
 * <p>The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.</p>

<p>In the so-called GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples in a very efficient way, which enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).</p>

 <p>In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. For that purpose, use MuTect2 instead.</p>

 <p>Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers.</p>
 *
 * <h3>How HaplotypeCaller works</h3>
 *
 * <br />
 * <h4>1. Define active regions </h4>
 *
 * <p>The program determines which regions of the genome it needs to operate on, based on the presence of significant
 * evidence for variation.</p>
 *
 * <br />
 * <h4>2. Determine haplotypes by assembly of the active region </h4>
 *
 * <p>For each ActiveRegion, the program builds a De Bruijn-like graph to reassemble the ActiveRegion, and identifies
 * what are the possible haplotypes present in the data. The program then realigns each haplotype against the reference
 * haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites. </p>
 *
 * <br />
 * <h4>3. Determine likelihoods of the haplotypes given the read data </h4>
 *
 * <p>For each ActiveRegion, the program performs a pairwise alignment of each read against each haplotype using the
 * PairHMM algorithm. This produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are
 * then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.   </p>
 *
 * <br />
 * <h4>4. Assign sample genotypes </h4>
 *
 * <p>For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the
 * read data to calculate the likelihoods of each genotype per sample given the read data observed for that
 * sample. The most likely genotype is then assigned to the sample.    </p>
 *
 * <h3>Input</h3>
 * <p>
 * Input bam file(s) from which to make calls
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Either a VCF or gVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant
 * recalibration (best) or hard-filtering before use in downstream analyses. If using the reference-confidence model
 * workflow for cohort analysis, the output is a GVCF file that must first be run through GenotypeGVCFs and then
 * filtering before further analysis.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <p>These are example commands that show how to run HaplotypeCaller for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values shown here may not be the latest recommended; see the
 * Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Single-sample GVCF calling on DNAseq (for `-ERC GVCF` cohort analysis workflow)</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     --emitRefConfidence GVCF \
 *     [--dbsnp dbSNP.vcf] \
 *     [-L targets.interval_list] \
 *     -o output.raw.snps.indels.g.vcf
 * </pre>
 *
 * <h4>Single-sample GVCF calling on DNAseq with allele-specific annotations (for allele-specific cohort analysis workflow)</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     --emitRefConfidence GVCF \
 *     [--dbsnp dbSNP.vcf] \
 *     [-L targets.interval_list] \
 *     -G Standard -G AS_Standard \
 *     -o output.raw.snps.indels.AS.g.vcf
 * </pre>
 *
 * <h4>Variant-only calling on DNAseq</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam [-I sample2.bam ...] \
 *     [--dbsnp dbSNP.vcf] \
 *     [-stand_call_conf 30] \
 *     [-stand_emit_conf 10] \
 *     [-L targets.interval_list] \
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h4>Variant-only calling on RNAseq</h4>
 * <pre>
 *   java -jar GenomeAnalysisTK.jar \
 *     -R reference.fasta \
 *     -T HaplotypeCaller \
 *     -I sample1.bam \
 *     [--dbsnp dbSNP.vcf] \
 *     -stand_call_conf 20 \
 *     -stand_emit_conf 20 \
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>We have not yet fully tested the interaction between the GVCF-based calling or the multisample calling and the
 * RNAseq-specific functionalities. Use those in combination at your own risk.</li>
 * <li>Many users have reported issues running HaplotypeCaller with the -nct argument, so we recommend using Queue to
 * parallelize HaplotypeCaller instead of multithreading.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle almost any ploidy (except very high ploidies in large pooled experiments); the ploidy can be specified using the -ploidy argument for non-diploid organisms.</p>
 *
 * <h3>Additional Notes</h3>
 * <ul>
 *     <li>When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).</li>
 *     <li>When running in `-ERC GVCF` or `-ERC BP_RESOLUTION` modes, the emitting and calling confidence thresholds
 *     are automatically set to 0. This cannot be overridden by the command line. The thresholds can be set manually
 *     to the desired levels in the next step of the workflow (GenotypeGVCFs)</li>
 * </ul>
 *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ActiveRegionTraversalParameters(extension=100, maxRegion=300)
@ReadFilters({HCMappingQualityFilter.class})
@Downsample(by= DownsampleType.BY_SAMPLE, toCoverage=500)
public class HaplotypeCaller extends ActiveRegionWalker<List<VariantContext>, Integer> implements AnnotatorCompatible, NanoSchedulable {
    // -----------------------------------------------------------------------------------------------
    // general haplotype caller arguments
    // -----------------------------------------------------------------------------------------------

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Hidden
    @Advanced
    @Argument(fullName="likelihoodCalculationEngine",shortName="likelihoodEngine",
            doc="What likelihood calculation engine to use to calculate the relative likelihood of reads vs haplotypes",required=false)
    protected ReadLikelihoodCalculationEngine.Implementation likelihoodEngineImplementation = ReadLikelihoodCalculationEngine.Implementation.PairHMM;

    @Hidden
    @Advanced
    @Argument(fullName="heterogeneousKmerSizeResolution",shortName="hksr",doc="How to solve heterogeneous kmer situations using the fast method",required=false)
    protected HeterogeneousKmerSizeResolution heterogeneousKmerSizeResolution = HeterogeneousKmerSizeResolution.COMBO_MIN;

    private HaplotypeBAMWriter haplotypeBAMWriter;

    /**
     * rsIDs from this file are used to populate the ID column of the output. Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    private double log10GlobalReadMismappingRate;

    /**
     * Active region trimmer reference.
     */
    @ArgumentCollection
    protected ActiveRegionTrimmer trimmer = new ActiveRegionTrimmer();

    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Input(fullName="comp", shortName = "comp", doc="Comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * Which annotations to add to the output VCF file. The single value 'none' removes the default annotations.
     * See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<>();

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the
     * -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other
     * options. When HaplotypeCaller is run with -ERC GVCF or -ERC BP_RESOLUTION, some annotations are excluded from the
     * output by default because they will only be meaningful once they have been recalculated by GenotypeGVCFs. As
     * of version 3.3 this concerns ChromosomeCounts, FisherStrand, StrandOddsRatio and QualByDepth.
     *
     */
    @Advanced
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<>();

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{StandardAnnotation.class.getSimpleName(), StandardHCAnnotation.class.getSimpleName() }));

    @ArgumentCollection
    private HaplotypeCallerArgumentCollection HCAC = new HaplotypeCallerArgumentCollection();

    @ArgumentCollection
    private LikelihoodEngineArgumentCollection LEAC = new LikelihoodEngineArgumentCollection();

    /**
     * You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This
     * is especially useful if your samples are all in the same file but you need to run them individually through HC
     * in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     */
    @Argument(fullName="sample_name", shortName = "sn", doc="Name of single sample to use from a multi-sample bam", required=false)
    protected String sampleNameToUse = null;

    @ArgumentCollection
    private ReadThreadingAssemblerArgumentCollection RTAC = new ReadThreadingAssemblerArgumentCollection();

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------


    /**
     * When HC is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar genotype quality (GQ) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     *
     * This argument allows you to set the GQ boundaries. HC expects a list of multiple GQ threshold values. To pass
     * multiple values, you provide them one by one with the argument, as in `-GQB 10 -GQB 20 -GQB 30` and so on. Note
     * that GQ values are capped at 99 in the GATK.
     */
    @Advanced
    @Argument(fullName="GVCFGQBands", shortName="GQB", doc="GQ thresholds for reference confidence bands", required = false)
    protected List<Integer> GVCFGQBands = new ArrayList<Integer>(70) {{
        for (int i=1; i<=60; ++i) add(i);
        add(70); add(80); add(90); add(99);
    }};

    /**
     * This parameter determines the maximum size of an indel considered as potentially segregating in the
     * reference model.  It is used to eliminate reads from being indel informative at a site, and determines
     * by that mechanism the certainty in the reference base.  Conceptually, setting this parameter to
     * X means that each informative read is consistent with any indel of size < X being present at a specific
     * position in the genome, given its alignment to the reference.
     */
    @Advanced
    @Argument(fullName="indelSizeToEliminateInRefModel", shortName="ERCIS", doc="The size of an indel to check for in the reference model", required = false)
    protected int indelSizeToEliminateInRefModel = 10;

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * Bases with a quality below this threshold will not be used for calling.
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public byte MIN_BASE_QUALTY_SCORE = 10;


    /**
     * If this flag is provided, the HaplotypeCaller will include unmapped reads (that have chromosomal coordinates) in the assembly and calling
     * when these reads occur in the region being analyzed.  This situation can occur in paired end analyses, when one read in the read pair
     * gets mapped but its mate is too divergent. In that case, the mate will be marked as unmapped and placed next to the first read, assigned to the same
     * contig and alignment start.  If this flag is provided, the HaplotypeCaller will see such reads, and may make use of them in assembly and calling, where possible.
     */
    @Hidden
    @Argument(fullName="includeUmappedReads", shortName="unmapped", doc="Include unmapped reads with chromosomal coordinates", required = false)
    protected boolean includeUnmappedReads = false;

    @Advanced
    @Argument(fullName="useAllelesTrigger", shortName="allelesTrigger", doc = "Use additional trigger on variants found in an external alleles file", required=false)
    protected boolean USE_ALLELES_TRIGGER = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName="doNotRunPhysicalPhasing", shortName="doNotRunPhysicalPhasing", doc="Disable physical phasing", required = false)
    protected boolean doNotRunPhysicalPhasing = false;

    // -----------------------------------------------------------------------------------------------
    // arguments for debugging / developing the haplotype caller
    // -----------------------------------------------------------------------------------------------

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use reads from this read group when making calls (but use all reads to build the assembly)", required = false)
    protected String keepRG = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName="justDetermineActiveRegions", shortName="justDetermineActiveRegions", doc = "Just determine ActiveRegions, don't perform assembly or calling", required=false)
    protected boolean justDetermineActiveRegions = false;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName="dontGenotype", shortName="dontGenotype", doc = "Perform assembly but do not genotype variants", required=false)
    protected boolean dontGenotype = false;

    @Advanced
    @Argument(fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", doc="Do not analyze soft clipped bases in the reads", required = false)
    protected boolean dontUseSoftClippedBases = false;

    @Hidden
    @Argument(fullName="captureAssemblyFailureBAM", shortName="captureAssemblyFailureBAM", doc="Write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", required = false)
    protected boolean captureAssemblyFailureBAM = false;

    // Parameters to control read error correction
    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName="errorCorrectReads", shortName="errorCorrectReads", doc = "Use an exploratory algorithm to error correct the kmers used during assembly", required=false)
    protected boolean errorCorrectReads = false;

    /**
     * When calculating the likelihood of variants, we can try to correct for PCR errors that cause indel artifacts.
     * The correction is based on the reference context, and acts specifically around repetitive sequences that tend
     * to cause PCR errors). The variant likelihoods are penalized in increasing scale as the context around a
     * putative indel is more repetitive (e.g. long homopolymer). The correction can be disabling by specifying
     * '-pcrModel NONE'; in that case the default base insertion/deletion qualities will be used (or taken from the
     * read if generated through the BaseRecalibrator). <b>VERY IMPORTANT: when using PCR-free sequencing data we
     * definitely recommend setting this argument to NONE</b>.
     */
    @Advanced
    @Argument(fullName = "pcr_indel_model", shortName = "pcrModel", doc = "The PCR indel model to use", required = false)
    public PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL pcrErrorModel = PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE;

    // -----------------------------------------------------------------------------------------------
    // done with Haplotype caller parameters
    // -----------------------------------------------------------------------------------------------

    // the UG engines
    private UnifiedGenotypingEngine activeRegionEvaluationGenotyperEngine = null;

    // the assembly engine
    private LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    protected CachingIndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 500;

    /**
     * When downsampling, level the coverage of the reads in each sample to no more than maxReadsInRegionPerSample reads,
     * not reducing coverage at any read start to less than minReadsPerAlignmentStart
     */
    @Argument(fullName = "maxReadsInRegionPerSample", shortName = "maxReadsInRegionPerSample", doc="Maximum reads in an active region", required = false)
    protected int maxReadsInRegionPerSample = 10000;

    @Argument(fullName = "minReadsPerAlignmentStart", shortName = "minReadsPerAlignStart", doc="Minimum number of reads sharing the same alignment start for each genomic location in an active region", required = false)
    protected int minReadsPerAlignmentStart = 10;

    private byte MIN_TAIL_QUALITY;
    private static final byte MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION = 6;

    /**
     * Minimum (exclusive) average number of high quality bases per soft-clip to consider that a set of soft-clips is a
     * high quality set.
     */
    private static final double AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD = 6.0;

    /**
     * Maximum-mininum confidence on a variant to exist to consider the position as a potential variant harbouring locus
     * when looking for active regions.
     */
    private static final double MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY = 4.0;

    /**
     * Minimum ploidy assumed when looking for loci that may harbour variation to identify active regions.
     * <p>
     * By default we take the putative ploidy provided by the user, but this turned out to be too insensitive
     * for low ploidy, notoriously with haploid samples. Therefore we impose this minimum.
     * </p>
     */
    private static final int MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY = 2;


    /**
     * Reads with length lower than this number, after clipping off overhands outside the active region,
     * won't be considered for genotyping.
     */
    private final static int READ_LENGTH_FILTER_THRESHOLD = 10;

    /**
     * Reads with mapping quality lower than this number won't be considered for genotyping, even if they have
     * passed earlier filters (e.g. the walkers' own min MQ filter).
     */
    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;

    private SampleList samplesList;

    private final static Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private final static Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    ReferenceConfidenceModel referenceConfidenceModel = null;


    ////////////////////////////////////////////////////////////////////////
    //// Deprecated Arguments					                        ////
    //// Keeping them here is meant to provide informative error messages //
    //// when an argument has been put out of service		            ////
    ////////////////////////////////////////////////////////////////////////
    /**
     * @deprecated
     * Deprecated: 2015-04-01, J.White
     * mergeVariantsViaLD = false made final
     */
    @Hidden
    @Deprecated
    @Argument(fullName="mergeVariantsViaLD", shortName="mergeVariantsViaLD", doc="DEPRECATED; This argument is no longer used in GATK versions 3.4 and newer. Please see the online documentation for the latest usage recommendations.", required = false)
    static final boolean mergeVariantsViaLD = false;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        super.initialize();

        if (HCAC.genotypeArgs.samplePloidy != HomoSapiensConstants.DEFAULT_PLOIDY && !doNotRunPhysicalPhasing) {
            doNotRunPhysicalPhasing = true;
            logger.info("Currently, physical phasing is not available when ploidy is different than " + HomoSapiensConstants.DEFAULT_PLOIDY + "; therefore it won't be performed");
        }

        if (dontGenotype && emitReferenceConfidence())
            throw new UserException("You cannot request gVCF output and 'do not genotype' at the same time");

        if ( emitReferenceConfidence() ) {

            if (HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES)
                throw new UserException.BadArgumentValue("ERC/gt_mode","you cannot request reference confidence output and GENOTYPE_GIVEN_ALLELES at the same time");

            HCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = -0.0;
            HCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = -0.0;

            // also, we don't need to output several of the annotations
            annotationsToExclude.add("ChromosomeCounts");
            annotationsToExclude.add("FisherStrand");
            annotationsToExclude.add("StrandOddsRatio");
            annotationsToExclude.add("QualByDepth");

            // but we definitely want certain other ones
            annotationsToUse.add("StrandBiasBySample");
            logger.info("Standard Emitting and Calling confidence set to 0.0 for reference-model confidence output");
            if (!HCAC.annotateAllSitesWithPLs)
                logger.info("All sites annotated with PLs forced to true for reference-model confidence output");
            HCAC.annotateAllSitesWithPLs = true;
        } else if ( ! doNotRunPhysicalPhasing ) {
            doNotRunPhysicalPhasing = true;
            logger.info("Disabling physical phasing, which is supported only for reference-model confidence output");
        }

        final GenomeAnalysisEngine toolkit = getToolkit();
        samplesList = toolkit.getReadSampleList();
        Set<String> sampleSet = SampleListUtils.asSet(samplesList);

        if (sampleNameToUse != null) {
            if (!sampleSet.contains(sampleNameToUse))
                throw new UserException.BadArgumentValue("sample_name", "Specified name does not exist in input bam files");
            if (sampleSet.size() == 1) {
                //No reason to incur performance penalty associated with filtering if they specified the name of the only sample
                sampleNameToUse = null;
            } else {
                samplesList = new IndexedSampleList(sampleNameToUse);
                sampleSet = SampleListUtils.asSet(samplesList);
            }
        }


        // create a UAC but with the exactCallsLog = null, so we only output the log for the HC caller itself, if requested
        final UnifiedArgumentCollection simpleUAC = HCAC.cloneTo(UnifiedArgumentCollection.class);
        simpleUAC.outputMode = OutputMode.EMIT_VARIANTS_ONLY;
        simpleUAC.genotypingOutputMode = GenotypingOutputMode.DISCOVERY;
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, HCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = Math.min(MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY, HCAC.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.CONTAMINATION_FRACTION = 0.0;
        simpleUAC.CONTAMINATION_FRACTION_FILE = null;
        simpleUAC.exactCallsLog = null;
        // Seems that at least with some test data we can lose genuine haploid variation if we use
        // UGs engine with ploidy == 1
        simpleUAC.genotypeArgs.samplePloidy = Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, HCAC.genotypeArgs.samplePloidy);

        activeRegionEvaluationGenotyperEngine = new UnifiedGenotypingEngine(simpleUAC,
                FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(),simpleUAC,logger), toolkit);
        activeRegionEvaluationGenotyperEngine.setLogger(logger);

        if( HCAC.CONTAMINATION_FRACTION_FILE != null )
            HCAC.setSampleContamination(AlleleBiasedDownsamplingUtils.loadContaminationFile(HCAC.CONTAMINATION_FRACTION_FILE, HCAC.CONTAMINATION_FRACTION, sampleSet, logger));

        if( HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES && RTAC.consensusMode )
            throw new UserException("HaplotypeCaller cannot be run in both GENOTYPE_GIVEN_ALLELES mode and in consensus mode at the same time. Please choose one or the other.");

        final GenomeLocParser genomeLocParser = toolkit.getGenomeLocParser();

        genotypingEngine = new HaplotypeCallerGenotypingEngine(HCAC, samplesList, genomeLocParser, FixedAFCalculatorProvider.createThreadSafeProvider(getToolkit(), HCAC,logger), !doNotRunPhysicalPhasing);
        // initialize the output VCF header
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(annotationGroupsToUse, annotationsToUse, annotationsToExclude, this, getToolkit());

        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        //initialize the annotations (this is particularly important to turn off RankSumTest dithering in integration tests)
        //do this before we write the header because SnpEff adds to header lines
        annotationEngine.invokeAnnotationInitializationMethods(headerInfo);

        headerInfo.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());
        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());
        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.DOWNSAMPLED_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        if ( ! doNotRunPhysicalPhasing ) {
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
            headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));
        }

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        initializeReferenceConfidenceModel(samplesList, headerInfo);

        vcfWriter.writeHeader(new VCFHeader(headerInfo, sampleSet));

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
        assemblyEngine.setDebug(HCAC.DEBUG);
        assemblyEngine.setDebugGraphTransformations(RTAC.debugGraphTransformations);
        assemblyEngine.setAllowCyclesInKmerGraphToGeneratePaths(RTAC.allowCyclesInKmerGraphToGeneratePaths);
        assemblyEngine.setRecoverDanglingBranches(!RTAC.doNotRecoverDanglingBranches);
        assemblyEngine.setMinDanglingBranchLength(RTAC.minDanglingBranchLength);
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

        genotypingEngine.setCrossHaplotypeEventMerger(variantMerger);

        genotypingEngine.setAnnotationEngine(annotationEngine);

        if ( HCAC.bamWriter != null ) {
            // we currently do not support multi-threaded BAM writing, so exception out
            if ( getToolkit().getTotalNumberOfThreads() > 1 )
                throw new UserException.BadArgumentValue("bamout", "Currently cannot emit a BAM file from the HaplotypeCaller in multi-threaded mode.");
            haplotypeBAMWriter = HaplotypeBAMWriter.create(HCAC.bamWriterType, HCAC.bamWriter, getToolkit().getSAMFileHeader());
        }

        trimmer.initialize(getToolkit().getGenomeLocParser(), HCAC.DEBUG,
                HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES,emitReferenceConfidence());
    }

    private void initializeReferenceConfidenceModel(final SampleList samples, final Set<VCFHeaderLine> headerInfo) {
        referenceConfidenceModel = new ReferenceConfidenceModel(getToolkit().getGenomeLocParser(), samples, getToolkit().getSAMFileHeader(), indelSizeToEliminateInRefModel);
        if ( emitReferenceConfidence() ) {
            if ( samples.sampleCount() != 1 )
                throw new UserException.BadArgumentValue("emitRefConfidence", "Can only be used in single sample mode currently. Use the sample_name argument to run on a single sample out of a multi-sample BAM file.");
            headerInfo.addAll(referenceConfidenceModel.getVCFHeaderLines());
            if ( HCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
                // A kluge to enforce the use of this indexing strategy - must set the gVCF indexing values if not a using a gVCF output file .
                // An output gVCF file automatically sets the indexing values because it has the .g.vcf extension.
                if (!GATKVCFUtils.usingGVCFIndexingArguments(getToolkit().getArguments().variant_index_type, getToolkit().getArguments().variant_index_parameter) && !isGVCF()) {
                    throw new UserException.GVCFIndexException(GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
                }

                try {
                    vcfWriter = new GVCFWriter(vcfWriter, GVCFGQBands, HCAC.genotypeArgs.samplePloidy);
                } catch ( IllegalArgumentException e ) {
                    throw new UserException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
                }
            }
        }
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        switch (likelihoodEngineImplementation) {
            case PairHMM:
                return new PairHMMLikelihoodCalculationEngine( (byte) LEAC.gcpHMM, LEAC.pairHMM, LEAC.pairHMMSub, LEAC.alwaysLoadVectorLoglessPairHMMLib, log10GlobalReadMismappingRate, LEAC.noFpga, pcrErrorModel );
            case GraphBased:
                return new GraphBasedLikelihoodCalculationEngine( (byte) LEAC.gcpHMM,log10GlobalReadMismappingRate, heterogeneousKmerSizeResolution, HCAC.DEBUG, RTAC.debugGraphTransformations);
            case Random:
                return new RandomLikelihoodCalculationEngine();
            default:
                //Note: we do not include in the error message list as it is of no grand public interest.
                throw new UserException("Unsupported likelihood calculation engine '" + likelihoodCalculationEngine +
                        "'. Please use one of the following instead: 'PairHMM' or 'GraphBased'.");
        }
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
        if ( includeUnmappedReads )
            throw new UserException.BadArgumentValue("includeUnmappedReads", "is not yet functional");
        else
            return EnumSet.of(
                    ActiveRegionReadState.PRIMARY,
                    ActiveRegionReadState.NONPRIMARY,
                    ActiveRegionReadState.EXTENDED);
    }

    @Override
    @Ensures({"result.isActiveProb >= 0.0", "result.isActiveProb <= 1.0"})
    public ActivityProfileState isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {

        if( HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vcFromAllelesRod = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(tracker, ref.getLocus(), false, logger, HCAC.alleles);
            if( vcFromAllelesRod != null ) {
                return new ActivityProfileState(ref.getLocus(), 1.0);
            }
        }

        if( USE_ALLELES_TRIGGER ) {
            return new ActivityProfileState( ref.getLocus(), tracker.getValues(HCAC.alleles, ref.getLocus()).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null || context.getBasePileup().isEmpty() )
            // if we don't have any data, just abort early
            return new ActivityProfileState(ref.getLocus(), 0.0);

        final int ploidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied
        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        final GenotypingModel genotypingModel = genotypingEngine.getGenotypingModel();
        for( final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet() ) {
            final String sampleName = sample.getKey();
            // The ploidy here is not dictated by the sample but by the simple genotyping-engine used to determine whether regions are active or not.
            final int activeRegionDetectionHackishSamplePloidy = activeRegionEvaluationGenotyperEngine.getConfiguration().genotypeArgs.samplePloidy;
            final double[] genotypeLikelihoods = referenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny(activeRegionDetectionHackishSamplePloidy,sample.getValue().getBasePileup(), ref.getBase(), MIN_BASE_QUALTY_SCORE, averageHQSoftClips).genotypeLikelihoods;
            genotypes.add( new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make() );
        }

        final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE , FAKE_ALT_ALLELE);
        final double isActiveProb;

        if (genotypes.size() == 1) {
            // Faster implementation avoiding the costly and over complicated Exact AFCalculator machinery:
            // This is the case when doing GVCF output.
            isActiveProb = activeRegionEvaluationGenotyperEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(genotypes.get(0).getLikelihoods().getAsVector());
        } else {
            final VariantCallContext vcOut = activeRegionEvaluationGenotyperEngine.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.SNP);
            isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb(vcOut.getPhredScaledQual());
        }
        return new ActivityProfileState( ref.getLocus(), isActiveProb, averageHQSoftClips.mean() > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean() );
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

        if (sampleNameToUse != null)
            removeReadsFromAllSamplesExcept(sampleNameToUse, originalActiveRegion);

        if( !originalActiveRegion.isActive() )
            // Not active so nothing to do!
            return referenceModelForNoVariation(originalActiveRegion, true);

        final List<VariantContext> givenAlleles = new ArrayList<>();
        if( HCAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            for ( final VariantContext vc : metaDataTracker.getValues(HCAC.alleles) ) {
                if ( vc.isNotFiltered() ) {
                    givenAlleles.add(vc); // do something with these VCs during GGA mode
                }
            }
            // No alleles found in this region so nothing to do!
            if ( givenAlleles.isEmpty() ) { return referenceModelForNoVariation(originalActiveRegion, true); }
        } else {
            // No reads here so nothing to do!
            if( originalActiveRegion.size() == 0 ) { return referenceModelForNoVariation(originalActiveRegion, true); }
        }

        // run the local assembler, getting back a collection of information on how we should proceed
        final AssemblyResultSet untrimmedAssemblyResult = assembleReads(originalActiveRegion, givenAlleles);

        final TreeSet<VariantContext> allVariationEvents = untrimmedAssemblyResult.getVariationEvents();
        // TODO - line bellow might be unnecessary : it might be that assemblyResult will always have those alleles anyway
        // TODO - so check and remove if that is the case:
        allVariationEvents.addAll(givenAlleles);

        final ActiveRegionTrimmer.Result trimmingResult = trimmer.trim(originalActiveRegion,allVariationEvents);

        if (!trimmingResult.isVariationPresent() && !HCAC.disableOptimizations)
            return referenceModelForNoVariation(originalActiveRegion,false);

        final AssemblyResultSet assemblyResult =
                trimmingResult.needsTrimming() ? untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        final ActiveRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.TRIMMMED, originalActiveRegion.getReads(), regionForGenotyping.getReads());
        }

        // filter out reads from genotyping which fail mapping quality based criteria
        //TODO - why don't do this before any assembly is done? Why not just once at the beginning of this method
        //TODO - on the originalActiveRegion?
        //TODO - if you move this up you might have to consider to change referenceModelForNoVariation
        //TODO - that does also filter reads.
        final Collection<GATKSAMRecord> filteredReads = filterNonPassingReads( regionForGenotyping );

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReads(DroppedReadsTracker.Reason.FILTERED, filteredReads);
        }

        final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );

        // abort early if something is out of the acceptable range
        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( ! assemblyResult.isVariationPresent() && ! HCAC.disableOptimizations)
            return referenceModelForNoVariation(originalActiveRegion, false);

        // For sure this is not true if gVCF is on.
        if (dontGenotype) return NO_CALLS; // user requested we not proceed


        // TODO is this ever true at this point??? perhaps GGA. Need to check.
        if( regionForGenotyping.size() == 0 && ! HCAC.disableOptimizations) {
            // no reads remain after filtering so nothing else to do!
            return referenceModelForNoVariation(originalActiveRegion, false);
        }

        // evaluate each sample's reads against all haplotypes
        //logger.info("Computing read likelihoods with " + assemblyResult.regionForGenotyping.size() + " reads");
        final List<Haplotype> haplotypes = assemblyResult.getHaplotypeList();
        final Map<String,List<GATKSAMRecord>> reads = splitReadsBySample( regionForGenotyping.getReads() );

        // Calculate the likelihoods: CPU intensive part.
        final ReadLikelihoods<Haplotype> readLikelihoods =
                likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult,samplesList,reads);

        // Realign reads to their best haplotype.
        final Map<GATKSAMRecord,GATKSAMRecord> readRealignments = realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc());

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(
                    DroppedReadsTracker.Reason.REALIGNMENT_FAILURE,
                    regionForGenotyping.getReads(),
                    readRealignments.values());
        }

        readLikelihoods.changeReads(readRealignments);

        // Note: we used to subset down at this point to only the "best" haplotypes in all samples for genotyping, but there
        //  was a bad interaction between that selection and the marginalization that happens over each event when computing
        //  GLs.  In particular, for samples that are heterozygous non-reference (B/C) the marginalization for B treats the
        //  haplotype containing C as reference (and vice versa).  Now this is fine if all possible haplotypes are included
        //  in the genotyping, but we lose information if we select down to a few haplotypes.  [EB]

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.assignGenotypeLikelihoods(
                haplotypes,
                readLikelihoods,
                perSampleFilteredReadList,
                assemblyResult.getFullReferenceWithPadding(),
                assemblyResult.getPaddedReferenceLoc(),
                regionForGenotyping.getLocation(),
                getToolkit().getGenomeLocParser(),
                metaDataTracker,
                (RTAC.consensusMode ? Collections.<VariantContext>emptyList() : givenAlleles),
                emitReferenceConfidence());

        if ( HCAC.bamWriter != null ) {
            final Set<Haplotype> calledHaplotypeSet = new HashSet<>(calledHaplotypes.getCalledHaplotypes());
            if (HCAC.disableOptimizations)
                calledHaplotypeSet.add(assemblyResult.getReferenceHaplotype());
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    assemblyResult.getPaddedReferenceLoc(),
                    haplotypes,
                    calledHaplotypeSet,
                    readLikelihoods);

            if ( HCAC.emitDroppedReads ) {
                haplotypeBAMWriter.writeDroppedReads();
            }
        }

        if( HCAC.DEBUG ) { logger.info("----------------------------------------------------------------------------------"); }


        if ( emitReferenceConfidence() ) {
            if ( !containsCalls(calledHaplotypes) ) {
                // no called all of the potential haplotypes
                return referenceModelForNoVariation(originalActiveRegion, false);
            } else {
                final List<VariantContext> result = new LinkedList<>();
                // output left-flanking non-variant section:
                if (trimmingResult.hasLeftFlankingRegion())
                    result.addAll(referenceModelForNoVariation(trimmingResult.nonVariantLeftFlankRegion(),false));
                // output variant containing region.
                result.addAll(referenceConfidenceModel.calculateRefConfidence(assemblyResult.getReferenceHaplotype(),
                        calledHaplotypes.getCalledHaplotypes(), assemblyResult.getPaddedReferenceLoc(), regionForGenotyping,
                        readLikelihoods, genotypingEngine.getPloidyModel(), genotypingEngine.getGenotypingModel(), calledHaplotypes.getCalls()));
                // output right-flanking non-variant section:
                if (trimmingResult.hasRightFlankingRegion())
                    result.addAll(referenceModelForNoVariation(trimmingResult.nonVariantRightFlankRegion(),false));
                return result;
            }
        } else
            return calledHaplotypes.getCalls();
    }

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
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

    private boolean containsCalls(final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes) {
        final List<VariantContext> calls = calledHaplotypes.getCalls();
        if (calls.isEmpty()) return false;
        for (final VariantContext call : calls)
            for (final Genotype genotype : call.getGenotypes())
                if (genotype.isCalled())
                    return true;
        return false;
    }

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
        if( HCAC.DEBUG ) { logger.info("Assembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:    (with overlap region = " + activeRegion.getExtendedLoc() + ")"); }

        final byte[] fullReferenceWithPadding = activeRegion.getActiveRegionReference(referenceReader, REFERENCE_PADDING);
        final GenomeLoc paddedReferenceLoc = getPaddedLoc(activeRegion);
        final Haplotype referenceHaplotype = createReferenceHaplotype(activeRegion, paddedReferenceLoc);

        // Create ReadErrorCorrector object if requested - will be used within assembly engine.
        ReadErrorCorrector readErrorCorrector = null;
        if (errorCorrectReads)
            readErrorCorrector = new ReadErrorCorrector(RTAC.kmerLengthForReadErrorCorrection, MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION, RTAC.minObservationsForKmerToBeSolid, HCAC.DEBUG, fullReferenceWithPadding);

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
     * Create an ref model result (ref model or no calls depending on mode) for an active region without any variation
     * (not is active, or assembled to just ref)
     *
     * @param region the region to return a no-variation result
     * @param needsToBeFinalized should the region be finalized before computing the ref model (should be false if already done)
     * @return a list of variant contexts (can be empty) to emit for this ref region
     */
    private List<VariantContext> referenceModelForNoVariation(final ActiveRegion region, final boolean needsToBeFinalized) {
        if ( emitReferenceConfidence() ) {
            //TODO - why the activeRegion cannot manage its own one-time finalization and filtering?
            //TODO - perhaps we can remove the last parameter of this method and the three lines bellow?
            if ( needsToBeFinalized )
                finalizeActiveRegion(region);
            filterNonPassingReads(region);

            final GenomeLoc paddedLoc = region.getExtendedLoc();
            final Haplotype refHaplotype = createReferenceHaplotype(region, paddedLoc);
            final List<Haplotype> haplotypes = Collections.singletonList(refHaplotype);
            return referenceConfidenceModel.calculateRefConfidence(refHaplotype, haplotypes,
                    paddedLoc, region, createDummyStratifiedReadMap(refHaplotype, samplesList, region),
                    genotypingEngine.getPloidyModel(), genotypingEngine.getGenotypingModel(), Collections.<VariantContext>emptyList());
        } else
            return NO_CALLS;
    }

    /**
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     * @param refHaplotype a non-null reference haplotype
     * @param samples a list of all samples
     * @param region the active region containing reads
     * @return a map from sample -> PerReadAlleleLikelihoodMap that maps each read to ref
     */
    public static ReadLikelihoods<Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                          final SampleList samples,
                                                                          final ActiveRegion region) {
        return new ReadLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                splitReadsBySample(samples, region.getReads()));
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
        genotypingEngine.printFinalMaxNumPLValuesWarning();
        if ( HCAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) ((GVCFWriter)vcfWriter).close(false); // GROSS -- engine forces us to close our own VCF writer since we wrapped it
        referenceConfidenceModel.close();
        //TODO remove the need to call close here for debugging, the likelihood output stream should be managed
        //TODO (open & close) at the walker, not the engine.
        likelihoodCalculationEngine.close();
        logger.info("Ran local assembly on " + result + " active regions");
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // private helper functions
    //
    //---------------------------------------------------------------------------------------------------------------

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if (activeRegion.isFinalized()) return;

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
                clippedRead = ReadClipper.hardClipToRegion( clippedRead, activeRegion.getExtendedLoc().getStart(), activeRegion.getExtendedLoc().getStop() );
                if( activeRegion.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 ) {
                    //logger.info("Keeping read " + clippedRead + " start " + clippedRead.getAlignmentStart() + " end " + clippedRead.getAlignmentEnd());
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.

        final List<GATKSAMRecord> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);

        if ( HCAC.bamWriter != null && HCAC.emitDroppedReads ) {
            haplotypeBAMWriter.addDroppedReadsFromDelta(DroppedReadsTracker.Reason.DOWNSAMPLED, activeRegion.getReads(), downsampledReads);
        }

        // handle overlapping read pairs from the same fragment
        cleanOverlappingReadPairs(downsampledReads);

        activeRegion.clearReads();
        activeRegion.addAll(downsampledReads);
        activeRegion.setFinalized(true);
    }

    private Set<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion ) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getReadLength() < READ_LENGTH_FILTER_THRESHOLD || rec.getMappingQuality() < READ_QUALITY_FILTER_THRESHOLD || BadMateFilter.hasBadMate(rec) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getExtendedLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getExtendedLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getExtendedLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getExtendedLoc().getContig(), padLeft, padRight);
    }

    private Map<String, List<GATKSAMRecord>> splitReadsBySample( final Collection<GATKSAMRecord> reads ) {
        return splitReadsBySample(samplesList, reads);
    }

    public static Map<String, List<GATKSAMRecord>> splitReadsBySample( final SampleList samplesList, final Collection<GATKSAMRecord> reads ) {
        final Map<String, List<GATKSAMRecord>> returnMap = new HashMap<>();
        final int sampleCount = samplesList.sampleCount();
        for (int i = 0; i < sampleCount; i++)
            returnMap.put(samplesList.sampleAt(i), new ArrayList<GATKSAMRecord>());

        for( final GATKSAMRecord read : reads )
            returnMap.get(read.getReadGroup().getSample()).add(read);

        return returnMap;
    }

    /**
     * Are we emitting a reference confidence in some form, or not?
     *
     * @return true if HC must emit reference confidence.
     */
    public boolean emitReferenceConfidence() {
        return HCAC.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private void cleanOverlappingReadPairs(final List<GATKSAMRecord> reads) {
        for ( final List<GATKSAMRecord> perSampleReadList : splitReadsBySample(reads).values() ) {
            final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create(perSampleReadList);
            for ( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() )
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
        }
    }

    private void removeReadsFromAllSamplesExcept(final String targetSample, final ActiveRegion activeRegion) {
        final Set<GATKSAMRecord> readsToRemove = new LinkedHashSet<>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( !rec.getReadGroup().getSample().equals(targetSample) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );

    }

    /**
     * Is writing to an output GVCF file?
     *
     * @return true if the VCF output file has a .g.vcf or .g.vcf.gz extension or if no output file
     */
    private boolean isGVCF() {
        final File file = ((VariantContextWriterStub) vcfWriter).getOutputFile();
        if ( file == null ){
            return true;
        } else {
            final String fileName = file.getName();
            return ( fileName.endsWith("." + GATKVCFUtils.GVCF_EXT) || fileName.endsWith("." + GATKVCFUtils.GVCF_GZ_EXT) );
        }
    }
}
