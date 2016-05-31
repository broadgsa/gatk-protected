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

package org.broadinstitute.gatk.tools.walkers.cancer.contamination;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedGenotypingEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.VariantCallContext;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.*;

import java.io.*;
import java.util.*;

/**
 * Estimate cross-sample contamination
 *
 * <p>This tool determine the percent contamination of an input bam by sample, by lane, or in aggregate across all the input reads.</p>
 *
 * <h3>Usage examples</h3>
 * <p>These are example commands that show how to run ContEst for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values and/or resources shown here may not be the latest recommended; see the Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Contamination estimation using a VCF containing the normal sample's genotypes (as might be derived from a genotyping array)</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar \
 *     -T ContEst \
 *     -R reference.fasta \
 *     -I tumor.bam \
 *     --genotypes normalGenotypes.vcf \
 *     --popFile populationAlleleFrequencies.vcf \
 *     -L populationSites.interval_list
 *     [-L targets.interval_list] \
 *     -isr INTERSECTION \
 *     -o output.txt
 * </pre>
 *
 * <br />
 * <h4>Contamination estimation using the normal BAM for genotyping on-the-fly</h4>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar \
 *     -T ContEst \
 *     -R reference.fasta \
 *     -I:eval tumor.bam \
 *     -I:genotype normal.bam \
 *     --popFile populationAlleleFrequencies.vcf \
 *     -L populationSites.interval_list
 *     [-L targets.interval_list] \
 *     -isr INTERSECTION \
 *     -o output.txt
 * </pre>
 *
 *<h3>Output</h3>
 * A text file containing estimated percent contamination, as well as error bars on this estimate.
 *
 * <h3>Notes</h3>
 * Multiple modes are supported simultaneously, e.g. contamination by sample and readgroup can be computed in the same run.
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Allows(value = {DataSource.READS, DataSource.REFERENCE})
@Requires(value = {DataSource.READS, DataSource.REFERENCE}, referenceMetaData = @RMD(name = "genotypes", type = VariantContext.class))
@By(DataSource.READS)
public class ContEst extends RodWalker<Map<String, Map<String, ContaminationStats>>, ContaminationResults> {

    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // Some constants we use
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    /** what type of run stats would we like: */
    public enum ContaminationRunType {
        SAMPLE,    // calculate contamination for each sample
        READGROUP, // for each read group
        META        // for all inputs as a single source
    }
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // inputs
    // ------------------------------------------------------------------------------------------------------------------------------------------------------

    // the genotypes ROD; this contains information about the genotypes from our sample
    @Input(fullName="genotypes", shortName = "genotypes", doc="the genotype information for our sample", required=false)
    public RodBinding<VariantContext> genotypes;

    // the population information; the allele frequencies for each position in known populations
    @Input(fullName="popfile", shortName = "pf", doc="the variant file containing information about the population allele frequencies", required=true)
    public RodBinding<VariantContext> pop;

    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // outputs and args
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    @Output
    PrintStream out; // the general output of the tool

    @Argument(fullName = "min_qscore", required = false, doc = "threshold for minimum base quality score")
    public int MIN_QSCORE = 20;

    @Argument(fullName = "min_mapq", required = false, doc = "threshold for minimum mapping quality score")
    public int MIN_MAPQ = 20;

    @Argument(fullName = "trim_fraction", doc = "at most, what fraction of sites should be trimmed based on BETA_THRESHOLD", required = false)
    public double TRIM_FRACTION = 0.01;

    @Argument(fullName = "beta_threshold", doc = "threshold for p(f>=0.5) to trim", required = false)
    public double BETA_THRESHOLD = 0.95;

    @Argument(shortName = "llc", fullName = "lane_level_contamination", doc = "set to META (default), SAMPLE or READGROUP to produce per-bam, per-sample or per-lane estimates", required = false)
    private Set<ContaminationRunType> laneStats = null;

    @Argument(shortName = "sn", fullName = "sample_name", doc = "The sample name; used to extract the correct genotypes from mutli-sample truth vcfs", required = false)
    private String sampleName = "unknown";

    @Argument(shortName = "pc", fullName = "precision", doc = "the degree of precision to which the contamination tool should estimate (e.g. the bin size)", required = false)
    private double precision = 0.1;

    @Argument(shortName = "br", fullName = "base_report", doc = "Where to write a full report about the loci we processed", required = false)
    public PrintStream baseReport = null;

    @Argument(shortName = "lf", fullName = "likelihood_file", doc = "write the likelihood values to the specified location", required = false)
    public PrintStream likelihoodFile = null;

    @Argument(shortName = "vs", fullName = "verify_sample", doc = "should we verify that the sample name is in the genotypes file?", required = false)
    public boolean verifySample = false;

    @Argument(shortName = "mbc", fullName = "minimum_base_count", doc = "what minimum number of bases do we need to see to call contamination in a lane / sample?", required = false)
    public Integer minBaseCount = 500;

    @Argument(shortName = "population", fullName = "population", doc = "evaluate contamination for just a single contamination population", required = false)
    public String population = "CEU";

    @Argument(shortName = "gm", fullName = "genotype_mode", doc = "which approach should we take to getting the genotypes (only in array-free mode)", required = false)
    public SeqGenotypeMode genotypeMode = SeqGenotypeMode.HARD_THRESHOLD;

    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // hidden arguments
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    @Hidden
    @Argument(fullName = "trim_interval", doc = "progressively trim from 0 to TRIM_FRACTION by this interval", required = false)
    public double TRIM_INTERVAL = 0;

    @Hidden
    @Argument(fullName = "min_site_depth", required = false, doc = "minimum depth at a site to consider in calculation")
    public int MIN_SITE_DEPTH = 0;

    @Hidden
    @Argument(fullName = "fixed_epsilon_qscore", required = false, doc = "use a constant epsilon (phred scale) for calculation")
    public Byte FIXED_EPSILON = null;

    @Hidden
    @Argument(fullName = "min_genotype_depth", required = false, doc = "what minimum depth is required to call a site in seq genotype mode")
    public int MIN_GENOTYPE_DEPTH_FOR_SEQ = 50;

    @Hidden
    @Argument(fullName = "min_genotype_ratio", required = false, doc = "the ratio of alt to other bases to call a site a hom non-ref variant")
    public double MIN_GENOTYPE_RATIO = 0.80;

    @Hidden
    @Argument(fullName = "min_genotype_llh", required = false, doc = "the min log likelihood for UG to call a genotype")
    public double MIN_UG_LOG_LIKELIHOOD = 5;
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // global variables to the walker
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    private static final Map<Integer,Allele> alleles = new HashMap<Integer,Allele>();                   // the set of alleles we work with
    private boolean verifiedSampleName = false;                                    // have we yet verified the sample name?
    private final Map<String, ContaminationRunType> contaminationNames = new LinkedHashMap<String, ContaminationRunType>();       // a list, containing the contamination names, be it read groups or bam file names
    private static String[] ALL_POPULATIONS = new String[]{"ALL", "CHD", "LWK", "CHB", "CEU", "MXL", "GIH", "MKK", "TSI", "CLM", "GBR", "ASW", "YRI", "IBS", "FIN", "PUR", "JPT", "CHS"};
    private String[] populationsToEvaluate;

    // variables involved in the array-free mode
    private boolean useSequencingGenotypes = false; // if false we're using the sequencing geneotypes; otherwise we require array genotypes
    public static final String EVAL_BAM_TAG = "eval";
    public static final String GENOTYPE_BAM_TAG = "genotype";
    String evalSample = null;
    String genotypeSample = null;


    // counts for each of the possible combinations
    int totalSites = 0;
    int countPopulationSites = 0;
    int countGenotypeNonHomVar = 0;
    int countGenotypeHomVar = 0;
    int countPassCoverage = 0;
    int countResults = 0;

    public enum SeqGenotypeMode { HARD_THRESHOLD, UNIFIED_GENOTYPER }
    // create our list of allele characters for conversion
    static {
        alleles.put(0,Allele.create((byte) 'A'));
        alleles.put(1,Allele.create((byte) 'C'));
        alleles.put(2,Allele.create((byte) 'G'));
        alleles.put(3,Allele.create((byte) 'T'));
    }

    // a bunch of setup to initialize the walker
    public void initialize() {
        // set the genotypes source - figure out what to do if we're not using arrays
        if (genotypes == null || !genotypes.isBound()) {
            logger.info("Running in sequencing mode");
            useSequencingGenotypes = true;
            // if were not using arrays, we need to figure out what samples are what
            for(SAMReaderID id : getToolkit().getReadsDataSource().getReaderIDs()) {
                if (id.getTags().getPositionalTags().isEmpty())
                    throw new UserException.BadInput("BAMs must be tagged with " + GENOTYPE_BAM_TAG + " and " + EVAL_BAM_TAG + " when running in array-free mode. Please see the ContEst documentation for more details");

                // now sort out what tags go with what bam
                for (String tag : id.getTags().getPositionalTags()) {
                    if (GENOTYPE_BAM_TAG.equalsIgnoreCase(tag)) {
                        try {
                            if (getToolkit().getReadsDataSource().getHeader(id).getReadGroups().isEmpty())
                                throw new RuntimeException("No Read Groups found for Genotyping BAM -- Read Groups are Required in sequencing genotype mode!");
                            genotypeSample = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();
                        } catch (NullPointerException npe) {
                            throw new UserException.BadInput("Unable to fetch read group from the bam files tagged with " + GENOTYPE_BAM_TAG);
                        }
                    } else if (EVAL_BAM_TAG.equalsIgnoreCase(tag)) {
                        try {
                            if (getToolkit().getReadsDataSource().getHeader(id).getReadGroups().isEmpty())
                                throw new RuntimeException("No Read Groups found for Genotyping BAM -- Read Groups are Required in sequencing genotype mode!");
                            evalSample = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();
                        } catch (NullPointerException npe) {
                            throw new UserException.BadInput("Unable to fetch read group from the bam files tagged with " + EVAL_BAM_TAG);
                        }
                    } else {
                        throw new UserException.BadInput("Unable to process " + tag + " tag, it's not either of the two accepted values: " + GENOTYPE_BAM_TAG + " or " + EVAL_BAM_TAG);
                    }
                }
            }
            if (evalSample == null || genotypeSample == null)
                throw new UserException.BadInput("You must provide both a " + GENOTYPE_BAM_TAG + " tagged bam and a " + EVAL_BAM_TAG + " tagged bam file.  Please see the ContEst documentation");

        } else {
            logger.info("Running in array mode");
        }
        if (laneStats == null) {
            laneStats = new HashSet<ContaminationRunType>();
            laneStats.add(ContaminationRunType.META);
        }

        for (ContaminationRunType type : laneStats) {
            if (type == ContaminationRunType.READGROUP) {
                for (SAMReadGroupRecord name : getToolkit().getSAMFileHeader().getReadGroups())
                    this.contaminationNames.put(name.getId(),ContaminationRunType.READGROUP);
            } else if (type == ContaminationRunType.SAMPLE) {
                for (SAMReadGroupRecord  name : getToolkit().getSAMFileHeader().getReadGroups())
                    this.contaminationNames.put(name.getSample(),ContaminationRunType.SAMPLE);
            } else if (type == ContaminationRunType.META)
                this.contaminationNames.put("META",ContaminationRunType.META);
            else
                throw new IllegalArgumentException("Unknown type name " + laneStats);
        }
        if (baseReport != null)
            baseReport.println("lane\tchrom\tposition\trs_id\tref\tfreq_major_allele\tfreq_minor_allele\tgeli_gt\tmaf\tmajor_allele_counts\tminor_allele_counts\ta_counts\tc_counts\tg_counts\tt_counts");

        this.populationsToEvaluate = (population == null || "EVERY".equals(population)) ? ALL_POPULATIONS : new String[]{population};

    }
    /**
     * our map function, which emits a contamination stats for each of the subgroups (lanes, samples, etc) that we encounter
     *
     * @param tracker the reference meta data tracker, from which we get the array truth data
     * @param ref     the reference information at this position
     * @param context the read context, where we get the alignment data
     * @return a mapping of our subgroup name to contamination estimate
     */
    @Override
    public Map<String, Map<String, ContaminationStats>> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        totalSites++;
        if (tracker == null) return null;
        if (context == null) return null;

        VariantContext popVC = tracker.getFirstValue(pop);
        byte referenceBase = ref.getBase();
        if (popVC == null) return null;
        countPopulationSites++;
        Genotype genotype = getGenotype(tracker,context,ref,useSequencingGenotypes);

        // only use homozygous sites
        if (genotype == null || !genotype.isHomVar()) {
            countGenotypeNonHomVar++;
            return null;
        } else {
            countGenotypeHomVar++;
        }


        // only use non-reference sites
        byte myBase = genotype.getAllele(0).getBases()[0];

        String rsNumber = "";

        // our map of contamination results
        Map<String, Map<String, ContaminationStats>> stats = new HashMap<String, Map<String, ContaminationStats>>();

        // get the base pileup.  This is only really required when we have both a genotyping and EVAL_BAM_TAG tagged bams
        // becuase we only want contamination estimates drawn from the eval tagged bam
        ReadBackedPileup defaultPile;
        if (this.useSequencingGenotypes)
            defaultPile = context.getBasePileup().getPileupForSample(evalSample);
        else
            defaultPile = context.getBasePileup();

        // if we're by-lane, get those stats
        for (Map.Entry<String, ContaminationRunType> namePair : contaminationNames.entrySet()) {
            ReadBackedPileup pile;
            if (namePair.getValue() == ContaminationRunType.READGROUP)
                pile = defaultPile.getPileupForReadGroup(namePair.getKey());
            else if (namePair.getValue() == ContaminationRunType.META)
                pile = defaultPile;
            else if (namePair.getValue() == ContaminationRunType.SAMPLE)
                pile = defaultPile.getPileupForSample(namePair.getKey());
            else
                throw new IllegalStateException("Unknown state, contamination type = " + laneStats + " is unsupported");
            if (pile != null) {

                ReadBackedPileup filteredPile =
                        pile.getBaseAndMappingFilteredPileup(MIN_QSCORE, MIN_MAPQ);

                byte[] bases = filteredPile.getBases();

                // restrict to sites that have greater than our required total depth
                if (bases.length < MIN_SITE_DEPTH) {
                    continue;
                } else {
                    countPassCoverage++;
                }

                byte[] quals;
                if (FIXED_EPSILON == null) {
                    quals = filteredPile.getQuals();
                } else {
                    quals = new byte[bases.length];
                    Arrays.fill(quals, FIXED_EPSILON);
                }

                Map<String, ContaminationStats> results =
                        calcStats(referenceBase,
                                bases,
                                quals,
                                myBase,
                                rsNumber,
                                popVC,
                                baseReport,
                                context.getLocation(),
                                precision,
                                namePair.getKey(),
                                populationsToEvaluate);

                if (!results.isEmpty()) {
                    countResults++;
                    stats.put(namePair.getKey(), results);
                }
            }
        }
        // return our collected stats
        return stats;
    }

    /**
     * get the genotype for the sample at the current position
     * @param tracker the reference meta data (RODs)
     * @param context the reads
     * @param referenceContext the reference information
     * @param useSeq are we using sequencing to get our genotypes
     * @return a genotype call, which could be null
     */
    private Genotype getGenotype(RefMetaDataTracker tracker, AlignmentContext context, ReferenceContext referenceContext, boolean useSeq) {
        if (!useSeq) {
            Genotype g = getGenotypeFromArray(tracker, this.genotypes,this.verifiedSampleName,this.verifySample,this.sampleName);
            if (g != null) this.verifiedSampleName = true;
            return g;
        } else {
            return getGenotypeFromSeq(
                    context,
                    referenceContext,
                    this.alleles,
                    this.genotypeMode,
                    this.MIN_GENOTYPE_RATIO,
                    this.MIN_GENOTYPE_DEPTH_FOR_SEQ,
                    this.MIN_UG_LOG_LIKELIHOOD,
                    this.genotypeSample,
                    this.sampleName,
                    this.getToolkit());
        }
    }
    
     static Genotype getGenotypeFromSeq(AlignmentContext context,
                                               ReferenceContext referenceContext, 
                                               Map<Integer, Allele> alleles, 
                                               SeqGenotypeMode genotypeMode, 
                                               double minGenotypeRatio,
                                               int minGenotypingDepth,
                                               double minGenotypingLOD,
                                               String genotypingSample,
                                               String sampleName,
                                               GenomeAnalysisEngine toolKit) {
        ReadBackedPileup pileup = context.getBasePileup().getPileupForSample(genotypingSample);
        if (pileup == null || pileup.isEmpty()) return null;

        // which genotyping mode are we using
        if (genotypeMode == SeqGenotypeMode.HARD_THRESHOLD) {
            if (sum(pileup.getBaseCounts()) < minGenotypingDepth) return null;
            int[] bases = pileup.getBaseCounts();
            int mx = maxPos(bases);
            int allGenotypes = sum(bases);
            String refBase = String.valueOf((char)referenceContext.getBase());
            if (bases[mx] / (float)allGenotypes >= minGenotypeRatio && !refBase.equals(alleles.get(mx).getBaseString())) {
                List<Allele> al = new ArrayList<Allele>();
                al.add(alleles.get(mx));
                GenotypeBuilder builder = new GenotypeBuilder(sampleName, al);
                return builder.make();
            }
        } else if (genotypeMode == SeqGenotypeMode.UNIFIED_GENOTYPER) {
            UnifiedArgumentCollection basicUAC = new UnifiedArgumentCollection();
            UnifiedGenotypingEngine engine = new UnifiedGenotypingEngine(basicUAC, FixedAFCalculatorProvider.createThreadSafeProvider(toolKit, basicUAC, logger),toolKit);
            AlignmentContext contextSubset = new AlignmentContext(context.getLocation(),pileup,0,true);
            List<VariantCallContext> callContexts = engine.calculateLikelihoodsAndGenotypes(null, referenceContext, contextSubset);
            if (callContexts != null && callContexts.size() == 1)
                for (Genotype g : callContexts.get(0).getGenotypes()){
                    if (g.isCalled() && g.isHomVar() && g.getLog10PError() > minGenotypingLOD)
                        return g;
                }
        }
        else {
            throw new GATKException("Unknown genotyping mode, being an enum this really shouldn't be seen ever.");
        }
        return null;
    }

    // utils
    private static int sum(int[] a) {int sm = 0; for (int i : a) {sm = sm + i;} return sm;}
    private static int maxPos(int[] a) {int mx = 0; for (int i = 0;i < a.length; i++) {if (a[i] > a[mx]) mx = i;} return mx;}
    
    private static Genotype getGenotypeFromArray(RefMetaDataTracker tracker, RodBinding<VariantContext> genotypes, boolean verifiedSampleName, boolean verifySample, String sampleName) {
        // get the truthForSample and the hapmap information for this site; if either are null we can't move forward
        Collection<VariantContext> truths = tracker.getValues(genotypes);
        if (truths == null || truths.isEmpty()) return null;

        VariantContext truthForSample = truths.iterator().next();

        // verify that the sample name exists in the input genotype file
        if (!verifiedSampleName && verifySample) {
            if (!truthForSample.getSampleNames().contains(sampleName))
                throw new UserException.BadInput("The sample name was set to " + sampleName + " but this sample isn't in your genotypes file.  Please Verify your sample name");
            verifiedSampleName = true;
        }

        GenotypesContext gt = truthForSample.getGenotypes();

        // if we are supposed to verify the sample name, AND the sample doesn't exist in the genotypes -- skip this site
        if (verifySample && !gt.containsSample(sampleName)) return null;

        // if the sample doesn't exist in genotypes AND there is more than one sample in the genotypes file -- skip this site
        if (!gt.containsSample(sampleName) && gt.size() != 1) return null;

        // if there is more than one sample in the genotypes file, get it by name.  Otherwise just get the sole sample genotype
        return gt.size() != 1 ? gt.get(sampleName) : gt.get(0);
    }


    private static class PopulationFrequencyInfo {
        private byte majorAllele;
        private byte minorAllele;
        private double minorAlleleFrequency;

        private PopulationFrequencyInfo(byte majorAllele, byte minorAllele, double minorAlleleFrequency) {
            this.majorAllele = majorAllele;
            this.minorAllele = minorAllele;
            this.minorAlleleFrequency = minorAlleleFrequency;
        }

        public byte getMajorAllele() {
            return majorAllele;
        }

        public byte getMinorAllele() {
            return minorAllele;
        }

        public double getMinorAlleleFrequency() {
            return minorAlleleFrequency;
        }
    }

    private static PopulationFrequencyInfo parsePopulationFrequencyInfo(VariantContext variantContext, String population) {
        PopulationFrequencyInfo info = null;

        List<String> values = (List<String>) variantContext.getAttribute(population);

        if (values != null) {
            byte majorAllele = 0;
            byte minorAllele = 0;
            double maf = -1;

            for (String str : values) {
                // strip off the curly braces and trim whitespace
                if (str.startsWith("{")) str = str.substring(1, str.length());
                if (str.contains("}")) str = str.substring(0, str.indexOf("}"));
                str = str.trim();
                String spl[] = str.split("=");

                byte allele = (byte) spl[0].trim().charAt(0);
                double af = Double.valueOf(spl[1].trim());

                if (af <= 0.5 && minorAllele == 0) {
                    minorAllele = allele;
                    maf = af;
                } else {
                    majorAllele = allele;
                }

            }

            info = new PopulationFrequencyInfo(majorAllele, minorAllele, maf);
        }
        return info;
    }


    /**
     * Calculate the contamination values per division, be it lane, meta, sample, etc
     * @param referenceBase the reference base
     * @param bases the bases seen
     * @param quals and the bases qual values
     * @param myAllele the allele we have (our hom var genotype allele)
     * @param rsNumber the dbsnp number if available
     * @param popVC the population variant context from hapmap
     * @param baseReport if we're writing a base report, write it here
     * @param loc our location
     * @param precision the percision we're aiming for
     * @param lane the lane name information
     * @param pops our pops to run over
     * @return a mapping of each target population to their estimated contamination
     */
    private static Map<String, ContaminationStats> calcStats(byte referenceBase,
                                                             byte[] bases,
                                                             byte[] quals,
                                                             byte myAllele,
                                                             String rsNumber,
                                                             VariantContext popVC,
                                                             PrintStream baseReport,
                                                             GenomeLoc loc,
                                                             Double precision,
                                                             String lane,
                                                             String[] pops) {
        int[] alts = new int[4];
        int total = 0;
        // get the depth ratio we are aiming for
        for (byte base : bases) {
            if (base == 'A' || base == 'a') alts[0]++;
            if (base == 'C' || base == 'c') alts[1]++;
            if (base == 'G' || base == 'g') alts[2]++;
            if (base == 'T' || base == 't') alts[3]++;
            total++;
        }

        Map<String, ContaminationStats> ret = new HashMap<String, ContaminationStats>();

        for (String pop : pops) {
            PopulationFrequencyInfo info = parsePopulationFrequencyInfo(popVC, pop);
            if ( info == null )
                throw new RuntimeException("No population frequency annotation for " + pop + " in " + popVC.toString());

            double alleleFreq = info.getMinorAlleleFrequency();
            if (alleleFreq > 0.5) {
                throw new RuntimeException("Minor allele frequency is greater than 0.5, this is an error; we saw AF of " + alleleFreq);
            }

            int majorCounts = alts[getBaseIndex(info.getMajorAllele())];
            int minorCounts = alts[getBaseIndex(info.getMinorAllele())];
            int otherCounts = total - majorCounts - minorCounts;


            // only use sites where this is the minor allele
            if (myAllele == info.minorAllele) {

                if (pops.length == 1) {
                    if (baseReport != null) {
                        baseReport.print(
                                StringUtil.join("\t",
                                        lane,
                                        loc.getContig(),
                                        "" + loc.getStart(),
                                        rsNumber,
                                        "" + (char) referenceBase,
                                        "" + (char) info.getMajorAllele(),
                                        "" + (char) info.getMinorAllele(),
                                        "" + (char) info.getMinorAllele() + "" + (char) info.getMinorAllele(),
                                        String.format("%1.4f", alleleFreq),
                                        "" + majorCounts,
                                        "" + minorCounts));

                        for (long cnt : alts)
                            baseReport.print("\t" + cnt);
                        baseReport.println();
                    }
                }

                ContaminationEstimate est = new ContaminationEstimate(precision, alleleFreq, bases, quals, info.getMinorAllele(), info.getMajorAllele(), pop, loc);
                ret.put(pop, new ContaminationStats(loc, 1, alleleFreq, minorCounts, majorCounts, otherCounts, alts, est));

            }

        }
        return ret;
    }

    private static int getBaseIndex(byte base) {
        if (base == 'A' || base == 'a') return 0;
        if (base == 'C' || base == 'c') return 1;
        if (base == 'G' || base == 'g') return 2;
        if (base == 'T' || base == 't') return 3;
        return -1;
    }

    // create a ContaminationResults to store the run information
    @Override
    public ContaminationResults reduceInit() {
        return new ContaminationResults(precision);
    }


    @Override
    public ContaminationResults reduce(Map<String, Map<String, ContaminationStats>> value, ContaminationResults sum) {
        if (value != null)
            sum.add(value);
        return sum;
    }

    /**
     * on traversal done, output all the stats to the appropriate files
     *
     * @param result the results of our contamination estimate
     */
    public void onTraversalDone(ContaminationResults result) {

        // filter out lanes / samples that don't have the minBaseCount
        Map<String, Map<String, ContaminationStats>> cleanedMap = new HashMap<String, Map<String, ContaminationStats>>();
        for (Map.Entry<String, Map<String, ContaminationStats>> entry : result.getStats().entrySet()) {

            Map<String, ContaminationStats> newMap = new HashMap<String, ContaminationStats>();

            Map<String, ContaminationStats> statMap = entry.getValue();
            for (String popKey : statMap.keySet()) {
                ContaminationStats stat = statMap.get(popKey);
                if (stat.getBasesMatching() + stat.getBasesMismatching() >= minBaseCount) newMap.put(popKey, stat);
            }


            if (!newMap.isEmpty())
                cleanedMap.put(entry.getKey(), newMap);
            else
                out.println("Warning: We're throwing out lane " + entry.getKey() + " since it has fewer than " + minBaseCount +
                        " read bases at genotyped positions");
        }

        // output results at the end, based on the input parameters
        result.setStats(cleanedMap);
        result.outputReport(precision, out, TRIM_FRACTION, TRIM_INTERVAL, BETA_THRESHOLD);
        if (likelihoodFile != null) result.writeCurves(likelihoodFile);
        logger.info("Total sites:  " + totalSites);
        logger.info("Population informed sites:  " + countPopulationSites);
        logger.info("Non homozygous variant sites: " + countGenotypeNonHomVar);
        logger.info("Homozygous variant sites: " + countGenotypeHomVar);
        logger.info("Passed coverage: " + countPassCoverage);
        logger.info("Results: " + countResults);
    }
}
