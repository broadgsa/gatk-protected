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

package org.broadinstitute.gatk.tools.walkers.phasing;

import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.MappingQualityZeroFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.HasGenomeLocation;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.*;
import java.util.*;

import static org.broadinstitute.gatk.engine.GATKVCFUtils.getVCFHeadersFromRods;

/**
 * Annotate physical phasing information
 *
 * <p>This tool identifies haplotypes based on the overlap between reads and uses this information to generate physical
 * phasing information for variants within these haplotypes.</p>
 *
 * <p>It operates by walking along all variant ROD loci, caching a user-defined window of VariantContext sites, and
 * then finishes phasing them when they go out of range (using upstream and downstream reads). The underlying algorithm
 * is based on building up 2^n local haplotypes, where n is the number of heterozygous SNPs in the local region we
 * expected to find phase-informative reads (and assumes a maximum value of maxPhaseSites, a user parameter). Then,
 * these 2^n haplotypes are used to determine, with sufficient certainty (the assigned PQ score), to which haplotype
 * the alleles of a genotype at a particular locus belong (denoted by the HP tag).</p>
 *
 * <p>
 * Performs physical phasing of SNP calls, based on sequencing reads.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * VCF file of SNP calls, BAM file of sequence reads.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Phased VCF file using HP tags to link alleles at (possibly non-consecutive) genotypes of the same sample.
 * </p>
 * <h4>Example</h4>
 * <pre>
 *     GT:GQ:HP    0/1:99:17690409-1,17690409-2
 *     GT:GQ:HP    0/1:99:17690409-2,17690409-1:1258.14
 * </pre>
 * <p>The second site's alternate allele (1) is on the same physical haplotype as the first site's reference allele (0),
 * and vice versa [second site's 0 goes with first site's 1]. This is based on the fact that the HP pairs line up in
 * reverse order between these two genotypes.</p>
 * <p>In an old notation that RBP used to output in much earlier versions, the genotypes would have been: 0/1 and 1|0,
 * respectively. This was changed because depending on the case it caused ambiguity, incompleteness, and possible
 * inconsistency with trio-based phasing. In contrast, the HP tag is much more explicitl for linking alleles, especially
 * if the genotypes are non-consecutive.</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 *    java -jar GenomeAnalysisTK.jar \
 *      -T ReadBackedPhasing \
 *      -R reference.fasta \
 *      -I reads.bam \
 *      --variant SNPs.vcf \
 *      -L SNPs.vcf \
 *      -o phased_SNPs.vcf \
 *      --phaseQualityThresh 20.0
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>The current implementation works for diploid SNPs, and will transparently (but properly) ignore other sites.</p>
 *
 * @author Menachem Fromer
 * @since July 2010
 */
@Allows(value = {DataSource.READS, DataSource.REFERENCE})
@Requires(value = {DataSource.READS, DataSource.REFERENCE})
@By(DataSource.READS)

// Filter out all reads with zero mapping quality
@ReadFilters({MappingQualityZeroFilter.class})

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class ReadBackedPhasing extends RodWalker<PhasingStatsAndOutput, PhasingStats> {
    @Argument(fullName="debug", shortName="debug", doc="If specified, print out very verbose debug information (if -l DEBUG is also specified)", required = false)
    protected boolean DEBUG = false;
    /**
     * The VCF file we are phasing variants from.
     *
     * All heterozygous variants found in this VCF file will be phased, where possible
     */
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc = "File to which variants should be written")
    protected VariantContextWriter writer = null;

    @Argument(fullName = "cacheWindowSize", shortName = "cacheWindow", doc = "The window size (in bases) to cache variant sites and their reads for the phasing procedure", required = false)
    protected Integer cacheWindow = 20000;

    @Argument(fullName = "maxPhaseSites", shortName = "maxSites", doc = "The maximum number of successive heterozygous sites permitted to be used by the phasing algorithm", required = false)
    protected Integer maxPhaseSites = 10; // 2^10 == 10^3 diploid haplotypes

    @Argument(fullName = "phaseQualityThresh", shortName = "phaseThresh", doc = "The minimum phasing quality score required to output phasing", required = false)
    protected Double phaseQualityThresh = 20.0; // PQ = 20.0 <=> P(error) = 10^(-20/10) = 0.01, P(correct) = 0.99

    @Hidden
    @Argument(fullName = "variantStatsFilePrefix", shortName = "variantStats", doc = "The prefix of the VCF/phasing statistics files [For DEBUGGING purposes only - DO NOT USE!]", required = false)
    protected String variantStatsFilePrefix = null;
    private PhasingQualityStatsWriter statsWriter = null;

    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for phasing", required = false)
    public int MIN_BASE_QUALITY_SCORE = 17;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for phasing", required = false)
    public int MIN_MAPPING_QUALITY_SCORE = 20;

    @Argument(fullName = "sampleToPhase", shortName = "sampleToPhase", doc = "Only include these samples when phasing", required = false)
    protected Set<String> samplesToPhase = null;

    @Hidden
    @Argument(fullName = "permitNoSampleOverlap", shortName = "permitNoSampleOverlap", doc = "Don't exit (just WARN) when the VCF and BAMs do not overlap in samples", required = false)
    private boolean permitNoSampleOverlap = false;

    private GenomeLoc mostDownstreamLocusReached = null;

    private LinkedList<VariantAndReads> unphasedSiteQueue = null;
    private CloneableIteratorLinkedList<UnfinishedVariantAndReads> partiallyPhasedSites = null; // the phased VCs to be emitted, and the alignment bases at these positions

    private static PreciseNonNegativeDouble ZERO = new PreciseNonNegativeDouble(0.0);

    // In order to detect phase inconsistencies:
    private static final double FRACTION_OF_MEAN_PQ_CHANGES = 0.1; // If the PQ decreases by this fraction of the mean PQ changes (thus far), then this read is inconsistent with previous reads
    private static final double MAX_FRACTION_OF_INCONSISTENT_READS = 0.1; // If there are more than this fraction of inconsistent reads, then flag this site

    @Argument(fullName = "enableMergePhasedSegregatingPolymorphismsToMNP", shortName = "enableMergeToMNP", doc = "Merge consecutive phased sites into MNP records", required = false)
    protected boolean enableMergePhasedSegregatingPolymorphismsToMNP = false;

    @Argument(fullName = "maxGenomicDistanceForMNP", shortName = "maxDistMNP", doc = "The maximum reference-genome distance between consecutive heterozygous sites to permit merging phased VCF records into a MNP record", required = false)
    protected int maxGenomicDistanceForMNP = 1;

    @Hidden
    @Argument(fullName = "outputMultipleBaseCountsFile", shortName = "outputMultipleBaseCountsFile", doc = "File to output cases where a single read has multiple bases at the same position [For DEBUGGING purposes only - DO NOT USE!]", required = false)
    protected File outputMultipleBaseCountsFile = null;
    private MultipleBaseCountsWriter outputMultipleBaseCountsWriter = null;

    public void initialize() {
        if (maxPhaseSites <= 2)
            maxPhaseSites = 2; // by definition, must phase a site relative to previous site [thus, 2 in total]

        /*
         Since we cap each base quality (BQ) by its read's mapping quality (MQ) [in Read.updateBaseAndQuality()], then:
         if minBQ > minMQ, then we require that MQ be >= minBQ as well.
         [Otherwise, we end up capping BQ by MQ only AFTER we tried removing bases with BQ < minBQ, which is WRONG!]

         To do this properly, we set: minMQ = max(minMQ, minBQ)
         */
        MIN_MAPPING_QUALITY_SCORE = Math.max(MIN_MAPPING_QUALITY_SCORE, MIN_BASE_QUALITY_SCORE);

        unphasedSiteQueue = new LinkedList<VariantAndReads>();
        partiallyPhasedSites = new CloneableIteratorLinkedList<UnfinishedVariantAndReads>();

        initializeVcfWriter();

        if (variantStatsFilePrefix != null)
            statsWriter = new PhasingQualityStatsWriter(variantStatsFilePrefix);

        if (outputMultipleBaseCountsFile != null)
            outputMultipleBaseCountsWriter = new MultipleBaseCountsWriter(outputMultipleBaseCountsFile);
    }

    private void initializeVcfWriter() {
        // Wrapper VCFWriters will take ownership of inner writers iff: inner writer != origWriter [which wasn't created here]
        VariantContextWriter origWriter = writer;

        if (enableMergePhasedSegregatingPolymorphismsToMNP)
            writer = new MergeSegregatingAlternateAllelesVCFWriter(writer, getToolkit().getGenomeLocParser(), getToolkit().getArguments().referenceFile, maxGenomicDistanceForMNP, logger, writer != origWriter);

        /* Due to discardIrrelevantPhasedSites(), the startDistance spanned by [partiallyPhasedSites.peek(), unphasedSiteQueue.peek()] is <= cacheWindow
           Due to processQueue(), the startDistance spanned by [unphasedSiteQueue.peek(), mostDownstreamLocusReached] is <= cacheWindow
           Hence, the startDistance between: partiallyPhasedSites.peek() --> mostDownstreamLocusReached is <= 2 * cacheWindow

           Therefore, can write the filtered records located at mostDownstreamLocusReached (if any) to SortingVCFWriter, even though partiallyPhasedSites.peek() has not yet been written.

           But, NOTE that map() is careful to pass out a list of records to be written that FIRST includes any records discarded due to having reached mostDownstreamLocusReached,
           and only THEN records located at mostDownstreamLocusReached.  The opposite order in map() would violate the startDistance limits imposed when contracting SortingVCFWriter with (2 * cacheWindow).
         */
        writer = VariantContextWriterFactory.sortOnTheFly(writer, 2 * cacheWindow, writer != origWriter);

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        // Phasing-specific INFO fields:
        hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.PHASE_QUALITY_KEY, true));
        hInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.RBP_HAPLOTYPE_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RBP_INCONSISTENT_KEY));

        // todo -- fix samplesToPhase
        String trackName = variantCollection.variants.getName();
        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));
        Set<String> vcfSamples = new TreeSet<String>(samplesToPhase == null ? rodNameToHeader.get(trackName).getGenotypeSamples() : samplesToPhase);
        writer.writeHeader(new VCFHeader(hInfo, vcfSamples));

        Set<String> readSamples = ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        readSamples.retainAll(vcfSamples);
        if (readSamples.isEmpty()) {
            String noPhaseString = "No common samples in VCF and BAM headers" + (samplesToPhase == null ? "" : " (limited to sampleToPhase parameters)") + ", so nothing could possibly be phased!";
            if (permitNoSampleOverlap)
                logger.warn(noPhaseString);
            else
                throw new UserException(noPhaseString);
        }
    }

    public PhasingStats reduceInit() {
        return new PhasingStats();
    }

    /**
     * For each site of interest, cache the current site and then use the cache to phase all sites
     * for which "sufficient" information has already been observed.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public PhasingStatsAndOutput map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        mostDownstreamLocusReached = ref.getLocus();
        if (DEBUG) logger.debug("map() at: " + mostDownstreamLocusReached);

        PhasingStats phaseStats = new PhasingStats();
        List<VariantContext> unprocessedList = new LinkedList<VariantContext>();

        for (VariantContext vc : tracker.getValues(variantCollection.variants, context.getLocation())) {
            if (samplesToPhase != null) vc = reduceVCToSamples(vc, samplesToPhase);

            if (ReadBackedPhasing.processVariantInPhasing(vc)) {
                VariantAndReads vr = new VariantAndReads(vc, context);
                unphasedSiteQueue.add(vr);

                if (DEBUG)
                    logger.debug("Added variant to queue = " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant));
            }
            else {
                unprocessedList.add(vc); // Finished with the unprocessed variant, and writer can enforce sorting on-the-fly

                if (DEBUG)
                    logger.debug("Unprocessed variant = " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));
            }

            int numReads = context.getBasePileup().getNumberOfElements();
            PhasingStats addInPhaseStats = new PhasingStats(numReads, 1);
            phaseStats.addIn(addInPhaseStats);
        }

        List<VariantContext> completedList = processQueue(phaseStats, false);
        completedList.addAll(unprocessedList); // add unprocessedList on to the END of completedList so that the processQueue() results, which are necessarily more upstream, are first!

        return new PhasingStatsAndOutput(phaseStats, completedList);
    }

    private static final Set<String> KEYS_TO_KEEP_IN_REDUCED_VCF = new HashSet<>(Arrays.asList(VCFConstants.PHASE_QUALITY_KEY));

    private VariantContext reduceVCToSamples(VariantContext vc, Set<String> samplesToPhase) {
//        for ( String sample : samplesToPhase )
//            logger.debug(String.format("  Sample %s has genotype %s, het = %s", sample, vc.getGenotype(sample), vc.getGenotype(sample).isHet() ));
        VariantContext subvc = vc.subContextFromSamples(samplesToPhase);
//        logger.debug("original VC = " + vc);
//        logger.debug("sub      VC = " + subvc);
        return GATKVariantContextUtils.pruneVariantContext(subvc, KEYS_TO_KEEP_IN_REDUCED_VCF);
    }

    // Phase all "waiting" genotypes in the unphasedSiteQueue, but only if we have sufficient downstream genotypes with which to phase them
    private List<VariantContext> processQueue(PhasingStats phaseStats, boolean processAll) {
        List<VariantContext> oldPhasedList = new LinkedList<VariantContext>();

        VariantAndReads prevVr = null;
        while (!unphasedSiteQueue.isEmpty()) {
            if (!processAll) { // otherwise, phase until the end of unphasedSiteQueue
                VariantContext nextToPhaseVc = unphasedSiteQueue.peek().variant;
                if (startDistancesAreInWindowRange(mostDownstreamLocusReached, GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), nextToPhaseVc))) {
                    /* mostDownstreamLocusReached is still not far enough ahead of nextToPhaseVc to have all phasing information for nextToPhaseVc
                     (note that we ASSUME that the VCF is ordered by <contig,locus>).
                      Note that this will always leave at least one entry (the last one), since mostDownstreamLocusReached is in range of itself.
                    */
                    break;
                }
                // Already saw all variant positions within cacheWindow startDistance ahead of vc (on its contig)
            }
            // Update partiallyPhasedSites before it's used in phaseSite:
            oldPhasedList.addAll(discardIrrelevantPhasedSites());
            if (DEBUG) logger.debug("oldPhasedList(1st) = " + toStringVCL(oldPhasedList));

            final VariantAndReads vr = unphasedSiteQueue.remove();

            // should try to merge variants if they are SNPs that are within the minimum merging distance from each other
            final boolean shouldTryToMerge = enableMergePhasedSegregatingPolymorphismsToMNP &&
                    prevVr != null &&
                    prevVr.variant.isSNP() && vr.variant.isSNP() &&
                    vr.variant.getStart() - prevVr.variant.getStart() <= maxGenomicDistanceForMNP;

            // if should try to merge, find if there are any reads that contain both SNPs
            boolean commonReads = false;
            if ( shouldTryToMerge ) {
                for ( final String readName : vr.variantReadNames) {
                    if (prevVr.variantReadNames.contains(readName)) {
                        commonReads = true;
                        break;
                    }
                }
                if ( DEBUG && !commonReads )
                    logger.debug("No common reads with previous variant for " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant));
            }

            if (DEBUG)
                logger.debug("Performing phasing for " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant));

            // phase the variant site, cannot phase if trying to merge variants that do not have any common reads
            phaseSite(vr, phaseStats, !(shouldTryToMerge && !commonReads));

            // save previous variant reads for next iteration
            prevVr = vr;
        }

        // Update partiallyPhasedSites after phaseSite is done:
        oldPhasedList.addAll(discardIrrelevantPhasedSites());
        if (DEBUG) logger.debug("oldPhasedList(2nd) = " + toStringVCL(oldPhasedList));

        if (outputMultipleBaseCountsWriter != null)
            outputMultipleBaseCountsWriter.outputMultipleBaseCounts();

        return oldPhasedList;
    }

    // Flush out sites with (possibly) phased genotypes, if those sites are no longer needed to phase other downstream sites
    private List<VariantContext> discardIrrelevantPhasedSites() {
        List<VariantContext> vcList = new LinkedList<VariantContext>();

        GenomeLoc nextToPhaseLoc = null;
        if (!unphasedSiteQueue.isEmpty())
            nextToPhaseLoc = GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), unphasedSiteQueue.peek().variant);

        while (!partiallyPhasedSites.isEmpty()) {
            if (nextToPhaseLoc != null) { // otherwise, unphasedSiteQueue.isEmpty(), and therefore no need to keep any of the "past"
                UnfinishedVariantAndReads partPhasedVr = partiallyPhasedSites.peek();

                if (startDistancesAreInWindowRange(partPhasedVr.unfinishedVariant.getLocation(), nextToPhaseLoc))
                    // nextToPhaseLoc is still not far enough ahead of partPhasedVr to exclude partPhasedVr from calculations
                    break;
            }
            UnfinishedVariantAndReads uvr = partiallyPhasedSites.remove();
            vcList.add(uvr.unfinishedVariant.toVariantContext());
        }

        return vcList;
    }

    /* Phase vc (removed head of unphasedSiteQueue) using all VariantContext objects in
       partiallyPhasedSites, and all in unphasedSiteQueue that are within cacheWindow startDistance ahead of vc (on its contig).

       ASSUMES: All VariantContexts in unphasedSiteQueue are in positions downstream of vc (head of queue).
     */

    /**
     * Phase the variant site relative to the previous site
     *
     * @param vr            A variant and the reads for each sample at that site:
     * @param phaseStats    Summary statistics about phasing rates for each sample
     * @param canPhase      Can phase variant site relative to the previous site
     */
    private void phaseSite(final VariantAndReads vr, final PhasingStats phaseStats, final boolean canPhase) {
        VariantContext vc = vr.variant;
        logger.debug("Will phase vc = " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));

        UnfinishedVariantAndReads uvr = new UnfinishedVariantAndReads(vr);
        UnfinishedVariantContext uvc = uvr.unfinishedVariant;

        // Perform per-sample phasing:
        GenotypesContext sampGenotypes = vc.getGenotypes();
        Map<String, PhaseCounts> samplePhaseStats = new TreeMap<String, PhaseCounts>();
        for (final Genotype gt : sampGenotypes) {
            String samp = gt.getSampleName();

            if (DEBUG) logger.debug("sample = " + samp);
            if (isUnfilteredCalledDiploidGenotype(gt)) {
                if (gt.isHet()) { // Attempt to phase this het genotype relative to *SOME* previous het genotype:

                    // Create the list of all het genotypes preceding this one (and in the phasing window as contained in partiallyPhasedSites):
                    List<GenotypeAndReadBases> prevHetGenotypes = new LinkedList<GenotypeAndReadBases>();
                    CloneableIteratorLinkedList.CloneableIterator<UnfinishedVariantAndReads> phasedIt = partiallyPhasedSites.iterator();
                    while (phasedIt.hasNext()) {
                        UnfinishedVariantAndReads phasedVr = phasedIt.next();
                        Genotype prevGt = phasedVr.unfinishedVariant.getGenotype(samp);
                        if (prevGt != null && isUnfilteredCalledDiploidGenotype(prevGt) && prevGt.isHet()) {
                            GenotypeAndReadBases grb = new GenotypeAndReadBases(prevGt, phasedVr.sampleReadBases.get(samp), phasedVr.unfinishedVariant.getLocation());
                            prevHetGenotypes.add(grb);
                            if (DEBUG) logger.debug("Using UPSTREAM het site = " + grb.loc);
                        }
                    }

                    SNPallelePair allelePair = new SNPallelePair(gt);
                    if (DEBUG) logger.debug("Want to phase TOP vs. BOTTOM for: " + "\n" + allelePair);

                    boolean phasedCurGenotypeRelativeToPrevious = false;
                    if ( canPhase ) {
                        for (int goBackFromEndOfPrevHets = 0; goBackFromEndOfPrevHets < prevHetGenotypes.size(); goBackFromEndOfPrevHets++) {
                            PhasingWindow phaseWindow = new PhasingWindow(vr, samp, prevHetGenotypes, goBackFromEndOfPrevHets);

                            PhaseResult pr = phaseSampleAtSite(phaseWindow);
                            phasedCurGenotypeRelativeToPrevious = passesPhasingThreshold(pr.phaseQuality);

                            if (pr.phasingContainsInconsistencies) {
                                if (DEBUG)
                                    logger.debug("MORE than " + (MAX_FRACTION_OF_INCONSISTENT_READS * 100) + "% of the reads are inconsistent for phasing of " + GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));
                                uvc.setPhasingInconsistent();
                            }

                            if (phasedCurGenotypeRelativeToPrevious) {
                                Genotype prevHetGenotype = phaseWindow.phaseRelativeToGenotype();
                                SNPallelePair prevAllelePair = new SNPallelePair(prevHetGenotype);
                                if (!prevHetGenotype.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY))
                                    throw new ReviewedGATKException("Internal error: missing haplotype markings for previous genotype, even though we put it there...");
                                String[] prevPairNames = (String[]) prevHetGenotype.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY);

                                String[] curPairNames = ensurePhasing(allelePair, prevAllelePair, prevPairNames, pr.haplotype);
                                Genotype phasedGt = new GenotypeBuilder(gt)
                                        .alleles(allelePair.getAllelesAsList())
                                        .attribute(VCFConstants.PHASE_QUALITY_KEY, pr.phaseQuality)
                                        .attribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY, curPairNames)
                                        .make();
                                uvc.setGenotype(samp, phasedGt);

                                if (DEBUG) {
                                    logger.debug("PREVIOUS CHROMOSOME NAMES: Top= " + prevPairNames[0] + ", Bot= " + prevPairNames[1]);
                                    logger.debug("PREVIOUS CHROMOSOMES:\n" + prevAllelePair + "\n");

                                    logger.debug("CURRENT CHROMOSOME NAMES: Top= " + curPairNames[0] + ", Bot= " + curPairNames[1]);
                                    logger.debug("CURRENT CHROMOSOMES:\n" + allelePair + "\n");
                                    logger.debug("\n");
                                }
                            }

                            if (statsWriter != null) {
                                GenomeLoc prevLoc = null;
                                int curIndex = 0;
                                for (GenotypeAndReadBases grb : prevHetGenotypes) {
                                    if (curIndex == prevHetGenotypes.size() - 1 - goBackFromEndOfPrevHets) {
                                        prevLoc = grb.loc;
                                        break;
                                    }
                                    ++curIndex;
                                }
                                statsWriter.addStat(samp, GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc), startDistance(prevLoc, vc), pr.phaseQuality, phaseWindow.readsAtHetSites.size(), phaseWindow.hetGenotypes.length);
                            }

                            PhaseCounts sampPhaseCounts = samplePhaseStats.get(samp);
                            if (sampPhaseCounts == null) {
                                sampPhaseCounts = new PhaseCounts();
                                samplePhaseStats.put(samp, sampPhaseCounts);
                            }
                            sampPhaseCounts.numTestedSites++;

                            if (pr.phasingContainsInconsistencies) {
                                if (phasedCurGenotypeRelativeToPrevious)
                                    sampPhaseCounts.numInconsistentSitesPhased++;
                                else
                                    sampPhaseCounts.numInconsistentSitesNotPhased++;
                            }

                            if (phasedCurGenotypeRelativeToPrevious)
                                sampPhaseCounts.numPhased++;

                            // Phased current relative to *SOME* previous het genotype, so break out of loop:
                            if (phasedCurGenotypeRelativeToPrevious)
                                break;
                        }
                    }

                    if (!phasedCurGenotypeRelativeToPrevious) { // Either no previous hets, or unable to phase relative to any previous het:
                        String locStr = Integer.toString(GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc).getStart());

                        Genotype startNewHaplotypeGt = new GenotypeBuilder(gt)
                                .attribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY, new String[]{locStr + "-1", locStr + "-2"})
                                .make();

                        uvc.setGenotype(samp, startNewHaplotypeGt);
                    }
                }
            }
        }

        partiallyPhasedSites.add(uvr); // only add it in now, since don't want it to be there during phasing
        phaseStats.addIn(new PhasingStats(samplePhaseStats));
    }

    public boolean passesPhasingThreshold(double PQ) {
        return PQ >= phaseQualityThresh;
    }

    // A genotype and the base pileup that supports it
    private static class GenotypeAndReadBases {
        public Genotype genotype;
        public ReadBasesAtPosition readBases;
        public GenomeLoc loc;

        public GenotypeAndReadBases(Genotype genotype, ReadBasesAtPosition readBases, GenomeLoc loc) {
            this.genotype = genotype;
            this.readBases = readBases;
            this.loc = loc;
        }
    }

    // Object to represent the local window of het genotypes for which haplotypes are being scored and ranked
    private class PhasingWindow {
        private Genotype[] hetGenotypes = null;

        private int phaseRelativeToIndex = -1;
        private int phasingSiteIndex = -1;

        private Map<String, PhasingRead> readsAtHetSites = null;

        public Genotype phaseRelativeToGenotype() {
            return hetGenotypes[phaseRelativeToIndex];
        }

        // ASSUMES that: isUnfilteredCalledDiploidGenotype(vrGt) && vrGt.isHet() [vrGt = vr.variant.getGenotype(sample)]

        public PhasingWindow(VariantAndReads vr, String sample, List<GenotypeAndReadBases> prevHetGenotypes, int goBackFromEndOfPrevHets) {
            if (prevHetGenotypes.isEmpty() || goBackFromEndOfPrevHets >= prevHetGenotypes.size()) // no previous sites against which to phase
                throw new ReviewedGATKException("Should never get empty set of previous sites to phase against");

            // Include these previously phased sites in the phasing computation:
            List<GenotypeAndReadBases> listHetGenotypes = new LinkedList<GenotypeAndReadBases>(prevHetGenotypes);

            phaseRelativeToIndex = listHetGenotypes.size() - 1 - goBackFromEndOfPrevHets;
            phasingSiteIndex = listHetGenotypes.size();

            // Add the (het) position to be phased [at phasingSiteIndex]:
            GenomeLoc phaseLocus = GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant);
            GenotypeAndReadBases grbPhase = new GenotypeAndReadBases(vr.variant.getGenotype(sample), vr.sampleReadBases.get(sample), phaseLocus);
            listHetGenotypes.add(grbPhase);
            if (DEBUG) logger.debug("PHASING het site = " + grbPhase.loc + " [phasingSiteIndex = " + phasingSiteIndex + "]");

            // Include as-of-yet unphased sites in the phasing computation:
            for (VariantAndReads nextVr : unphasedSiteQueue) {
                if (!startDistancesAreInWindowRange(vr.variant, nextVr.variant)) //nextVr too far ahead of the range used for phasing vc
                    break;
                Genotype gt = nextVr.variant.getGenotype(sample);
                if (gt != null && isUnfilteredCalledDiploidGenotype(gt) && gt.isHet()) {
                    GenotypeAndReadBases grb = new GenotypeAndReadBases(gt, nextVr.sampleReadBases.get(sample), GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), nextVr.variant));
                    listHetGenotypes.add(grb);
                    if (DEBUG) logger.debug("Using DOWNSTREAM het site = " + grb.loc);
                }
            }

            // First, assemble the "sub-reads" from the COMPLETE WINDOW-BASED SET of heterozygous positions for this sample:
            buildReadsAtHetSites(listHetGenotypes, sample, grbPhase.loc);

            // Remove extraneous reads (those that do not "connect" the two core phasing sites):
            Set<String> onlyKeepReads = removeExtraneousReads(listHetGenotypes.size());

            // Dynamically modify the window to only include sites which have a non-empty set of reads:
            listHetGenotypes = removeExtraneousSites(listHetGenotypes);

            // In any case, must still trim the window size to be "feasible"
            // [**NOTE**: May want to do this to try maximize the preservation of paths from phaseRelativeToIndex to phasingSiteIndex]:
            if (listHetGenotypes.size() > maxPhaseSites) {
                listHetGenotypes = trimWindow(listHetGenotypes, sample, phaseLocus);

                // Can now remove any extra reads (and then sites):
                buildReadsAtHetSites(listHetGenotypes, onlyKeepReads);
                onlyKeepReads = removeExtraneousReads(listHetGenotypes.size());
                listHetGenotypes = removeExtraneousSites(listHetGenotypes);
            }

            // Lastly, assemble the "sub-reads" from the FINAL SET of heterozygous positions for this sample:
            buildReadsAtHetSites(listHetGenotypes, onlyKeepReads);

            // Copy to a fixed-size array:
            if (DEBUG) logger.debug("FINAL phasing window of " + listHetGenotypes.size() + " sites:\n" + toStringGRL(listHetGenotypes));
            hetGenotypes = new Genotype[listHetGenotypes.size()];
            int index = 0;
            for (GenotypeAndReadBases copyGrb : listHetGenotypes)
                hetGenotypes[index++] = copyGrb.genotype;
        }

        // Build the read sub-sequences at the het genomic positions:
        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, String sample, GenomeLoc phasingLoc) {
            buildReadsAtHetSites(listHetGenotypes, sample, phasingLoc, null);
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, Set<String> onlyKeepReads) {
            buildReadsAtHetSites(listHetGenotypes, null, null, onlyKeepReads);
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, String sample, GenomeLoc phasingLoc, Set<String> onlyKeepReads) {
            readsAtHetSites = new HashMap<String, PhasingRead>();

            int index = 0;
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                ReadBasesAtPosition readBases = grb.readBases;
                if (readBases != null) {
                    for (ReadBase rb : readBases) {
                        String readName = rb.readName;
                        if (onlyKeepReads != null && !onlyKeepReads.contains(readName)) // if onlyKeepReads exists, ignore reads not in onlyKeepReads
                            continue;

                        PhasingRead rd = readsAtHetSites.get(readName);
                        if (rd == null) {
                            rd = new PhasingRead(listHetGenotypes.size(), rb.mappingQual);
                            readsAtHetSites.put(readName, rd);
                        }
                        else if (outputMultipleBaseCountsWriter != null && rd.getBase(index) != null // rd already has a base at index
                                && sample != null && phasingLoc != null) {
                            outputMultipleBaseCountsWriter.setMultipleBases(new SampleReadLocus(sample, readName, grb.loc), phasingLoc, rd.getBase(index), rb.base);
                        }

                        // Arbitrarily updates to the last base observed for this sample and read (rb.base):
                        rd.updateBaseAndQuality(index, rb.base, rb.baseQual);
                    }
                }
                index++;
            }
            if (DEBUG) logger.debug("Number of sites in window = " + index);

            if (DEBUG && logger.isDebugEnabled()) {
                logger.debug("ALL READS [phasingSiteIndex = " + phasingSiteIndex + "]:");
                for (Map.Entry<String, PhasingRead> nameToReads : readsAtHetSites.entrySet()) {
                    String rdName = nameToReads.getKey();
                    PhasingRead rd = nameToReads.getValue();
                    logger.debug(rd + "\t" + rdName);
                }
            }
        }

        // Object to represent a pair of genomic sites, and all reads overlapping those 2 sites (though possibly others)
        private class EdgeToReads {
            private TreeMap<PhasingGraphEdge, List<String>> edgeReads;

            public EdgeToReads() {
                this.edgeReads = new TreeMap<PhasingGraphEdge, List<String>>(); // implemented GraphEdge.compareTo()
            }

            public void addRead(PhasingGraphEdge e, String readName) {
                List<String> reads = edgeReads.get(e);
                if (reads == null) {
                    reads = new LinkedList<String>();
                    edgeReads.put(e, reads);
                }
                reads.add(readName);
            }

            public List<String> getReads(PhasingGraphEdge e) {
                return edgeReads.get(e);
            }
        }

        private class IntegerSet implements Iterable<Integer> {
            private Set<Integer> list;

            public IntegerSet(Set<Integer> list) {
                this.list = list;
            }

            public boolean contains(int i) {
                return list.contains(i);
            }

            public Iterator<Integer> iterator() {
                return list.iterator();
            }

            public String toString() {
                StringBuilder sb = new StringBuilder();
                for (int i : this) {
                    sb.append(i + ", ");
                }
                return sb.toString();
            }
        }

        // Remove any reads that add no "connections" (PhasingGraphEdge) between pairs of het sites:
        public Set<String> removeExtraneousReads(int numHetSites) {
            PhasingGraph readGraph = new PhasingGraph(numHetSites);
            EdgeToReads edgeToReads = new EdgeToReads();
            Set<Integer> sitesWithEdges = new TreeSet<Integer>();

            for (Map.Entry<String, PhasingRead> nameToReads : readsAtHetSites.entrySet()) {
                String rdName = nameToReads.getKey();
                PhasingRead rd = nameToReads.getValue();

                int[] siteInds = rd.getNonNullIndices();
                // Connect each pair of non-null sites in rd:
                for (int i = 0; i < siteInds.length; i++) {
                    for (int j = i + 1; j < siteInds.length; j++) {
                        PhasingGraphEdge e = new PhasingGraphEdge(siteInds[i], siteInds[j]);
                        if (DEBUG) logger.debug("Read = " + rdName + " is adding edge: " + e);
                        readGraph.addEdge(e);

                        edgeToReads.addRead(e, rdName);

                        sitesWithEdges.add(e.getV1());
                        sitesWithEdges.add(e.getV2());
                    }
                }
            }
            if (DEBUG) logger.debug("Read graph:\n" + readGraph);
            Set<String> keepReads = new HashSet<String>();

            /* Check which Reads are involved in acyclic paths from phaseRelativeToIndex to (phasingSiteIndex):

               In detail:
               Every Read links EACH pair of sites for which it contains bases.  Then, each such edge is added to a "site connectivity graph".
               A read provides non-trivial bias toward the final haplotype decision if it participates in a path from prev ---> cur.  This is tested by
               considering each edge that the read contributes.  For edge e=(v1,v2), if there exists a path from prev ---> v1 [that doesn't include v2] and
               cur ---> v2 [that doesn't include v1], then there is a path from prev ---> cur that uses e, hence making the read significant.
               By excluding each vertex's edges and then calculating connected components, we are able to make the determination, for example,
               if a path exists from prev ---> v1 that excludes v2.

               Furthermore, if the path DOES use other edges that exist solely due to the read, then that's fine, since adding in the read will give those edges as well.
               And, if the path uses edges from other reads, then keeping all other reads that contribute those edges
               [which will happen since those edges are also in paths from prev ---> cur] is sufficient for this path to exist.

               NOTE:
               If we would use NON-UNIFORM priors for the various haplotypes consistent with a margnialized haplotype, then this calculation would not be correct, since the equivalence of:
               1. The read affects the final marginal haplotype posterior probability (for general mapping and base quality values).
               2. The read has edges involved in a path from prev ---> cur.
               DEPENDS STRONGLY on the fact that all haplotypes have the same EXACT prior.

               This is due to the following:
               [We denote:
               R = set of all reads
               r = a single read
               "AA + CC" = AA on top chromosome, CC on bottom chromosome]

               Note that since there are only two haplotype possibilities:
               P(AA + CC | R) + P(AC + CA | R) = 1

               Now, if we assume that all haplotypes consistent with AA + CC have the same prior probability [P(AA + CC | R)], then:
               P(AA + CC | R)
               = P(AAAA + CCCC | R) + ... + P(AACC + CCAA | R)
               = [P(AAAA + CCCC , R) + ... + P(AACC + CCAA , R)] / P(R)
               \propto P(AAAA + CCCC , R) + ... + P(AACC + CCAA , R)
               = P(R | AAAA + CCCC)*P(AAAA + CCCC) + ... + P(R | AACC + CCAA)*P(AACC + CCAA)
               = P(AA + CC | R) * [P(R | AAAA + CCCC) + ... + P(R | AACC + CCAA)]
               
               Since we assume independence between reads given a particular haplotype [P(R | AAAA + CCCC) = \prod_r P(r | AAAA + CCCC)],
               a new read r affects P(AA + CC | R) by multiplying each of the terms in the sum by, e.g., P(r | AAAA + CCCC).
               Therefore, if these values do not affect the ratio of:
               (I) [P(R | AAAA + CCCC) + ... + P(R | AACC + CCAA)] / [P(R | ACAA + CACC) + ... + P(R | ACCC + CAAA)]
               then they do not affect the value of:
               (II) P(AA + CC | R) / P(AC + CA | R)   [which uniquely defines their values, since they sum to 1]

               And, the P(r | AAAA + CCCC), ..., P(r | ACCC + CAAA) do not affect ratio (I) iff r's edges do not take part in a path from prev to cur in combination with the other reads in R.
             */
            int prev = phaseRelativeToIndex;
            int cur = phasingSiteIndex;

            if (!readGraph.getConnectedComponents().inSameSet(prev, cur)) { // There is NO path between cur and prev
                if (DEBUG)
                    logger.debug("NO READ PATH between PHASE site [" + cur + "] and UPSTREAM site [" + prev + "]");
                readsAtHetSites.clear();
                return keepReads;
            }

            /* Check the connected components of prev and cur when removing each individual vertex's edges:
               [Total run-time: for each vertex, calculate connected components after removing it's edges: O(V * E)]
             */
            IntegerSet[] removedSiteSameCCAsPrev = new IntegerSet[numHetSites];
            IntegerSet[] removedSiteSameCCAsCur = new IntegerSet[numHetSites];
            for (int i : sitesWithEdges) {
                if (DEBUG) logger.debug("Calculating CC after removing edges of site: " + i);

                // Remove all edges incident to i and see which positions have paths to prev and cur:
                Collection<PhasingGraphEdge> removedEdges = readGraph.removeAllIncidentEdges(i);

                // Run-time for efficiently calculating connected components using DisjointSet: O(E)
                DisjointSet ccAfterRemove = readGraph.getConnectedComponents();
                removedSiteSameCCAsPrev[i] = new IntegerSet(ccAfterRemove.inSameSetAs(prev, sitesWithEdges));
                removedSiteSameCCAsCur[i] = new IntegerSet(ccAfterRemove.inSameSetAs(cur, sitesWithEdges));

                if (DEBUG) logger.debug("Same CC as previous [" + prev + "]: " + removedSiteSameCCAsPrev[i]);
                if (DEBUG) logger.debug("Same CC as current  [" + cur + "]: " + removedSiteSameCCAsCur[i]);

                // Add the removed edges back in:
                readGraph.addEdges(removedEdges);
            }

            for (PhasingGraphEdge e : readGraph) {
                if (DEBUG) logger.debug("Testing the path-connectivity of Edge: " + e);

                /* Edge e={v1,v2} contributes a path between prev and cur for testRead iff:
                   testRead[v1] != null, testRead[v2] != null, and there is a path from prev ---> v1 -> v2 ---> cur  [or vice versa].
                   Note that the path from prev ---> v1 will NOT contain v2, since we removed all of v2's edges,
                   and the path from v2 ---> cur will NOT contain v1.
                 */
                boolean prevTo2and1ToCur = removedSiteSameCCAsPrev[e.getV1()].contains(e.getV2()) && removedSiteSameCCAsCur[e.getV2()].contains(e.getV1());
                boolean prevTo1and2ToCur = removedSiteSameCCAsPrev[e.getV2()].contains(e.getV1()) && removedSiteSameCCAsCur[e.getV1()].contains(e.getV2());

                if (prevTo2and1ToCur || prevTo1and2ToCur) {
                    for (String readName : edgeToReads.getReads(e)) {
                        keepReads.add(readName);

                        if (DEBUG && logger.isDebugEnabled()) {
                            if (prevTo2and1ToCur)
                                logger.debug("Keep read " + readName + " due to path: " + prev + " ---> " + e.getV2() + " -> " + e.getV1() + " ---> " + cur);
                            else
                                logger.debug("Keep read " + readName + " due to path: " + prev + " ---> " + e.getV1() + " -> " + e.getV2() + " ---> " + cur);
                        }
                    }
                }
            }

            // Retain only the reads that contain an edge in a path connecting prev and cur:
            Iterator<Map.Entry<String, PhasingRead>> readIt = readsAtHetSites.entrySet().iterator();
            while (readIt.hasNext()) {
                Map.Entry<String, PhasingRead> nameToReads = readIt.next();
                String rdName = nameToReads.getKey();
                if (!keepReads.contains(rdName)) {
                    readIt.remove();
                    if (DEBUG) logger.debug("Removing extraneous read: " + rdName);
                }
            }

            return keepReads;
        }

        // Remove all het sites that have no reads (which may occur if all of the reads supporting the original call don't contain an additional het site and were thus removed above):
        private List<GenotypeAndReadBases> removeExtraneousSites(List<GenotypeAndReadBases> listHetGenotypes) {
            Set<Integer> sitesWithReads = new HashSet<Integer>();
            for (Map.Entry<String, PhasingRead> nameToReads : readsAtHetSites.entrySet()) {
                PhasingRead rd = nameToReads.getValue();
                for (int i : rd.getNonNullIndices())
                    sitesWithReads.add(i);
            }

            // Remove all sites that have no read bases:
            List<GenotypeAndReadBases> keepHetSites = new LinkedList<GenotypeAndReadBases>();
            int index = 0;
            int numPrecedingPhaseRelativeToSiteRemoved = 0;
            int numPrecedingPhasingSiteRemoved = 0;
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                boolean keepSite = sitesWithReads.contains(index);
                if (DEBUG && logger.isDebugEnabled() && !keepSite)
                    logger.debug("Removing read-less site " + grb.loc);

                if (keepSite || index == phasingSiteIndex || index == phaseRelativeToIndex) {
                    keepHetSites.add(grb);
                    if (!keepSite)
                        if (DEBUG)
                            logger.debug("Although current or previous sites have no relevant reads, continuing empty attempt to phase them [for sake of program flow]...");
                }
                else {
                    if (index <= phaseRelativeToIndex)
                        numPrecedingPhaseRelativeToSiteRemoved++;
                    if (index <= phasingSiteIndex)
                        numPrecedingPhasingSiteRemoved++;
                }

                index++;
            }

            phaseRelativeToIndex -= numPrecedingPhaseRelativeToSiteRemoved;
            phasingSiteIndex -= numPrecedingPhasingSiteRemoved;
            return keepHetSites;
        }

        /* Auxilary object to sort candidate het sites with which to phase the index site,
           where sorting is performed based on distance to the index site
           (since presumably closer sites will have greater numbers of overlapping reads)
         */
        private class SortSitesBySumOfDist implements Comparator<Integer> {
            private Vector<GenotypeAndReadBases> grb;

            public SortSitesBySumOfDist(List<GenotypeAndReadBases> listHetGenotypes) {
                grb = new Vector<GenotypeAndReadBases>(listHetGenotypes);
            }

            public int compare(Integer i1, Integer i2) {
                int d1 = calcGenomicDist(i1);
                int d2 = calcGenomicDist(i2);

                if (d1 != d2)
                    return d1 - d2;

                int id1 = calcIndexDist(i1);
                int id2 = calcIndexDist(i2);
                if (id1 != id2)
                    return id1 - id2;

                return i1 - i2;
            }

            private int calcGenomicDist(int i) {
                int d1 = grb.get(i).loc.distance(grb.get(phaseRelativeToIndex).loc);
                int d2 = grb.get(i).loc.distance(grb.get(phasingSiteIndex).loc);

                return d1 + d2;
            }

            private int calcIndexDist(int i) {
                int d1 = Math.abs(i - phaseRelativeToIndex);
                int d2 = Math.abs(i - phasingSiteIndex);

                return d1 + d2;
            }
        }

        // Create a "phasing window" of het sites to use for phasing the index site, but limiting to only maxPhaseSites het sites to incorporate [as specified by the user]
        private List<GenotypeAndReadBases> trimWindow(List<GenotypeAndReadBases> listHetGenotypes, String sample, GenomeLoc phaseLocus) {
            if (DEBUG)
                logger.warn("Trying to phase sample " + sample + " at locus " + phaseLocus + " within a window of " + cacheWindow + " bases yields " + listHetGenotypes.size() + " heterozygous sites to phase:\n" + toStringGRL(listHetGenotypes));

            Set<Integer> scoreAllIndices = new TreeSet<Integer>(new SortSitesBySumOfDist(listHetGenotypes));
            for (int i = 0; i < listHetGenotypes.size(); ++i) {
                if (i != phaseRelativeToIndex && i != phasingSiteIndex)
                    scoreAllIndices.add(i);
            }

            Set<Integer> keepIndices = new TreeSet<Integer>();
            // always keep these two indices:
            keepIndices.add(phaseRelativeToIndex);
            keepIndices.add(phasingSiteIndex);
            for (int addInd : scoreAllIndices) {
                if (keepIndices.size() >= maxPhaseSites)
                    break;
                else // keepIndices.size() < maxPhaseSites
                    keepIndices.add(addInd);
            }

            List<GenotypeAndReadBases> newListHetGenotypes = new LinkedList<GenotypeAndReadBases>();
            int newPhaseRelativeToIndex = -1;
            int newPhasingSiteIndex = -1;
            int oldIndex = 0;
            int newIndex = 0;
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                if (keepIndices.contains(oldIndex)) {
                    newListHetGenotypes.add(grb);

                    if (oldIndex == phaseRelativeToIndex)
                        newPhaseRelativeToIndex = newIndex;
                    if (oldIndex == phasingSiteIndex)
                        newPhasingSiteIndex = newIndex;

                    ++newIndex;
                }
                ++oldIndex;
            }

            phaseRelativeToIndex = newPhaseRelativeToIndex;
            phasingSiteIndex = newPhasingSiteIndex;
            listHetGenotypes = newListHetGenotypes;
            if (DEBUG)
                logger.warn("NAIVELY REDUCED to " + listHetGenotypes.size() + " sites:\n" + toStringGRL(listHetGenotypes));

            return listHetGenotypes;
        }
    }

    // Phase a particular sample's het genotype using a constructed PhasingWindow:
    private PhaseResult phaseSampleAtSite(PhasingWindow phaseWindow) {
        /* Will map a phase and its "complement" to a single representative phase,
          and marginalizeAsNewTable() marginalizes to 2 positions [starting at the previous position, and then the current position]:
        */
        int[] marginalizeInds = {phaseWindow.phaseRelativeToIndex, phaseWindow.phasingSiteIndex};
        HaplotypeTableCreator tabCreator = new TableCreatorOfHaplotypeAndComplementForDiploidAlleles(phaseWindow.hetGenotypes, marginalizeInds);
        PhasingTable sampleHaps = tabCreator.getNewTable();

        if (DEBUG && logger.isDebugEnabled()) {
            logger.debug("Number of USED reads [connecting the two positions to be phased] at sites: " + phaseWindow.readsAtHetSites.size());
            logger.debug("USED READS:");
            for (Map.Entry<String, PhasingRead> nameToReads : phaseWindow.readsAtHetSites.entrySet()) {
                String rdName = nameToReads.getKey();
                PhasingRead rd = nameToReads.getValue();
                logger.debug(rd + "\t" + rdName);
            }
        }

        // Update the phasing table based on each of the sub-reads for this sample:
        MaxHaplotypeAndQuality prevMaxHapAndQual = null;

        int numHighQualityIterations = 0;
        int numInconsistentIterations = 0;

        double totalAbsPQchange = 0;
        int numPQchangesObserved = 0;

        for (Map.Entry<String, PhasingRead> nameToReads : phaseWindow.readsAtHetSites.entrySet()) {
            PhasingRead rd = nameToReads.getValue();
            if (DEBUG) logger.debug("\nrd = " + rd + "\tname = " + nameToReads.getKey());

            for (PhasingTable.PhasingTableEntry pte : sampleHaps) {
                PhasingScore score = rd.matchHaplotypeClassScore(pte.getHaplotypeClass());
                pte.getScore().integrateReadScore(score);
                if (DEBUG) logger.debug("score(" + rd + ", " + pte.getHaplotypeClass() + ") = " + score);
            }

            // Check the current best haplotype assignment and compare it to the previous one:
            MaxHaplotypeAndQuality curMaxHapAndQual = new MaxHaplotypeAndQuality(sampleHaps, false);
            if (DEBUG)
                logger.debug("CUR MAX hap:\t" + curMaxHapAndQual.maxEntry.getHaplotypeClass() + "\tcurPhaseQuality:\t" + curMaxHapAndQual.phaseQuality);
            if (prevMaxHapAndQual != null) {
                double changeInPQ = prevMaxHapAndQual.phaseQuality - curMaxHapAndQual.phaseQuality;

                if (passesPhasingThreshold(prevMaxHapAndQual.phaseQuality)) {
                    numHighQualityIterations++;
                    if (!curMaxHapAndQual.hasSameRepresentativeHaplotype(prevMaxHapAndQual) || // switched phase
                            (numPQchangesObserved > 0 && changeInPQ > FRACTION_OF_MEAN_PQ_CHANGES * (totalAbsPQchange / numPQchangesObserved))) { // a "significant" decrease in PQ
                        if (DEBUG) logger.debug("Inconsistent read found!");
                        numInconsistentIterations++;
                    }
                }

                totalAbsPQchange += Math.abs(changeInPQ);
                numPQchangesObserved++;
            }
            prevMaxHapAndQual = curMaxHapAndQual;
        }

        if (DEBUG) logger.debug("\nPhasing table [AFTER CALCULATION]:\n" + sampleHaps + "\n");
        MaxHaplotypeAndQuality maxHapQual = new MaxHaplotypeAndQuality(sampleHaps, DEBUG);
        double posteriorProb = maxHapQual.maxEntry.getScore().getValue();

        if (DEBUG)
            logger.debug("MAX hap:\t" + maxHapQual.maxEntry.getHaplotypeClass() + "\tposteriorProb:\t" + posteriorProb + "\tphaseQuality:\t" + maxHapQual.phaseQuality);
        if (DEBUG)
            logger.debug("Number of used reads " + phaseWindow.readsAtHetSites.size() + "; number of high PQ iterations " + numHighQualityIterations + "; number of inconsistencies " + numInconsistentIterations);

        boolean phasingContainsInconsistencies = false;
        if (numInconsistentIterations / (double) numHighQualityIterations > MAX_FRACTION_OF_INCONSISTENT_READS)
            phasingContainsInconsistencies = true;

        return new PhaseResult(maxHapQual.getRepresentative(), maxHapQual.phaseQuality, phasingContainsInconsistencies);
    }

    // Object represents the maximum-scoring haplotype and its corresponding quality score
    private static class MaxHaplotypeAndQuality {
        public PhasingTable.PhasingTableEntry maxEntry;
        public double phaseQuality;

        public MaxHaplotypeAndQuality(PhasingTable hapTable, boolean printDebug) {
            // Marginalize each haplotype to its first 2 positions:
            hapTable = HaplotypeTableCreator.marginalizeAsNewTable(hapTable);
            if (printDebug)
                logger.debug("\nPhasing table [AFTER MAPPING]:\n" + hapTable + "\n");

            calculateMaxHapAndPhasingQuality(hapTable, printDebug);
        }

        // Calculates maxEntry and its PQ (within table hapTable):
        private void calculateMaxHapAndPhasingQuality(PhasingTable hapTable, boolean printDebug) {
            hapTable.normalizeScores();
            if (printDebug)
                logger.debug("\nPhasing table [AFTER NORMALIZATION]:\n" + hapTable + "\n");

            // Determine the phase at this position:
            this.maxEntry = hapTable.maxEntry();

            // convert posteriorProb to PHRED scale, but do NOT cap the quality as in QualityUtils.trueProbToQual(posteriorProb):
            PreciseNonNegativeDouble sumErrorProbs = new PreciseNonNegativeDouble(ZERO);
            for (PhasingTable.PhasingTableEntry pte : hapTable) {
                if (pte != maxEntry)
                    sumErrorProbs.plusEqual(pte.getScore());
            }
            this.phaseQuality = -10.0 * (sumErrorProbs.getLog10Value());
        }

        // Comparator that compares if 2 haplotypes map back to the same "representative" haplotype (accounts for reverse complementarity)
        public boolean hasSameRepresentativeHaplotype(MaxHaplotypeAndQuality that) {
            return this.getRepresentative().equals(that.getRepresentative());
        }

        private Haplotype getRepresentative() {
            return maxEntry.getHaplotypeClass().getRepresentative();
        }
    }

    /*
        Ensure that curAllelePair is phased relative to prevAllelePair as specified by hap.
     */

    public static String[] ensurePhasing(SNPallelePair curAllelePair, SNPallelePair prevAllelePair, String[] prevPairNames, Haplotype hap) {
        if (hap.size() < 2)
            throw new ReviewedGATKException("LOGICAL ERROR: Only considering haplotypes of length > 2!");

        String[] curPairNames = prevPairNames;

        byte prevBase = hap.getBase(0); // The 1st base in the haplotype
        byte curBase = hap.getBase(1);  // The 2nd base in the haplotype

        boolean chosePrevTopChrom = prevAllelePair.matchesTopBase(prevBase);
        boolean choseCurTopChrom = curAllelePair.matchesTopBase(curBase);
        if (chosePrevTopChrom != choseCurTopChrom) {
            //curAllelePair.swapAlleles();

            /* Instead of swapping the alleles (as we used to above),
               we swap the haplotype names to fit the unswapped alleles as they are ordered in the Genotype:
            */
            curPairNames = new String[]{prevPairNames[1], prevPairNames[0]};
        }

        return curPairNames;
    }

    private boolean startDistancesAreInWindowRange(VariantContext vc1, VariantContext vc2) {
        return startDistancesAreInWindowRange(GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc1), GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc2));
    }

    private boolean startDistancesAreInWindowRange(GenomeLoc loc1, GenomeLoc loc2) {
        return loc1.distance(loc2) <= cacheWindow; // distance() checks: loc1.onSameContig(loc2)
    }

    private int startDistance(GenomeLoc gl1, VariantContext vc2) {
        return gl1.distance(GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc2));
    }

    public PhasingStats reduce(PhasingStatsAndOutput statsAndList, PhasingStats stats) {
        if (statsAndList != null) {
            writeVcList(statsAndList.output);
            stats.addIn(statsAndList.ps);
        }
        return stats;
    }

    /**
     * Phase anything left in the cached unphasedSiteQueue, and report the number of reads and VariantContexts processed.
     *
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(PhasingStats result) {
        List<VariantContext> finalList = processQueue(result, true); // process all remaining data
        writeVcList(finalList);
        writer.close();

        if (statsWriter != null)
            statsWriter.close();

        if (outputMultipleBaseCountsWriter != null)
            outputMultipleBaseCountsWriter.close();

        System.out.println("Coverage over ALL samples:");
        System.out.println("Number of reads observed: " + result.getNumReads());
        System.out.println("Number of variant sites observed: " + result.getNumVarSites());
        System.out.println("Average coverage: " + ((double) result.getNumReads() / result.getNumVarSites()));

        System.out.println("\n--- Phasing summary [minimal haplotype quality (PQ): " + phaseQualityThresh + ", maxPhaseSites: " + maxPhaseSites + ", cacheWindow: " + cacheWindow + "] ---");
        for (Map.Entry<String, PhaseCounts> sampPhaseCountEntry : result.getPhaseCounts()) {
            PhaseCounts pc = sampPhaseCountEntry.getValue();
            System.out.print("Sample: " + sampPhaseCountEntry.getKey() + "\tSites tested: " + pc.numTestedSites + "\tSites phased: " + pc.numPhased);
            System.out.println("\tPhase-inconsistent sites: " + (pc.numInconsistentSitesPhased + pc.numInconsistentSitesNotPhased) + " [phased: " + pc.numInconsistentSitesPhased + ", unphased:" + pc.numInconsistentSitesNotPhased + "]");
        }
        System.out.println("");
    }

    private void writeVcList(List<VariantContext> varContList) {
        for (VariantContext vc : varContList)
            writeVCF(vc);
    }

    private void writeVCF(VariantContext vc) {
        if (samplesToPhase == null || vc.isNotFiltered())
            //if ( samplesToPhase == null || (vc.isVariant() && vc.isNotFiltered())) // if we are only operating on specific samples, don't write out all sites, just those where the VC is variant
            writer.add(vc);
    }

    public static boolean processVariantInPhasing(VariantContext vc) {
        return vc.isNotFiltered() && ((vc.isSNP() && vc.isBiallelic()) || !vc.isVariant()); // we can handle the non-variant case as well
        //return isUnfilteredBiallelicSNP(vc);
    }


    /*
      Inner classes:
    */

    // A variant and the reads for each sample at that site:
    private class VariantAndReads {
        public VariantContext variant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;
        public Set<String> variantReadNames;

        public VariantAndReads(VariantContext variant, HashMap<String, ReadBasesAtPosition> sampleReadBases) {
            this.variant = variant;
            this.sampleReadBases = sampleReadBases;
        }

        public VariantAndReads(VariantContext variant, AlignmentContext alignment) {
            this.variant = variant;
            this.sampleReadBases = new HashMap<String, ReadBasesAtPosition>();

            if (alignment != null) {
                ReadBackedPileup pileup = alignment.getBasePileup();
                if (pileup != null) {
                    // filter the read-base pileup based on min base and mapping qualities:
                    pileup = pileup.getBaseAndMappingFilteredPileup(MIN_BASE_QUALITY_SCORE, MIN_MAPPING_QUALITY_SCORE);
                    if (pileup != null) {
                        for (final String sample : pileup.getSamples()) {
                            ReadBackedPileup samplePileup = pileup.getPileupForSample(sample);
                            ReadBasesAtPosition readBases = new ReadBasesAtPosition();
                            for (PileupElement p : samplePileup) {
                                if (!p.isDeletion()) // IGNORE deletions for now
                                    readBases.putReadBase(p);
                            }
                            sampleReadBases.put(sample, readBases);
                        }
                    }

                    // if merging SNPs, save the read names overlapping the variant
                    if (enableMergePhasedSegregatingPolymorphismsToMNP && variant.isSNP()) {
                        variantReadNames = new HashSet<>();
                        for ( final GATKSAMRecord read : pileup.getReads() ) {
                            // get the SNP position in the read
                            Pair<Integer, Boolean> pair = ReadUtils.getReadCoordinateForReferenceCoordinate(read, variant.getStart());

                            // get the reads containing the SNP
                            for (final Allele altAllele : variant.getAlternateAlleles()) {
                                if (read.getReadBases()[pair.first] == altAllele.getBases()[0]) {
                                    variantReadNames.add(read.getReadName());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Object to represent a variant that has yet to be phased, along with its underlying base pileups:
    private class UnfinishedVariantAndReads {
        public UnfinishedVariantContext unfinishedVariant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;

        public UnfinishedVariantAndReads(VariantAndReads vr) {
            this.unfinishedVariant = new UnfinishedVariantContext(vr.variant);
            this.sampleReadBases = vr.sampleReadBases;
        }
    }

    // COULD replace with MutableVariantContext if it worked [didn't throw exceptions when trying to call its set() methods]...

    private class UnfinishedVariantContext implements HasGenomeLocation {
        private String name;
        private String contig;
        private int start;
        private int stop;
        private Collection<Allele> alleles;
        private Map<String,Genotype> genotypes;
        private double log10PError;
        private Set<String> filters;
        private Map<String, Object> attributes;
        private String id;

        public UnfinishedVariantContext(VariantContext vc) {
            this.name = vc.getSource();
            this.id = vc.getID();
            this.contig = vc.getChr();
            this.start = vc.getStart();
            this.stop = vc.getEnd();
            this.alleles = vc.getAlleles();

            this.genotypes = new HashMap<String, Genotype>();
            for ( final Genotype g : vc.getGenotypes() ) {
                this.genotypes.put(g.getSampleName(), g);
            }

            this.log10PError = vc.getLog10PError();
            this.filters = vc.filtersWereApplied() ? vc.getFilters() : null;
            this.attributes = new HashMap<String, Object>(vc.getAttributes());
        }

        public VariantContext toVariantContext() {
            GenotypesContext gc = GenotypesContext.copy(this.genotypes.values());
            return new VariantContextBuilder(name, contig, start, stop, alleles).id(id)
                    .genotypes(gc).log10PError(log10PError).filters(filters).attributes(attributes).make();
        }

        public GenomeLoc getLocation() {
            return getToolkit().getGenomeLocParser().createGenomeLoc(contig, start, stop);
        }

        public Genotype getGenotype(String sample) {
            return genotypes.get(sample);
        }

        public void setGenotype(String sample, Genotype newGt) {
            this.genotypes.put(sample, newGt);
        }

        public void setPhasingInconsistent() {
            attributes.put(GATKVCFConstants.RBP_INCONSISTENT_KEY, true);
        }
    }

    private static String toStringGRL(List<GenotypeAndReadBases> grbList) {
        boolean first = true;
        StringBuilder sb = new StringBuilder();
        for (GenotypeAndReadBases grb : grbList) {
            if (first)
                first = false;
            else
                sb.append(" -- ");

            sb.append(grb.loc);
        }
        return sb.toString();
    }

    private String toStringVCL(List<VariantContext> vcList) {
        boolean first = true;
        StringBuilder sb = new StringBuilder();
        for (VariantContext vc : vcList) {
            if (first)
                first = false;
            else
                sb.append(" -- ");

            sb.append(GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));
        }
        return sb.toString();
    }

//
// THIS IMPLEMENTATION WILL FAIL WHEN NOT DEALING WITH SNP Alleles (e.g., MNP or INDEL), SINCE THEN THE Allele.getBases()
// FUNCTION WILL RETURN VARIABLE-LENGTH Byte ARRAYS.  IN THAT CASE, BaseArray/Haplotype/Read WILL NEED TO BE REPLACED WITH
// AN ArrayList OF Allele [OR SIMILAR OBJECT], and WON'T USE: getSingleBase(alleleI)
//

    /* Creates table of all 2^n local haplotypes,
       where n is the number of heterozygous SNPs in the local region we expected to find phase-informative reads
     */
    private static abstract class HaplotypeTableCreator {
        protected Genotype[] genotypes;

        public HaplotypeTableCreator(Genotype[] hetGenotypes) {
            this.genotypes = hetGenotypes;
        }

        abstract public PhasingTable getNewTable();

        protected List<Haplotype> getAllHaplotypes() {
            int numSites = genotypes.length;
            int[] genotypeCards = new int[numSites];
            for (int i = 0; i < numSites; i++)
                genotypeCards[i] = genotypes[i].getPloidy();

            LinkedList<Haplotype> allHaps = new LinkedList<Haplotype>();
            CardinalityCounter alleleCounter = new CardinalityCounter(genotypeCards);
            for (int[] alleleInds : alleleCounter) {
                byte[] hapBases = new byte[numSites];
                for (int i = 0; i < numSites; i++) {
                    Allele alleleI = genotypes[i].getAllele(alleleInds[i]);
                    hapBases[i] = SNPallelePair.getSingleBase(alleleI);
                }
                allHaps.add(new Haplotype(hapBases));
            }
            return allHaps;
        }

        /* For phasing site X relative to site X-1, we sum the probabilities over all haplotypes of the phases of [X-1, X].
           That is, we aggregate probability mass over all haplotypes consistent with a particular phase at the [X-1, X] pair.
         */
        public static PhasingTable marginalizeAsNewTable(PhasingTable table) {
            TreeMap<Haplotype, PreciseNonNegativeDouble> hapMap = new TreeMap<Haplotype, PreciseNonNegativeDouble>();
            for (PhasingTable.PhasingTableEntry pte : table) {
                Haplotype rep = pte.getHaplotypeClass().getRepresentative();
                PreciseNonNegativeDouble score = hapMap.get(rep);
                if (score == null) {
                    score = new PreciseNonNegativeDouble(ZERO);
                    hapMap.put(rep, score);
                }
                score.plusEqual(pte.getScore());
            }

            PhasingTable margTable = new PhasingTable();
            for (Map.Entry<Haplotype, PreciseNonNegativeDouble> hapClassAndScore : hapMap.entrySet()) {
                Haplotype rep = hapClassAndScore.getKey();
                ArrayList<Haplotype> hapList = new ArrayList<Haplotype>();
                hapList.add(rep);

                HaplotypeClass hc = new HaplotypeClass(hapList, rep);
                margTable.addEntry(hc, hapClassAndScore.getValue());
            }
            return margTable;
        }
    }

    // Implementation for diploid alleles (thus assuming 2^n haplotypes):
    private static class TableCreatorOfHaplotypeAndComplementForDiploidAlleles extends HaplotypeTableCreator {
        private SNPallelePair[] SNPallelePairs;
        Set<Integer> marginalizeInds;

        public TableCreatorOfHaplotypeAndComplementForDiploidAlleles(Genotype[] hetGenotypes, int[] marginalizeInds) {
            super(hetGenotypes);

            this.SNPallelePairs = new SNPallelePair[genotypes.length];
            for (int i = 0; i < genotypes.length; i++)
                SNPallelePairs[i] = new SNPallelePair(genotypes[i]);

            this.marginalizeInds = new TreeSet<Integer>();
            for (int mind : marginalizeInds)
                this.marginalizeInds.add(mind);
        }

        public PhasingTable getNewTable() {
            int startIndex = marginalizeInds.iterator().next();

            PhasingTable table = new PhasingTable();
            for (Haplotype hap : getAllHaplotypes()) {
                if (SNPallelePairs[startIndex].matchesTopBase(hap.getBase(startIndex))) {
                    /* hap is the "representative" haplotype [DEFINED here to be
                      the one with the top base at the startIndex position.
                      NOTE that it is CRITICAL that this definition be consistent with the representative sub-haplotypes defined below!]
                    */
                    ArrayList<Haplotype> hapList = new ArrayList<Haplotype>();
                    hapList.add(hap);
                    hapList.add(complement(hap));

                    Haplotype rep = hap.subHaplotype(marginalizeInds);
                    double hapClassPrior = getHaplotypeRepresentativePrior(rep); // Note that prior is ONLY a function of the representative haplotype

                    HaplotypeClass hapClass = new HaplotypeClass(hapList, rep);
                    table.addEntry(hapClass, hapClassPrior);
                }
            }
            return table;
        }

        // Can change later to weight the representative Haplotypes differently:

        private double getHaplotypeRepresentativePrior(Haplotype rep) {
            return 1.0;
        }

        /* Since assuming biallelic genotypes, we use this to map a haplotype to the corresponding haplotype,
           where the other allele is chosen at each site
         */
        private Haplotype complement(Haplotype hap) {
            int numSites = SNPallelePairs.length;
            if (hap.size() != numSites)
                throw new ReviewedGATKException("INTERNAL ERROR: hap.size() != numSites");

            // Take the other base at EACH position of the Haplotype:
            byte[] complementBases = new byte[numSites];
            for (int i = 0; i < numSites; i++)
                complementBases[i] = SNPallelePairs[i].getOtherBase(hap.getBase(i));

            return new Haplotype(complementBases);
        }
    }

    // Table to represent the list of all haplotypes and their scores:
    private static class PhasingTable implements Iterable<PhasingTable.PhasingTableEntry> {
        private LinkedList<PhasingTableEntry> table;

        public PhasingTable() {
            this.table = new LinkedList<PhasingTableEntry>();
        }

        public PhasingTableEntry addEntry(HaplotypeClass haplotypeClass, PreciseNonNegativeDouble initialScore) {
            PhasingTableEntry pte = new PhasingTableEntry(haplotypeClass, new PhasingScore(initialScore));
            table.add(pte);
            return pte;
        }

        public PhasingTableEntry addEntry(HaplotypeClass haplotypeClass, double initialScore) {
            return addEntry(haplotypeClass, new PreciseNonNegativeDouble(initialScore));
        }

        public Iterator<PhasingTableEntry> iterator() {
            return table.iterator();
        }

        public boolean isEmpty() {
            return table.isEmpty();
        }

        public PhasingTableEntry maxEntry() {
            if (table.isEmpty())
                return null;

            PhasingTableEntry maxPte = null;
            for (PhasingTableEntry pte : table) {
                if (maxPte == null || pte.getScore().gt(maxPte.getScore())) {
                    maxPte = pte;
                }
            }
            return maxPte;
        }

        // Normalize all the scores of the phasing table by their sum total:
        public void normalizeScores() {
            PreciseNonNegativeDouble normalizeBy = new PreciseNonNegativeDouble(ZERO);
            for (PhasingTableEntry pte : table)
                normalizeBy.plusEqual(pte.getScore());

            if (!normalizeBy.equals(ZERO)) { // prevent precision problems
                for (PhasingTableEntry pte : table)
                    pte.getScore().divEqual(normalizeBy);
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("-------------------\n");
            for (PhasingTableEntry pte : this) {
                sb.append("Haplotypes:\t" + pte.getHaplotypeClass() + "\tScore:\t" + pte.getScore() + "\n");
            }
            sb.append("-------------------\n");
            return sb.toString();
        }

        // An entry in the phasing table for a particular set of equivalent haplotypes (e.g., a haplotype and its "complement" -- see above)
        public static class PhasingTableEntry implements Comparable<PhasingTableEntry> {
            private HaplotypeClass haplotypeClass;
            private PhasingScore score;

            public PhasingTableEntry(HaplotypeClass haplotypeClass, PhasingScore score) {
                this.haplotypeClass = haplotypeClass;
                this.score = score;
            }

            public HaplotypeClass getHaplotypeClass() {
                return haplotypeClass;
            }

            public PhasingScore getScore() {
                return score;
            }

            public int compareTo(PhasingTableEntry that) {
                return this.getScore().compareTo(that.getScore());
            }
        }
    }

    private static class PhaseResult {
        public Haplotype haplotype;
        public double phaseQuality;
        public boolean phasingContainsInconsistencies;

        public PhaseResult(Haplotype haplotype, double phaseQuality, boolean phasingContainsInconsistencies) {
            this.haplotype = haplotype;
            this.phaseQuality = phaseQuality;
            this.phasingContainsInconsistencies = phasingContainsInconsistencies;
        }
    }

    public static boolean isUnfilteredBiallelicSNP(VariantContext vc) {
        return (vc.isNotFiltered() && vc.isSNP() && vc.isBiallelic());
    }

    public static boolean isUnfilteredCalledDiploidGenotype(Genotype gt) {
        return (! gt.isFiltered() && gt.isCalled() && gt.getPloidy() == 2);
    }

    // Class to output verbose information on instances where a single read has multiple bases at the same position (e.g., from paired-end overlap with a base error):
    private class MultipleBaseCountsWriter {
        private BufferedWriter writer = null;
        private TreeMap<SampleReadLocus, MultipleBaseCounts> multipleBaseCounts = null;

        public MultipleBaseCountsWriter(File outputMultipleBaseCountsFile) {
            FileOutputStream output;
            try {
                output = new FileOutputStream(outputMultipleBaseCountsFile);
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Unable to create multiple base count file at location: " + outputMultipleBaseCountsFile);
            }
            this.writer = new BufferedWriter(new OutputStreamWriter(output));

            this.multipleBaseCounts = new TreeMap<SampleReadLocus, MultipleBaseCounts>(); // implemented SampleReadLocus.compareTo()
        }

        public void setMultipleBases(SampleReadLocus srl, GenomeLoc phasingLoc, byte prevBase, byte newBase) {
            MultipleBaseCounts mbc = multipleBaseCounts.get(srl);
            if (mbc == null) {
                mbc = new MultipleBaseCounts(phasingLoc);
                mbc.incrementBaseCount(prevBase); // only now, do we know to note this
                multipleBaseCounts.put(srl, mbc);
            }
            if (mbc.samePhasingLocAs(phasingLoc)) // otherwise, don't want to count these multiple base counts again
                mbc.incrementBaseCount(newBase);

        }

        public void outputMultipleBaseCounts() {
            GenomeLoc nextToPhaseLoc = null;
            if (!unphasedSiteQueue.isEmpty())
                nextToPhaseLoc = GATKVariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), unphasedSiteQueue.peek().variant);

            outputMultipleBaseCounts(nextToPhaseLoc);
        }

        private void outputMultipleBaseCounts(GenomeLoc nextToPhaseLoc) {
            try {
                Iterator<Map.Entry<SampleReadLocus, MultipleBaseCounts>> multBaseCountIt = multipleBaseCounts.entrySet().iterator();
                while (multBaseCountIt.hasNext()) {
                    Map.Entry<SampleReadLocus, MultipleBaseCounts> sampleReadLocBaseCountsEntry = multBaseCountIt.next();
                    SampleReadLocus srl = sampleReadLocBaseCountsEntry.getKey();
                    if (nextToPhaseLoc == null || !startDistancesAreInWindowRange(srl.getLocus(), nextToPhaseLoc)) {
                        // Done with entry, so print it and remove it from map:
                        writer.write(srl + "\t" + sampleReadLocBaseCountsEntry.getValue() + "\n");
                        multBaseCountIt.remove();
                    }
                }
                writer.flush();
            } catch (IOException e) {
                throw new RuntimeException("Unable to write to outputMultipleBaseCountsFile", e);
            }
        }

        public void close() {
            outputMultipleBaseCounts(null);

            try {
                writer.flush();
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Unable to close outputMultipleBaseCountsFile");
            }
        }
    }
}


class PhasingScore extends PreciseNonNegativeDouble {
    public PhasingScore(double score) {
        super(score);
    }

    public PhasingScore(PreciseNonNegativeDouble val) {
        super(val);
    }

    public PhasingScore integrateReadScore(PhasingScore score) {
        timesEqual(score);
        return this;
    }
}

class HaplotypeClass implements Iterable<Haplotype> {
    private ArrayList<Haplotype> haps;
    private Haplotype rep;

    public HaplotypeClass(ArrayList<Haplotype> haps, Haplotype rep) {
        this.haps = haps;
        this.rep = rep;
    }

    public Iterator<Haplotype> iterator() {
        return haps.iterator();
    }

    public Haplotype getRepresentative() {
        return rep;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        boolean isFirst = true;
        for (Haplotype h : haps) {
            if (isFirst)
                isFirst = false;
            else
                sb.append(" + ");

            sb.append(h);
        }
        sb.append(" [").append(rep).append("]");
        return sb.toString();
    }
}

// Summary statistics about phasing rates, for each sample
class PhasingStats {
    private int numReads;
    private int numVarSites;

    // Map of: sample -> PhaseCounts:
    private Map<String, PhaseCounts> samplePhaseStats;

    public PhasingStats() {
        this(new TreeMap<String, PhaseCounts>());
    }

    public PhasingStats(int numReads, int numVarSites) {
        this.numReads = numReads;
        this.numVarSites = numVarSites;
        this.samplePhaseStats = new TreeMap<String, PhaseCounts>();
    }

    public PhasingStats(Map<String, PhaseCounts> samplePhaseStats) {
        this.numReads = 0;
        this.numVarSites = 0;
        this.samplePhaseStats = samplePhaseStats;
    }

    public void addIn(PhasingStats other) {
        this.numReads += other.numReads;
        this.numVarSites += other.numVarSites;

        for (Map.Entry<String, PhaseCounts> sampPhaseEntry : other.samplePhaseStats.entrySet()) {
            String sample = sampPhaseEntry.getKey();
            PhaseCounts otherCounts = sampPhaseEntry.getValue();
            PhaseCounts thisCounts = this.samplePhaseStats.get(sample);
            if (thisCounts == null) {
                thisCounts = new PhaseCounts();
                this.samplePhaseStats.put(sample, thisCounts);
            }
            thisCounts.addIn(otherCounts);
        }
    }

    public int getNumReads() {
        return numReads;
    }

    public int getNumVarSites() {
        return numVarSites;
    }

    public Collection<Map.Entry<String, PhaseCounts>> getPhaseCounts() {
        return samplePhaseStats.entrySet();
    }
}

class PhaseCounts {
    public int numTestedSites; // number of het sites directly succeeding het sites
    public int numInconsistentSitesPhased;
    public int numInconsistentSitesNotPhased;
    public int numPhased;

    public PhaseCounts() {
        this.numTestedSites = 0;
        this.numInconsistentSitesPhased = 0;
        this.numInconsistentSitesNotPhased = 0;
        this.numPhased = 0;
    }

    public void addIn(PhaseCounts other) {
        this.numTestedSites += other.numTestedSites;
        this.numInconsistentSitesPhased += other.numInconsistentSitesPhased;
        this.numInconsistentSitesNotPhased += other.numInconsistentSitesNotPhased;
        this.numPhased += other.numPhased;
    }
}

class PhasingStatsAndOutput {
    public PhasingStats ps;
    public List<VariantContext> output;

    public PhasingStatsAndOutput(PhasingStats ps, List<VariantContext> output) {
        this.ps = ps;
        this.output = output;
    }
}

class PhasingQualityStatsWriter {
    private String variantStatsFilePrefix;
    private HashMap<String, BufferedWriter> sampleToStatsWriter = new HashMap<String, BufferedWriter>();

    public PhasingQualityStatsWriter(String variantStatsFilePrefix) {
        this.variantStatsFilePrefix = variantStatsFilePrefix;
    }

    public void addStat(String sample, GenomeLoc locus, int startDistanceFromPrevious, double phasingQuality, int numReads, int windowSize) {
        BufferedWriter sampWriter = sampleToStatsWriter.get(sample);
        if (sampWriter == null) {
            String fileName = variantStatsFilePrefix + "." + sample + ".locus_distance_PQ_numReads_windowSize.txt";

            FileOutputStream output;
            try {
                output = new FileOutputStream(fileName);
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Unable to create phasing quality stats file at location: " + fileName);
            }
            sampWriter = new BufferedWriter(new OutputStreamWriter(output));
            sampleToStatsWriter.put(sample, sampWriter);
        }
        try {
            sampWriter.write(locus + "\t" + startDistanceFromPrevious + "\t" + phasingQuality + "\t" + numReads + "\t" + windowSize + "\n");
            sampWriter.flush();
        } catch (IOException e) {
            throw new RuntimeException("Unable to write to per-sample phasing quality stats file", e);
        }
    }

    public void close() {
        for (Map.Entry<String, BufferedWriter> sampWriterEntry : sampleToStatsWriter.entrySet()) {
            BufferedWriter sampWriter = sampWriterEntry.getValue();
            try {
                sampWriter.flush();
                sampWriter.close();
            } catch (IOException e) {
                throw new RuntimeException("Unable to close per-sample phasing quality stats file");
            }
        }
    }
}

class SampleReadLocus implements Comparable<SampleReadLocus> {
    private String sample;
    private String read;
    private GenomeLoc locus;

    public SampleReadLocus(String sample, String read, GenomeLoc locus) {
        this.sample = sample;
        this.read = read;
        this.locus = locus;
    }

    public GenomeLoc getLocus() {
        return locus;
    }

    public int compareTo(SampleReadLocus that) {
        int comp = this.sample.compareTo(that.sample);
        if (comp != 0)
            return comp;

        comp = this.read.compareTo(that.read);
        if (comp != 0)
            return comp;

        return this.locus.compareTo(that.locus);
    }

    public String toString() {
        return "Sample " + sample + ", read " + read + ", locus " + locus;
    }
}

class MultipleBaseCounts {
    private Map<Integer, Integer> baseCounts;
    private GenomeLoc phasingLocus;

    public MultipleBaseCounts(GenomeLoc phasingLoc) {
        this.baseCounts = new HashMap<Integer, Integer>();
        this.phasingLocus = phasingLoc;
    }

    public boolean samePhasingLocAs(GenomeLoc loc) {
        return phasingLocus.equals(loc);
    }

    public void incrementBaseCount(byte base) {
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
        Integer cnt = baseCounts.get(baseIndex);
        if (cnt == null)
            cnt = 0;

        baseCounts.put(baseIndex, cnt + 1);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append("Base counts");
        for (Map.Entry<Integer, Integer> baseCountEntry : baseCounts.entrySet()) {
            byte base = BaseUtils.baseIndexToSimpleBase(baseCountEntry.getKey());
            int cnt = baseCountEntry.getValue();
            sb.append("\t" + (char) base + ": " + cnt);
        }

        return sb.toString();
    }
}
