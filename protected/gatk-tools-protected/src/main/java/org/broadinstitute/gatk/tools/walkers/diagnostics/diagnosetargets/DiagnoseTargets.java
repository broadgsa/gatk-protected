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

package org.broadinstitute.gatk.tools.walkers.diagnostics.diagnosetargets;

import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.io.PrintStream;
import java.util.*;

/**
 * Analyze coverage distribution and validate read mates per interval and per sample
 *
 * <p>
 * This tool is useful for diagnosing regions with bad coverage, mapping, or read mate pairs. It analyzes each sample independently and aggregates results over intervals of interest. Low-coverage regions can be identified by using e.g. FindCoveredIntervals with the -uncovered argument.
 * </p>
 * <h3>Input</h3>
 * <ul>
 * <li>A reference file</li>
 * <li>one or more input BAMs</li>
 * <li>One or more intervals</li>
 * </ul>
 * <h3>Output</h3>
 * <p>
 * A modified VCF detailing each interval by sample and information for each interval according to the thresholds used.
 * Interval information includes GC Content, average interval depth, callable status among others. If you use the
 * --missing option, you can get as a second output a intervals file with the loci that have missing data.
 * This file can then be used as input to QualifyMissingIntervals for full qualification and interpretation of why
 * the data is missing.
 * </p>
 * <h3>Usage example</h3>
 * <pre>
 *    java -jar GenomeAnalysisTK.jar
 *              -T DiagnoseTargets \
 *              -R reference.fasta \
 *              -I sample1.bam \
 *              -I sample2.bam \
 *              -I sample3.bam \
 *              -L intervals.interval_list \
 *              -o output.vcf
 *  </pre>
 *
 * @author Mauricio Carneiro, Roger Zurawicki
 * @since 5/8/12
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@By(value = DataSource.READS)
@PartitionBy(PartitionType.INTERVAL)
@Downsample(by = DownsampleType.NONE)
public class DiagnoseTargets extends LocusWalker<Long, Long> {

    @Output(doc = "File to which interval statistics should be written")
    private VariantContextWriter vcfWriter = null;

    @ArgumentCollection
    private ThresHolder thresholds = new ThresHolder();

    private Map<GenomeLoc, IntervalStratification> intervalMap = null;          // maps each interval => statistics
    private PeekableIterator<GenomeLoc> intervalListIterator;                   // an iterator to go over all the intervals provided as we traverse the genome
    private Set<String> samples = null;                                         // all the samples being processed
    private static final Allele SYMBOLIC_ALLELE = Allele.create("<DT>", false); // avoid creating the symbolic allele multiple times
    private static final Allele UNCOVERED_ALLELE = Allele.create("A", true);    // avoid creating the 'fake' ref allele for uncovered intervals multiple times
    private static final int INITIAL_HASH_SIZE = 50;                            // enough room for potential overlapping intervals plus recently finished intervals

    @Override
    public void initialize() {
        super.initialize();

        if (getToolkit().getIntervals() == null || getToolkit().getIntervals().isEmpty())
            throw new UserException("This tool only works if you provide one or more intervals (use the -L argument). If you want to run whole genome, use -T DepthOfCoverage instead.");

        intervalMap = new LinkedHashMap<>(INITIAL_HASH_SIZE);
        intervalListIterator = new PeekableIterator<>(getToolkit().getIntervals().iterator());

        // get all of the unique sample names for the VCF Header
        samples = ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        vcfWriter.writeHeader(new VCFHeader(getHeaderInfo(), samples));

        // pre load all the statistics classes because it is costly to operate on the JVM and we only want to do it once.
        loadAllPlugins(thresholds);
    }

    @Override
    public Long map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();

        // process and remove any intervals in the map that are don't overlap the current locus anymore
        // and add all new intervals that may overlap this reference locus
        addNewOverlappingIntervals(refLocus);
        outputFinishedIntervals(refLocus, ref.getBase());

        // at this point, all intervals in intervalMap overlap with this locus, so update all of them
        for (IntervalStratification intervalStratification : intervalMap.values())
            intervalStratification.addLocus(context, ref);

        return 1L;
    }

    @Override
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Not sure what we are going to do here
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return a long
     */
    @Override
    public Long reduce(Long value, Long sum) {
        return sum + value;
    }

    /**
     * Process all remaining intervals
     *
     * @param result number of loci processed by the walker
     */
    @Override
    public void onTraversalDone(final Long result) {
        for (GenomeLoc interval : intervalMap.keySet())
            outputStatsToVCF(intervalMap.get(interval), UNCOVERED_ALLELE);

        GenomeLoc interval = intervalListIterator.peek();
        while (interval != null) {
            outputStatsToVCF(createIntervalStatistic(interval), UNCOVERED_ALLELE);
            intervalListIterator.next();
            interval = intervalListIterator.peek();
        }

        if (thresholds.missingTargets != null) {
            thresholds.missingTargets.close();
        }
    }

    /**
     * Outputs all intervals that are behind the current reference locus
     *
     * @param refLocus the current reference locus
     * @param refBase  the reference allele
     */
    private void outputFinishedIntervals(final GenomeLoc refLocus, final byte refBase) {
        // output any intervals that were finished
        final List<GenomeLoc> toRemove = new LinkedList<>();
        for (GenomeLoc key : intervalMap.keySet()) {
            if (key.isBefore(refLocus)) {
                final IntervalStratification intervalStats = intervalMap.get(key);
                outputStatsToVCF(intervalStats, Allele.create(refBase, true));
                if (hasMissingLoci(intervalStats)) {
                    outputMissingInterval(intervalStats);
                }
                toRemove.add(key);
            }
        }
        for (GenomeLoc key : toRemove) {
            intervalMap.remove(key);
        }
    }

    /**
     * Adds all intervals that overlap the current reference locus to the intervalMap
     *
     * @param refLocus the current reference locus
     */
    private void addNewOverlappingIntervals(final GenomeLoc refLocus) {
        GenomeLoc interval = intervalListIterator.peek();
        while (interval != null && !interval.isPast(refLocus)) {
            intervalMap.put(interval, createIntervalStatistic(interval));
            intervalListIterator.next();
            interval = intervalListIterator.peek();
        }
    }

    /**
     * Takes the interval, finds it in the stash, prints it to the VCF
     *
     * @param stats     The statistics of the interval
     * @param refAllele the reference allele
     */
    private void outputStatsToVCF(final IntervalStratification stats, final Allele refAllele) {
        GenomeLoc interval = stats.getInterval();

        final List<Allele> alleles = new ArrayList<>();
        final Map<String, Object> attributes = new HashMap<>();
        final ArrayList<Genotype> genotypes = new ArrayList<>();

        for (String sample : samples) {
            final GenotypeBuilder gb = new GenotypeBuilder(sample);

            SampleStratification sampleStat = stats.getSampleStatistics(sample);
            gb.attribute(GATKVCFConstants.AVG_INTERVAL_DP_BY_SAMPLE_KEY, sampleStat.averageCoverage(interval.size()));
            gb.attribute(GATKVCFConstants.LOW_COVERAGE_LOCI, sampleStat.getNLowCoveredLoci());
            gb.attribute(GATKVCFConstants.ZERO_COVERAGE_LOCI, sampleStat.getNUncoveredLoci());
            gb.filters(statusToStrings(stats.getSampleStatistics(sample).callableStatuses(), false));

            genotypes.add(gb.make());
        }
        alleles.add(refAllele);
        alleles.add(SYMBOLIC_ALLELE);
        VariantContextBuilder vcb = new VariantContextBuilder("DiagnoseTargets", interval.getContig(), interval.getStart(), interval.getStop(), alleles);

        vcb = vcb.log10PError(VariantContext.NO_LOG10_PERROR);
        vcb.filters(new LinkedHashSet<>(statusToStrings(stats.callableStatuses(), true)));

        attributes.put(VCFConstants.END_KEY, interval.getStop());
        attributes.put(GATKVCFConstants.AVG_INTERVAL_DP_KEY, stats.averageCoverage(interval.size()));
        attributes.put(GATKVCFConstants.INTERVAL_GC_CONTENT_KEY, stats.gcContent());

        vcb = vcb.attributes(attributes);
        vcb = vcb.genotypes(genotypes);

        vcfWriter.add(vcb.make());
    }

    private boolean hasMissingStatuses(AbstractStratification stats) {
        return !stats.callableStatuses().isEmpty();
    }

    private boolean hasMissingLoci(final IntervalStratification stats) {
        return thresholds.missingTargets != null && hasMissingStatuses(stats);
    }

    private void outputMissingInterval(final IntervalStratification stats) {
        final GenomeLoc interval = stats.getInterval();
        final boolean missing[] = new boolean[interval.size()];
        Arrays.fill(missing, true);
        for (AbstractStratification sample : stats.getElements()) {
            if (hasMissingStatuses(sample)) {
                int pos = 0;
                for (AbstractStratification locus : sample.getElements()) {
                    if (locus.callableStatuses().isEmpty()) {
                        missing[pos] = false;
                    }
                    pos++;
                }
            }
        }
        int start = -1;
        boolean insideMissing = false;
        for (int i = 0; i < missing.length; i++) {
            if (missing[i] && !insideMissing) {
                start = interval.getStart() + i;
                insideMissing = true;
            } else if (!missing[i] && insideMissing) {
                final int stop = interval.getStart() + i - 1;
                outputMissingInterval(interval.getContig(), start, stop);
                insideMissing = false;
            }
        }
        if (insideMissing) {
            outputMissingInterval(interval.getContig(), start, interval.getStop());
        }
    }

    private void outputMissingInterval(final String contig, final int start, final int stop) {
        final PrintStream out = thresholds.missingTargets;
            out.println(String.format("%s:%d-%d", contig, start, stop));
    }

    /**
     * Function that process a set of statuses into strings
     *
     * @param statuses the set of statuses to be converted
     * @return a matching set of strings
     */
    private List<String> statusToStrings(Iterable<CallableStatus> statuses, final boolean isInfoField) {
        List<String> output = new LinkedList<>();

        for (CallableStatus status : statuses)
            if ( isInfoField || status != CallableStatus.PASS )
                output.add(status.name());

        return output;
    }

    private IntervalStratification createIntervalStatistic(GenomeLoc interval) {
        return new IntervalStratification(samples, interval, thresholds);
    }

    protected static void loadAllPlugins(final ThresHolder thresholds) {
        for (Class<?> stat : new PluginManager<LocusMetric>(LocusMetric.class).getPlugins()) {
            try {
                final LocusMetric stats = (LocusMetric) stat.newInstance();
                stats.initialize(thresholds);
                thresholds.locusMetricList.add(stats);
            } catch (Exception e) {
                throw new DynamicClassResolutionException(stat, e);
            }
        }

        for (Class<?> stat : new PluginManager<SampleMetric>(SampleMetric.class).getPlugins()) {
            try {
                final SampleMetric stats = (SampleMetric) stat.newInstance();
                stats.initialize(thresholds);
                thresholds.sampleMetricList.add(stats);
            } catch (Exception e) {
                throw new DynamicClassResolutionException(stat, e);
            }
        }

        for (Class<?> stat : new PluginManager<IntervalMetric>(IntervalMetric.class).getPlugins()) {
            try {
                final IntervalMetric stats = (IntervalMetric) stat.newInstance();
                stats.initialize(thresholds);
                thresholds.intervalMetricList.add(stats);
            } catch (Exception e) {
                throw new DynamicClassResolutionException(stat, e);
            }
        }
    }

    /**
     * Gets the header lines for the VCF writer
     *
     * @return A set of VCF header lines
     */
    private static Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerLines = new HashSet<>();

        // INFO fields for overall data
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AVG_INTERVAL_DP_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.INTERVAL_GC_CONTENT_KEY));
        headerLines.add(new VCFInfoHeaderLine("Diagnose Targets", 0, VCFHeaderLineType.Flag, "DiagnoseTargets mode"));

        // FORMAT fields for each genotype
        headerLines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.AVG_INTERVAL_DP_BY_SAMPLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.LOW_COVERAGE_LOCI));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ZERO_COVERAGE_LOCI));

        // FILTER fields
        for (CallableStatus stat : CallableStatus.values())
            headerLines.add(new VCFFilterHeaderLine(stat.name(), stat.description));

        return headerLines;
    }

}
