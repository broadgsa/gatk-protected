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

package org.broadinstitute.sting.gatk.walkers.diagnostics.diagnosetargets;

import net.sf.picard.util.PeekableIterator;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.*;

import java.util.*;

/**
 * Analyzes coverage distribution and validates read mates for a given interval and sample.
 * <p/>
 * <p>
 * Used to diagnose regions with bad coverage, mapping, or read mating. Analyzes each sample independently in addition
 * to interval wide analysis.
 * </p>
 * <p/>
 * <p/>
 * <h3>Input</h3>
 * <p>
 * <ul>
 * <li>A reference file</li>
 * <li>one or more input BAMs</li>
 * <li>One or more intervals</li>
 * </ul>
 * </p>
 * <p/>
 * <h3>Output</h3>
 * <p>
 * A modified VCF detailing each interval by sample
 * </p>
 * <p/>
 * <h3>Examples</h3>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *              -T DiagnoseTargets \
 *              -R reference.fasta \
 *              -o output.vcf \
 *              -I sample1.bam \
 *              -I sample2.bam \
 *              -I sample3.bam \
 *              -L intervals.interval_list
 *  </pre>
 *
 * @author Mauricio Carneiro, Roger Zurawicki
 * @since 5/8/12
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@By(value = DataSource.READS)
@PartitionBy(PartitionType.INTERVAL)
public class DiagnoseTargets extends LocusWalker<Long, Long> {

    private static final String AVG_INTERVAL_DP_KEY = "IDP";

    @Output(doc = "File to which interval statistics should be written")
    private VariantContextWriter vcfWriter = null;

    @ArgumentCollection
    private ThresHolder thresholds = new ThresHolder();

    private Map<GenomeLoc, IntervalStratification> intervalMap = null;              // maps each interval => statistics
    private PeekableIterator<GenomeLoc> intervalListIterator;                   // an iterator to go over all the intervals provided as we traverse the genome
    private Set<String> samples = null;                                         // all the samples being processed
    private static final Allele SYMBOLIC_ALLELE = Allele.create("<DT>", false); // avoid creating the symbolic allele multiple times
    private static final Allele UNCOVERED_ALLELE = Allele.create("A", true);    // avoid creating the 'fake' ref allele for uncovered intervals multiple times

    private static final int INITIAL_HASH_SIZE = 500000;

    @Override
    public void initialize() {
        super.initialize();

        if (getToolkit().getIntervals() == null || getToolkit().getIntervals().isEmpty())
            throw new UserException("This tool only works if you provide one or more intervals (use the -L argument). If you want to run whole genome, use -T DepthOfCoverage instead.");

        intervalMap = new HashMap<GenomeLoc, IntervalStratification>(INITIAL_HASH_SIZE);
        intervalListIterator = new PeekableIterator<GenomeLoc>(getToolkit().getIntervals().iterator());

        // get all of the unique sample names for the VCF Header
        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        vcfWriter.writeHeader(new VCFHeader(getHeaderInfo(), samples));

        // pre load all the statistics classes because it is costly to operate on the JVM and we only want to do it once.
        loadAllPlugins(thresholds);
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();

        // process and remove any intervals in the map that are don't overlap the current locus anymore
        // and add all new intervals that may overlap this reference locus
        outputFinishedIntervals(refLocus, ref.getBase());
        addNewOverlappingIntervals(refLocus);

        // at this point, all intervals in intervalMap overlap with this locus, so update all of them
        for (IntervalStratification intervalStratification : intervalMap.values())
            intervalStratification.addLocus(context);

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
    public void onTraversalDone(Long result) {
        for (GenomeLoc interval : intervalMap.keySet())
            outputStatsToVCF(intervalMap.get(interval), UNCOVERED_ALLELE);

        GenomeLoc interval = intervalListIterator.peek();
        while (interval != null) {
            outputStatsToVCF(createIntervalStatistic(interval), UNCOVERED_ALLELE);
            intervalListIterator.next();
            interval = intervalListIterator.peek();
        }
    }

    /**
     * Outputs all intervals that are behind the current reference locus
     *
     * @param refLocus the current reference locus
     * @param refBase  the reference allele
     */
    private void outputFinishedIntervals(final GenomeLoc refLocus, final byte refBase) {
        GenomeLoc interval = intervalListIterator.peek();

        // output empty statistics for uncovered intervals
        while (interval != null && interval.isBefore(refLocus)) {
            final IntervalStratification stats = intervalMap.get(interval);
            outputStatsToVCF(stats != null ? stats : createIntervalStatistic(interval), UNCOVERED_ALLELE);
            if (stats != null) intervalMap.remove(interval);
            intervalListIterator.next();
            interval = intervalListIterator.peek();
        }

        // remove any potential leftover interval in intervalMap (this will only happen when we have overlapping intervals)
        for (GenomeLoc key : intervalMap.keySet()) {
            if (key.isBefore(refLocus)) {
                outputStatsToVCF(intervalMap.get(key), Allele.create(refBase, true));
                intervalMap.remove(key);
            }
        }
    }

    /**
     * Adds all intervals that overlap the current reference locus to the intervalMap
     *
     * @param refLocus the current reference locus
     */
    private void addNewOverlappingIntervals(GenomeLoc refLocus) {
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
    private void outputStatsToVCF(IntervalStratification stats, Allele refAllele) {
        GenomeLoc interval = stats.getInterval();


        List<Allele> alleles = new ArrayList<Allele>();
        Map<String, Object> attributes = new HashMap<String, Object>();
        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();

        alleles.add(refAllele);
        alleles.add(SYMBOLIC_ALLELE);
        VariantContextBuilder vcb = new VariantContextBuilder("DiagnoseTargets", interval.getContig(), interval.getStart(), interval.getStop(), alleles);

        vcb = vcb.log10PError(VariantContext.NO_LOG10_PERROR);
        vcb.filters(new LinkedHashSet<String>(statusToStrings(stats.callableStatuses(), true)));

        attributes.put(VCFConstants.END_KEY, interval.getStop());
        attributes.put(AVG_INTERVAL_DP_KEY, stats.averageCoverage(interval.size()));

        vcb = vcb.attributes(attributes);
        for (String sample : samples) {
            final GenotypeBuilder gb = new GenotypeBuilder(sample);

            SampleStratification sampleStat = stats.getSampleStatistics(sample);
            gb.attribute(AVG_INTERVAL_DP_KEY, sampleStat.averageCoverage(interval.size()));

            gb.filters(statusToStrings(stats.getSampleStatistics(sample).callableStatuses(), false));

            genotypes.add(gb.make());
        }
        vcb = vcb.genotypes(genotypes);

        vcfWriter.add(vcb.make());
    }

    /**
     * Function that process a set of statuses into strings
     *
     * @param statuses the set of statuses to be converted
     * @return a matching set of strings
     */
    private List<String> statusToStrings(Iterable<CallableStatus> statuses, final boolean isInfoField) {
        List<String> output = new LinkedList<String>();

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
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

        // INFO fields for overall data
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        headerLines.add(new VCFInfoHeaderLine(AVG_INTERVAL_DP_KEY, 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a loci divided by interval size."));
        headerLines.add(new VCFInfoHeaderLine("Diagnose Targets", 0, VCFHeaderLineType.Flag, "DiagnoseTargets mode"));

        // FORMAT fields for each genotype
        headerLines.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
        headerLines.add(new VCFFormatHeaderLine(AVG_INTERVAL_DP_KEY, 1, VCFHeaderLineType.Float, "Average sample depth across the interval. Sum of the sample specific depth in all loci divided by interval size."));

        // FILTER fields
        for (CallableStatus stat : CallableStatus.values())
            headerLines.add(new VCFFilterHeaderLine(stat.name(), stat.description));

        return headerLines;
    }

}
