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

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import net.sf.picard.util.PeekableIterator;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

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
 * <h2>Input</h2>
 * <p>
 * <ul>
 * <li>A reference file</li>
 * <li>one or more input BAMs</li>
 * <li>One or more intervals</li>
 * </ul>
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * A modified VCF detailing each interval by sample
 * </p>
 * <p/>
 * <h2>Examples</h2>
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

    @Output(doc = "File to which variants should be written", required = true)
    private VariantContextWriter vcfWriter = null;

    @Argument(fullName = "minimum_base_quality", shortName = "BQ", doc = "The minimum Base Quality that is considered for calls", required = false)
    private int minimumBaseQuality = 20;

    @Argument(fullName = "minimum_mapping_quality", shortName = "MQ", doc = "The minimum read mapping quality considered for calls", required = false)
    private int minimumMappingQuality = 20;

    @Argument(fullName = "minimum_coverage", shortName = "min", doc = "The minimum allowable coverage, used for calling LOW_COVERAGE", required = false)
    private int minimumCoverage = 5;

    @Argument(fullName = "maximum_coverage", shortName = "max", doc = "The maximum allowable coverage, used for calling EXCESSIVE_COVERAGE", required = false)
    private int maximumCoverage = 700;

    @Argument(fullName = "minimum_median_depth", shortName = "med", doc = "The minimum allowable median coverage, used for calling LOW_MEDIAN_DEPTH", required = false)
    private int minMedianDepth = 10;

    @Argument(fullName = "maximum_insert_size", shortName = "ins", doc = "The maximum allowed distance between a read and its mate", required = false)
    private int maxInsertSize = 500;

    @Argument(fullName = "voting_status_threshold", shortName = "stV", doc = "The needed percentage of samples containing a call for the interval to adopt the call ", required = false)
    private double votePercentage = 0.50;

    @Argument(fullName = "low_median_depth_status_threshold", shortName = "stMED", doc = "The percentage of the loci needed for calling LOW_MEDIAN_DEPTH", required = false)
    private double lowMedianDepthPercentage = 0.20;

    @Argument(fullName = "bad_mate_status_threshold", shortName = "stBM", doc = "The percentage of the loci needed for calling BAD_MATE", required = false)
    private double badMateStatusThreshold = 0.50;

    @Argument(fullName = "coverage_status_threshold", shortName = "stC", doc = "The percentage of the loci needed for calling LOW_COVERAGE and COVERAGE_GAPS", required = false)
    private double coverageStatusThreshold = 0.20;

    @Argument(fullName = "excessive_coverage_status_threshold", shortName = "stXC", doc = "The percentage of the loci needed for calling EXCESSIVE_COVERAGE", required = false)
    private double excessiveCoverageThreshold = 0.20;

    @Argument(fullName = "quality_status_threshold", shortName = "stQ", doc = "The percentage of the loci needed for calling POOR_QUALITY", required = false)
    private double qualityStatusThreshold = 0.50;

    @Argument(fullName = "print_debug_log", shortName = "dl", doc = "Used only for debugging the walker. Prints extra info to screen", required = false)
    private boolean debug = false;

    private HashMap<GenomeLoc, IntervalStatistics> intervalMap = null;                                                  // maps each interval => statistics
    private PeekableIterator<GenomeLoc> intervalListIterator;                                                           // an iterator to go over all the intervals provided as we traverse the genome
    private Set<String> samples = null;                                                                                 // all the samples being processed
    private final Allele SYMBOLIC_ALLELE = Allele.create("<DT>", false);                                                // avoid creating the symbolic allele multiple times
    private ThresHolder thresholds = null;

    @Override
    public void initialize() {
        super.initialize();

        if (getToolkit().getIntervals() == null)
            throw new UserException("This tool only works if you provide one or more intervals. ( Use the -L argument )");

        thresholds = new ThresHolder(minimumBaseQuality, minimumMappingQuality, minimumCoverage, maximumCoverage, minMedianDepth, maxInsertSize, votePercentage, lowMedianDepthPercentage, badMateStatusThreshold, coverageStatusThreshold, excessiveCoverageThreshold, qualityStatusThreshold);

        intervalMap = new HashMap<GenomeLoc, IntervalStatistics>();
        intervalListIterator = new PeekableIterator<GenomeLoc>(getToolkit().getIntervals().iterator());

        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());                                       // get all of the unique sample names for the VCF Header
        vcfWriter.writeHeader(new VCFHeader(ThresHolder.getHeaderInfo(), samples));                                                 // initialize the VCF header
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();

        removePastIntervals(refLocus, ref.getBase());                                                                   // process and remove any intervals in the map that are don't overlap the current locus anymore
        addNewOverlappingIntervals(refLocus);                                                                           // add all new intervals that may overlap this reference locus    

        for (IntervalStatistics intervalStatistics : intervalMap.values())
            intervalStatistics.addLocus(context, ref, thresholds);                                                      // Add current locus to stats

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
            outputStatsToVCF(intervalMap.get(interval), Allele.create("A", true));
    }

    private GenomeLoc getIntervalMapSpan() {
        GenomeLoc loc = null;
        for (GenomeLoc interval : intervalMap.keySet()) {
            if (loc == null) {
                loc = interval;
            } else
                loc = interval.union(loc);
        }

        return loc;
    }

    private GenomeLoc getFinishedIntervalSpan(GenomeLoc pos) {
        GenomeLoc loc = null;
        for (GenomeLoc interval : intervalMap.keySet()) {
            if (interval.isBefore(pos)) {
                if (loc == null)
                    loc = interval;
                else
                    loc = interval.union(loc);
            }
        }

        return loc;
    }

    /**
     * Removes all intervals that are behind the current reference locus from the intervalMap
     *
     * @param refLocus the current reference locus
     * @param refBase  the reference allele
     */
    private void removePastIntervals(GenomeLoc refLocus, byte refBase) {
        // if there are statistics to output/ check to see that we can output them in order
        if (getFinishedIntervalSpan(refLocus) != null &&
                getIntervalMapSpan().getStart() == getFinishedIntervalSpan(refLocus).getStart()) {

            for (GenomeLoc interval : intervalMap.keySet()) {
                if (interval.isBefore(refLocus)) {
                    outputStatsToVCF(intervalMap.get(interval), Allele.create(refBase, true));
                    intervalMap.remove(interval);
                }
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

        // skip any intervals with no coverage that we have passed
        while (interval != null && interval.isBefore(refLocus)) {
            intervalListIterator.next();                                                                                // discard the interval (we've already added it to the map)
            interval = intervalListIterator.peek();
        }

        // add any intervals that overlap this one
        while (interval != null && !interval.isPast(refLocus)) {
            intervalMap.put(interval, createIntervalStatistic(interval));
            intervalListIterator.next();                                                                                // discard the interval (we've already added it to the map)
            interval = intervalListIterator.peek();
        }
    }

    /**
     * Takes the interval, finds it in the stash, prints it to the VCF
     *
     * @param stats     The statistics of the interval
     * @param refAllele the reference allele
     */
    private void outputStatsToVCF(IntervalStatistics stats, Allele refAllele) {
        GenomeLoc interval = stats.getInterval();


        List<Allele> alleles = new ArrayList<Allele>();
        Map<String, Object> attributes = new HashMap<String, Object>();
        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();

        alleles.add(refAllele);
        alleles.add(SYMBOLIC_ALLELE);
        VariantContextBuilder vcb = new VariantContextBuilder("DiagnoseTargets", interval.getContig(), interval.getStart(), interval.getStop(), alleles);

        vcb = vcb.log10PError(VariantContext.NO_LOG10_PERROR);                                                          // QUAL field makes no sense in our VCF
        vcb.filters(new HashSet<String>(statusesToStrings(stats.callableStatuses(thresholds), true)));

        attributes.put(VCFConstants.END_KEY, interval.getStop());
        attributes.put(ThresHolder.AVG_INTERVAL_DP_KEY, stats.averageCoverage());

        vcb = vcb.attributes(attributes);
        if (debug) {
            System.out.printf("Output -- Interval: %s, Coverage: %.2f%n", stats.getInterval(), stats.averageCoverage());
        }
        for (String sample : samples) {
            final GenotypeBuilder gb = new GenotypeBuilder(sample);

            SampleStatistics sampleStat = stats.getSample(sample);
            gb.attribute(ThresHolder.AVG_INTERVAL_DP_KEY, sampleStat.averageCoverage());
            gb.attribute("Q1", sampleStat.getQuantileDepth(0.25));
            gb.attribute("MED", sampleStat.getQuantileDepth(0.50));
            gb.attribute("Q3", sampleStat.getQuantileDepth(0.75));

            if (debug) {
                System.out.printf("Found %d bad mates out of %d reads %n", sampleStat.getnBadMates(), sampleStat.getnReads());
            }
            gb.filters(statusesToStrings(stats.getSample(sample).getCallableStatuses(thresholds), false));

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
    private List<String> statusesToStrings(Set<CallableStatus> statuses, final boolean includePASS) {
        List<String> output = new ArrayList<String>(statuses.size());

        for (CallableStatus status : statuses)
            if ( includePASS || status != CallableStatus.PASS ) // adding pass => results in a filter for genotypes
                output.add(status.name());

        return output;
    }

    private IntervalStatistics createIntervalStatistic(GenomeLoc interval) {
        return new IntervalStatistics(samples, interval);
    }
}
