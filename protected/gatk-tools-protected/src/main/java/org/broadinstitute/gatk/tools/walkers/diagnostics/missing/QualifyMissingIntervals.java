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

package org.broadinstitute.gatk.tools.walkers.diagnostics.missing;

import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Gather;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportGatherer;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.List;

/**
 * Collect quality metrics for a set of intervals
 *
 * <p>This tool collects the following metrics:</p>
 * <ul>
 *     <li>Average Base Quality</li>
 *     <li>Average Mapping Quality</li>
 *     <li>Average Depth</li>
 *     <li>GC Content</li>
 *     <li>Position in the target (Integer.MIN_VALUE if no overlap)</li>
 *     <li>Length of the overlapping target (zero if no overlap)</li>
 *     <li>Coding Sequence / Intron (optional)</li>
 *     <li>Length of the uncovered interval</li>
 * </ul>
 *
 * <p>It is meant to be run on a set of intervals that have been identified as problematic in earlier stages of quality control and are considered "missing" from the sequence dataset.</p>
 *
 * <h3>Input</h3>
 * <p>
 *  A reference file (for GC content), the input bam file (for base and mapping quality calculation), the missing intervals (in the -L), the baits/targets used to sequence (in the -targets) and a bed file with the coding sequence intervals of the genome (in the -cds).
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 *  GC content, distance from the end of the target, coding sequence intersection, mapping and base quality averages and average depth per "missing" interval.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T QualifyMissingIntervals \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -o output.grp \
 *   -L input.intervals \
 *   -cds cds.intervals \
 *   -targets targets.intervals
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@By(DataSource.REFERENCE)
@PartitionBy(PartitionType.INTERVAL)
public final class QualifyMissingIntervals extends LocusWalker<Metrics, Metrics> implements NanoSchedulable {
    /**
     * A single GATKReport table with the qualifications on why the intervals passed by the -L argument were missing.
     */
    @Gather(GATKReportGatherer.class)
    @Output
    protected PrintStream out;

    /**
     * List of targets used in the experiment. This file will be used to calculate the distance your missing
     * intervals are to the targets (usually exons). Typically this is your hybrid selection targets file
     * (e.g. Agilent exome target list)
     */
    @Argument(shortName = "targets", required = true)
    public String targetsFile;

    /**
     * List of baits to distinguish untargeted intervals from those that are targeted but not covered
     */
    @Argument(shortName = "baits", required = false)
    public String baitsFile = null;

    /**
     * This value will be used to determine whether or not an interval had too high or too low GC content to be
     * sequenced. This is only applied if there was not enough data in the interval.
     */
    @Argument(doc = "upper and lower bound for an interval to be considered high/low GC content",
            shortName = "gc", required = false)
    public double gcThreshold = 0.3;

    /**
     * The coverage of a missing interval may determine whether or not an interval is sequenceable. A low coverage will
     * trigger gc content, mapping, base qualities and other checks to figure out why this interval was deemed
     * unsequenceable.
     */
    @Argument(doc = "minimum coverage to be considered sequenceable",
            shortName = "cov", required = false)
    public int coverageThreshold = 20;

    /**
     * An average mapping quality above this value will determine the interval to be mappable.
     */
    @Argument(doc = "minimum mapping quality for it to be considered usable",
            shortName = "mmq", required = false)
    public byte mappingThreshold = 20;

    /**
     * An average base quality above this value will rule out the possibility of context specific problems with the
     * sequencer.
     */
    @Argument(doc = "minimum base quality for it to be considered usable",
            shortName = "mbq", required = false)
    public byte qualThreshold = 20;

    /**
     * Intervals that are too small generate biased analysis. For example an interval of size 1 will have GC content
     * 1 or 0. To avoid misinterpreting small intervals, all intervals below this threshold will be ignored in the
     * interpretation.
     */
    @Argument(doc = "minimum interval length to be considered",
            shortName = "size", required = false)
    public byte intervalSizeThreshold = 10;

    enum Interpretation {
        UNKNOWN,
        UNMAPPABLE,
        UNSEQUENCEABLE,
        GCCONTENT,
        NO_DATA,
        SMALL_INTERVAL
    }

    GATKReport simpleReport;
    GenomeLocSortedSet targets;
    GenomeLocSortedSet baits;

    public boolean isReduceByInterval() {
        return true;
    }

    public void initialize() {
        // if cds file is not provided, just use the targets file (no harm done)
        if (baitsFile == null)
            baitsFile = targetsFile;

        simpleReport = GATKReport.newSimpleReport("QualifyMissingIntervals", "INTERVAL", "GC", "BQ", "MQ", "DP", "POS_IN_TARGET", "TARGET_SIZE", "BAITED", "MISSING_SIZE", "INTERPRETATION");
        final GenomeLocParser parser = getToolkit().getGenomeLocParser();
        targets = new GenomeLocSortedSet(parser, IntervalUtils.intervalFileToList(parser, targetsFile));
        baits = new GenomeLocSortedSet(parser, IntervalUtils.intervalFileToList(parser, baitsFile));
    }

    public Metrics reduceInit() {
        return new Metrics();
    }

    public Metrics map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        final Metrics metrics = new Metrics();
        final byte baseIndex = ref.getBase();
        final ReadBackedPileup pileup = context.getBasePileup();
        final int nBases = pileup.getNumberOfElements();

        double baseQual = 0.0;
        for (byte qual : pileup.getQuals()) {
            baseQual += qual;
        }
        double mapQual = 0.0;
        for (int qual : pileup.getMappingQuals()) {
            mapQual += qual;
        }

        metrics.baseQual(baseQual);
        metrics.mapQual(mapQual);
        metrics.gccontent(baseIndex == 'C' || baseIndex == 'G' ? 1.0 : 0.0);
        metrics.reads(nBases);
        metrics.refs(1);

        return metrics;
    }

    @Override
    public Metrics reduce(Metrics value, Metrics sum) {
        return sum.combine(value);
    }

    public void onTraversalDone(List<Pair<GenomeLoc, Metrics>> results) {
        for (Pair<GenomeLoc, Metrics> r : results) {
            final GenomeLoc interval = r.getFirst();
            final Metrics metrics = r.getSecond();
            final List<GenomeLoc> overlappingIntervals = targets.getOverlapping(interval);

            simpleReport.addRow(
                    interval.toString(),
                    metrics.gccontent(),
                    metrics.baseQual(),
                    metrics.mapQual(),
                    metrics.depth(),
                    getPositionInTarget(interval, overlappingIntervals),
                    getTargetSize(overlappingIntervals),
                    baits.overlaps(interval),
                    interval.size(),
                    interpret(metrics, interval)
            );
        }
        simpleReport.print(out);
        out.close();
    }

    protected static int getPositionInTarget(final GenomeLoc interval, final List<GenomeLoc> targets) {
        if (targets.size() > 0) {
            final GenomeLoc target = targets.get(0);

            // interval is larger on both ends than the target -- return the maximum distance to either side as a negative number. (min of 2 negative numbers)
            if (interval.getStart() < target.getStart() && interval.getStop() > target.getStop())
                return Math.min(target.getStart() - interval.getStart(), target.getStop() - interval.getStop());

            // interval is a left overlap -- return a negative number representing the distance between the two starts
            else if (interval.getStart() < target.getStart())
                return interval.getStart() - target.getStart();

            // interval is a right overlap -- return a negative number representing the distance between the two stops
            else if (interval.getStop() > target.getStop())
                return target.getStop() - interval.getStop();

            // interval is fully contained -- return the smallest distance to the edge of the target (left or right) as a positive number
            return Math.min(interval.getStart() - target.getStart(), target.getStop() - interval.getStop());
        }
        // if there is no overlapping interval, return int min value.
        return Integer.MIN_VALUE;
    }

    private int getTargetSize(final List<GenomeLoc> overlappingIntervals) {
        return overlappingIntervals.size() > 0 ? overlappingIntervals.get(0).size() : -1;
    }

    String interpret(final Metrics metrics, final GenomeLoc interval) {
        if (interval.size() < intervalSizeThreshold) {
            return Interpretation.SMALL_INTERVAL.toString();
        }
        else if (metrics.depth() == 0.0) {
            return Interpretation.NO_DATA.toString();
        }
        return trim(checkMappability(metrics) + checkGCContent(metrics) + checkContext(metrics));
    }

    String checkMappability(Metrics metrics) {
        return metrics.depth() >= coverageThreshold && metrics.mapQual() < mappingThreshold ?
                Interpretation.UNMAPPABLE + ", " : "";
    }

    String checkGCContent(Metrics metrics) {
        return metrics.depth() < coverageThreshold && (metrics.gccontent() < gcThreshold || metrics.gccontent() > 1.0-gcThreshold) ?
                Interpretation.GCCONTENT + ", " : "";
    }

    String checkContext(Metrics metrics) {
        return metrics.depth() < coverageThreshold && metrics.baseQual() < qualThreshold ?
                Interpretation.UNSEQUENCEABLE + ", " : "";
    }

    String trim (String s) {
        if (s.isEmpty())
            return Interpretation.UNKNOWN.toString();

        s = s.trim();
        if (s.endsWith(","))
            s = s.substring(0, s.length() - 1);
        return s;
    }




}
