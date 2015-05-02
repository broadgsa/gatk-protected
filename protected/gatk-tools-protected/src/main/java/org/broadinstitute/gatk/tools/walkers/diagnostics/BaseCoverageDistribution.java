/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.diagnostics;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Evaluate coverage distribution per base
 *
 * <p>
 * This tool reports the distribution of coverage per base. It includes reads with deletions in the counts unless
 * otherwise specified. Quality filters can be applied before the coverage is calculated.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * The BAM file and an optional interval list
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A GATK Report with the coverage distribution per base
 *
 * <p/>
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T BaseCoverageDistribution \
 *   -I myData.bam \
 *   -L intervals.list \
 *   -fd \
 *   -o report.grp
 * </pre>
 *
 * @author carneiro
 * @since 1/27/13
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class BaseCoverageDistribution extends LocusWalker<ArrayList<Integer>, Map<Integer, ArrayList<Long>>> {
    /**
     * The name of the file to output the GATK Report table. See the FAQs for more information on the GATK Report format.
     */
    @Output(doc = "Output filename")
    private PrintStream out;

    /**
     * Whether or not a deletion should be counted towards the coverage of a site
     */
    @Argument(required = false, shortName="del", fullName = "include_deletions", doc ="Include reads with deletions")
    private boolean includeDeletions = true;

    /**
     * Whether or not to apply quality filters before calculating coverage distribution. Filtering will use the
     * minimum_mapping_quality and minimum_base_quality parameters below.
     */
    @Argument(required = false, shortName="fd", fullName = "filtered_distribution", doc ="Apply quality filters")
    private boolean calculateFilteredDistribution = false;

    /**
     * The minimum mapping quality a read must have to be counted towards the filtered coverage of a site
     */
    @Argument(required = false, shortName="mmq", fullName = "minimum_mapping_quality", doc ="Minimum read mapping quality of a read to pass filters")
    private byte minMappingQuality = 20;

    /**
     * The minimum base quality a base must have to be counted towards the filtered coverage of a site
     */
    @Argument(required = false, shortName="mbq", fullName = "minimum_base_quality", doc ="Minimum base quality to pass filters")
    private byte minBaseQuality = 17;

    private GenomeLoc previousLocus = null;
    private long uncoveredBases = 0L;
    private final LinkedList<GenomeLoc> intervalList = new LinkedList<GenomeLoc>();

    @Override
    public boolean includeReadsWithDeletionAtLoci() {
        return includeDeletions;
    }

    @Override
    public void initialize() {
        if (getToolkit().getIntervals() != null)
            intervalList.addAll(getToolkit().getIntervals()); // if the user provided intervals, keep track of them for uncovered bases calculation
    }

    @Override
    public ArrayList<Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ArrayList<Integer> result = new ArrayList<Integer>(2);
        GenomeLoc currentLocus = ref.getLocus();
        tallyUncoveredBases(currentLocus);
        previousLocus = currentLocus;
        result.add(context.getBasePileup().getReads().size()); // I want the reads instead of the base pileup because I want to count deletions.
        if (calculateFilteredDistribution)
            result.add(context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality).getReads().size()); // filtered pileup
        else {
            result.add(result.get(0)); // repeat the same value as the unfiltered pileup if filters are not on
        }
        return result;
    }

    @Override
    public Map<Integer, ArrayList<Long>> reduceInit() {
        return new HashMap<Integer, ArrayList<Long>>(10000);
    }

    @Override
    public  Map<Integer, ArrayList<Long>> reduce(ArrayList<Integer> value, Map<Integer, ArrayList<Long>> sum) {
        final int unfilteredCoverage = value.get(0);
        final int filteredCoverage = value.get(1);
        incrementSumArray(sum, unfilteredCoverage, 0);
        incrementSumArray(sum, filteredCoverage, 1);
        return sum;
    }

    @Override
    public void onTraversalDone(Map<Integer, ArrayList<Long>> result) {
        tallyUncoveredBasesTillEndOfTraversal();
        GATKReport report;

        if (calculateFilteredDistribution) {
            report = GATKReport.newSimpleReport("BaseCoverageDistribution", "Coverage", "Count", "Filtered");
        } else {
            report = GATKReport.newSimpleReport("BaseCoverageDistribution", "Coverage", "Count");
            report.addRow(0, uncoveredBases); // preemptively add the uncovered bases row (since they'll never exist in the Map)
        }

        for (Map.Entry<Integer, ArrayList<Long>> entry : result.entrySet()) {
            final ArrayList<Long> values = entry.getValue();
            final int coverage = entry.getKey();
            if (calculateFilteredDistribution) {
                if (coverage == 0) { // special case for the uncovered bases. The filtered pileups may have an entry, but the unfiltered ones won't.
                    report.addRow(coverage, uncoveredBases, uncoveredBases + values.get(1));
                } else {
                    report.addRow(coverage, values.get(0), values.get(1));
                }
            } else {
                report.addRow(coverage, values.get(0));
            }
        }
        // In case the filtered distribution never had a pileup filtered down to zero coverage, output the overall uncovered bases for both
        if (calculateFilteredDistribution && !result.containsKey(0)) {
            report.addRow(0, uncoveredBases, uncoveredBases);
        }
        report.print(out);
    }

    /**
     * Initializes the ArrayList if needed. Returns the initialized element (or previously initialized)
     * this method is used directly by the incrementSumArray.
     *
     * @param sum      the map
     * @param coverage the key to the map to extract the array list
     * @return if the ArrayList exists, return it. Otherwise, initialize it with 0 counters.
     */
    private ArrayList<Long> initializeSumArray(final Map<Integer, ArrayList<Long>> sum, final int coverage) {
        ArrayList<Long> curr = sum.get(coverage);
        if (curr == null) {
            curr = new ArrayList<Long>(2);
            curr.add(0L); // number of bases with this unfiltered coverage
            curr.add(0L); // number of bases with this filtered coverage
            sum.put(coverage, curr);
        }
        return curr;
    }

    /**
     * Increments the counter for the given arrayindex (type of coverage : filtered or unfiltered) initializing if necessary
     *
     * @param sum        the hash
     * @param coverage   the hash key
     * @param arrayIndex which distribution to increment, 0 for unfiltered, 1 for filtered.
     */
    private void incrementSumArray(final Map<Integer, ArrayList<Long>> sum, final int coverage, final int arrayIndex) {
        final ArrayList<Long> currentTally = initializeSumArray(sum, coverage);
        currentTally.set(arrayIndex, currentTally.get(arrayIndex) + 1);
    }

    /**
     * Counts all the uncovered loci after the end of traversal.
     *
     * - Modifies the global variable uncoveredBases
     * - Uses global variables: intervalList and previousLocus
     *
     * takes into account that the traversal may have been due over a set of intervals, or over the whole genome.
     */
    private void tallyUncoveredBasesTillEndOfTraversal() {
        GenomeLocParser parser = getToolkit().getGenomeLocParser();
        GenomeLoc lastLocus;
        if (intervalList.isEmpty()) { // whole genome, add up all contigs past previousLocus
            final int lastContigIndex = getToolkit().getSAMFileHeader().getSequenceDictionary().size() - 1;
            final int lastContigLength = getToolkit().getSAMFileHeader().getSequence(lastContigIndex).getSequenceLength();
            final String lastContigName = getToolkit().getSAMFileHeader().getSequence(lastContigIndex).getSequenceName();
            lastLocus = parser.createGenomeLoc(lastContigName, lastContigIndex, lastContigLength, lastContigLength);
        } else {
            GenomeLoc lastInterval = intervalList.getLast();
            lastLocus = parser.createGenomeLoc(lastInterval.getContig(), lastInterval.getContigIndex(), lastInterval.getStop(), lastInterval.getStop());
        }
        tallyUncoveredBases(lastLocus);
    }

    /**
     * Counts all the uncovered loci that have been skipped since the last visited locus. This method allows coverage
     * tools to run with @By(DataSource.READS) instead of @By(DataSource.REFERENCE), while still accurately calculating
     * uncovered bases
     *
     * //todo -- make this a generic capability of Coverage and DiagnoseTargets
     *
     * - Modifies the global variable uncoveredBases
     * - Uses global variables: intervalList and previousLocus
     *
     * takes into account that the traversal may have been due over a set of intervals, or over the whole genome.
     *
     * @param currentLocus the locus we are visiting right now
     */
    private void tallyUncoveredBases(GenomeLoc currentLocus) {
        long distance = 0;
        if (previousLocus == null) { // first base visited
            GenomeLocParser parser = getToolkit().getGenomeLocParser();
            if (intervalList.isEmpty()) { // if this is whole genome (no intervals requested), add what we missed.
                final GenomeLoc zeroLoc = parser.createGenomeLoc(getToolkit().getSAMFileHeader().getSequence(0).getSequenceName(), 0, 1, 1);
                distance += currentLocus.distanceAcrossContigs(zeroLoc, getToolkit().getSAMFileHeader());
            } else { // if we are running on an interval list, add all intervals before the current locus to the uncovered bases counter
                while (!intervalList.peek().containsP(currentLocus)) {
                    GenomeLoc interval = intervalList.removeFirst();
                    distance += interval.size();
                }
                distance += currentLocus.getStart() - intervalList.peek().getStart(); // now this is the interval that contains the current locus. Discount the bases from the beginning.
            }
        } else {
            final GenomeLoc previousInterval = intervalList.peekFirst();  // peekFirst returns null if interval list is empty (WGS).
            distance = currentLocus.distanceAcrossContigs(previousLocus, getToolkit().getSAMFileHeader()) - 1;
            if (previousInterval != null && !previousInterval.containsP(currentLocus)) {
                intervalList.removeFirst(); // we're done with the previous interval
                final GenomeLoc currentInterval = intervalList.peekFirst();
                distance -= currentInterval.distanceAcrossContigs(previousInterval, getToolkit().getSAMFileHeader()) - 1;
            }
        }

        uncoveredBases += distance;
    }
}