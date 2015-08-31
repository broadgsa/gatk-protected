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
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.*;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

/**
 * Define intervals to target for local realignment
 *
 * <p>
 * The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases
 * is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion
 * or deletion (indels) in the individual's genome with respect to the reference genome.  Such alignment artifacts result in many bases mismatching
 * the reference near the misalignment, which are easily mistaken as SNPs.  Moreover, since read mapping algorithms operate on each read independently,
 * it is impossible to place reads on the reference genome such that mismatches are minimized across all reads.  Consequently, even when some reads are
 * correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel,
 * also requiring realignment.  Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus
 * indel suitable for standard variant discovery approaches.  Unlike most mappers, this tool uses the full alignment context to determine whether an
 * appropriate alternate reference (i.e. indel) exists.
 * </p>
 *     <p>There are 2 steps to the realignment process:</p>
 *     <ol>
 *     <li>Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)</li>
 *     <li>Running the realigner over those intervals (see the IndelRealigner tool)</li>
 *     </ol>
 *     <p>
 * <p>
 *     For more details, see <a href="http://www.broadinstitute.org/gatk/guide/article?id=38">the indel realignment method documentation</a>.
 * </p>
 *
 * <h3>Inputs</h3>
 * <p>
 * One or more aligned BAM files and optionally, one or more lists of known indels.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A list of target intervals to pass to the IndelRealigner.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T RealignerTargetCreator \
 *   -R reference.fasta \
 *   -I input.bam \
 *   --known indels.vcf \
 *   -o forIndelRealigner.intervals
 * </pre>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>The input BAM(s), reference, and known indel file(s) should be the same ones to be used for the IndelRealigner step.</li>
 *     <li>When multiple potential indels are found by the tool in the same general region, the tool will choose the most likely
 *          one for realignment to the exclusion of the others.  This is a known limitation of the tool.</li>
 *     <li>Because reads produced from the 454 technology inherently contain false indels, the realigner will not work with them
 * (or with reads from similar technologies).</li>
 *     <li>This tool also ignores MQ0 reads and reads with consecutive indel operators in the CIGAR string.</li>
 * </ul>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
@ReadFilters({MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class, BadMateFilter.class, Platform454Filter.class})
@Reference(window=@Window(start=-1,stop=50))
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@By(DataSource.REFERENCE)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
public class RealignerTargetCreator extends RodWalker<RealignerTargetCreator.Event, RealignerTargetCreator.EventPair> implements TreeReducible<RealignerTargetCreator.EventPair> {

    /**
     * The target intervals for realignment.
     */
    @Output
    protected File out;

    /**
     * Any number of VCF files representing known SNPs and/or indels.  Could be e.g. dbSNP and/or official 1000 Genomes indel calls.
     * SNPs in these files will be ignored unless the --mismatchFraction argument is used.
     */
    @Input(fullName="known", shortName = "known", doc="Input VCF file with known indels", required=false)
    public List<RodBinding<VariantContext>> known = Collections.emptyList();

    /**
     * Any two SNP calls and/or high entropy positions are considered clustered when they occur no more than this many basepairs apart. Must be > 1.
     */
    @Argument(fullName="windowSize", shortName="window", doc="window size for calculating entropy or SNP clusters", required=false)
    protected int windowSize = 10;

    /**
     * To disable this behavior, set this value to <= 0 or > 1.  This feature is really only necessary when using an ungapped aligner
     * (e.g. MAQ in the case of single-end read data) and should be used in conjunction with '--model USE_SW' in the IndelRealigner.
     */
    @Argument(fullName="mismatchFraction", shortName="mismatch", doc="fraction of base qualities needing to mismatch for a position to have high entropy", required=false)
    protected double mismatchThreshold = 0.0;

    @Argument(fullName="minReadsAtLocus", shortName="minReads", doc="minimum reads at a locus to enable using the entropy calculation", required=false)
    protected int minReadsAtLocus = 4;

    /**
     * Because the realignment algorithm is N^2, allowing too large an interval might take too long to completely realign.
     */
    @Argument(fullName="maxIntervalSize", shortName="maxInterval", doc="maximum interval size; any intervals larger than this value will be dropped", required=false)
    protected int maxIntervalSize = 500;


    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }


    private boolean lookForMismatchEntropy;
    
    public void initialize() {
        if ( windowSize < 2 )
            throw new UserException.BadArgumentValue("windowSize", "Window Size must be an integer greater than 1");

        lookForMismatchEntropy = mismatchThreshold > 0.0 && mismatchThreshold <= 1.0;
    }

    public Event map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        boolean hasIndel = false;
        boolean hasInsertion = false;
        boolean hasPointEvent = false;

        int furthestStopPos = -1;

        // look at the rods for indels or SNPs
        if ( tracker != null ) {
            for ( VariantContext vc : tracker.getValues(known) ) {
                switch ( vc.getType() ) {
                    case INDEL:
                        hasIndel = true;
                        if ( vc.isSimpleInsertion() )
                            hasInsertion = true;
                        break;
                    case SNP:
                        hasPointEvent = true;
                        break;
                    case MIXED:
                        hasPointEvent = true;
                        hasIndel = true;
                        if ( vc.isSimpleInsertion() )
                            hasInsertion = true;
                        break;
                    default:
                        break;
                }
                if ( hasIndel )
                    furthestStopPos = vc.getEnd();
            }
        }

        // look at the normal context to get deletions and positions with high entropy
        final ReadBackedPileup pileup = context.getBasePileup();

        int mismatchQualities = 0, totalQualities = 0;
        final byte refBase = ref.getBase();
        for ( PileupElement p : pileup ) {

            // check the ends of the reads to see how far they extend
            furthestStopPos = Math.max(furthestStopPos, p.getRead().getAlignmentEnd());

            // is it a deletion or insertion?
            if ( p.isDeletion() || p.isBeforeInsertion() ) {
                hasIndel = true;
                if ( p.isBeforeInsertion() )
                    hasInsertion = true;
            }

            // look for mismatches
            else if ( lookForMismatchEntropy ) {
                if ( p.getBase() != refBase )
                    mismatchQualities += p.getQual();
                totalQualities += p.getQual();
            }
        }

        // make sure we're supposed to look for high entropy
        if ( lookForMismatchEntropy &&
                pileup.getNumberOfElements() >= minReadsAtLocus &&
                (double)mismatchQualities / (double)totalQualities >= mismatchThreshold )
            hasPointEvent = true;

        // return null if no event occurred
        if ( !hasIndel && !hasPointEvent )
            return null;

        // return null if we didn't find any usable reads/rods associated with the event
        if ( furthestStopPos == -1 )
            return null;

        GenomeLoc eventLoc = context.getLocation();
        if ( hasInsertion )
            eventLoc =  getToolkit().getGenomeLocParser().createGenomeLoc(eventLoc.getContig(), eventLoc.getStart(), eventLoc.getStart()+1);

        EVENT_TYPE eventType = (hasIndel ? (hasPointEvent ? EVENT_TYPE.BOTH : EVENT_TYPE.INDEL_EVENT) : EVENT_TYPE.POINT_EVENT);

        return new Event(eventLoc, furthestStopPos, eventType);
    }

    public void onTraversalDone(EventPair sum) {
        if ( sum.left != null && sum.left.isReportableEvent() )
            sum.intervals.add(sum.left.getLoc());
        if ( sum.right != null && sum.right.isReportableEvent() )
            sum.intervals.add(sum.right.getLoc());

        if ( FilenameUtils.getExtension(out.getName()).equals("interval_list") ) {
            final SAMFileHeader masterSequenceDictionaryHeader = new SAMFileHeader();
            masterSequenceDictionaryHeader.setSequenceDictionary(getToolkit().getMasterSequenceDictionary());
            final IntervalList intervalList = new IntervalList(masterSequenceDictionaryHeader);
            for ( GenomeLoc loc : sum.intervals ) {
                intervalList.add(new Interval(loc.getContig(), loc.getStart(), loc.getStop()));
            }
            intervalList.write(out);
        } else {
            try ( BufferedWriter bufferedWriter = IOUtil.openFileForBufferedWriting(out) ) {
                for ( GenomeLoc loc : sum.intervals ) {
                    bufferedWriter.write(loc.toString());
                    bufferedWriter.newLine();
                }
            } catch (final IOException e) {
                throw new GATKException("Error writing out intervals to file: " + out.getAbsolutePath(), e);
            }
        }
    }

    public EventPair reduceInit() {
        return new EventPair(null, null);
    }

    public EventPair treeReduce(EventPair lhs, EventPair rhs) {
        EventPair result;

        if ( lhs.left == null ) {
            result = rhs;
        } else if ( rhs.left == null ) {
            result = lhs;
        } else if ( lhs.right == null ) {
            if ( rhs.right == null ) {
                if ( canBeMerged(lhs.left, rhs.left) )
                    result = new EventPair(mergeEvents(lhs.left, rhs.left), null, lhs.intervals, rhs.intervals);
                else
                    result = new EventPair(lhs.left, rhs.left, lhs.intervals, rhs.intervals);
            } else {
                if ( canBeMerged(lhs.left, rhs.left) )
                    result = new EventPair(mergeEvents(lhs.left, rhs.left), rhs.right, lhs.intervals, rhs.intervals);
                else {
                    if ( rhs.left.isReportableEvent() )
                        rhs.intervals.add(rhs.left.getLoc());
                    result = new EventPair(lhs.left, rhs.right, lhs.intervals, rhs.intervals);
                }
            }
        } else if ( rhs.right == null ) {
            if ( canBeMerged(lhs.right, rhs.left) )
                result = new EventPair(lhs.left, mergeEvents(lhs.right, rhs.left), lhs.intervals, rhs.intervals);
            else {
                if ( lhs.right.isReportableEvent() )
                    lhs.intervals.add(lhs.right.getLoc());
                result = new EventPair(lhs.left, rhs.left, lhs.intervals, rhs.intervals);
            }
        } else {
            if ( canBeMerged(lhs.right, rhs.left) ) {
                Event merge = mergeEvents(lhs.right, rhs.left);
                if ( merge.isReportableEvent() )
                    lhs.intervals.add(merge.getLoc());
            } else {
                if ( lhs.right.isReportableEvent() )
                    lhs.intervals.add(lhs.right.getLoc());
                if ( rhs.left.isReportableEvent() )
                    rhs.intervals.add(rhs.left.getLoc());
            }

            result = new EventPair(lhs.left, rhs.right, lhs.intervals, rhs.intervals);
        }

        return result;
    }

    public EventPair reduce(Event value, EventPair sum) {
        if ( value == null ) {
            ; // do nothing
        } else if ( sum.left == null ) {
            sum.left = value;
        } else if ( sum.right == null ) {
            if ( canBeMerged(sum.left, value) )
                sum.left = mergeEvents(sum.left, value);
            else
                sum.right = value;
        } else {
            if ( canBeMerged(sum.right, value) )
                sum.right = mergeEvents(sum.right, value);
            else {
                if ( sum.right.isReportableEvent() )
                    sum.intervals.add(sum.right.getLoc());
                sum.right = value;
            }
        }

        return sum;
    }

    static private boolean canBeMerged(Event left, Event right) {
        return left.loc.getContigIndex() == right.loc.getContigIndex() && left.furthestStopPos >= right.loc.getStart();
    }

    @com.google.java.contract.Requires({"left != null", "right != null"})
    static private Event mergeEvents(Event left, Event right) {
        left.merge(right);
        return left;
    }

    private enum EVENT_TYPE { POINT_EVENT, INDEL_EVENT, BOTH }

    static class EventPair {
        public Event left, right;
        public TreeSet<GenomeLoc> intervals = new TreeSet<GenomeLoc>();

        public EventPair(Event left, Event right) {
            this.left = left;
            this.right = right;
        }

        public EventPair(Event left, Event right, TreeSet<GenomeLoc> set1, TreeSet<GenomeLoc> set2) {
            this.left = left;
            this.right = right;
            intervals.addAll(set1);
            intervals.addAll(set2);
        }
    }

    class Event {
        public int furthestStopPos;

        private GenomeLoc loc;
        private int eventStartPos;
        private int eventStopPos;
        private EVENT_TYPE type;
        private ArrayList<Integer> pointEvents = new ArrayList<Integer>();

        public Event(GenomeLoc loc, int furthestStopPos, EVENT_TYPE type) {
            this.loc = loc;
            this.furthestStopPos = furthestStopPos;
            this.type = type;

            if ( type == EVENT_TYPE.INDEL_EVENT || type == EVENT_TYPE.BOTH ) {
                eventStartPos = loc.getStart();
                eventStopPos = loc.getStop();
            } else {
                eventStartPos = -1;
                eventStopPos = -1;
            }

            if ( type == EVENT_TYPE.POINT_EVENT || type == EVENT_TYPE.BOTH ) {
                pointEvents.add(loc.getStart());
            }
        }

        public void merge(Event e) {

            // merges only get called for events with certain types
            if ( e.type == EVENT_TYPE.INDEL_EVENT || e.type == EVENT_TYPE.BOTH ) {
                if ( eventStartPos == -1 )
                    eventStartPos = e.eventStartPos;
                eventStopPos = e.eventStopPos;
                furthestStopPos = e.furthestStopPos;
            }

            if ( e.type == EVENT_TYPE.POINT_EVENT || e.type == EVENT_TYPE.BOTH ) {
                int newPosition = e.pointEvents.get(0);
                if ( pointEvents.size() > 0 ) {
                    int lastPosition = pointEvents.get(pointEvents.size()-1);
                    if ( newPosition - lastPosition < windowSize ) {
                        eventStopPos = Math.max(eventStopPos, newPosition);
                        furthestStopPos = e.furthestStopPos;

                        if ( eventStartPos == -1 )
                            eventStartPos = lastPosition;
                        else
                            eventStartPos = Math.min(eventStartPos, lastPosition);
                    } else if ( eventStartPos == -1 && e.eventStartPos != -1 ) {
                        eventStartPos = e.eventStartPos;
                        eventStopPos = e.eventStopPos;
                        furthestStopPos = e.furthestStopPos;
                    }
                }
                pointEvents.add(newPosition);
            }
        }

        public boolean isReportableEvent() {
            return getToolkit().getGenomeLocParser().isValidGenomeLoc(loc.getContig(), eventStartPos, eventStopPos, true) && eventStopPos >= 0 && eventStopPos - eventStartPos < maxIntervalSize;
        }

        public GenomeLoc getLoc() {
            return getToolkit().getGenomeLocParser().createGenomeLoc(loc.getContig(), eventStartPos, eventStopPos);
        }
    }
}