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

import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.*;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * A locally resorting, mate fixing sam file writer that supports an idiom where reads are only moved around if
 * the ISIZE of the pair is < X and reads are not allowed to move any more than Y bp from their original positions.
 *
 * To understand this data structure, let's begin by asking -- when are we certain we know the position of read R added
 * to the writer and its mate M given that R has been added to the writer (but M may not be), their ISIZE in R, at the
 * moment that a read K is added to the writer, under the constraints X and Y?  Complex I know.  First, because
 * reads cannot move more than Y bp in either direction, we know that R originated at most R.pos + Y bp from its
 * current position.  Also, we know that K is at most K.pos + Y bp from it's original position.  If R is maximally
 * shifted to the right, and K shifted to the left, then they could at most move 2Y together.  So if the distance
 * between R and K > 2Y, we know that there are no reads left in the original stream that could be moved before R.
 *
 * Now, we also need to be certain if we have a mate pair M, that won't emit R before we can incorporate any move of
 * M into the mate pair info R.  There are two cases to consider here:
 *
 * If ISIZE > X, we know that we won't move M when we see it, so we can safely emit R knowing that
 * M is fixed in place.
 *
 * If ISIZE <= X, M might be moved, and it we have to wait until we see M in the stream to know it's position.
 * So R must be buffered until either M arrives, or we see a read K that's more than 2Y units past the original position
 * of M.
 *
 * So the worst-case memory consumption here is proportional to the number of reads
 * occurring between R and M + 2 Y, and so is proportional to the depth of the data and X and Y.
 *
 * This leads to the following simple algorithm:
 *
 * addAlignment(newRead):
 *   addReadToListOfReads(newRead)
 *   update mate pair of newRead if present in list of reads
 *
 *   for ( read in list of reads [in order of increasing read.pos] ):
 *     if read.pos < newRead.pos - 2Y && (read.isize >= X || read.matePos < newRead.pos - 2 * Y):
 *        emit read and remove from list of reads
 *     else:
 *        break
 *
 * @author depristo, ebanks
 * @version 0.2
 */
public class ConstrainedMateFixingManager {

    protected static final Logger logger = Logger.getLogger(ConstrainedMateFixingManager.class);
    private static final boolean DEBUG = false;

    /** How often do we check whether we want to emit reads? */
    protected final static int EMIT_FREQUENCY = 1000;

    /**
     * How much could a single read move in position from its original position?
     */
    final int MAX_POS_MOVE_ALLOWED;

    /**
     * How many reads should we store in memory before flushing the queue?
     */
    final int MAX_RECORDS_IN_MEMORY;

    /** how we order our SAM records */
    private final SAMRecordComparator comparer = new SAMRecordCoordinateComparator();

    /** The place where we ultimately write out our records */
    final SAMFileWriter writer;

    /**
     * what is the maximum isize of a pair of reads that can move?  Reads with isize > this value
     * are assumes to not be allowed to move in the incoming read stream.
     */
    final int maxInsertSizeForMovingReadPairs;
    final int initialCapacity = 5000;

    final GenomeLocParser genomeLocParser;
    private GenomeLoc lastLocFlushed = null;

    int counter = 0;

    /** read.name -> records */
    HashMap<String, SAMRecordHashObject> forMateMatching = new HashMap<String, SAMRecordHashObject>();
    PriorityQueue<SAMRecord> waitingReads = new PriorityQueue<SAMRecord>(initialCapacity, comparer);

    private SAMRecord remove(PriorityQueue<SAMRecord> queue) {
        SAMRecord first = queue.poll();
        if (first == null)
            throw new UserException("Error caching SAM record -- priority queue is empty, and yet there was an attempt to poll it -- which is usually caused by malformed SAM/BAM files in which multiple identical copies of a read are present.");
        return first;
    }

    private static class SAMRecordHashObject {
        public SAMRecord record;
        public boolean wasModified;

        public SAMRecordHashObject(SAMRecord record, boolean wasModified) {
            this.record = record;
            this.wasModified = wasModified;
        }
    }


    //private SimpleTimer timer = new SimpleTimer("ConstrainedWriter");
    //private long PROGRESS_PRINT_FREQUENCY = 10 * 1000;             // in milliseconds
    //private long lastProgressPrintTime = -1;                       // When was the last time we printed progress log?


    /**
     *
     * @param writer                                 actual writer
     * @param genomeLocParser                        the GenomeLocParser object
     * @param maxInsertSizeForMovingReadPairs        max insert size allowed for moving pairs
     * @param maxMoveAllowed                         max positional move allowed for any read
     * @param maxRecordsInMemory                     max records to keep in memory
     */
    public ConstrainedMateFixingManager(final SAMFileWriter writer,
                                        final GenomeLocParser genomeLocParser,
                                        final int maxInsertSizeForMovingReadPairs,
                                        final int maxMoveAllowed,
                                        final int maxRecordsInMemory) {
        this.writer = writer;
        this.genomeLocParser = genomeLocParser;
        this.maxInsertSizeForMovingReadPairs = maxInsertSizeForMovingReadPairs;
        this.MAX_POS_MOVE_ALLOWED = maxMoveAllowed;
        this.MAX_RECORDS_IN_MEMORY = maxRecordsInMemory;

        //timer.start();
        //lastProgressPrintTime = timer.currentTime();
    }

    public int getNReadsInQueue() { return waitingReads.size(); }

    /**
     * For testing purposes only
     *
     * @return the list of reads currently in the queue
     */
    protected List<SAMRecord> getReadsInQueueForTesting() {
        return new ArrayList<SAMRecord>(waitingReads);
    }

    public boolean canMoveReads(GenomeLoc earliestPosition) {
        if ( DEBUG ) logger.info("Refusing to realign? " + earliestPosition + " vs. " + lastLocFlushed);

        return lastLocFlushed == null ||
                lastLocFlushed.compareContigs(earliestPosition) != 0 ||
                lastLocFlushed.distance(earliestPosition) > maxInsertSizeForMovingReadPairs;
    }

    private boolean noReadCanMoveBefore(int pos, SAMRecord addedRead) {
        return pos + 2 * MAX_POS_MOVE_ALLOWED < addedRead.getAlignmentStart();
    }

    public void addRead(SAMRecord newRead, boolean readWasModified) {
        addRead(newRead, readWasModified, true);
    }

    public void addReads(List<GATKSAMRecord> newReads, Set<GATKSAMRecord> modifiedReads) {
        for ( GATKSAMRecord newRead : newReads )
            addRead(newRead, modifiedReads.contains(newRead), false);
    }

    protected void addRead(SAMRecord newRead, boolean readWasModified, boolean canFlush) {
        if ( DEBUG ) logger.info("New read pos " + newRead.getAlignmentStart() + " OP = " + newRead.getAttribute("OP") + " " + readWasModified);

        //final long curTime = timer.currentTime();
        //if ( curTime - lastProgressPrintTime > PROGRESS_PRINT_FREQUENCY ) {
        //    lastProgressPrintTime = curTime;
        //    System.out.println("WaitingReads.size = " + waitingReads.size() + ", forMateMatching.size = " + forMateMatching.size());
        //}

        // if the new read is on a different contig or we have too many reads, then we need to flush the queue and clear the map
        boolean tooManyReads = getNReadsInQueue() >= MAX_RECORDS_IN_MEMORY;
        if ( (canFlush && tooManyReads) || (getNReadsInQueue() > 0 && !waitingReads.peek().getReferenceIndex().equals(newRead.getReferenceIndex())) ) {
            if ( DEBUG ) logger.warn("Flushing queue on " + (tooManyReads ? "too many reads" : ("move to new contig: " + newRead.getReferenceName() + " from " + waitingReads.peek().getReferenceName())) + " at " + newRead.getAlignmentStart());

            while ( getNReadsInQueue() > 1 ) {
                // emit to disk
                writeRead(remove(waitingReads));
            }

            SAMRecord lastRead = remove(waitingReads);
            lastLocFlushed = (lastRead.getReferenceIndex() == -1) ? null : genomeLocParser.createGenomeLoc(lastRead);
            writeRead(lastRead);

            if ( !tooManyReads )
                forMateMatching.clear();
            else
                purgeUnmodifiedMates();
        }

        // fix mates, as needed
        // Since setMateInfo can move reads, we potentially need to remove the mate, and requeue
        // it to ensure proper sorting
        if ( isMateFixableRead(newRead) ) {
            SAMRecordHashObject mate = forMateMatching.get(newRead.getReadName());
            if ( mate != null ) {
                // 1. Frustratingly, Picard's setMateInfo() method unaligns (by setting the reference contig
                // to '*') read pairs when both of their flags have the unmapped bit set.  This is problematic
                // when trying to emit reads in coordinate order because all of a sudden we have reads in the
                // middle of the bam file that now belong at the end - and any mapped reads that get emitted
                // after them trigger an exception in the writer.  For our purposes, because we shouldn't be
                // moving read pairs when they are both unmapped anyways, we'll just not run fix mates on them.
                // 2. Furthermore, when reads get mapped to the junction of two chromosomes (e.g. MT since it
                // is actually circular DNA), their unmapped bit is set, but they are given legitimate coordinates.
                // The Picard code will come in and move the read all the way back to its mate (which can be
                // arbitrarily far away).  However, we do still want to move legitimately unmapped reads whose
                // mates are mapped, so the compromise will be that if the mate is still in the queue then we'll
                // move the read and otherwise we won't.
                boolean doNotFixMates = newRead.getReadUnmappedFlag() && (mate.record.getReadUnmappedFlag() || !waitingReads.contains(mate.record));
                if ( !doNotFixMates ) {

                    boolean reQueueMate = mate.record.getReadUnmappedFlag() && ! newRead.getReadUnmappedFlag();
                    if ( reQueueMate ) {
                        // the mate was unmapped, but newRead was mapped, so the mate may have been moved
                        // to be next-to newRead, so needs to be reinserted into the waitingReads queue
                        // note -- this must be called before the setMateInfo call below
                        if ( ! waitingReads.remove(mate.record) )
                            // we must have hit a region with too much depth and flushed the queue
                            reQueueMate = false;
                    }

                    // we've already seen our mate -- set the mate info and remove it from the map;
                    // add/update the mate cigar if appropriate
                    SamPairUtil.setMateInfo(mate.record, newRead, true);
                    if ( reQueueMate ) waitingReads.add(mate.record);
                }

                forMateMatching.remove(newRead.getReadName());
            } else if ( pairedReadIsMovable(newRead) ) {
                forMateMatching.put(newRead.getReadName(), new SAMRecordHashObject(newRead, readWasModified));
            }
        }

        waitingReads.add(newRead);

        if ( ++counter % EMIT_FREQUENCY == 0 ) {
            while ( ! waitingReads.isEmpty() ) { // there's something in the queue
                SAMRecord read = waitingReads.peek();

                if ( noReadCanMoveBefore(read.getAlignmentStart(), newRead) &&
                        (!pairedReadIsMovable(read)                               // we won't try to move such a read
                           || noReadCanMoveBefore(read.getMateAlignmentStart(), newRead ) ) ) { // we're already past where the mate started

                    // remove reads from the map that we have emitted -- useful for case where the mate never showed up
                    if ( !read.getNotPrimaryAlignmentFlag() )
                        forMateMatching.remove(read.getReadName());

                    if ( DEBUG )
                        logger.warn(String.format("EMIT!  At %d: read %s at %d with isize %d, mate start %d, op = %s",
                                newRead.getAlignmentStart(), read.getReadName(), read.getAlignmentStart(),
                                read.getInferredInsertSize(), read.getMateAlignmentStart(), read.getAttribute("OP")));
                    // emit to disk
                    writeRead(remove(waitingReads));
                } else {
                    if ( DEBUG )
                        logger.warn(String.format("At %d: read %s at %d with isize %d couldn't be emited, mate start %d",
                                newRead.getAlignmentStart(), read.getReadName(), read.getAlignmentStart(), read.getInferredInsertSize(), read.getMateAlignmentStart()));
                    break;
                }
            }

            if ( DEBUG ) logger.warn(String.format("At %d: Done with emit cycle", newRead.getAlignmentStart()));
        }
    }

    private void writeRead(SAMRecord read) {
        try {
            if ( writer != null )
                writer.addAlignment(read);
        } catch (IllegalArgumentException e) {
            throw new UserException("If the maximum allowable reads in memory is too small, it may cause reads to be written out of order when trying to write the BAM; please see the --maxReadsInMemory argument for details.  " + e.getMessage(), e);
        }
    }

    /**
     * Is the given read one for which we can fix its mate?
     *
     * @param read  the read
     * @return true if we could fix its mate, false otherwise
     */
    protected boolean isMateFixableRead(final SAMRecord read) {
        return read.getReadPairedFlag() && !read.isSecondaryOrSupplementary();
    }

    /**
     * @param read  the read
     * @return true if the read shouldn't be moved given the constraints of this SAMFileWriter
     */
    public boolean iSizeTooBigToMove(SAMRecord read) {
        return iSizeTooBigToMove(read, maxInsertSizeForMovingReadPairs);               // we won't try to move such a read
    }

    public static boolean iSizeTooBigToMove(SAMRecord read, int maxInsertSizeForMovingReadPairs) {
        return ( read.getReadPairedFlag() && ! read.getMateUnmappedFlag() && !read.getReferenceName().equals(read.getMateReferenceName()) ) // maps to different chromosomes
                || Math.abs(read.getInferredInsertSize()) > maxInsertSizeForMovingReadPairs;     // we won't try to move such a read
    }

    private void purgeUnmodifiedMates() {
        HashMap<String, SAMRecordHashObject> forMateMatchingCleaned = new HashMap<String, SAMRecordHashObject>();
        for ( Map.Entry<String, SAMRecordHashObject> entry : forMateMatching.entrySet() ) {
            if ( entry.getValue().wasModified )
                forMateMatchingCleaned.put(entry.getKey(), entry.getValue());
        }

        forMateMatching.clear(); // explicitly clear the memory
        forMateMatching = forMateMatchingCleaned;
    }

    private boolean pairedReadIsMovable(SAMRecord read) {
        return read.getReadPairedFlag()                                          // we're a paired read
                && (!read.getReadUnmappedFlag() || !read.getMateUnmappedFlag())  // at least one read is mapped
                && !iSizeTooBigToMove(read);                                     // insert size isn't too big

    }

    public void close() {
        // write out all of the remaining reads
        while ( ! waitingReads.isEmpty() ) { // there's something in the queue
            writeRead(remove(waitingReads));
        }
    }
}
