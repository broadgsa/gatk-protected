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

package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.gatk.downsampling.ReservoirDownsampler;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 8/3/11
 * Time: 2:24 PM
 */
public class SlidingWindow {

    // Sliding Window data
    final private TreeSet<GATKSAMRecord> readsInWindow;
    final private LinkedList<HeaderElement> windowHeader;
    protected int contextSize;                                                                                          // the largest context size (between mismatches and indels)
    protected String contig;
    protected int contigIndex;
    protected SAMFileHeader samHeader;
    protected GATKSAMReadGroupRecord readGroupAttribute;
    protected int downsampleCoverage;

    // Running consensus data
    protected SyntheticRead runningConsensus;
    protected int consensusCounter;
    protected String consensusReadName;

    // Filtered Data Consensus data
    protected SyntheticRead filteredDataConsensus;
    protected int filteredDataConsensusCounter;
    protected String filteredDataReadName;


    // Additional parameters
    protected double MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT;                                                        // proportion has to be greater than this value to trigger variant region due to mismatches
    protected double MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;                                                      // proportion has to be greater than this value to trigger variant region due to deletions
    protected int MIN_BASE_QUAL_TO_COUNT;                                                                               // qual has to be greater than or equal to this value
    protected int MIN_MAPPING_QUALITY;

    protected ReduceReads.DownsampleStrategy downsampleStrategy;
    private boolean hasIndelQualities;

    private boolean allowPolyploidReductionInGeneral;

    private static CompressionStash emptyRegions = new CompressionStash();

    /**
     * The types of synthetic reads to use in the finalizeAndAdd method
     */
    private enum ConsensusType {
        CONSENSUS,
        FILTERED,
        BOTH
    }

    public int getStopLocation() {
        return getStopLocation(windowHeader);
    }

    private int getStopLocation(LinkedList<HeaderElement> header) {
        return getStartLocation(header) + header.size() - 1;
    }

    public String getContig() {
        return contig;
    }

    public int getContigIndex() {
        return contigIndex;
    }

    public int getStartLocation(LinkedList<HeaderElement> header) {
        return header.isEmpty() ? -1 : header.peek().getLocation();
    }

    // for testing only
    protected SlidingWindow(final String contig, final int contigIndex, final int startLocation) {
        this.contig = contig;
        this.contigIndex = contigIndex;

        contextSize = 10;

        this.windowHeader = new LinkedList<HeaderElement>();
        windowHeader.addFirst(new HeaderElement(startLocation));
        this.readsInWindow = new TreeSet<GATKSAMRecord>();
    }

    public SlidingWindow(String contig, int contigIndex, int contextSize, SAMFileHeader samHeader, GATKSAMReadGroupRecord readGroupAttribute, int windowNumber, final double minAltProportionToTriggerVariant, final double minIndelProportionToTriggerVariant, int minBaseQual, int minMappingQuality, int downsampleCoverage, final ReduceReads.DownsampleStrategy downsampleStrategy, boolean hasIndelQualities, boolean allowPolyploidReduction) {
        this.contextSize = contextSize;
        this.downsampleCoverage = downsampleCoverage;

        this.MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT = minAltProportionToTriggerVariant;
        this.MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT = minIndelProportionToTriggerVariant;
        this.MIN_BASE_QUAL_TO_COUNT = minBaseQual;
        this.MIN_MAPPING_QUALITY = minMappingQuality;

        this.windowHeader = new LinkedList<HeaderElement>();
        this.readsInWindow = new TreeSet<GATKSAMRecord>(new Comparator<GATKSAMRecord>() {
            @Override
            public int compare(GATKSAMRecord read1, GATKSAMRecord read2) {
                final int difference = read1.getSoftEnd() - read2.getSoftEnd();
                return difference != 0 ? difference : read1.getReadName().compareTo(read2.getReadName());
            }
        });

        this.contig = contig;
        this.contigIndex = contigIndex;
        this.samHeader = samHeader;
        this.readGroupAttribute = readGroupAttribute;

        this.consensusCounter = 0;
        this.consensusReadName = "Consensus-" + windowNumber + "-";

        this.filteredDataConsensusCounter = 0;
        this.filteredDataReadName = "Filtered-" + windowNumber + "-";

        this.runningConsensus = null;
        this.filteredDataConsensus = null;
        
        this.downsampleStrategy = downsampleStrategy;
        this.hasIndelQualities = hasIndelQualities;

        this.allowPolyploidReductionInGeneral = allowPolyploidReduction;
    }

    /**
     * Add a read to the sliding window and slides the window accordingly.
     * 
     * Reads are assumed to be in order, therefore, when a read is added the sliding window can
     * assume that no more reads will affect read.getUnclippedStart() - contextSizeMismatches. The window
     * slides forward to that position and returns all reads that may have been finalized in the
     * sliding process.
     *
     * @param read the read
     * @return a non-null list of reads (in the CompressionStash) that have been finished by sliding the window.
     */
    @Requires({"read != null"})
    @Ensures("result != null")
    public CompressionStash addRead(GATKSAMRecord read) {
        addToHeader(windowHeader, read);                                                                                // update the window header counts
        readsInWindow.add(read);                                                                                        // add read to sliding reads
        return slideWindow(read.getUnclippedStart());
    }

    /**
     * Returns the next complete (or incomplete if closeLastRegion is true) variant region between 'from' (inclusive) and 'to' (exclusive)
     * but converted to global coordinates.
     *
     * @param from         beginning window header index of the search window (inclusive) in local (to the windowHeader) coordinates
     * @param to           end window header index of the search window (exclusive) in local (to the windowHeader) coordinates
     * @param variantSite  boolean array with true marking variant regions
     * @param closeLastRegion  if the last index is variant (so it's an incomplete region), should we close (and return as an interval) the location or ignore it?
     * @return null if nothing is variant, start/stop if there is a complete variant region, start/-1 if there is an incomplete variant region.  All coordinates returned are global.
     */
    @Requires({"from >= 0", "from <= to", "to <= variantSite.length"})
    private FinishedGenomeLoc findNextVariantRegion(int from, int to, boolean[] variantSite, boolean closeLastRegion) {
        boolean foundStart = false;
        final int windowHeaderStart = getStartLocation(windowHeader);
        int variantRegionStartIndex = 0;
        for (int i=from; i<to; i++) {
            if (variantSite[i] && !foundStart) {
                variantRegionStartIndex = i;
                foundStart = true;
            }
            else if(!variantSite[i] && foundStart) {
                return(new FinishedGenomeLoc(contig, contigIndex, windowHeaderStart + variantRegionStartIndex, windowHeaderStart + i - 1, true));
            }
        }
        final int refStart = windowHeaderStart + variantRegionStartIndex;
        final int refStop  = windowHeaderStart + to - 1;
        return (foundStart && closeLastRegion) ? new FinishedGenomeLoc(contig, contigIndex, refStart, refStop, true) : null;
    }

    /**
     * Creates a list with all the complete and incomplete variant regions within 'from' (inclusive) and 'to' (exclusive)
     *
     * @param from         beginning window header index of the search window (inclusive) in local (to the windowHeader) coordinates
     * @param to           end window header index of the search window (exclusive) in local (to the windowHeader) coordinates
     * @param variantSite  boolean array with true marking variant regions
     * @return a list with start/stops of variant regions following findNextVariantRegion description in global coordinates
     */
    @Requires({"from >= 0", "from <= to", "to <= variantSite.length"})
    @Ensures("result != null")
    protected CompressionStash findVariantRegions(int from, int to, boolean[] variantSite, boolean closeLastRegion) {
        final int windowHeaderStart = getStartLocation(windowHeader);

        CompressionStash regions = new CompressionStash();
        int index = from;
        while(index < to) {
            // returns results in global coordinates
            FinishedGenomeLoc result = findNextVariantRegion(index, to, variantSite, closeLastRegion);
            if (result == null)
                break;

            regions.add(result);
            if (!result.isFinished())
                break;

            index = result.getStop() - windowHeaderStart + 1; // go back to local coordinates
        }
        return regions;
    }

    /**
     * Determines if the window can be slid given the new incoming read.
     *
     * We check from the start of the window to the (unclipped) start of the new incoming read if there
     * is any variant.
     * If there are variant sites, we check if it's time to close the variant region.
     *
     * @param incomingReadUnclippedStart the incoming read's start position. Must be the unclipped start!
     * @return all reads that have fallen to the left of the sliding window after the slide
     */
    protected CompressionStash slideWindow(final int incomingReadUnclippedStart) {
        final int windowHeaderStartLocation = getStartLocation(windowHeader);
        CompressionStash regions = emptyRegions;
        boolean forceClose = true;

        if (incomingReadUnclippedStart - contextSize > windowHeaderStartLocation) {
            markSites(incomingReadUnclippedStart);
            int readStartHeaderIndex = incomingReadUnclippedStart - windowHeaderStartLocation;
            int breakpoint = Math.max(readStartHeaderIndex - contextSize - 1, 0);                                       // this is the limit of what we can close/send to consensus (non-inclusive)

            regions = findVariantRegions(0, breakpoint, markedSites.getVariantSiteBitSet(), !forceClose);
        }

        while (!readsInWindow.isEmpty() && readsInWindow.first().getSoftEnd() < windowHeaderStartLocation) {
                readsInWindow.pollFirst();
        }

        return regions;
    }


    protected final class MarkedSites {

        private boolean[] siteIsVariant = new boolean[0];
        private int startLocation = 0;

        public MarkedSites() {}

        public boolean[] getVariantSiteBitSet() { return siteIsVariant; }

        protected int getStartLocation() { return startLocation; }

        /**
         * Updates the variant site bitset given the new startlocation and size of the region to mark.
         *
         * @param newStartLocation   the new start location of the bitset
         * @param sizeOfRegion       the new size of the region to be represented
         *
         * @return the end position (newStartLocation + index) of the region marked by this method; the calling method is responsible for the remainder.
         */
        public int updateRegion(final int newStartLocation, final int sizeOfRegion) {
            int lastPositionMarked = sizeOfRegion;

            // if this is the first time we set the array and we can't reuse anything, just create a new array from scratch
            if ( newStartLocation >= this.startLocation + siteIsVariant.length || newStartLocation < this.startLocation ) {
                siteIsVariant = new boolean[sizeOfRegion];
                lastPositionMarked = 0;
            }
            // if the dimensions change, copy what we can and continue
            else if ( newStartLocation != this.startLocation || sizeOfRegion != siteIsVariant.length ) {
                final boolean[] tempArray = new boolean[sizeOfRegion];
                final int differenceInStartPositions = newStartLocation - this.startLocation;
                lastPositionMarked = Math.min(siteIsVariant.length - differenceInStartPositions, sizeOfRegion);
                System.arraycopy(siteIsVariant, differenceInStartPositions, tempArray, 0, lastPositionMarked);
                siteIsVariant = null;   // explicitly allow garbage collection
                siteIsVariant = tempArray;
            }

            this.startLocation = newStartLocation;

            return lastPositionMarked + newStartLocation;
        }
    }

    private final MarkedSites markedSites = new MarkedSites();

    /**
     * returns an array marked with variant and non-variant regions (it uses
     * markVariantRegions to make the marks)
     *
     * @param stop check the window from start to stop (not-inclusive)
     */
    protected void markSites(final int stop) {

        final int windowHeaderStartLocation = getStartLocation(windowHeader);
        final int sizeOfMarkedRegion = stop - windowHeaderStartLocation + contextSize + 1;

        // copy over as many bits as we can from the previous calculation.  Note that we can't trust the
        // last (contextSize - 1) worth of bits because we may not have actually looked at variant regions there.
        final int lastPositionMarked = markedSites.updateRegion(windowHeaderStartLocation, sizeOfMarkedRegion) - contextSize - 1;
        final int locationToProcess = Math.min(lastPositionMarked, stop - contextSize);

        // update the iterator to the correct position
        Iterator<HeaderElement> headerElementIterator = windowHeader.iterator();
        for (int i = windowHeaderStartLocation; i < locationToProcess; i++) {
            if (headerElementIterator.hasNext())
                headerElementIterator.next();
        }

        // process a contextSize worth of region from scratch in case there's a variant there
        for (int i = locationToProcess; i < stop; i++) {
            if (headerElementIterator.hasNext()) {
                HeaderElement headerElement = headerElementIterator.next();

                if (headerElement.isVariant(MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT, MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT))
                    markVariantRegion(markedSites, i - windowHeaderStartLocation);

            } else
                break;
        }
    }

    /**
     * Marks the sites around the variant site (as true)
     *
     * @param markedSites         the boolean array to bear the marks
     * @param variantSiteLocation the location where a variant site was found
     */
    protected void markVariantRegion(final MarkedSites markedSites, final int variantSiteLocation) {
        int from = (variantSiteLocation < contextSize) ? 0 : variantSiteLocation - contextSize;
        int to = (variantSiteLocation + contextSize + 1 > markedSites.getVariantSiteBitSet().length) ? markedSites.getVariantSiteBitSet().length : variantSiteLocation + contextSize + 1;
        for (int i = from; i < to; i++)
            markedSites.getVariantSiteBitSet()[i] = true;
    }

    /**
     * Adds bases to the running consensus or filtered data accordingly
     * 
     * If adding a sequence with gaps, it will finalize multiple consensus reads and keep the last running consensus
     *
     * @param header  the window header
     * @param start   the first header index to add to consensus
     * @param end     the first header index NOT TO add to consensus
     * @param isNegativeStrand  should the synthetic read be represented as being on the negative strand?
     * @return a non-null list of consensus reads generated by this call. Empty list if no consensus was generated.
     */
    @Requires({"start >= 0 && (end >= start || end == 0)"})
    @Ensures("result != null")
    protected List<GATKSAMRecord> addToSyntheticReads(LinkedList<HeaderElement> header, int start, int end, boolean isNegativeStrand) {
        LinkedList<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        if (start < end) {
            ListIterator<HeaderElement> headerElementIterator = header.listIterator(start);

            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException(String.format("Requested to add to synthetic reads a region that contains no header element at index: %d  - %d / %d", start, header.size(), end));

            HeaderElement headerElement = headerElementIterator.next();

            if (headerElement.hasConsensusData()) {
                reads.addAll(finalizeAndAdd(ConsensusType.FILTERED));

                int endOfConsensus = findNextNonConsensusElement(header, start, end);
                addToRunningConsensus(header, start, endOfConsensus, isNegativeStrand);

                if (endOfConsensus <= start)
                    throw new ReviewedStingException(String.format("next start is <= current start: (%d <= %d)", endOfConsensus, start));

                reads.addAll(addToSyntheticReads(header, endOfConsensus, end, isNegativeStrand));
            } else if (headerElement.hasFilteredData()) {
                reads.addAll(finalizeAndAdd(ConsensusType.CONSENSUS));

                int endOfFilteredData = findNextNonFilteredDataElement(header, start, end);
                reads.addAll(addToFilteredData(header, start, endOfFilteredData, isNegativeStrand));

                if (endOfFilteredData <= start)
                    throw new ReviewedStingException(String.format("next start is <= current start: (%d <= %d)", endOfFilteredData, start));

                reads.addAll(addToSyntheticReads(header, endOfFilteredData, end, isNegativeStrand));
            } else if (headerElement.isEmpty()) {
                reads.addAll(finalizeAndAdd(ConsensusType.BOTH));

                int endOfEmptyData = findNextNonEmptyElement(header, start, end);

                if (endOfEmptyData <= start)
                    throw new ReviewedStingException(String.format("next start is <= current start: (%d <= %d)", endOfEmptyData, start));

                reads.addAll(addToSyntheticReads(header, endOfEmptyData, end, isNegativeStrand));
            } else
                throw new ReviewedStingException(String.format("Header Element %d is neither Consensus, Data or Empty. Something is wrong.", start));

        }

        return reads;
    }

    /**
     * Finalizes one or more synthetic reads.
     *
     * @param type the synthetic reads you want to close
     * @return a possibly null list of GATKSAMRecords generated by finalizing the synthetic reads
     */
    private List<GATKSAMRecord> finalizeAndAdd(ConsensusType type) {
        GATKSAMRecord read = null;
        List<GATKSAMRecord> list = new LinkedList<GATKSAMRecord>();

        switch (type) {
            case CONSENSUS:
                read = finalizeRunningConsensus();
                break;
            case FILTERED:
                read = finalizeFilteredDataConsensus();
                break;
            case BOTH:
                read = finalizeRunningConsensus();
                if (read != null) list.add(read);
                read = finalizeFilteredDataConsensus();
        }
        if (read != null)
            list.add(read);

        return list;
    }

    /**
     * Looks for the next position without consensus data
     *
     * @param start beginning of the filtered region
     * @param upTo  limit to search for another consensus element
     * @return next position in local coordinates (relative to the windowHeader) with consensus data; otherwise, the start position
     */
    private int findNextNonConsensusElement(LinkedList<HeaderElement> header, int start, int upTo) {
        Iterator<HeaderElement> headerElementIterator = header.listIterator(start);
        int index = start;
        while (index < upTo) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("There are no more header elements in this window");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasConsensusData())
                break;
            index++;
        }
        return index;
    }

    /**
     * Looks for the next position without filtered data
     *
     * @param start beginning of the region
     * @param upTo  limit to search for
     * @return next position in local coordinates (relative to the windowHeader) with no filtered data; otherwise, the start position
     */
    private int findNextNonFilteredDataElement(LinkedList<HeaderElement> header, int start, int upTo) {
        Iterator<HeaderElement> headerElementIterator = header.listIterator(start);
        int index = start;
        while (index < upTo) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("There are no more header elements in this window");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasFilteredData() || headerElement.hasConsensusData())
                break;
            index++;
        }
        return index;
    }

    /**
     * Looks for the next non-empty header element
     *
     * @param start beginning of the region
     * @param upTo  limit to search for
     * @return next position in local coordinates (relative to the windowHeader) with non-empty element; otherwise, the start position
     */
    private int findNextNonEmptyElement(LinkedList<HeaderElement> header, int start, int upTo) {
        ListIterator<HeaderElement> headerElementIterator = header.listIterator(start);
        int index = start;
        while (index < upTo) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("There are no more header elements in this window");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.isEmpty())
                break;
            index++;
        }
        return index;
    }


    /**
     * Adds bases to the filtered data synthetic read.
     * 
     * Different from the addToConsensus method, this method assumes a contiguous sequence of filteredData bases.
     *
     * @param header  the window header
     * @param start   the first header index to add to consensus
     * @param end     the first header index NOT TO add to consensus
     * @param isNegativeStrand  should the synthetic read be represented as being on the negative strand?
     * @return a non-null list of GATKSAMRecords representing finalized filtered consensus data. Empty list if no consensus was generated.
     */
    @Requires({"start >= 0 && (end >= start || end == 0)"})
    @Ensures("result != null")
    private List<GATKSAMRecord> addToFilteredData(LinkedList<HeaderElement> header, int start, int end, boolean isNegativeStrand) {
        List<GATKSAMRecord> result = new ArrayList<GATKSAMRecord>(0);

        if (filteredDataConsensus == null)
            filteredDataConsensus = new SyntheticRead(samHeader, readGroupAttribute, contig, contigIndex, filteredDataReadName + filteredDataConsensusCounter++, header.get(start).getLocation(), GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG, hasIndelQualities, isNegativeStrand);

        ListIterator<HeaderElement> headerElementIterator = header.listIterator(start);
        for (int index = start; index < end; index++) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("Requested to create a filtered data synthetic read from " + start + " to " + end + " but " + index + " does not exist");

            HeaderElement headerElement = headerElementIterator.next();
            if (headerElement.hasConsensusData())
                throw new ReviewedStingException("Found consensus data inside region to add to filtered data.");

            if (!headerElement.hasFilteredData())
                throw new ReviewedStingException("No filtered data in " + index);

            if ( filteredDataConsensus.getRefStart() + filteredDataConsensus.size() != headerElement.getLocation() ) {
                result.add(finalizeFilteredDataConsensus());
                filteredDataConsensus = new SyntheticRead(samHeader, readGroupAttribute, contig, contigIndex, filteredDataReadName + filteredDataConsensusCounter++, headerElement.getLocation(), GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG, hasIndelQualities, isNegativeStrand);
            }

            genericAddBaseToConsensus(filteredDataConsensus, headerElement.getFilteredBaseCounts(), headerElement.getRMS());
        }

        return result;
    }

    /**
     * Adds bases to the filtered data synthetic read.
     * 
     * Different from the addToConsensus method, this method assumes a contiguous sequence of filteredData
     * bases.
     *
     * @param header  the window header
     * @param start the first header index to add to consensus
     * @param end   the first header index NOT TO add to consensus
     * @param isNegativeStrand  should the synthetic read be represented as being on the negative strand?
     */
    @Requires({"start >= 0 && (end >= start || end == 0)"})
    private void addToRunningConsensus(LinkedList<HeaderElement> header, int start, int end, boolean isNegativeStrand) {
        if (runningConsensus == null)
            runningConsensus = new SyntheticRead(samHeader, readGroupAttribute, contig, contigIndex, consensusReadName + consensusCounter++, header.get(start).getLocation(), GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG, hasIndelQualities, isNegativeStrand);

        Iterator<HeaderElement> headerElementIterator = header.listIterator(start);
        for (int index = start; index < end; index++) {
            if (!headerElementIterator.hasNext())
                throw new ReviewedStingException("Requested to create a running consensus synthetic read from " + start + " to " + end + " but " + index + " does not exist");

            HeaderElement headerElement = headerElementIterator.next();
            if (!headerElement.hasConsensusData())
                throw new ReviewedStingException("No CONSENSUS data in " + index);

            genericAddBaseToConsensus(runningConsensus, headerElement.getConsensusBaseCounts(), headerElement.getRMS());
        }
    }

    /**
     * Generic accessor to add base and qualities to a synthetic read
     *
     * @param syntheticRead the synthetic read to add to
     * @param baseCounts    the base counts object in the header element
     * @param rms           the rms mapping quality in the header element
     */
    private void genericAddBaseToConsensus(SyntheticRead syntheticRead, BaseAndQualsCounts baseCounts, double rms) {
        final BaseIndex base = baseCounts.baseIndexWithMostProbability();
        byte count = (byte) Math.min(baseCounts.countOfBase(base), Byte.MAX_VALUE);
        byte qual = baseCounts.averageQualsOfBase(base);
        byte insQual = baseCounts.averageInsertionQualsOfBase(base);
        byte delQual = baseCounts.averageDeletionQualsOfBase(base);
        syntheticRead.add(base, count, qual, insQual, delQual, rms);
    }

    /**
     * Method to compress a variant region and return the associated reduced reads
     *
     * @param start   the first window header index in the variant region (inclusive)
     * @param stop    the last window header index of the variant region (inclusive)
     * @param disallowPolyploidReductionAtThisPosition       should we disallow polyploid (het) compression here?
     * @return a non-null list of all reads contained in the variant region
     */
    @Requires({"start >= 0 && (stop >= start || stop == 0)"})
    @Ensures("result != null")
    protected List<GATKSAMRecord> compressVariantRegion(final int start, final int stop, final boolean disallowPolyploidReductionAtThisPosition) {
        List<GATKSAMRecord> allReads = new LinkedList<GATKSAMRecord>();

        // Try to compress into a polyploid consensus
        int nVariantPositions = 0;
        int hetRefPosition = -1;
        boolean canCompress = true;
        Object[] header = windowHeader.toArray();

        // foundEvent will remain false if we don't allow polyploid reduction
        if ( allowPolyploidReductionInGeneral && !disallowPolyploidReductionAtThisPosition ) {
            for (int i = start; i<=stop; i++) {

                int nAlleles = ((HeaderElement) header[i]).getNumberOfAlleles(MIN_ALT_BASE_PROPORTION_TO_TRIGGER_VARIANT);

                // we will only work on diploid cases because we just don't want to handle/test other scenarios
                if ( nAlleles > 2 ) {
                    canCompress = false;
                    break;
                } else if ( nAlleles == 2 ) {
                    nVariantPositions++;

                    // make sure that there is only 1 site in the variant region that contains more than one allele
                    if ( nVariantPositions == 1 ) {
                        hetRefPosition = i;
                    } else if ( nVariantPositions > 1 ) {
                        canCompress = false;
                        break;
                    }
                }
            }
        }

        // Try to compress the variant region; note that using the hetRefPosition protects us from trying to compress
        // variant regions that are created by insertions (since we can't confirm here that they represent the same allele)
        if ( canCompress && hetRefPosition != -1 ) {
            allReads = createPolyploidConsensus(start, stop, ((HeaderElement) header[hetRefPosition]).getLocation());
        }

        // Return all reads that overlap the variant region and remove them from the window header entirely
        // also remove all reads preceding the variant region (since they will be output as consensus right after compression
        else {
            final int refStart = windowHeader.get(start).getLocation();
            final int refStop = windowHeader.get(stop).getLocation();

            LinkedList<GATKSAMRecord> toRemove = new LinkedList<GATKSAMRecord>();
            for (GATKSAMRecord read : readsInWindow) {
                if (read.getSoftStart() <= refStop) {
                    if (read.getAlignmentEnd() >= refStart) {
                        allReads.add(read);
                        removeFromHeader(windowHeader, read);
                    }
                    toRemove.add(read);
                }
            }
            removeReadsFromWindow(toRemove);
        }
        return allReads;
    }

    /**
     * Finalizes a variant region, any adjacent synthetic reads.
     *
     * @param start   the first window header index in the variant region (inclusive)
     * @param stop    the last window header index of the variant region (inclusive)
     * @param disallowPolyploidReductionAtThisPosition       should we disallow polyploid (het) compression here?
     * @return a non-null list of all reads contained in the variant region plus any adjacent synthetic reads
     */
    @Requires({"start >= 0 && (stop >= start || stop == 0)"})
    @Ensures("result != null")
    protected List<GATKSAMRecord> closeVariantRegion(final int start, final int stop, final boolean disallowPolyploidReductionAtThisPosition) {
        List<GATKSAMRecord> allReads = compressVariantRegion(start, stop, disallowPolyploidReductionAtThisPosition);

        List<GATKSAMRecord> result = (downsampleCoverage > 0) ? downsampleVariantRegion(allReads) : allReads;
        result.addAll(addToSyntheticReads(windowHeader, 0, stop, false));
        result.addAll(finalizeAndAdd(ConsensusType.BOTH));

        return result; // finalized reads will be downsampled if necessary
    }

    public Set<GATKSAMRecord> closeVariantRegions(CompressionStash regions) {
        TreeSet<GATKSAMRecord> allReads = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        if (!regions.isEmpty()) {
            int lastStop = -1;
            int windowHeaderStart = getStartLocation(windowHeader);

            for (GenomeLoc region : regions) {
                if (((FinishedGenomeLoc)region).isFinished() && region.getContig() == contig && region.getStart() >= windowHeaderStart && region.getStop() < windowHeaderStart + windowHeader.size()) {
                    int start = region.getStart() - windowHeaderStart;
                    int stop = region.getStop() - windowHeaderStart;

                    allReads.addAll(closeVariantRegion(start, stop, regions.size() > 1)); // todo -- add condition here dependent on dbSNP track
                    lastStop = stop;
                }
            }

            // clean up the window header elements up until the end of the variant region.
            // note that we keep the last element of the region in the event that the following element has a read that starts with insertion.
            if ( lastStop >= 0 ) {
                for (int i = 0; i < lastStop; i++)
                    windowHeader.remove();
                final HeaderElement lastOfRegion = windowHeader.remove();
                if ( lastOfRegion.hasInsertionToTheRight() )
                    windowHeader.addFirst(new HeaderElement(lastOfRegion.getLocation(), lastOfRegion.numInsertionsToTheRight()));
            }
        }
        return allReads;
    }

    /**
     * Downsamples a variant region to the downsample coverage of the sliding window.
     *
     * It will use the downsampling strategy defined by the SlidingWindow
     *
     * @param allReads  a non-null list of reads to select from (all reads that cover the window)
     * @return a non-null list of reads selected by the downsampler to cover the window to at least the desired coverage
     */
    @Requires({"allReads != null"})
    @Ensures("result != null")
    protected List<GATKSAMRecord> downsampleVariantRegion(final List<GATKSAMRecord> allReads) {
        int nReads = allReads.size();
        if (nReads == 0)
            return allReads;

        if (downsampleCoverage >= nReads)
            return allReads;

        ReservoirDownsampler <GATKSAMRecord> downsampler = new ReservoirDownsampler<GATKSAMRecord>(downsampleCoverage);
        downsampler.submit(allReads);
        return downsampler.consumeFinalizedItems();
    }


    /**
     * Properly closes a Sliding Window, finalizing all consensus and variant
     * regions that still exist regardless of being able to fulfill the
     * context size requirement in the end.
     *
     * @return A non-null set/list of all reads generated
     */
    @Ensures("result != null")
    public Pair<Set<GATKSAMRecord>, CompressionStash> close() {
        // mark variant regions
        Set<GATKSAMRecord> finalizedReads = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        CompressionStash regions = new CompressionStash();
        boolean forceCloseUnfinishedRegions = true;

        if (!windowHeader.isEmpty()) {
            markSites(getStopLocation(windowHeader) + 1);
            regions = findVariantRegions(0, windowHeader.size(), markedSites.getVariantSiteBitSet(), forceCloseUnfinishedRegions);
            finalizedReads = closeVariantRegions(regions);

            if (!windowHeader.isEmpty()) {
                finalizedReads.addAll(addToSyntheticReads(windowHeader, 0, windowHeader.size(), false));
                finalizedReads.addAll(finalizeAndAdd(ConsensusType.BOTH));                                              // if it ended in running consensus, finish it up
            }
        }

        return new Pair<Set<GATKSAMRecord>, CompressionStash>(finalizedReads, regions);
    }

    /**
     * generates the SAM record for the running consensus read and resets it (to null)
     *
     * @return the read contained in the running consensus or null
     */
    protected GATKSAMRecord finalizeRunningConsensus() {
        GATKSAMRecord finalizedRead = null;
        if (runningConsensus != null) {
            if (runningConsensus.size() > 0)
                finalizedRead = runningConsensus.close();
            else
                consensusCounter--;

            runningConsensus = null;
        }
        return finalizedRead;
    }

    /**
     * generates the SAM record for the filtered data consensus and resets it (to null)
     *
     * @return the read contained in the running consensus or null
     */
    protected GATKSAMRecord finalizeFilteredDataConsensus() {
        GATKSAMRecord finalizedRead = null;
        if (filteredDataConsensus != null) {
            if (filteredDataConsensus.size() > 0)
                finalizedRead = filteredDataConsensus.close();
            else
                filteredDataConsensusCounter--;

            filteredDataConsensus = null;
        }
        return finalizedRead;
    }

    /**
     * Finalizes a variant region, any adjacent synthetic reads.
     *
     * @param start   the first window header index in the variant region (inclusive)
     * @param stop    the last window header index of the variant region (inclusive)
     * @param hetRefPosition    reference position (in global coordinates) of the het site
     * @return a non-null list of all reads contained in the variant region as a polyploid consensus
     */
    @Requires({"start >= 0 && (stop >= start || stop == 0)"})
    @Ensures("result != null")
    private List<GATKSAMRecord> createPolyploidConsensus(final int start, final int stop, final int hetRefPosition) {
        // we will create two (positive strand, negative strand) headers for each contig
        List<LinkedList<HeaderElement>> headersPosStrand = new ArrayList<LinkedList<HeaderElement>>();
        List<LinkedList<HeaderElement>> headersNegStrand = new ArrayList<LinkedList<HeaderElement>>();
        List<GATKSAMRecord> hetReads = new LinkedList<GATKSAMRecord>();
        Map<Byte, Integer> haplotypeHeaderMap = new HashMap<Byte, Integer>(2);
        int currentHaplotype = 0;
        int refStart = windowHeader.get(start).getLocation();
        int refStop = windowHeader.get(stop).getLocation();
        List<GATKSAMRecord> toRemove = new LinkedList<GATKSAMRecord>();
        for (GATKSAMRecord read : readsInWindow) {
            int haplotype;

            // check if the read is either before or inside the variant region
            if (read.getSoftStart() <= refStop) {
                // check if the read is inside the variant region
                if (read.getMappingQuality() >= MIN_MAPPING_QUALITY && read.getSoftEnd() >= refStart) {
                    // check if the read contains the het site
                    if (read.getSoftStart() <= hetRefPosition && read.getSoftEnd() >= hetRefPosition) {
                        int readPos = ReadUtils.getReadCoordinateForReferenceCoordinate(read, hetRefPosition, ReadUtils.ClippingTail.LEFT_TAIL);
                        // TODO -- THIS IS A HUGE BUG AS IT WILL NOT WORK FOR DELETIONS; see commented out unit test
                        byte base = read.getReadBases()[readPos];
                        byte qual = read.getBaseQualities(EventType.BASE_SUBSTITUTION)[readPos];

                        // check if base passes the filters!
                        if (qual >= MIN_BASE_QUAL_TO_COUNT) {
                            // check which haplotype this read represents and take the index of it from the list of headers
                            if (haplotypeHeaderMap.containsKey(base)) {
                                haplotype = haplotypeHeaderMap.get(base);
                            }
                            // create new lists if this haplotype has not been seen yet
                            else {
                                haplotype = currentHaplotype;
                                haplotypeHeaderMap.put(base, currentHaplotype);
                                headersPosStrand.add(new LinkedList<HeaderElement>());
                                headersNegStrand.add(new LinkedList<HeaderElement>());
                                currentHaplotype++;
                            }
                            LinkedList<HeaderElement> header = read.getReadNegativeStrandFlag() ? headersNegStrand.get(haplotype) : headersPosStrand.get(haplotype);
                            // add to the polyploid header
                            addToHeader(header, read);
                            // remove from the standard header so that we don't double count it
                            removeFromHeader(windowHeader, read);
                        }
                    }
                }

                // we remove all reads before and inside the variant region from the window
                toRemove.add(read);
            }
        }

        for (LinkedList<HeaderElement> header : headersPosStrand) {
            if (header.size() > 0)
                hetReads.addAll(addToSyntheticReads(header, 0, header.size(), false));
            if (runningConsensus != null)
                hetReads.add(finalizeRunningConsensus());
        }
        for (LinkedList<HeaderElement> header : headersNegStrand) {
            if (header.size() > 0)
                hetReads.addAll(addToSyntheticReads(header, 0, header.size(), true));
            if (runningConsensus != null)
                hetReads.add(finalizeRunningConsensus());
        }

        removeReadsFromWindow(toRemove);

        return hetReads;
    }


    private void addToHeader(LinkedList<HeaderElement> header, GATKSAMRecord read) {
        updateHeaderCounts(header, read, false);
    }

    private void removeFromHeader(LinkedList<HeaderElement> header, GATKSAMRecord read) {
        updateHeaderCounts(header, read, true);
    }


    /**
     * Updates the sliding window's header counts with the incoming read bases, insertions
     * and deletions.
     *
     * @param header the sliding window header to use
     * @param read the incoming read to be added to the sliding window
     * @param removeRead if we are removing the read from the header or adding
     */
    private void updateHeaderCounts(final LinkedList<HeaderElement> header, final GATKSAMRecord read, final boolean removeRead) {
        byte[] bases = read.getReadBases();
        byte[] quals = read.getBaseQualities();
        byte[] insQuals = read.getExistingBaseInsertionQualities();
        byte[] delQuals = read.getExistingBaseDeletionQualities();
        int readStart = read.getSoftStart();
        int readEnd = read.getSoftEnd();
        Cigar cigar = read.getCigar();

        int readBaseIndex = 0;
        int startLocation = getStartLocation(header);
        int locationIndex = startLocation < 0 ? 0 : readStart - startLocation;
        int stopLocation = getStopLocation(header);

        if (removeRead && locationIndex < 0)
            throw new ReviewedStingException("read is behind the Sliding Window. read: " + read + " start " + read.getUnclippedStart() + "," + read.getUnclippedEnd() + " cigar: " + read.getCigarString() + " window: " + startLocation + "," + stopLocation);

        if (!removeRead) {                                                                                              // we only need to create new header elements if we are adding the read, not when we're removing it
            if (locationIndex < 0) {                                                                                    // Do we need to add extra elements before the start of the header? -- this may happen if the previous read was clipped and this alignment starts before the beginning of the window
                for (int i = 1; i <= -locationIndex; i++)
                    header.addFirst(new HeaderElement(startLocation - i));

                startLocation = readStart;                                                               // update start location accordingly
                locationIndex = 0;
            }

            if (stopLocation < readEnd) {                                                                // Do we need to add extra elements to the header?
                int elementsToAdd = (stopLocation < 0) ? readEnd - readStart + 1 : readEnd - stopLocation;
                while (elementsToAdd-- > 0)
                    header.addLast(new HeaderElement(readEnd - elementsToAdd));
            }

            // Special case for leading insertions before the beginning of the sliding read
            if (ReadUtils.readStartsWithInsertion(read).getFirst() && (readStart == startLocation || startLocation < 0)) {
                header.addFirst(new HeaderElement(readStart - 1));                                 // create a new first element to the window header with no bases added
                locationIndex = 1;                                                                                      // This allows the first element (I) to look at locationIndex - 1 in the subsequent switch and do the right thing.
            }
        }

        Iterator<HeaderElement> headerElementIterator = header.listIterator(locationIndex);
        HeaderElement headerElement;
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            switch (cigarElement.getOperator()) {
                case H:
                    break;
                case I:
                    if (removeRead && locationIndex == 0) {                                                             // special case, if we are removing a read that starts in insertion and we don't have the previous header element anymore, don't worry about it.
                        break;
                    }

                    headerElement = header.get(locationIndex - 1);                                                // insertions are added to the base to the left (previous element)

                    if (removeRead) {
                        headerElement.removeInsertionToTheRight();
                    }
                    else {
                        headerElement.addInsertionToTheRight();
                    }
                    readBaseIndex += cigarElement.getLength();
                    break;                                                                                              // just ignore the insertions at the beginning of the read
                case D:
                    int nDeletions = cigarElement.getLength();
                    while (nDeletions-- > 0) {                                                                          // deletions are added to the baseCounts with the read mapping quality as it's quality score
                        headerElement = headerElementIterator.next();
                        byte mq = (byte) read.getMappingQuality();
                        if (removeRead)
                            headerElement.removeBase((byte) 'D', mq, mq, mq, mq, MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, false);
                        else
                            headerElement.addBase((byte) 'D', mq, mq, mq, mq, MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, false);

                        locationIndex++;
                    }
                    break;
                case S:
                case M:
                case P:
                case EQ:
                case X:
                    int nBasesToAdd = cigarElement.getLength();
                    while (nBasesToAdd-- > 0) {
                        headerElement = headerElementIterator.next();
                        byte insertionQuality = insQuals == null ? -1 : insQuals[readBaseIndex];                        // if the read doesn't have indel qualities, use -1 (doesn't matter the value because it won't be used for anything)
                        byte deletionQuality = delQuals == null ? -1 : delQuals[readBaseIndex];
                        if (removeRead)
                            headerElement.removeBase(bases[readBaseIndex], quals[readBaseIndex], insertionQuality, deletionQuality, read.getMappingQuality(), MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, cigarElement.getOperator() == CigarOperator.S);
                        else
                            headerElement.addBase(bases[readBaseIndex], quals[readBaseIndex], insertionQuality, deletionQuality, read.getMappingQuality(), MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, cigarElement.getOperator() == CigarOperator.S);

                        readBaseIndex++;
                        locationIndex++;
                    }
                    break;
            }
        }
    }

    private void removeReadsFromWindow (List<GATKSAMRecord> readsToRemove) {
        for (GATKSAMRecord read : readsToRemove) {
            readsInWindow.remove(read);
        }
    }
}

