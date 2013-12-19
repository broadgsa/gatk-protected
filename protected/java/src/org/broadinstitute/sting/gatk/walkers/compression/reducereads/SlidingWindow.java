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
import it.unimi.dsi.fastutil.bytes.Byte2IntArrayMap;
import it.unimi.dsi.fastutil.bytes.Byte2IntMap;
import it.unimi.dsi.fastutil.objects.*;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.gatk.downsampling.ReservoirDownsampler;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.UnvalidatingGenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
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
    final protected PriorityQueue<GATKSAMRecord> readsInWindow;
    final protected LinkedList<HeaderElement> windowHeader;
    protected int contextSize;                                                                                          // the largest context size (between mismatches and indels)
    protected String contig;
    protected int contigIndex;
    protected SAMFileHeader samHeader;
    protected GATKSAMReadGroupRecord readGroupAttribute;
    protected int downsampleCoverage;

    // Running consensus data
    protected int consensusCounter;
    protected String consensusReadName;

    // Filtered Data Consensus data
    protected int filteredDataConsensusCounter;
    protected String filteredDataReadName;

    // Additional parameters
    protected double MIN_ALT_PVALUE_TO_TRIGGER_VARIANT;                                                                 // pvalue has to be greater than this value to trigger variant region due to mismatches
    protected double MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT;                                                             // proportion has to be greater than this value to trigger variant region due to mismatches
    protected double MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT;                                                      // proportion has to be greater than this value to trigger variant region due to deletions
    protected int MIN_BASE_QUAL_TO_COUNT;                                                                               // qual has to be greater than or equal to this value
    protected int MIN_MAPPING_QUALITY;

    protected ReduceReads.DownsampleStrategy downsampleStrategy;
    private boolean hasIndelQualities;

    private static CompressionStash emptyRegions = new CompressionStash();

    /**
     * The types of synthetic reads
     */
    protected enum ConsensusType {
        POSITIVE_CONSENSUS,
        NEGATIVE_CONSENSUS,
        FILTERED
    }

    public int getStopLocation() {
        return getStopLocation(windowHeader);
    }

    private int getStopLocation(final LinkedList<HeaderElement> header) {
        return header.isEmpty() ? -1 : header.peekLast().getLocation();
    }

    public String getContig() {
        return contig;
    }

    public int getContigIndex() {
        return contigIndex;
    }

    public int getStartLocation(final LinkedList<HeaderElement> header) {
        return header.isEmpty() ? -1 : header.peek().getLocation();
    }

    // for testing only
    protected SlidingWindow(final String contig, final int contigIndex, final int startLocation) {
        this.contig = contig;
        this.contigIndex = contigIndex;

        contextSize = 10;

        this.windowHeader = new LinkedList<>();
        windowHeader.addFirst(new HeaderElement(startLocation));
        this.readsInWindow = new PriorityQueue<>(100, new Comparator<GATKSAMRecord>() {
            @Override
            public int compare(GATKSAMRecord read1, GATKSAMRecord read2) {
                return read1.getSoftEnd() - read2.getSoftEnd();
            }
        });
    }

    public SlidingWindow(final String contig, final int contigIndex, final int contextSize, final SAMFileHeader samHeader,
                         final GATKSAMReadGroupRecord readGroupAttribute, final int windowNumber,
                         final double minAltPValueToTriggerVariant, final double minAltProportionToTriggerVariant, final double minIndelProportionToTriggerVariant,
                         final int minBaseQual, final int minMappingQuality, final int downsampleCoverage,
                         final ReduceReads.DownsampleStrategy downsampleStrategy, final boolean hasIndelQualities) {
        this.contextSize = contextSize;
        this.downsampleCoverage = downsampleCoverage;

        this.MIN_ALT_PVALUE_TO_TRIGGER_VARIANT = minAltPValueToTriggerVariant;
        this.MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT = minAltProportionToTriggerVariant;
        this.MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT = minIndelProportionToTriggerVariant;
        this.MIN_BASE_QUAL_TO_COUNT = minBaseQual;
        this.MIN_MAPPING_QUALITY = minMappingQuality;

        this.windowHeader = new LinkedList<>();
        this.readsInWindow = new PriorityQueue<>(1000, new Comparator<GATKSAMRecord>() {
            @Override
            public int compare(GATKSAMRecord read1, GATKSAMRecord read2) {
                return read1.getSoftEnd() - read2.getSoftEnd();
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

        this.downsampleStrategy = downsampleStrategy;
        this.hasIndelQualities = hasIndelQualities;
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
        // no need to track low mapping quality reads
        if ( read.getMappingQuality() >= MIN_MAPPING_QUALITY )
            readsInWindow.add(read);                                                                                    // add read to sliding reads
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

        while (!readsInWindow.isEmpty() && readsInWindow.peek().getSoftEnd() < windowHeaderStartLocation) {
            readsInWindow.poll();
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
     * returns the MarkedSites object so that it can be tested after adding data to the Sliding Window
     *
     * @return the Marked Sites object used by this Sliding Window
     */
    protected MarkedSites getMarkedSitesForTesting() { return markedSites; }

    /**
     * returns an array marked with variant and non-variant regions (it uses markVariantRegion to make the marks)
     *
     * @param stop check the window from start to stop (not-inclusive); given in global coordinates
     */
    protected void markSites(final int stop) {

        final int windowHeaderStartLocation = getStartLocation(windowHeader);
        final int sizeOfMarkedRegion = stop - windowHeaderStartLocation + contextSize + 1;

        // copy over as many bits as we can from the previous calculation.  Note that we can't trust the
        // last (contextSize - 1) worth of bits because we may not have actually looked at variant regions there.
        final int lastPositionMarked = markedSites.updateRegion(windowHeaderStartLocation, sizeOfMarkedRegion) - contextSize - 1;
        final int locationToProcess = Math.max(windowHeaderStartLocation, Math.min(lastPositionMarked, stop - contextSize));

        final ListIterator<HeaderElement> headerElementIterator = windowHeader.listIterator(locationToProcess - windowHeaderStartLocation);

        // process a contextSize worth of region from scratch in case there's a variant there
        for (int i = locationToProcess; i < stop; i++) {
            if (headerElementIterator.hasNext()) {
                HeaderElement headerElement = headerElementIterator.next();

                if (headerElement.isVariant(MIN_ALT_PVALUE_TO_TRIGGER_VARIANT, MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT, MIN_INDEL_BASE_PROPORTION_TO_TRIGGER_VARIANT))
                    markVariantRegion(i - windowHeaderStartLocation);

            } else
                break;
        }
    }

    /**
     * Marks the sites around the variant site (as true)
     *
     * @param variantSiteLocation the location where a variant site was found
     */
    protected void markVariantRegion(final int variantSiteLocation) {
        int from = (variantSiteLocation < contextSize) ? 0 : variantSiteLocation - contextSize;
        int to = (variantSiteLocation + contextSize + 1 > markedSites.getVariantSiteBitSet().length) ? markedSites.getVariantSiteBitSet().length - 1 : variantSiteLocation + contextSize;
        markRegionAs(from, to, true);
    }

    /**
     * Marks the sites around the variant site (as true)
     *
     * @param from              the start index (inclusive) to mark
     * @param to                the end index (inclusive) to mark
     * @param isVariant         mark the region with this boolean value
     */
    private void markRegionAs(final int from, final int to, final boolean isVariant) {
        for (int i = from; i <= to; i++)
            markedSites.getVariantSiteBitSet()[i] = isVariant;
    }

    /**
     * Adds bases to the running consensus
     * 
     * If adding a sequence with gaps, it will finalize multiple consensus reads and keep the last running consensus
     *
     * @param header  the header to use
     * @param start   the first header index to add to consensus
     * @param end     the first header index NOT TO add to consensus
     * @param consensusType the consensus type to use
     * @return a non-null list of consensus reads generated by this call. Empty list if no consensus was generated.
     */
    @Requires({"start >= 0 && (end >= start || end == 0)"})
    @Ensures("result != null")
    protected ObjectArrayList<GATKSAMRecord> addToSyntheticReads(final LinkedList<HeaderElement> header, final int start, final int end, final ConsensusType consensusType) {
        final ObjectArrayList<GATKSAMRecord> reads = new ObjectArrayList<>();

        SyntheticRead consensus = null;
        final ListIterator<HeaderElement> headerElementIterator = header.listIterator(start);
        boolean wasInConsensus = false;

        for ( int currentPosition = start; currentPosition < end; currentPosition++ ) {

            if ( ! headerElementIterator.hasNext() )
                throw new IllegalStateException(String.format("Requested to add to synthetic reads a region that contains no header element at index: %d  - %d / %d", start, windowHeader.size(), end));
            final HeaderElement headerElement = headerElementIterator.next();

            if ( headerElement.hasConsensusData(consensusType) ) {
                wasInConsensus = true;

                // add to running consensus
                if ( consensus == null )
                    consensus = createNewConsensus(consensusType, headerElement.getLocation());

                genericAddBaseToConsensus(consensus, headerElement.getBaseCounts(consensusType));

            } else {

                // add any outstanding consensus data
                if ( wasInConsensus ) {
                    reads.addAll(finalizeAndAdd(consensus, consensusType));
                    consensus = null;
                }

                wasInConsensus = false;
            }
        }

        // add any outstanding consensus data
        reads.addAll(finalizeAndAdd(consensus, consensusType));

        return reads;
    }

    private SyntheticRead createNewConsensus(final ConsensusType consensusType, final int start) {
        if ( consensusType == ConsensusType.FILTERED )
            return new SyntheticRead(samHeader, readGroupAttribute, contig, contigIndex, filteredDataReadName + filteredDataConsensusCounter++, start, hasIndelQualities, SyntheticRead.StrandType.STRANDLESS);
        return new SyntheticRead(samHeader, readGroupAttribute, contig, contigIndex, consensusReadName + consensusCounter++, start, hasIndelQualities, consensusType == ConsensusType.POSITIVE_CONSENSUS ? SyntheticRead.StrandType.POSITIVE : SyntheticRead.StrandType.NEGATIVE);
    }

    /**
     * Finalizes a synthetic read.
     *
     * @param consensus the consensus to finalize
     * @param type the synthetic reads you want to close
     * @return a possibly empty list of GATKSAMRecords generated by finalizing the synthetic reads
     */
    private ObjectArrayList<GATKSAMRecord> finalizeAndAdd(final SyntheticRead consensus, final ConsensusType type) {

        final ObjectArrayList<GATKSAMRecord> list = new ObjectArrayList<>();

        final GATKSAMRecord read;
        if ( type == ConsensusType.FILTERED )
            read = finalizeFilteredDataConsensus(consensus);
        else
            read = finalizeRunningConsensus(consensus);

        if ( read != null )
            list.add(read);

        return list;
    }

    /**
     * Generic accessor to add base and qualities to a synthetic read
     *
     * @param syntheticRead the synthetic read to add to
     * @param baseCounts    the base counts object in the header element
     */
    private void genericAddBaseToConsensus(final SyntheticRead syntheticRead, final BaseAndQualsCounts baseCounts) {
        final BaseIndex base = baseCounts.baseIndexWithMostProbability();
        final int count = baseCounts.countOfBase(base);
        final byte qual = baseCounts.averageQualsOfBase(base);
        final byte insQual = baseCounts.averageInsertionQualsOfBase(base);
        final byte delQual = baseCounts.averageDeletionQualsOfBase(base);
        syntheticRead.add(base, count, qual, insQual, delQual, baseCounts.getRMS());
    }

    /**
     * Method to compress a variant region and return the associated reduced reads
     *
     * @param start   the first window header index in the variant region (inclusive)
     * @param stop    the last window header index of the variant region (inclusive)
     * @param knownSnpPositions  the set of known SNPs used to determine whether to allow polyploid consensus creation here; can be null (to allow polyploid consensus anywhere)
     * @return a non-null object representing all reads contained in the variant region
     */
    @Requires({"start >= 0 && (stop >= start || stop == 0)"})
    @Ensures("result != null")
    protected CloseVariantRegionResult compressVariantRegion(final int start, final int stop, final ObjectSortedSet<GenomeLoc> knownSnpPositions) {
        final CloseVariantRegionResult allReads = new CloseVariantRegionResult(stop);

        // Try to compress into a polyploid consensus
        // Optimization: don't bother if there are no known SNPs here
        final int hetRefPosition = (knownSnpPositions != null && knownSnpPositions.isEmpty()) ? -1 : findSinglePolyploidCompressiblePosition(start, stop);

        // Note that using the hetRefPosition protects us from trying to compress variant regions that are created by
        //   insertions (which we don't want because we can't confirm that they represent the same allele).
        // Also, we only allow polyploid consensus creation at known sites if provided.
        if ( hetRefPosition != -1 && matchesKnownPosition(windowHeader.get(hetRefPosition).getLocation(), knownSnpPositions) ) {
            // try to create the polyploid consensus
            allReads.reads.addAll(createPolyploidConsensus(hetRefPosition));
            allReads.stopPerformed = hetRefPosition;  // we stopped at the het position
        }
        // if we can't create a polyploid consensus here, return all reads that overlap the variant region and remove them
        // from the window header entirely; also remove all reads preceding the variant region (since they will be output
        // as consensus right after compression)
        else {
            final int refStart = windowHeader.get(start).getLocation();
            final int refStop = windowHeader.get(stop).getLocation();

            final ObjectList<GATKSAMRecord> toRemoveFromWindow = new ObjectArrayList<>();
            final ObjectList<GATKSAMRecord> toEmit = new ObjectArrayList<>();
            for ( final GATKSAMRecord read : readsInWindow ) {
                if ( read.getSoftStart() <= refStop ) {
                    if ( read.getAlignmentEnd() >= refStart ) {
                        toEmit.add(read);
                        removeFromHeader(windowHeader, read);
                    }
                    toRemoveFromWindow.add(read);
                }
            }

            // remove all used reads
            for ( final GATKSAMRecord read : toRemoveFromWindow )
                readsInWindow.remove(read);

            // down-sample the unreduced reads if needed
            allReads.reads.addAll(downsampleCoverage > 0 ? downsampleVariantRegion(toEmit) : toEmit);
        }

        return allReads;
    }

    /**
     * Determines whether the given position match one of the known sites
     *
     * @param targetPosition     the position of the het site
     * @param knownSnpPositions  the set of known SNPs used to determine whether to allow polyploid consensus creation here; can be null (to allow polyploid consensus anywhere)
     * @return true if the targetPosition matches a known SNP position, false otherwise
     */
    @Requires({"targetPosition >= 1 && knownSnpPositions != null"})
    protected boolean matchesKnownPosition(final int targetPosition, final ObjectSortedSet<GenomeLoc> knownSnpPositions) {
        final GenomeLoc targetLoc = new UnvalidatingGenomeLoc(contig, contigIndex, targetPosition, targetPosition);
        return knownSnpPositions == null || knownSnpPositions.contains(targetLoc);
    }

    /*
     * Finds the het variant position located within start and stop (inclusive) if one exists.
     *
     * @param start   the first header index in the region to check (inclusive)
     * @param stop    the last header index of the region to check (inclusive)
     * @return the window header index of the single het position or -1 if either none or more than one exists
     */
    @Requires("start >= 0 && (stop >= start || stop == 0)")
    protected int findSinglePolyploidCompressiblePosition(final int start, final int stop) {
        int hetRefPosition = -1;

        for ( int i = start; i <= stop; i++ ) {

            final int nAlleles = windowHeader.get(i).getNumberOfBaseAlleles(MIN_ALT_PVALUE_TO_TRIGGER_VARIANT, MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT);

            // we will only work on diploid non-indel cases because we just don't want to handle/test other scenarios
            if ( nAlleles > 2 || nAlleles == -1 )
                return -1;

            if ( nAlleles == 2 ) {

                // make sure that there is only 1 site in the region that contains more than one allele
                if ( hetRefPosition != -1 )
                    return -1;

                hetRefPosition = i;
            }
        }

        return hetRefPosition;
    }

    /*
     * Checks whether there's a position in the header with a significant number of softclips or a variant.
     *
     * @param header          the window header to examine
     * @param positionToSkip  the global position to skip in the examination (use negative number if you don't want to make use of this argument)
     * @return true if there exists a position with significant softclips, false otherwise
     */
    @Requires("header != null")
    protected boolean hasPositionWithSignificantSoftclipsOrVariant(final List<HeaderElement> header, final int positionToSkip) {

        for ( final HeaderElement headerElement : header ) {

            if ( headerElement.getLocation() == positionToSkip )
                continue;

            if ( headerElement.hasSignificantSoftclips(MIN_ALT_PVALUE_TO_TRIGGER_VARIANT, MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT) ||
                 headerElement.getNumberOfBaseAlleles(MIN_ALT_PVALUE_TO_TRIGGER_VARIANT, MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT) != 1 )
                return true;
        }

        return false;
    }

    /**
     * Finalizes a variant region, any adjacent synthetic reads.
     *
     * @param start   the first window header index in the variant region (inclusive)
     * @param stop    the last window header index of the variant region (inclusive)
     * @param knownSnpPositions  the set of known SNPs used to determine whether to allow polyploid consensus creation here; can be null (to allow polyploid consensus anywhere)
     * @return a non-null object representing all reads contained in the variant region plus any adjacent synthetic reads
     */
    @Requires({"start >= 0 && (stop >= start || stop == 0)"})
    @Ensures("result != null")
    protected CloseVariantRegionResult closeVariantRegion(final int start, final int stop, final ObjectSortedSet<GenomeLoc> knownSnpPositions) {
        final CloseVariantRegionResult allReads = compressVariantRegion(start, stop, knownSnpPositions);
        allReads.reads.addAll(addAllSyntheticReadTypes(0, allReads.stopPerformed + 1));
        return allReads;
    }

    /**
     * Adds reads for all possible strands (positive, negative, filtered) from the global windowHeader object
     *
     * @param start   the start position (inclusive)
     * @param end     the end position (exclusive)
     * @return non-null but possibly empty array list with reduced reads
     */
    private ObjectArrayList<GATKSAMRecord> addAllSyntheticReadTypes(final int start, final int end) {
        final ObjectArrayList<GATKSAMRecord> reads = new ObjectArrayList<>();
        reads.addAll(addToSyntheticReads(windowHeader, start, end, ConsensusType.POSITIVE_CONSENSUS));
        reads.addAll(addToSyntheticReads(windowHeader, start, end, ConsensusType.NEGATIVE_CONSENSUS));
        reads.addAll(addToSyntheticReads(windowHeader, start, end, ConsensusType.FILTERED));
        return reads;
    }

    /*
     * @see #closeVariantRegions(CompressionStash, ObjectSortedSet<GenomeLoc>, boolean) with forceCloseFullRegions set to false
     */
    public ObjectSet<GATKSAMRecord> closeVariantRegions(final CompressionStash regions, final ObjectSortedSet<GenomeLoc> knownSnpPositions) {
        return closeVariantRegions(regions, knownSnpPositions, false);
    }

    private static final class CloseVariantRegionResult {
        final private ObjectList<GATKSAMRecord> reads = new ObjectArrayList<>();
        private int stopPerformed;

        public CloseVariantRegionResult(final int stopPerformed) { this.stopPerformed = stopPerformed; }
    }

    /*
     * Finalizes the list of regions requested (and any regions preceding them)
     *
     * @param regions            the list of regions to finalize
     * @param knownSnpPositions  the set of known SNP positions; can be null (to allow polyploid consensus anywhere)
     * @param forceCloseFullRegions if true, requires this method to make sure all regions are fully closed; otherwise, we may decide not to close up to the very end (e.g. during het compression)
     * @return a non-null set of reduced reads representing the finalized regions
     */
    public ObjectSet<GATKSAMRecord> closeVariantRegions(final CompressionStash regions, final ObjectSortedSet<GenomeLoc> knownSnpPositions, final boolean forceCloseFullRegions) {
        final ObjectAVLTreeSet<GATKSAMRecord> allReads = new ObjectAVLTreeSet<>(new AlignmentStartWithNoTiesComparator());
        if ( !regions.isEmpty() ) {

            int windowHeaderStart = getStartLocation(windowHeader);
            HeaderElement lastCleanedElement = null;

            for ( final GenomeLoc region : regions ) {
                if (((FinishedGenomeLoc)region).isFinished() && region.getContig().equals(contig) && region.getStart() >= windowHeaderStart && region.getStop() < windowHeaderStart + windowHeader.size()) {
                    final int start = region.getStart() - windowHeaderStart;
                    int stop = region.getStop() - windowHeaderStart;

                    // make sure the bitset is complete given the region (it might not be in multi-sample mode)
                    if ( region.getStop() > markedSites.getStartLocation() + markedSites.getVariantSiteBitSet().length - 1 )
                        markSites(region.getStop());

                    CloseVariantRegionResult closeVariantRegionResult = closeVariantRegion(start, stop, knownSnpPositions);
                    allReads.addAll(closeVariantRegionResult.reads);

                    // check whether we didn't close the whole region that was requested
                    if ( stop > 0 && closeVariantRegionResult.stopPerformed < stop ) {
                        // we should update the variant sites bitset because the context size's worth of bases after the variant position are no longer "variant"
                        markRegionAs(closeVariantRegionResult.stopPerformed + 1, stop, false);

                        // if the calling method said that it didn't care then we are okay so update the stop
                        if ( !forceCloseFullRegions ) {
                            stop = closeVariantRegionResult.stopPerformed;
                        }
                        // otherwise, we need to forcibly push the stop that we originally requested
                        else {
                            while ( closeVariantRegionResult.stopPerformed < stop ) {
                                // first clean up used header elements so they don't get reused
                                for ( int i = 0; i <= closeVariantRegionResult.stopPerformed; i++ )
                                    windowHeader.remove();
                                stop -= (closeVariantRegionResult.stopPerformed + 1);

                                closeVariantRegionResult = closeVariantRegion(0, stop, knownSnpPositions);
                                allReads.addAll(closeVariantRegionResult.reads);
                            }
                        }
                    }

                    // We need to clean up the window header elements up until the end of the requested region so that they don't get used for future regions.
                    // Note that this cleanup used to happen outside the above for-loop, but that was causing an occasional doubling of the reduced reads
                    //  (in the case where there are multiple regions to close we'd reuse the reads for each region).
                    if ( stop >= 0 ) {
                        for ( int i = 0; i < stop; i++ )
                            windowHeader.remove();
                        lastCleanedElement = windowHeader.remove();
                        windowHeaderStart = getStartLocation(windowHeader);
                    }
                }
            }

            // we need to keep the last element of the last cleaned region in the event that the following element has a read that starts with an insertion.
            if ( lastCleanedElement != null && lastCleanedElement.hasInsertionToTheRight() )
                windowHeader.addFirst(new HeaderElement(lastCleanedElement.getLocation(), lastCleanedElement.numInsertionsToTheRight()));
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
    protected ObjectList<GATKSAMRecord> downsampleVariantRegion(final ObjectList<GATKSAMRecord> allReads) {
        int nReads = allReads.size();
        if (nReads == 0)
            return allReads;

        if (downsampleCoverage >= nReads)
            return allReads;

        ReservoirDownsampler <GATKSAMRecord> downsampler = new ReservoirDownsampler<>(downsampleCoverage);
        downsampler.submit(allReads);
        return new ObjectArrayList<>(downsampler.consumeFinalizedItems());
    }


    /**
     * Properly closes a Sliding Window, finalizing all consensus and variant
     * regions that still exist regardless of being able to fulfill the
     * context size requirement in the end.
     *
     * @param knownSnpPositions  the set of known SNP positions; can be null (to allow polyploid consensus anywhere)
     * @return A non-null set/list of all reads generated
     */
    @Ensures("result != null")
    public Pair<ObjectSet<GATKSAMRecord>, CompressionStash> close(final ObjectSortedSet<GenomeLoc> knownSnpPositions) {
        // mark variant regions
        ObjectSet<GATKSAMRecord> finalizedReads = new ObjectAVLTreeSet<>(new AlignmentStartWithNoTiesComparator());
        CompressionStash regions = new CompressionStash();

        if (!windowHeader.isEmpty()) {
            markSites(getStopLocation(windowHeader) + 1);
            regions = findVariantRegions(0, windowHeader.size(), markedSites.getVariantSiteBitSet(), true);
            finalizedReads = closeVariantRegions(regions, knownSnpPositions, true);

            if (!windowHeader.isEmpty())
                finalizedReads.addAll(addAllSyntheticReadTypes(0, windowHeader.size()));
        }

        return new Pair<>(finalizedReads, regions);
    }

    /**
     * generates the SAM record for the running consensus read and resets it (to null)
     *
     * @param runningConsensus the consensus to finalize
     * @return the read contained in the running consensus or null
     */
    protected GATKSAMRecord finalizeRunningConsensus(final SyntheticRead runningConsensus) {
        GATKSAMRecord finalizedRead = null;

        if ( runningConsensus != null ) {
            if ( runningConsensus.size() > 0 )
                finalizedRead = runningConsensus.close();
            else
                consensusCounter--;
        }

        return finalizedRead;
    }

    /**
     * generates the SAM record for the filtered data consensus and resets it (to null)
     *
     * @param filteredDataConsensus the consensus to finalize
     * @return the read contained in the running consensus or null
     */
    protected GATKSAMRecord finalizeFilteredDataConsensus(final SyntheticRead filteredDataConsensus) {
        GATKSAMRecord finalizedRead = null;
        if (filteredDataConsensus != null) {
            if (filteredDataConsensus.size() > 0)
                finalizedRead = filteredDataConsensus.close();
            else
                filteredDataConsensusCounter--;
        }
        return finalizedRead;
    }

    // define this so that we can use Java generics below
    private final static class HeaderElementList extends LinkedList<HeaderElement> {}

    private final static class SingleStrandConsensusData {
        final HeaderElementList consensus = new HeaderElementList();
        final ObjectList<GATKSAMRecord> reads = new ObjectArrayList<>();
    }

    /**
     * Finalizes a variant region - and any adjacent synthetic reads - for point mutations (indel sites are not
     * supported) with polyploid compression.
     *
     * @param hetRefPosition    window header index of the het site; MUST NOT BE AN INDEL SITE!
     * @return a non-null list of all reads contained in the variant region as a polyploid consensus
     */
    @Requires({"start >= 0 && (stop >= start || stop == 0)"})
    @Ensures({"result != null"})
    protected ObjectList<GATKSAMRecord> createPolyploidConsensus(final int hetRefPosition) {
        // we will create two (positive strand, negative strand) headers for each haplotype
        final SingleStrandConsensusData[] headersPosStrand = new SingleStrandConsensusData[2];
        final SingleStrandConsensusData[] headersNegStrand = new SingleStrandConsensusData[2];

        final int globalHetRefPosition = windowHeader.get(hetRefPosition).getLocation();

        // initialize the mapping from base (allele) to header
        final Byte2IntMap alleleHeaderMap = new Byte2IntArrayMap(2);
        alleleHeaderMap.defaultReturnValue(-1);
        for ( final BaseIndex allele : windowHeader.get(hetRefPosition).getAlleles(MIN_ALT_PVALUE_TO_TRIGGER_VARIANT, MIN_ALT_PROPORTION_TO_TRIGGER_VARIANT) ) {
            final int currentIndex = alleleHeaderMap.size();
            if ( currentIndex > 1 )
                throw new IllegalStateException("There are more than 2 alleles present when creating a diploid consensus");

            alleleHeaderMap.put(allele.b, currentIndex);
            headersPosStrand[currentIndex] = new SingleStrandConsensusData();
            headersNegStrand[currentIndex] = new SingleStrandConsensusData();
        }

        // sanity check that we saw 2 alleles
        if ( alleleHeaderMap.size() != 2 )
            throw new IllegalStateException("We expected to see 2 alleles when creating a diploid consensus but saw " + alleleHeaderMap.size());

        final ObjectList<GATKSAMRecord> readsToRemove = new ObjectArrayList<>();

        for ( final GATKSAMRecord read : readsInWindow ) {

            // if the read falls after the het position, just skip it for now (we'll get to it later)
            if ( read.getSoftStart() > globalHetRefPosition )
                continue;

            // remove all other reads from the read cache since we're going to use them here
            readsToRemove.add(read);

            // if the read falls before the het position or has low MQ, we don't need to look at it
            if ( read.getSoftEnd() < globalHetRefPosition || read.getMappingQuality() < MIN_MAPPING_QUALITY)
                continue;

            // remove all spanning reads from the consensus header since we're going to incorporate them into a consensus here instead
            removeFromHeader(windowHeader, read);

            // where on the read is the het position?
            final int readPosOfHet = ReadUtils.getReadCoordinateForReferenceCoordinate(read, globalHetRefPosition, ReadUtils.ClippingTail.LEFT_TAIL);

            // this is safe because indels are not supported
            final byte base = read.getReadBases()[readPosOfHet];

            // check which allele this read represents
            final int allele = alleleHeaderMap.get(base);

            // ignore the read if it represents a base that's not part of the consensus
            if ( allele != -1 ) {
                // add to the appropriate polyploid header
                final SingleStrandConsensusData header = read.getReadNegativeStrandFlag() ? headersNegStrand[allele] : headersPosStrand[allele];
                header.reads.add(read);
                addToHeader(header.consensus, read);
            }
        }

        for ( final GATKSAMRecord read : readsToRemove )
            readsInWindow.remove(read);

        // create the polyploid synthetic reads if we can
        final ObjectList<GATKSAMRecord> hetReads = new ObjectArrayList<>();

        // sanity check that no new "variant region" exists on just a single consensus strand due to softclips
        // or multi-allelic sites now that we've broken everything out into their component parts.  if one does
        // exist then we need to back out the consensus for that strand only.
        for ( final SingleStrandConsensusData header : headersPosStrand ) {
            if ( hasPositionWithSignificantSoftclipsOrVariant(header.consensus, globalHetRefPosition) )
                hetReads.addAll(header.reads);
            else
                finalizeHetConsensus(header.consensus, false, hetReads);
        }
        for ( final SingleStrandConsensusData header : headersNegStrand ) {
            if ( hasPositionWithSignificantSoftclipsOrVariant(header.consensus, globalHetRefPosition) )
                hetReads.addAll(header.reads);
            else
                finalizeHetConsensus(header.consensus, true, hetReads);
        }

        return hetReads;
    }

    /*
     * Finalizes a particular het consensus for the given header representation
     *
     * @param header            the list of header elements representing the header for the consensus
     * @param isNegativeStrand  does this header represent reads on the negative strand?
     * @param result            list in which to store results
     */
    protected void finalizeHetConsensus(final LinkedList<HeaderElement> header, final boolean isNegativeStrand, final ObjectList<GATKSAMRecord> result) {
        if ( header.size() > 0 ) {
            if ( isNegativeStrand )
                result.addAll(addToSyntheticReads(header, 0, header.size(), ConsensusType.NEGATIVE_CONSENSUS));
            else
                result.addAll(addToSyntheticReads(header, 0, header.size(), ConsensusType.POSITIVE_CONSENSUS));
        }
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
     * @param header      the sliding window header to use
     * @param read        the incoming read to be added to the sliding window
     * @param removeRead  if we are removing the read from the header or adding
     */
    protected void updateHeaderCounts(final LinkedList<HeaderElement> header, final GATKSAMRecord read, final boolean removeRead) {
        final int readStart = read.getSoftStart();
        final int headerStart = getStartLocation(header);
        int locationIndex = headerStart < 0 ? 0 : readStart - headerStart;

        if ( removeRead && locationIndex < 0 )
            throw new IllegalStateException("Provided read is behind the Sliding Window! Read = " + read + ", readStart = " + readStart + ", cigar = " + read.getCigarString() + ", window = " + headerStart + "-" + getStopLocation(header));

        // we only need to create new header elements if we are adding the read, not when we're removing it
        if ( !removeRead )
            locationIndex = createNewHeaderElements(header, read, locationIndex);

        actuallyUpdateHeaderForRead(header, read, removeRead, locationIndex);
    }

    /*
     * Creates new header elements if needed for the given read.
     *
     * @param header        the sliding window header to use
     * @param read          the incoming read to be added to the sliding window
     * @param startIndex    the start location index into the header for this read
     *
     * @return an updated index into the modified header
     */
    @Requires("header != null && read != null")
    protected int createNewHeaderElements(final LinkedList<HeaderElement> header, final GATKSAMRecord read, final int startIndex) {

        int headerStart = getStartLocation(header);
        int locationIndex = startIndex;

        // Do we need to add extra elements before the start of the header?  This could happen if the previous read was
        // clipped and this alignment starts before the beginning of the window
        final int readStart = read.getSoftStart();
        if ( startIndex < 0 ) {
            for ( int i = 1; i <= -startIndex; i++ )
                header.addFirst(new HeaderElement(headerStart - i));

            // update the start location accordingly
            headerStart = readStart;
            locationIndex = 0;
        }

        // Do we need to add extra elements to the end of the header?
        final int headerStop = getStopLocation(header);
        final int readEnd = read.getSoftEnd();
        if ( headerStop < readEnd ) {
            final int elementsToAdd = (headerStop < 0) ? readEnd - readStart + 1 : readEnd - headerStop;
            for ( int i = elementsToAdd - 1; i >= 0; i-- )
                header.addLast(new HeaderElement(readEnd - i));
        }

        // Special case for leading insertions before the beginning of the sliding read
        if ( (readStart == headerStart || headerStart < 0) && ReadUtils.readStartsWithInsertion(read.getCigar(), false) != null ) {
            // create a new first element to the window header with no bases added
            header.addFirst(new HeaderElement(readStart - 1));
            // this allows the first element (I) to look at locationIndex - 1 when we update the header and do the right thing
            locationIndex = 1;
        }

        return locationIndex;
    }

    /*
     * Actually updates the sliding window's header counts with the incoming read bases and quals (including insertion and deletion quals).
     *
     * @param header        the sliding window header to use
     * @param read          the incoming read to be added to the sliding window
     * @param removeRead    if we are removing the read from the header or adding
     * @param startIndex    the start location index into the header for this read
     */
    @Requires("header != null && read != null && startIndex >= 0")
    protected void actuallyUpdateHeaderForRead(final LinkedList<HeaderElement> header, final GATKSAMRecord read, final boolean removeRead, final int startIndex) {

        final Iterator<HeaderElement> headerElementIterator = header.listIterator(startIndex);
        final int mappingQuality = read.getMappingQuality();
        final boolean isNegativeStrand = read.getReadNegativeStrandFlag();

        // iterator variables
        int locationIndex = startIndex;
        int readBaseIndex = 0;
        HeaderElement headerElement;

        for ( final CigarElement cigarElement : read.getCigar().getCigarElements() ) {
            switch ( cigarElement.getOperator() ) {
                case H:
                    break;
                case I:
                    readBaseIndex += cigarElement.getLength();

                    // special case, if we don't have the previous header element anymore, don't worry about it.
                    if ( locationIndex == 0 )
                        break;

                    // insertions are added to the base to the left (previous element)
                    headerElement = header.get(locationIndex - 1);

                    if ( removeRead )
                        headerElement.removeInsertionToTheRight();
                    else
                        headerElement.addInsertionToTheRight();

                    break;
                case D:
                    // deletions are added to the baseCounts with the read mapping quality as its quality score
                    final int nDeletionBases = cigarElement.getLength();
                    final byte MQbyte = mappingQuality > Byte.MAX_VALUE ? Byte.MAX_VALUE : (byte)mappingQuality;
                    for ( int i = 0; i < nDeletionBases; i++ ) {
                        headerElement = headerElementIterator.next();
                        if (removeRead)
                            headerElement.removeBase(BaseUtils.Base.D.base, MQbyte, MQbyte, MQbyte, mappingQuality, MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, false, isNegativeStrand);
                        else
                            headerElement.addBase(BaseUtils.Base.D.base, MQbyte, MQbyte, MQbyte, mappingQuality, MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, false, isNegativeStrand);
                    }
                    locationIndex += nDeletionBases;
                    break;
                case S:
                case M:
                case P:
                case EQ:
                case X:
                    final int nBasesToAdd = cigarElement.getLength();
                    final boolean isSoftClip = cigarElement.getOperator() == CigarOperator.S;
                    final byte[] readBases = read.getReadBases();
                    final byte[] readQuals = read.getBaseQualities();
                    final boolean readHasIndelQuals = read.hasBaseIndelQualities();
                    final byte[] insertionQuals = readHasIndelQuals ? read.getBaseInsertionQualities() : null;
                    final byte[] deletionQuals = readHasIndelQuals ? read.getBaseDeletionQualities() : null;

                    for ( int i = 0; i < nBasesToAdd; i++ ) {
                        headerElement = headerElementIterator.next();
                        final byte insertionQuality = readHasIndelQuals ? insertionQuals[readBaseIndex] : -1;
                        final byte deletionQuality = readHasIndelQuals ? deletionQuals[readBaseIndex] : -1;

                        if ( removeRead )
                            headerElement.removeBase(readBases[readBaseIndex], readQuals[readBaseIndex], insertionQuality, deletionQuality, mappingQuality, MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, isSoftClip, isNegativeStrand);
                        else
                            headerElement.addBase(readBases[readBaseIndex], readQuals[readBaseIndex], insertionQuality, deletionQuality, mappingQuality, MIN_BASE_QUAL_TO_COUNT, MIN_MAPPING_QUALITY, isSoftClip, isNegativeStrand);

                        readBaseIndex++;
                    }
                    locationIndex += nBasesToAdd;
                    break;
                default:
                    break;
            }
        }
    }
}

