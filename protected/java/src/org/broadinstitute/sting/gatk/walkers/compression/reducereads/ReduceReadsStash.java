package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class implements a "read stash" that keeps reads always sorted in alignment order. Useful
 * for read walkers that alter the alignment information of the incoming reads, but need to
 * maintain the reads sorted for the reduce step. (e.g. ReduceReads)
 */

public class ReduceReadsStash {
    protected MultiSampleCompressor compressor;
    SortedSet<GATKSAMRecord> outOfOrderReads;

    /**
     * Creates a stash with the default sorting order (read alignment)
     * @param compressor the MultiSampleCompressor object to be used with this stash (for stash.close())
     */
    public ReduceReadsStash(MultiSampleCompressor compressor) {
        this.compressor = compressor;
        this.outOfOrderReads = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
    }

    /**
     * Get all reads before a given read (for processing)
     *
     * @param read the original read
     * @return all reads that have alignment start before the original read.
     */
    public List<GATKSAMRecord> getAllReadsBefore(GATKSAMRecord read) {
        List<GATKSAMRecord> result = new LinkedList<GATKSAMRecord>();
        GATKSAMRecord newHead = null;

        for (GATKSAMRecord stashedRead : outOfOrderReads) {
            if (ReadUtils.compareSAMRecords(stashedRead, read) <= 0)
                result.add(stashedRead);
            else {
                newHead = stashedRead;
                break;
            }
        }

        if (result.size()  > 0) {
            if (result.size() == outOfOrderReads.size())
                outOfOrderReads.clear();
            else
                outOfOrderReads = new TreeSet<GATKSAMRecord>(outOfOrderReads.tailSet(newHead));
        }

        return result;
    }

    /**
     * sends the read to the MultiSampleCompressor
     *
     * @param read the read to be compressed
     * @return any compressed reads that may have resulted from adding this read to the machinery (due to the sliding window)
     */
    public Iterable<GATKSAMRecord> compress(GATKSAMRecord read) {
        return compressor.addAlignment(read);
    }

    /**
     * Add a read to the stash
     *
     * @param read any read
     */
    public void add(GATKSAMRecord read) {
        outOfOrderReads.add(read);
    }

    /**
     * Close the stash, processing all remaining reads in order
     *
     * @return a list of all the reads produced by the SlidingWindow machinery)
     */
    public Iterable<GATKSAMRecord> close() {
        LinkedList<GATKSAMRecord> result = new LinkedList<GATKSAMRecord>();

        // compress all the stashed reads (in order)
        for (GATKSAMRecord read : outOfOrderReads)
            for (GATKSAMRecord compressedRead : compressor.addAlignment(read))
                result.add(compressedRead);

        // output any remaining reads from the compressor
        for (GATKSAMRecord read : compressor.close())
            result.add(read);

        return result;
    }

    /**
     * Useful debug functionality, outputs all elements in the stash
     */
    public void print() {
        int i = 1;
        System.out.println("Stash Contents:");
        for (GATKSAMRecord read : outOfOrderReads) 
            System.out.println(String.format("%3d: %s %d %d", i++, read.getCigarString(), read.getAlignmentStart(), read.getAlignmentEnd()));
        System.out.println();
    }

}