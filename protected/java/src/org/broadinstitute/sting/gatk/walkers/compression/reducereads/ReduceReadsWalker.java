/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.util.SequenceUtil;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocComparator;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * Reduces the BAM file using read based compression that keeps only essential information for variant calling
 * <p/>
 * <p>
 * This walker will generated reduced versions of the BAM files that still follow the BAM spec
 * and contain all the information necessary for the GSA variant calling pipeline. Some options
 * allow you to tune in how much compression you want to achieve. The default values have been
 * shown to reduce a typical whole exome BAM file 100x. The higher the coverage, the bigger the
 * savings in file size and performance of the downstream tools.
 * <p/>
 * <h2>Input</h2>
 * <p>
 * The BAM file to be compressed
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * The compressed (reduced) BAM file.
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ReduceReads \
 *   -I myData.bam \
 *   -o myData.reduced.bam
 * </pre>
 */

@PartitionBy(PartitionType.INTERVAL)
@ReadFilters({UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class, BadCigarFilter.class})
public class ReduceReadsWalker extends ReadWalker<LinkedList<GATKSAMRecord>, ReduceReadsStash> {

    @Output
    protected StingSAMFileWriter out;

    /**
     * The number of bases to keep around mismatches (potential variation)
     */
    @Argument(fullName = "context_size", shortName = "cs", doc = "", required = false)
    protected int contextSize = 10;

    /**
     * The minimum mapping quality to be considered for the consensus synthetic read. Reads that have
     * mapping quality below this threshold will not be counted towards consensus, but are still counted
     * towards variable regions.
     */
    @Argument(fullName = "minimum_mapping_quality", shortName = "minmap", doc = "", required = false)
    protected int minMappingQuality = 20;

    /**
     * The minimum base quality to be considered for the consensus synthetic read. Reads that have
     * base quality below this threshold will not be counted towards consensus, but are still counted
     * towards variable regions.
     */
    @Argument(fullName = "minimum_base_quality_to_consider", shortName = "minqual", doc = "", required = false)
    protected byte minBaseQual = 20;

    /**
     * Reads have notoriously low quality bases on the tails (left and right). Consecutive bases with quality
     * lower than this threshold will be hard clipped off before entering the reduce reads algorithm.
     */
    @Argument(fullName = "minimum_tail_qualities", shortName = "mintail", doc = "", required = false)
    protected byte minTailQuality = 2;

    /**
     * Do not simplify read (strip away all extra information of the read -- anything other than bases, quals
     * and read group).
     */
    @Argument(fullName = "dont_simplify_reads", shortName = "nosimplify", doc = "", required = false)
    protected boolean DONT_SIMPLIFY_READS = false;

    /**
     * Do not hard clip adaptor sequences. Note: You don't have to turn this on for reads that are not mate paired.
     * The program will behave correctly in those cases.
     */
    @Argument(fullName = "dont_hardclip_adaptor_sequences", shortName = "noclip_ad", doc = "", required = false)
    protected boolean DONT_CLIP_ADAPTOR_SEQUENCES = false;

    /**
     * Do not hard clip the low quality tails of the reads. This option overrides the argument of minimum tail
     * quality.
     */
    @Argument(fullName = "dont_hardclip_low_qual_tails", shortName = "noclip_tail", doc = "", required = false)
    protected boolean DONT_CLIP_LOW_QUAL_TAILS = false;

    /**
     * Do not use high quality soft-clipped bases. By default, ReduceReads will hard clip away any low quality soft clipped
     * base left by the aligner and use the high quality soft clipped bases in it's traversal algorithm to identify variant
     * regions. The minimum quality for soft clipped bases is the same as the minimum base quality to consider (minqual)
     */
    @Argument(fullName = "dont_use_softclipped_bases", shortName = "no_soft", doc = "", required = false)
    protected boolean DONT_USE_SOFTCLIPPED_BASES = false;

    /**
     * Do not compress read names. By default, ReduceReads will compress read names to numbers and guarantee 
     * uniqueness and reads with similar name will still have similar compressed names. Note: If you scatter/gather
     * there is no guarantee that read name uniqueness will be maintained -- in this case we recommend not compressing. 
     */
    @Argument(fullName = "dont_compress_read_names", shortName = "nocmp_names", doc = "", required = false)
    protected boolean DONT_COMPRESS_READ_NAMES = false;

    /**
     * Optionally hard clip all incoming reads to the desired intervals. The hard clips will happen exactly at the interval
     * border.
     */
    @Argument(fullName = "hard_clip_to_interval", shortName = "clip_int", doc = "", required = false)
    protected boolean HARD_CLIP_TO_INTERVAL = false;

    /**
     * Minimum proportion of mismatches in a site to trigger a variant region. Anything below this will be
     * considered consensus.
     */
    @Argument(fullName = "minimum_alt_proportion_to_trigger_variant", shortName = "minvar", doc = "", required = false)
    protected double minAltProportionToTriggerVariant = 0.05;

    /**
     * Minimum proportion of indels in a site to trigger a variant region. Anything below this will be
     * considered consensus.
     */
    @Argument(fullName = "minimum_del_proportion_to_trigger_variant", shortName = "mindel", doc = "", required = false)
    protected double minIndelProportionToTriggerVariant = 0.05;

    /**
     * Downsamples the coverage of a variable region approximately (guarantees the minimum to be equal to this).
     * A value of 0 turns downsampling off.
     */
    @Argument(fullName = "downsample_coverage", shortName = "ds", doc = "", required = false)
    protected int downsampleCoverage = 0;

    @Hidden
    @Argument(fullName = "", shortName = "dl", doc = "", required = false)
    protected int debugLevel = 0;

    @Hidden
    @Argument(fullName = "", shortName = "dr", doc = "", required = false)
    protected String debugRead = "";

    @Hidden
    @Argument(fullName = "downsample_strategy", shortName = "dm", doc = "", required = false)
    protected DownsampleStrategy downsampleStrategy = DownsampleStrategy.Normal;
    
    @Hidden 
    @Argument(fullName = "no_pg_tag", shortName = "npt", doc ="", required = false)
    private boolean NO_PG_TAG = false;

    public enum DownsampleStrategy {
        Normal,
        Adaptive
    }
    
    protected int totalReads = 0;
    int nCompressedReads = 0;

    HashMap<String, Long> readNameHash;                                     // This hash will keep the name of the original read the new compressed name (a number).
    Long nextReadNumber = 1L;                                               // The next number to use for the compressed read name.

    SortedSet<GenomeLoc> intervalList;
    
    private static final String PROGRAM_RECORD_NAME = "GATK ReduceReads";   // The name that will go in the @PG tag

    /**
     * Basic generic initialization of the readNameHash and the intervalList. Output initialization
     * is done at the reduceInit method
     */
    @Override
    public void initialize() {
        super.initialize();
        GenomeAnalysisEngine toolkit = getToolkit();
        readNameHash = new HashMap<String, Long>();                         // prepare the read name hash to keep track of what reads have had their read names compressed
        intervalList = new TreeSet<GenomeLoc>(new GenomeLocComparator());   // get the interval list from the engine. If no interval list was provided, the walker will work in WGS mode

        if (toolkit.getIntervals() != null)
            intervalList.addAll(toolkit.getIntervals());

        if (!NO_PG_TAG)
            Utils.setupWriter(out, toolkit, false, true, this, PROGRAM_RECORD_NAME);
        else
            out.setPresorted(false);
    }

    /**
     * Takes in a read and prepares it for the SlidingWindow machinery by performing the
     * following optional clipping operations:
     * 1. Hard clip adaptor sequences
     * 2. Hard clip low quality tails
     * 3. Hard clip all remaining soft clipped bases
     * 4. Hard clip read to the intervals in the interval list (this step may produce multiple reads)
     *
     * @param ref             default map parameter
     * @param read            default map parameter
     * @param metaDataTracker default map parameter
     * @return a linked list with all the reads produced by the clipping operations
     */
    @Override
    public LinkedList<GATKSAMRecord> map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        LinkedList<GATKSAMRecord> mappedReads;
        totalReads++;
        if (!debugRead.isEmpty() && read.getReadName().contains(debugRead))
            System.out.println("Found debug read!");

        if (debugLevel == 1)
            System.out.printf("\nOriginal: %s %s %d %d\n", read, read.getCigar(), read.getAlignmentStart(), read.getAlignmentEnd());

        // we write the actual alignment starts to their respectiv alignment shift tags in the temporary
        // attribute hash so we can determine later if we need to write down the alignment shift to the reduced BAM file
        read.setTemporaryAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT, read.getAlignmentStart());
        read.setTemporaryAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_END_SHIFT, read.getAlignmentEnd());
        
        if (!DONT_SIMPLIFY_READS)
            read.simplify();                                                                                            // Clear all unnecessary attributes
        if (!DONT_CLIP_ADAPTOR_SEQUENCES)
            read = ReadClipper.hardClipAdaptorSequence(read);                                                           // Strip away adaptor sequences, if any.
        if (!DONT_CLIP_LOW_QUAL_TAILS)
            read = ReadClipper.hardClipLowQualEnds(read, minTailQuality);                                               // Clip low quality tails
        if (!isWholeGenome()) {
            if (HARD_CLIP_TO_INTERVAL)
                mappedReads = hardClipReadToInterval(read);                                                             // Hard clip the remainder of the read to the desired interval
            else {
                mappedReads = new LinkedList<GATKSAMRecord>();
                mappedReads.add(read);
            }
        }
        else {
            mappedReads = new LinkedList<GATKSAMRecord>();
            if (!read.isEmpty())
                mappedReads.add(read);
        }

        if (!mappedReads.isEmpty() && !DONT_USE_SOFTCLIPPED_BASES) {
            LinkedList<GATKSAMRecord> tempList = new LinkedList<GATKSAMRecord>();
            for (GATKSAMRecord mRead : mappedReads) {
                GATKSAMRecord clippedRead = ReadClipper.hardClipLowQualitySoftClips(mRead, minBaseQual);
                if (!clippedRead.isEmpty())
                    tempList.add(clippedRead);
            }
            mappedReads = tempList;
        }

        if (debugLevel == 1)
            for (GATKSAMRecord mappedRead : mappedReads)
                System.out.printf("MAPPED: %s %d %d\n", mappedRead.getCigar(), mappedRead.getAlignmentStart(), mappedRead.getAlignmentEnd());

        return mappedReads;

    }

    /**
     * Initializes the ReduceReadsStash that keeps track of all reads that are waiting to
     * enter the SlidingWindow machinery. The stash makes sure reads are served in order
     * even though map() may generate reads that are only supposed to enter the machinery
     * in the future.
     *
     * @return the empty stash
     */
    @Override
    public ReduceReadsStash reduceInit() {
        return new ReduceReadsStash(new MultiSampleCompressor(getToolkit().getSAMFileHeader(), contextSize, downsampleCoverage, minMappingQuality, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, downsampleStrategy));
    }

    /**
     * Takes the list of reads produced by map(), adds them to the stash (which keeps them sorted) and process
     * all reads that come before the original read (the read that was passed to map) including the original
     * read. This is where we send reads, in order, to the SlidingWindow machinery.
     *
     * @param mappedReads the list of reads sent by map
     * @param stash       the stash that keeps the reads in order for processing
     * @return the stash with all reads that have not been processed yet
     */
    public ReduceReadsStash reduce(LinkedList<GATKSAMRecord> mappedReads, ReduceReadsStash stash) {
        if (debugLevel == 1)
            stash.print();

        boolean firstRead = true;
        for (GATKSAMRecord read : mappedReads) {
            boolean originalRead = firstRead && isOriginalRead(mappedReads, read);

            if (read.getReadLength() == 0)
                throw new ReviewedStingException("Empty read sent to reduce, this should never happen! " + read.getReadName() + " -- " + read.getCigar() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd());

            if (originalRead) {
                List<GATKSAMRecord> readsReady = new LinkedList<GATKSAMRecord>();
                readsReady.addAll(stash.getAllReadsBefore(read));
                readsReady.add(read);

                for (GATKSAMRecord readReady : readsReady) {
                    if (debugLevel == 1)
                        System.out.println("REDUCE: " + readReady.getCigar() + " " + readReady.getAlignmentStart() + " " + readReady.getAlignmentEnd());

                    for (GATKSAMRecord compressedRead : stash.compress(readReady))
                        outputRead(compressedRead);

                }
            } else
                stash.add(read);

            firstRead = false;
        }

        return stash;
    }

    /**
     * Now that now more reads will come, we process all the remaining reads in the stash, in order.
     *
     * @param stash the ReduceReadsStash with all unprocessed reads (from reduce)
     */
    @Override
    public void onTraversalDone(ReduceReadsStash stash) {

        // output any remaining reads in the compressor
        for (GATKSAMRecord read : stash.close())
            outputRead(read);
    }

    /**
     * Hard clips away all parts of the read that doesn't agree with the intervals selected.
     *
     * Note: If read overlaps more than one interval, it will be hard clipped to all
     * the intervals it overlaps with
     *
     * @param read the read to be hard clipped to the interval.
     * @return a shallow copy of the read hard clipped to the interval
     */
    private LinkedList<GATKSAMRecord> hardClipReadToInterval(GATKSAMRecord read) {
        LinkedList<GATKSAMRecord> clippedReads = new LinkedList<GATKSAMRecord>();

        GenomeLoc intervalOverlapped = null;       // marks the interval to which the original read overlapped (so we can cut all previous intervals from the list)

        boolean originalRead = true;               // false if this is the right tail of the original read
        boolean overlap;                           // keeps track of the interval that overlapped the original read
        boolean doneClipping;                      // triggers an early exit if we are done clipping this read

        if (isWholeGenome())
            clippedReads.add(read);                // if we don't have intervals (wgs) the read goes in unchanged

        for (GenomeLoc interval : intervalList) {

            if (read.isEmpty())                    // nothing to do with an empty read (could have been fully clipped before)
                break;

            GATKSAMRecord clippedRead = null;      // this will hold the read clipped to the interval to be added in the end of the switch

            switch (ReadUtils.getReadAndIntervalOverlapType(read, interval)) {
                case NO_OVERLAP_RIGHT:             // no reads on this interval, check the next interval if this is the original read
                    if (!originalRead)             // something went wrong if this is the tail of the read
                        throw new ReviewedStingException("tail of the read should never NO_OVERLAP_RIGHT the following interval. " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());
                    overlap = false;
                    doneClipping = false;
                    break;


                case NO_OVERLAP_HARDCLIPPED_RIGHT: // read used to overlap but got hard clipped and doesn't overlap anymore
                    if (originalRead) {
                        overlap = true;            // effectively, we have found the read's location and now we are going to try and match it's tail (which happens to be the entire read).
                        clippedRead = GATKSAMRecord.emptyRead(read);
                    } else
                        overlap = false;

                    doneClipping = false;
                    break;

                case NO_OVERLAP_CONTIG:            // read is in a different contig
                    if (originalRead) {                                                                                 // the original read can be in a bigger contig, but not on a smaller one.
                        if (read.getReferenceIndex() < interval.getContigIndex())
                            throw new ReviewedStingException("read is behind interval list. (contig) " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());
                        else {
                            overlap = false;
                            doneClipping = false;
                        }
                    }                                                                                                   // tail read CANNOT be in a different contig.
                    else {
                        if (read.getReferenceIndex() < interval.getContigIndex()) {
                            overlap = false;
                            doneClipping = true;
                        } else
                            throw new ReviewedStingException("Tail read is in bigger contig than interval traversal. " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());

                    }
                    break;

                case NO_OVERLAP_LEFT:
                    if (originalRead)                                                                                   // if this is the first read this should never happen.
                        throw new ReviewedStingException("original read cannot be behind the first interval. (position) " + read.getReadName() + " -- " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd() + " x " + interval.getLocation().toString());

                    overlap = false;
                    doneClipping = true;
                    break;

                case NO_OVERLAP_HARDCLIPPED_LEFT:                                                                       // read used to overlap but got hard clipped and doesn't overlap anymore
                    overlap = originalRead;                                                                             // if this is the original read, we should not advance the interval list, the original overlap was here.
                    doneClipping = true;
                    break;

                case OVERLAP_LEFT:                                                                                      // clip the left tail of the read
                    clippedRead = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, interval.getStart() - 1);

                    overlap = true;
                    doneClipping = true;
                    break;

                case OVERLAP_RIGHT:                                                                                     // clip the right tail of the read and try to match it to the next interval
                    clippedRead = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, interval.getStop() + 1);
                    read = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, interval.getStop());

                    overlap = true;
                    doneClipping = false;
                    break;

                case OVERLAP_LEFT_AND_RIGHT:                                                                            // clip both left and right ends of the read
                    clippedRead = ReadClipper.hardClipBothEndsByReferenceCoordinates(read, interval.getStart() - 1, interval.getStop() + 1);
                    read = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, interval.getStop());

                    overlap = true;
                    doneClipping = false;
                    break;

                case OVERLAP_CONTAINED:                                                                                 // don't do anything to the read
                    clippedRead = read;

                    overlap = true;
                    doneClipping = true;
                    break;

                default:
                    throw new ReviewedStingException("interval overlap returned an unknown / unhandled state. If new state was added to intervalOverlap, it should be handled by hardClipReadToInterval.");
            }

            if (overlap && originalRead)
                intervalOverlapped = interval;

            if (clippedRead != null) {
                originalRead = false;

                if (!clippedRead.isEmpty())
                    clippedReads.add(clippedRead); // if the read overlaps the interval entirely within a deletion, it will be entirely clipped off
            }

            if (doneClipping)
                break;
        }

        if (intervalOverlapped != null)
            intervalList = intervalList.tailSet(intervalOverlapped);

        return clippedReads;
    }

    /**
     * Compresses the read name and adds it to output BAM file (reduced BAM)
     * after performing some quality control
     *
     * @param read any read
     */
    private void outputRead(GATKSAMRecord read) {
        if (debugLevel == 2) {
            checkForHighMismatch(read);
            checkCigar(read);
        }

        if (read.isReducedRead())
            nCompressedReads++;
        else {
            int originalAlignmentStart =  (Integer) read.getTemporaryAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT);
            int originalAlignmentEnd   =  (Integer) read.getTemporaryAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_END_SHIFT);

            int startShift = originalAlignmentStart - read.getUnclippedStart();                                         // we annotate the shifts for better compression
            int endShift = read.getUnclippedEnd() - originalAlignmentEnd;                                               // we annotate the shifts for better compression

            if (startShift > 0)
                read.setAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT, startShift);               // If the read had any soft clips before getting chopped (variant region) annotate it's original alignment (start)
            if (endShift > 0)
                read.setAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_END_SHIFT, endShift);                   // If the read had any soft clips before getting chopped (variant region) annotate it's original alignment (end)
            
            totalReads++;
        }

        if (debugLevel == 1)
            System.out.println("BAM: " + read.getCigar() + " " + read.getAlignmentStart() + " " + read.getAlignmentEnd());

//        if (!DONT_USE_SOFTCLIPPED_BASES)
//            reSoftClipBases(read);

        if (!DONT_COMPRESS_READ_NAMES)
            compressReadName(read);

        out.addAlignment(read);
    }

    private void reSoftClipBases(GATKSAMRecord read) {
        Integer left = (Integer) read.getTemporaryAttribute("SL");
        Integer right = (Integer) read.getTemporaryAttribute("SR");
        if (left != null || right != null) {
            Cigar newCigar = new Cigar();
            for (CigarElement element : read.getCigar().getCigarElements()) {
                newCigar.add(new CigarElement(element.getLength(), element.getOperator()));
            }

            if (left != null) {
                newCigar = updateFirstSoftClipCigarElement(left, newCigar);
                read.setAlignmentStart(read.getAlignmentStart() + left);
            }

            if (right != null) {
                Cigar invertedCigar = invertCigar(newCigar);
                newCigar = invertCigar(updateFirstSoftClipCigarElement(right, invertedCigar));
            }
            read.setCigar(newCigar);
        }
    }

    /**
     * Facility routine to revert the first element of a Cigar string (skipping hard clips) into a soft-clip.
     * To be used on both ends if provided a flipped Cigar
     *
     * @param softClipSize  the length of the soft clipped element to add
     * @param originalCigar the original Cigar string
     * @return a new Cigar object with the soft clips added
     */
    private Cigar updateFirstSoftClipCigarElement (int softClipSize, Cigar originalCigar) {
        Cigar result = new Cigar();
        CigarElement leftElement = new CigarElement(softClipSize, CigarOperator.S);
        boolean updated = false;
        for (CigarElement element : originalCigar.getCigarElements()) {
            if (!updated && element.getOperator() == CigarOperator.M) {
                result.add(leftElement);
                int newLength = element.getLength() - softClipSize;
                if (newLength > 0)
                    result.add(new CigarElement(newLength, CigarOperator.M));
                updated = true;
            }
            else
                result.add(element);
        }
        return result;
    }

    /**
     * Given a cigar string, returns the inverted cigar string.
     *
     * @param cigar the original cigar
     * @return the inverted cigar
     */
    private Cigar invertCigar(Cigar cigar) {
        Stack<CigarElement> stack = new Stack<CigarElement>();
        for (CigarElement e : cigar.getCigarElements())
            stack.push(e);
        Cigar inverted = new Cigar();
        while (!stack.empty()) {
            inverted.add(stack.pop());
        }
        return inverted;
    }


    /**
     * Quality control procedure that checks if the consensus reads contains too many
     * mismatches with the reference. This should never happen and is a good trigger for
     * errors with the algorithm.
     *
     * @param read any read
     */
    private void checkForHighMismatch(GATKSAMRecord read) {
        final int start = read.getAlignmentStart();
        final int stop = read.getAlignmentEnd();
        final byte[] ref = getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(read.getReferenceName(), start, stop).getBases();
        final int nm = SequenceUtil.countMismatches(read, ref, start - 1);
        final int readLen = read.getReadLength();
        final double nmFraction = nm / (1.0 * readLen);
        if (nmFraction > 0.4 && readLen > 20 && read.getAttribute(GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG) != null && read.getReadName().startsWith("Consensus"))
            throw new ReviewedStingException("BUG: High mismatch fraction found in read " + read.getReadName() + " position: " + read.getReferenceName() + ":" + read.getAlignmentStart() + "-" + read.getAlignmentEnd());
    }

    private void checkCigar (GATKSAMRecord read) {
        if (read.getCigar().isValid(null, -1) != null) {
            throw new ReviewedStingException("BUG: cigar string is not valid: " + read.getCigarString());
        }

    }


    /**
     * Compresses the read name using the readNameHash if we have already compressed
     * this read name before.
     *
     * @param read any read
     */
    private void compressReadName(GATKSAMRecord read) {
        String name = read.getReadName();
        String compressedName = read.isReducedRead() ? "C" : "";
        if (readNameHash.containsKey(name))
            compressedName += readNameHash.get(name).toString();
        else {
            readNameHash.put(name, nextReadNumber);
            compressedName += nextReadNumber.toString();
            nextReadNumber++;
        }

        read.setReadName(compressedName);
    }

    /**
     * Returns true if the read is the original read that went through map().
     *
     * This is important to know so we can decide what reads to pull from the stash. Only reads that came before the original read should be pulled.
     *
     * @param list the list
     * @param read the read
     * @return Returns true if the read is the original read that went through map().
     */
    private boolean isOriginalRead(LinkedList<GATKSAMRecord> list, GATKSAMRecord read) {
        return isWholeGenome() || (list.getFirst().equals(read) && ReadUtils.getReadAndIntervalOverlapType(read, intervalList.first()) == ReadUtils.ReadAndIntervalOverlap.OVERLAP_CONTAINED);
    }

    /**
     * Checks whether or not the intervalList is empty, meaning we're running in WGS mode.
     *
     * @return whether or not we're running in WGS mode.
     */
    private boolean isWholeGenome() {
        return intervalList.isEmpty();
    }

}
