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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.util.SequenceUtil;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.BySampleSAMFileWriter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * Reduces the BAM file using read based compression that keeps only essential information for variant calling
 *
 * <p>
 * This walker will generated reduced versions of the BAM files that still follow the BAM spec
 * and contain all the information necessary for the GSA variant calling pipeline. Some options
 * allow you to tune in how much compression you want to achieve. The default values have been
 * shown to reduce a typical whole exome BAM file 100x. The higher the coverage, the bigger the
 * savings in file size and performance of the downstream tools.
 *
 * <h2>Input</h2>
 * <p>
 * The BAM file to be compressed
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * The compressed (reduced) BAM file.
 *
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

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.CONTIG)
@ReadFilters({UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class, BadCigarFilter.class})
@Downsample(by=DownsampleType.BY_SAMPLE, toCoverage=40)
public class ReduceReads extends ReadWalker<LinkedList<GATKSAMRecord>, ReduceReadsStash> {

    @Output
    private StingSAMFileWriter out = null;
    private SAMFileWriter writerToUse = null;

    /**
     * The number of bases to keep around mismatches (potential variation)
     */
    @Argument(fullName = "context_size", shortName = "cs", doc = "", required = false)
    private int contextSize = 10;

    /**
     * The minimum mapping quality to be considered for the consensus synthetic read. Reads that have
     * mapping quality below this threshold will not be counted towards consensus, but are still counted
     * towards variable regions.
     */
    @Argument(fullName = "minimum_mapping_quality", shortName = "minmap", doc = "", required = false)
    private int minMappingQuality = 20;

    /**
     * The minimum base quality to be considered for the consensus synthetic read. Reads that have
     * base quality below this threshold will not be counted towards consensus, but are still counted
     * towards variable regions.
     */
    @Argument(fullName = "minimum_base_quality_to_consider", shortName = "minqual", doc = "", required = false)
    private byte minBaseQual = 20;

    /**
     * Reads have notoriously low quality bases on the tails (left and right). Consecutive bases with quality
     * lower than this threshold will be hard clipped off before entering the reduce reads algorithm.
     */
    @Argument(fullName = "minimum_tail_qualities", shortName = "mintail", doc = "", required = false)
    private byte minTailQuality = 2;

    /**
     * Allow the experimental polyploid-based reduction capabilities of this tool
     */
    @Argument(fullName = "allow_polyploid_reduction", shortName = "polyploid", doc = "", required = false)
    private boolean USE_POLYPLOID_REDUCTION = false;

    /**
     * Do not simplify read (strip away all extra information of the read -- anything other than bases, quals
     * and read group).
     */
    @Argument(fullName = "dont_simplify_reads", shortName = "nosimplify", doc = "", required = false)
    private boolean DONT_SIMPLIFY_READS = false;

    /**
     * Do not hard clip adaptor sequences. Note: You don't have to turn this on for reads that are not mate paired.
     * The program will behave correctly in those cases.
     */
    @Argument(fullName = "dont_hardclip_adaptor_sequences", shortName = "noclip_ad", doc = "", required = false)
    private boolean DONT_CLIP_ADAPTOR_SEQUENCES = false;

    /**
     * Do not hard clip the low quality tails of the reads. This option overrides the argument of minimum tail
     * quality.
     */
    @Argument(fullName = "dont_hardclip_low_qual_tails", shortName = "noclip_tail", doc = "", required = false)
    private boolean DONT_CLIP_LOW_QUAL_TAILS = false;

    /**
     * Do not use high quality soft-clipped bases. By default, ReduceReads will hard clip away any low quality soft clipped
     * base left by the aligner and use the high quality soft clipped bases in it's traversal algorithm to identify variant
     * regions. The minimum quality for soft clipped bases is the same as the minimum base quality to consider (minqual)
     */
    @Argument(fullName = "dont_use_softclipped_bases", shortName = "no_soft", doc = "", required = false)
    private boolean DONT_USE_SOFTCLIPPED_BASES = false;

    /**
     * Do not compress read names. By default, ReduceReads will compress read names to numbers and guarantee 
     * uniqueness and reads with similar name will still have similar compressed names. Note: If you scatter/gather
     * there is no guarantee that read name uniqueness will be maintained -- in this case we recommend not compressing. 
     */
    @Argument(fullName = "dont_compress_read_names", shortName = "nocmp_names", doc = "", required = false)
    private boolean DONT_COMPRESS_READ_NAMES = false;

    /**
     * Optionally hard clip all incoming reads to the desired intervals. The hard clips will happen exactly at the interval
     * border.
     */
    @Argument(fullName = "hard_clip_to_interval", shortName = "clip_int", doc = "", required = false)
    private boolean HARD_CLIP_TO_INTERVAL = false;

    /**
     * Minimum proportion of mismatches in a site to trigger a variant region. Anything below this will be
     * considered consensus.
     */
    @Argument(fullName = "minimum_alt_proportion_to_trigger_variant", shortName = "minvar", doc = "", required = false)
    private double minAltProportionToTriggerVariant = 0.05;

    /**
     * Minimum proportion of indels in a site to trigger a variant region. Anything below this will be
     * considered consensus.
     */
    @Argument(fullName = "minimum_del_proportion_to_trigger_variant", shortName = "mindel", doc = "", required = false)
    private double minIndelProportionToTriggerVariant = 0.05;

    /**
     * Downsamples the coverage of a variable region approximately (guarantees the minimum to be equal to this).
     * A value of 0 turns downsampling off.
     */
    @Argument(fullName = "downsample_coverage", shortName = "ds", doc = "", required = false)
    private int downsampleCoverage = 250;

    @Hidden
    @Argument(fullName = "nwayout", shortName = "nw", doc = "", required = false)
    private boolean nwayout = false;

    @Hidden
    @Argument(fullName = "", shortName = "dl", doc = "", required = false)
    private int debugLevel = 0;

    @Hidden
    @Argument(fullName = "", shortName = "dr", doc = "", required = false)
    private String debugRead = "";

    @Hidden
    @Argument(fullName = "downsample_strategy", shortName = "dm", doc = "", required = false)
    private DownsampleStrategy downsampleStrategy = DownsampleStrategy.Normal;
    
    @Hidden 
    @Argument(fullName = "no_pg_tag", shortName = "npt", doc ="", required = false)
    private boolean NO_PG_TAG = false;

    public enum DownsampleStrategy {
        Normal,
        Adaptive
    }
    
    int nCompressedReads = 0;

    HashMap<String, Long> readNameHash;                                     // This hash will keep the name of the original read the new compressed name (a number).
    Long nextReadNumber = 1L;                                               // The next number to use for the compressed read name.

    SortedSet<GenomeLoc> intervalList;

    // IMPORTANT: DO NOT CHANGE THE VALUE OF THIS CONSTANT VARIABLE; IT IS NOW PERMANENTLY THE @PG NAME THAT EXTERNAL TOOLS LOOK FOR IN THE BAM HEADER
    public static final String PROGRAM_RECORD_NAME = "GATK ReduceReads";   // The name that will go in the @PG tag
    private static final String PROGRAM_FILENAME_EXTENSION = ".reduced.bam";

    /**
     * Basic generic initialization of the readNameHash and the intervalList. Output initialization
     * is done at the reduceInit method
     */
    @Override
    public void initialize() {
        super.initialize();
        GenomeAnalysisEngine toolkit = getToolkit();
        readNameHash = new HashMap<String, Long>();           // prepare the read name hash to keep track of what reads have had their read names compressed
        intervalList = new TreeSet<GenomeLoc>();              // get the interval list from the engine. If no interval list was provided, the walker will work in WGS mode

        if (toolkit.getIntervals() != null)
            intervalList.addAll(toolkit.getIntervals());


        final boolean preSorted = true;
        final boolean indexOnTheFly = true;
        final boolean keep_records = true;
        final SAMFileHeader.SortOrder sortOrder = SAMFileHeader.SortOrder.coordinate;
        if (nwayout) {
            SAMProgramRecord programRecord = NO_PG_TAG ? null : Utils.createProgramRecord(toolkit, this, PROGRAM_RECORD_NAME);
            writerToUse = new BySampleSAMFileWriter(toolkit, PROGRAM_FILENAME_EXTENSION, sortOrder, preSorted, indexOnTheFly, NO_PG_TAG, programRecord, true);
        }
        else {
            writerToUse = out;
            out.setPresorted(false);
            if (!NO_PG_TAG) {
                Utils.setupWriter(out, toolkit, toolkit.getSAMFileHeader(), !preSorted, keep_records, this, PROGRAM_RECORD_NAME);
            }
        }
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
    public LinkedList<GATKSAMRecord> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        LinkedList<GATKSAMRecord> mappedReads;
        if (!debugRead.isEmpty() && read.getReadName().contains(debugRead))
                System.out.println("Found debug read!");

        if (debugLevel == 1)
            System.out.printf("\nOriginal: %s %s %d %d\n", read, read.getCigar(), read.getAlignmentStart(), read.getAlignmentEnd());

        // we write the actual alignment starts to their respective alignment shift tags in the temporary
        // attribute hash so we can determine later if we need to write down the alignment shift to the reduced BAM file
        read.setTemporaryAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_START_SHIFT, read.getAlignmentStart());
        read.setTemporaryAttribute(GATKSAMRecord.REDUCED_READ_ORIGINAL_ALIGNMENT_END_SHIFT, read.getAlignmentEnd());

        // Check if the read goes beyond the boundaries of the chromosome, and hard clip those boundaries.
        int chromosomeLength = ref.getGenomeLocParser().getContigInfo(read.getReferenceName()).getSequenceLength();
        if (read.getSoftStart() < 0)
            read = ReadClipper.hardClipByReadCoordinates(read, 0, -read.getSoftStart());
        if (read.getSoftEnd() > chromosomeLength)
            read = ReadClipper.hardClipByReadCoordinates(read, chromosomeLength - read.getSoftStart() + 1, read.getReadLength() - 1);

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
        return new ReduceReadsStash(new MultiSampleCompressor(getToolkit().getSAMFileHeader(), contextSize, downsampleCoverage, minMappingQuality, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, downsampleStrategy, USE_POLYPLOID_REDUCTION));
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

        if (nwayout)
            writerToUse.close();
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
        }

        if (debugLevel == 1)
            System.out.println("BAM: " + read.getCigar() + " " + read.getAlignmentStart() + " " + read.getAlignmentEnd());

        if (!DONT_COMPRESS_READ_NAMES)
            compressReadName(read);

        writerToUse.addAlignment(read);
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
        final Long readNumber = readNameHash.get(name);
        if (readNumber != null) {
            compressedName += readNumber.toString();
        } else {
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
        return isWholeGenome() || list.getFirst().equals(read);
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
