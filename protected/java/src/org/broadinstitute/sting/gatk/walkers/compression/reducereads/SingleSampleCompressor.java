package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Set;
import java.util.TreeSet;

/**
 *
 * @author carneiro, depristo
 * @version 3.0
 */
public class SingleSampleCompressor {
    final private int contextSize;
    final private int downsampleCoverage;
    final private int minMappingQuality;
    final private double minAltProportionToTriggerVariant;
    final private double minIndelProportionToTriggerVariant;
    final private int minBaseQual;
    final private ReduceReads.DownsampleStrategy downsampleStrategy;
    final private int nContigs;
    final private boolean allowPolyploidReduction;

    private SlidingWindow slidingWindow;
    private int slidingWindowCounter;

    public static Pair<Set<GATKSAMRecord>, CompressionStash> emptyPair = new Pair<Set<GATKSAMRecord>,CompressionStash>(new TreeSet<GATKSAMRecord>(), new CompressionStash());

    public SingleSampleCompressor(final int contextSize,
                                  final int downsampleCoverage,
                                  final int minMappingQuality,
                                  final double minAltProportionToTriggerVariant,
                                  final double minIndelProportionToTriggerVariant,
                                  final int minBaseQual,
                                  final ReduceReads.DownsampleStrategy downsampleStrategy,
                                  final int nContigs,
                                  final boolean allowPolyploidReduction) {
        this.contextSize = contextSize;
        this.downsampleCoverage = downsampleCoverage;
        this.minMappingQuality = minMappingQuality;
        this.slidingWindowCounter = 0;
        this.minAltProportionToTriggerVariant = minAltProportionToTriggerVariant;
        this.minIndelProportionToTriggerVariant = minIndelProportionToTriggerVariant;
        this.minBaseQual = minBaseQual;
        this.downsampleStrategy = downsampleStrategy;
        this.nContigs = nContigs;
        this.allowPolyploidReduction = allowPolyploidReduction;
    }

    public Pair<Set<GATKSAMRecord>, CompressionStash> addAlignment( GATKSAMRecord read ) {
        Set<GATKSAMRecord> reads = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        CompressionStash stash = new CompressionStash();
        int readOriginalStart = read.getUnclippedStart();

        // create a new window if:
        if ((slidingWindow != null) &&
            ( ( read.getReferenceIndex() != slidingWindow.getContigIndex() ) ||        // this is a brand new contig
              (readOriginalStart - contextSize > slidingWindow.getStopLocation()))) {  // this read is too far away from the end of the current sliding window

            // close the current sliding window
            Pair<Set<GATKSAMRecord>, CompressionStash> readsAndStash = slidingWindow.close();
            reads = readsAndStash.getFirst();
            stash = readsAndStash.getSecond();
            slidingWindow = null;                                                      // so we create a new one on the next if
        }

        if ( slidingWindow == null) {                                                  // this is the first read
            slidingWindow = new SlidingWindow(read.getReferenceName(), read.getReferenceIndex(), contextSize, read.getHeader(), read.getReadGroup(), slidingWindowCounter, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, minMappingQuality, downsampleCoverage, downsampleStrategy, read.hasBaseIndelQualities(), nContigs, allowPolyploidReduction);
            slidingWindowCounter++;
        }

        stash.addAll(slidingWindow.addRead(read));
        return new Pair<Set<GATKSAMRecord>, CompressionStash>(reads, stash);
    }

    public Pair<Set<GATKSAMRecord>, CompressionStash> close() {
        return (slidingWindow != null) ? slidingWindow.close() : emptyPair;
    }

    public Set<GATKSAMRecord> closeVariantRegions(CompressionStash regions) {
        return slidingWindow == null ? Collections.<GATKSAMRecord>emptySet() : slidingWindow.closeVariantRegions(regions);
    }

}

