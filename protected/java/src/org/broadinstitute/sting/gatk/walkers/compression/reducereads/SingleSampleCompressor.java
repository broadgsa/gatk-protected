package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.TreeSet;

/**
 *
 * @author depristo
 * @version 0.1
 */
public class SingleSampleCompressor implements Compressor {
    protected static final Logger logger = Logger.getLogger(SingleSampleCompressor.class);

    protected final int contextSize;
    protected final int downsampleCoverage;
    protected int minMappingQuality;
    protected int slidingWindowCounter;

    protected final String sampleName;

    protected SlidingWindow slidingWindow;
    protected double minAltProportionToTriggerVariant;
    protected double minIndelProportionToTriggerVariant;
    protected int minBaseQual;

    protected ReduceReadsWalker.DownsampleStrategy downsampleStrategy;

    public SingleSampleCompressor(final String sampleName,
                                  final int contextSize,
                                  final int downsampleCoverage,
                                  final int minMappingQuality,
                                  final double minAltProportionToTriggerVariant,
                                  final double minIndelProportionToTriggerVariant,
                                  final int minBaseQual,
                                  final ReduceReadsWalker.DownsampleStrategy downsampleStrategy) {
        this.sampleName = sampleName;
        this.contextSize = contextSize;
        this.downsampleCoverage = downsampleCoverage;
        this.minMappingQuality = minMappingQuality;
        this.slidingWindowCounter = 0;
        this.minAltProportionToTriggerVariant = minAltProportionToTriggerVariant;
        this.minIndelProportionToTriggerVariant = minIndelProportionToTriggerVariant;
        this.minBaseQual = minBaseQual;
        this.downsampleStrategy = downsampleStrategy;
    }

    /**
     * @{inheritDoc}
     */
    @Override
    public Iterable<GATKSAMRecord> addAlignment( GATKSAMRecord read ) {
        TreeSet<GATKSAMRecord> result = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        int readOriginalStart = read.getUnclippedStart();

        // create a new window if:
        if ((slidingWindow != null) &&
            ( ( read.getReferenceIndex() != slidingWindow.getContigIndex() ) ||        // this is a brand new contig
              (readOriginalStart - contextSize > slidingWindow.getStopLocation()))) {  // this read is too far away from the end of the current sliding window

            // close the current sliding window
            result.addAll(slidingWindow.close());
            slidingWindow = null;                                                      // so we create a new one on the next if
        }

        if ( slidingWindow == null) {                                                  // this is the first read
            slidingWindow = new SlidingWindow(read.getReferenceName(), read.getReferenceIndex(), contextSize, read.getHeader(), read.getReadGroup(), slidingWindowCounter, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, minMappingQuality, downsampleCoverage, downsampleStrategy, read.hasBaseIndelQualities());
            slidingWindowCounter++;
        }

        result.addAll(slidingWindow.addRead(read));
        return result;
    }

    @Override
    public Iterable<GATKSAMRecord> close() {
        return (slidingWindow != null) ? slidingWindow.close() : new TreeSet<GATKSAMRecord>();
    }

}

