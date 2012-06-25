package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import net.sf.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentStartWithNoTiesComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

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
 *
 * @author depristo
 */
public class MultiSampleCompressor implements Compressor {
    protected static final Logger logger = Logger.getLogger(MultiSampleCompressor.class);

    protected Map<String, SingleSampleCompressor> compressorsPerSample = new HashMap<String, SingleSampleCompressor>();

    public MultiSampleCompressor(SAMFileHeader header,
                                 final int contextSize,
                                 final int downsampleCoverage,
                                 final int minMappingQuality,
                                 final double minAltProportionToTriggerVariant,
                                 final double minIndelProportionToTriggerVariant,
                                 final int minBaseQual,
                                 final ReduceReadsWalker.DownsampleStrategy downsampleStrategy) {
        for ( String name : SampleUtils.getSAMFileSamples(header) ) {
            compressorsPerSample.put(name,
                    new SingleSampleCompressor(name, contextSize, downsampleCoverage,
                                    minMappingQuality, minAltProportionToTriggerVariant, minIndelProportionToTriggerVariant, minBaseQual, downsampleStrategy));
        }
    }

    @Override
    public Iterable<GATKSAMRecord> addAlignment(GATKSAMRecord read) {
        String sample = read.getReadGroup().getSample();
        SingleSampleCompressor compressor = compressorsPerSample.get(sample);
        if ( compressor == null )
            throw new ReviewedStingException("No compressor for sample " + sample);
        return compressor.addAlignment(read);
    }

    @Override
    public Iterable<GATKSAMRecord> close() {
        SortedSet<GATKSAMRecord> reads = new TreeSet<GATKSAMRecord>(new AlignmentStartWithNoTiesComparator());
        for ( SingleSampleCompressor comp : compressorsPerSample.values() )
            for ( GATKSAMRecord read : comp.close() )
                reads.add(read);
        return reads;
    }
}
