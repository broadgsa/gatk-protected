/*
* Copyright (c) 2012 The Broad Institute
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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSingleSampleReadStream;
import org.broadinstitute.sting.utils.sam.ArtificialSingleSampleReadStreamAnalyzer;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SlidingWindowUnitTest extends BaseTest {

    //////////////////////////////////////////////////////////////////////////////////////
    //// This section tests the findVariantRegions() method and related functionality ////
    //////////////////////////////////////////////////////////////////////////////////////

    private static final int variantRegionLength = 1000;
    private static final int globalStartPosition = 1000000;
    private static final SimpleGenomeLoc loc90to95 = new SimpleGenomeLoc("1", 0, 1000090, 1000095, false);
    private static final SimpleGenomeLoc loc96to99 = new SimpleGenomeLoc("1", 0, 1000096, 1000099, false);
    private static final SimpleGenomeLoc loc100to110 = new SimpleGenomeLoc("1", 0, 1000100, 1000110, false);
    private static final SimpleGenomeLoc loc999 = new SimpleGenomeLoc("1", 0, 1000999, 1000999, false);

    private class FindVariantRegionsTest {
        public List<SimpleGenomeLoc> locs, expectedResult;
        public boolean[] variantRegionBitset;

        private FindVariantRegionsTest(final List<SimpleGenomeLoc> locs) {
            this.locs = locs;
            this.expectedResult = locs;
            variantRegionBitset = createBitset(locs);
        }

        private FindVariantRegionsTest(final List<SimpleGenomeLoc> locs, final List<SimpleGenomeLoc> expectedResult) {
            this.locs = locs;
            this.expectedResult = expectedResult;
            variantRegionBitset = createBitset(locs);
        }
    }

    private static boolean[] createBitset(final List<SimpleGenomeLoc> locs) {
        boolean[] variantRegionBitset = new boolean[variantRegionLength];
        for ( SimpleGenomeLoc loc : locs ) {
            final int stop = loc.getStop() - globalStartPosition;
            for ( int i = loc.getStart() - globalStartPosition; i <= stop; i++ )
                variantRegionBitset[i] = true;
        }
        return variantRegionBitset;
    }

    @DataProvider(name = "findVariantRegions")
    public Object[][] createFindVariantRegionsData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<SimpleGenomeLoc>asList(loc90to95))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<SimpleGenomeLoc>asList(loc90to95, loc100to110))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<SimpleGenomeLoc>asList(loc90to95, loc96to99, loc100to110), Arrays.<SimpleGenomeLoc>asList(new SimpleGenomeLoc("1", 0, 1000090, 1000110, false)))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<SimpleGenomeLoc>asList(loc90to95, loc999))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<SimpleGenomeLoc>asList(loc999))});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "findVariantRegions", enabled = true)
    public void testFindVariantRegions(FindVariantRegionsTest test) {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, globalStartPosition);
        final CompressionStash locs = slidingWindow.findVariantRegions(0, variantRegionLength, test.variantRegionBitset, true);
        int index = 0;
        for ( final SimpleGenomeLoc loc : locs ) {
            Assert.assertTrue(loc.equals(test.expectedResult.get(index++)));
        }
    }

    @Test(enabled = true)
    public void testNoClosingRegions() {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, globalStartPosition);
        final CompressionStash locs = slidingWindow.findVariantRegions(0, variantRegionLength, createBitset(Arrays.<SimpleGenomeLoc>asList(loc90to95, loc999)), false);
        Assert.assertEquals(locs.size(), 1);
        Assert.assertEquals(locs.iterator().next(), loc90to95);
    }













    /*

    private static class DownsamplingReadsIteratorTest extends TestDataProvider {
        private DownsamplingReadsIterator downsamplingIter;
        private int targetCoverage;
        private ArtificialSingleSampleReadStream stream;
        private ArtificialSingleSampleReadStreamAnalyzer streamAnalyzer;

        public DownsamplingReadsIteratorTest( ArtificialSingleSampleReadStream stream, int targetCoverage ) {
            super(DownsamplingReadsIteratorTest.class);

            this.stream = stream;
            this.targetCoverage = targetCoverage;

            setName(String.format("%s: targetCoverage=%d numContigs=%d stacksPerContig=%d readsPerStack=%d-%d distanceBetweenStacks=%d-%d readLength=%d-%d unmappedReads=%d",
                    getClass().getSimpleName(),
                    targetCoverage,
                    stream.getNumContigs(),
                    stream.getNumStacksPerContig(),
                    stream.getMinReadsPerStack(),
                    stream.getMaxReadsPerStack(),
                    stream.getMinDistanceBetweenStacks(),
                    stream.getMaxDistanceBetweenStacks(),
                    stream.getMinReadLength(),
                    stream.getMaxReadLength(),
                    stream.getNumUnmappedReads()));
        }

        public void run() {
            streamAnalyzer = new PositionallyDownsampledArtificialSingleSampleReadStreamAnalyzer(stream, targetCoverage);
            downsamplingIter = new DownsamplingReadsIterator(stream.getStingSAMIterator(), new SimplePositionalDownsampler<SAMRecord>(targetCoverage));

            streamAnalyzer.analyze(downsamplingIter);

            // Check whether the observed properties of the downsampled stream are what they should be
            streamAnalyzer.validate();

            // Allow memory used by this test to be reclaimed
            stream = null;
            streamAnalyzer = null;
            downsamplingIter = null;
        }
    }

    @DataProvider(name = "DownsamplingReadsIteratorTestDataProvider")
    public Object[][] createDownsamplingReadsIteratorTests() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(5, 1, 10000);
        String readGroupID = "testReadGroup";
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord(readGroupID);
        readGroup.setSample("testSample");
        header.addReadGroup(readGroup);

        // Values that don't vary across tests
        int targetCoverage = 10;
        int minReadLength = 50;
        int maxReadLength = 100;
        int minDistanceBetweenStacks = 1;
        int maxDistanceBetweenStacks = maxReadLength + 1;

        GenomeAnalysisEngine.resetRandomGenerator();

        // brute force testing!
        for ( int numContigs : Arrays.asList(1, 2, 5) ) {
            for ( int stacksPerContig : Arrays.asList(1, 2, 10) ) {
                for ( int minReadsPerStack : Arrays.asList(1, targetCoverage / 2, targetCoverage, targetCoverage - 1, targetCoverage + 1, targetCoverage * 2) ) {
                    for ( int maxReadsPerStack : Arrays.asList(1, targetCoverage / 2, targetCoverage, targetCoverage - 1, targetCoverage + 1, targetCoverage * 2) ) {
                        for ( int numUnmappedReads : Arrays.asList(0, 1, targetCoverage, targetCoverage * 2) ) {
                            // Only interested in sane read stream configurations here
                            if ( minReadsPerStack <= maxReadsPerStack ) {
                                new DownsamplingReadsIteratorTest(new ArtificialSingleSampleReadStream(header,
                                                                                                       readGroupID,
                                                                                                       numContigs,
                                                                                                       stacksPerContig,
                                                                                                       minReadsPerStack,
                                                                                                       maxReadsPerStack,
                                                                                                       minDistanceBetweenStacks,
                                                                                                       maxDistanceBetweenStacks,
                                                                                                       minReadLength,
                                                                                                       maxReadLength,
                                                                                                       numUnmappedReads),
                                                                  targetCoverage);
                            }
                        }
                    }
                }
            }
        }

        return DownsamplingReadsIteratorTest.getTests(DownsamplingReadsIteratorTest.class);
    }

    @Test(dataProvider = "DownsamplingReadsIteratorTestDataProvider")
    public void runDownsamplingReadsIteratorTest( DownsamplingReadsIteratorTest test ) {
        logger.warn("Running test: " + test);

        GenomeAnalysisEngine.resetRandomGenerator();
        test.run();
    }

    */
}
