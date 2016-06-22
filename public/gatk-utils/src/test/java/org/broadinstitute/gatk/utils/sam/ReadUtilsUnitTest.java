/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


public class ReadUtilsUnitTest extends BaseTest {
    private interface GetAdaptorFunc {
        public int getAdaptor(final GATKSAMRecord record);
    }

    @DataProvider(name = "AdaptorGetter")
    public Object[][] makeActiveRegionCutTests() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        tests.add( new Object[]{ new GetAdaptorFunc() {
            @Override public int getAdaptor(final GATKSAMRecord record) { return ReadUtils.getAdaptorBoundary(record); }
        }});

        tests.add( new Object[]{ new GetAdaptorFunc() {
            @Override public int getAdaptor(final GATKSAMRecord record) { return record.getAdaptorBoundary(); }
        }});

        return tests.toArray(new Object[][]{});
    }

    private GATKSAMRecord makeRead(final int fragmentSize, final int mateStart) {
        final byte[] bases = {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
        final byte[] quals = {30, 30, 30, 30, 30, 30, 30, 30};
        final String cigar = "8M";
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, cigar);
        read.setProperPairFlag(true);
        read.setReadPairedFlag(true);
        read.setMateAlignmentStart(mateStart);
        read.setInferredInsertSize(fragmentSize);
        return read;
    }

    @Test(dataProvider = "AdaptorGetter")
    public void testGetAdaptorBoundary(final GetAdaptorFunc get) {
        final int fragmentSize = 10;
        final int mateStart = 1000;
        final int BEFORE = mateStart - 2;
        final int AFTER = mateStart + 2;
        int myStart, boundary;
        GATKSAMRecord read;

        // Test case 1: positive strand, first read
        read = makeRead(fragmentSize, mateStart);
        myStart = BEFORE;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(false);
        read.setMateNegativeStrandFlag(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, myStart + fragmentSize + 1);

        // Test case 2: positive strand, second read
        read = makeRead(fragmentSize, mateStart);
        myStart = AFTER;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(false);
        read.setMateNegativeStrandFlag(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, myStart + fragmentSize + 1);

        // Test case 3: negative strand, second read
        read = makeRead(fragmentSize, mateStart);
        myStart = AFTER;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(true);
        read.setMateNegativeStrandFlag(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, mateStart - 1);

        // Test case 4: negative strand, first read
        read = makeRead(fragmentSize, mateStart);
        myStart = BEFORE;
        read.setAlignmentStart(myStart);
        read.setReadNegativeStrandFlag(true);
        read.setMateNegativeStrandFlag(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, mateStart - 1);

        // Test case 5: mate is mapped to another chromosome (test both strands)
        read = makeRead(fragmentSize, mateStart);
        read.setInferredInsertSize(0);
        read.setReadNegativeStrandFlag(true);
        read.setMateNegativeStrandFlag(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setReadNegativeStrandFlag(false);
        read.setMateNegativeStrandFlag(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setInferredInsertSize(10);

        // Test case 6: read is unmapped
        read = makeRead(fragmentSize, mateStart);
        read.setReadUnmappedFlag(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setReadUnmappedFlag(false);

        // Test case 7:  reads don't overlap and look like this:
        //    <--------|
        //                 |------>
        // first read:
        read = makeRead(fragmentSize, mateStart);
        myStart = 980;
        read.setAlignmentStart(myStart);
        read.setInferredInsertSize(20);
        read.setReadNegativeStrandFlag(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // second read:
        read = makeRead(fragmentSize, mateStart);
        myStart = 1000;
        read.setAlignmentStart(myStart);
        read.setInferredInsertSize(20);
        read.setMateAlignmentStart(980);
        read.setReadNegativeStrandFlag(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // Test case 8: read doesn't have proper pair flag set
        read = makeRead(fragmentSize, mateStart);
        read.setReadPairedFlag(true);
        read.setProperPairFlag(false);
        Assert.assertEquals(get.getAdaptor(read), ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // Test case 9: read and mate have same negative flag setting
        for ( final boolean negFlag: Arrays.asList(true, false) ) {
            read = makeRead(fragmentSize, mateStart);
            read.setAlignmentStart(BEFORE);
            read.setReadPairedFlag(true);
            read.setProperPairFlag(true);
            read.setReadNegativeStrandFlag(negFlag);
            read.setMateNegativeStrandFlag(!negFlag);
            Assert.assertTrue(get.getAdaptor(read) != ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY, "Get adaptor should have succeeded");

            read = makeRead(fragmentSize, mateStart);
            read.setAlignmentStart(BEFORE);
            read.setReadPairedFlag(true);
            read.setProperPairFlag(true);
            read.setReadNegativeStrandFlag(negFlag);
            read.setMateNegativeStrandFlag(negFlag);
            Assert.assertEquals(get.getAdaptor(read), ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY, "Get adaptor should have failed for reads with bad alignment orientation");
        }
    }

    @Test (enabled = true)
    public void testGetBasesReverseComplement() {
        int iterations = 1000;
        Random random = Utils.getRandomGenerator();
        while(iterations-- > 0) {
            final int l = random.nextInt(1000);
            GATKSAMRecord read = GATKSAMRecord.createRandomRead(l);
            byte [] original = read.getReadBases();
            byte [] reconverted = new byte[l];
            String revComp = ReadUtils.getBasesReverseComplement(read);
            for (int i=0; i<l; i++) {
                reconverted[l-1-i] = BaseUtils.getComplement((byte) revComp.charAt(i));
            }
            Assert.assertEquals(reconverted, original);
        }
    }

    @Test (enabled = true)
    public void testGetMaxReadLength() {
        for( final int minLength : Arrays.asList( 5, 30, 50 ) ) {
            for( final int maxLength : Arrays.asList( 50, 75, 100 ) ) {
                final List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
                for( int readLength = minLength; readLength <= maxLength; readLength++ ) {
                    reads.add( ReadUtils.createRandomRead( readLength ) );
                }
                Assert.assertEquals(ReadUtils.getMaxReadLength(reads), maxLength, "max length does not match");
            }
        }

        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        Assert.assertEquals(ReadUtils.getMaxReadLength(reads), 0, "Empty list should have max length of zero");
    }

    @Test (enabled = true)
    public void testReadWithNsRefIndexInDeletion() throws FileNotFoundException {

        final ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final int readLength = 76;

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 8975, readLength);
        read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
        read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
        read.setCigarString("3M414N1D73M");

        final int result = ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, 9392, ReadUtils.ClippingTail.LEFT_TAIL);
        Assert.assertEquals(result, 2);
    }

    @Test (enabled = true)
    public void testReadWithNsRefAfterDeletion() throws FileNotFoundException {

        final ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final int readLength = 76;

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 8975, readLength);
        read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
        read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
        read.setCigarString("3M414N1D73M");

        final int result = ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, 9393, ReadUtils.ClippingTail.LEFT_TAIL);
        Assert.assertEquals(result, 3);
    }

    @DataProvider(name = "HasWellDefinedFragmentSizeData")
    public Object[][] makeHasWellDefinedFragmentSizeData() throws Exception {
        final List<Object[]> tests = new LinkedList<Object[]>();

        // setup a basic read that will work
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 10, 10);
        read.setReadPairedFlag(true);
        read.setProperPairFlag(true);
        read.setReadUnmappedFlag(false);
        read.setMateUnmappedFlag(false);
        read.setAlignmentStart(100);
        read.setCigarString("50M");
        read.setMateAlignmentStart(130);
        read.setInferredInsertSize(80);
        read.setFirstOfPairFlag(true);
        read.setReadNegativeStrandFlag(false);
        read.setMateNegativeStrandFlag(true);

        tests.add( new Object[]{ "basic case", read.clone(), true });

        {
            final GATKSAMRecord bad1 = (GATKSAMRecord)read.clone();
            bad1.setReadPairedFlag(false);
            tests.add( new Object[]{ "not paired", bad1, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setProperPairFlag(false);
            // we currently don't require the proper pair flag to be set
            tests.add( new Object[]{ "not proper pair", bad, true });
//            tests.add( new Object[]{ "not proper pair", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setReadUnmappedFlag(true);
            tests.add( new Object[]{ "read is unmapped", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setMateUnmappedFlag(true);
            tests.add( new Object[]{ "mate is unmapped", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setMateNegativeStrandFlag(false);
            tests.add( new Object[]{ "read and mate both on positive strand", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setReadNegativeStrandFlag(true);
            tests.add( new Object[]{ "read and mate both on negative strand", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setInferredInsertSize(0);
            tests.add( new Object[]{ "insert size is 0", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setAlignmentStart(1000);
            tests.add( new Object[]{ "positve read starts after mate end", bad, false });
        }

        {
            final GATKSAMRecord bad = (GATKSAMRecord)read.clone();
            bad.setReadNegativeStrandFlag(true);
            bad.setMateNegativeStrandFlag(false);
            bad.setMateAlignmentStart(1000);
            tests.add( new Object[]{ "negative strand read ends before mate starts", bad, false });
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "HasWellDefinedFragmentSizeData")
    private void testHasWellDefinedFragmentSize(final String name, final GATKSAMRecord read, final boolean expected) {
        Assert.assertEquals(ReadUtils.hasWellDefinedFragmentSize(read), expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetSAMFileMissingHeader() {
        ReadUtils.getSAMFileSamples(null);
    }

    @Test(expectedExceptions = UserException.class)
    public void testGetSAMFileMissingReadGroupsSamples() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        final SAMReadGroupRecord samGroup = new SAMReadGroupRecord("id");
        final List<SAMReadGroupRecord> list = new ArrayList<>(Arrays.asList(samGroup));
        header.setReadGroups(list);
        ReadUtils.getSAMFileSamples(header);
    }
}
