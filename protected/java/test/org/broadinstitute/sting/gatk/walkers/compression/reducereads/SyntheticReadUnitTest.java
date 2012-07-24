package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;

public class SyntheticReadUnitTest extends BaseTest {
    final SAMFileHeader artificialSAMHeader = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1);
    final GATKSAMReadGroupRecord artificialGATKRG = new GATKSAMReadGroupRecord("synthetic");
    final String artificialContig = "1";
    final int artificialContigIndex = 0;
    final String artificialReadName = "synth";
    final int artificialRefStart = 1;
    final double artificialMappingQuality = 60;

    final Random random = new Random(8854875);


@Test
public void testBaseCounts() {
        BaseIndex [] bases = new BaseIndex[] {BaseIndex.A,BaseIndex.A,BaseIndex.A,BaseIndex.A};
        Byte[] quals = new Byte[] {20, 20, 20, 20 };

        TestRead [] testReads = new TestRead [] {
                new TestRead(bases, quals, new Byte[] {100, 100, 100, 101}, new byte [] {100, 0, 0, 1}),
                new TestRead(bases, quals, new Byte[] {1, 100, 100, 0},     new byte [] {1, 99, 99, -1}),
                new TestRead(bases, quals, new Byte[] {127, 100, 0, 1},     new byte [] {127, -27, -127, -126}),
                new TestRead(bases, quals, new Byte[] {1, 127, 51, 126},    new byte [] {1, 126, 50, 125})};

        for (TestRead testRead : testReads) {
            SyntheticRead syntheticRead = new SyntheticRead(Arrays.asList(testRead.getBases()), Arrays.asList(testRead.getCounts()), Arrays.asList(testRead.getQuals()), Arrays.asList(testRead.getInsQuals()), Arrays.asList(testRead.getDelQuals()), artificialMappingQuality, GATKSAMRecord.REDUCED_READ_CONSENSUS_TAG, artificialSAMHeader, artificialGATKRG, artificialContig, artificialContigIndex, artificialReadName, artificialRefStart, false);
            Assert.assertEquals(syntheticRead.convertBaseCounts(), testRead.getExpectedCounts());
        }
}

private class TestRead {
    BaseIndex[] bases;
    Byte[] quals;
    Byte[] insQuals;
    Byte[] delQuals;
    Byte[] counts;
    byte [] expectedCounts;

    private TestRead(BaseIndex[] bases, Byte[] quals, Byte[] counts, byte[] expectedCounts) {
        this.bases = bases;
        this.quals = quals;
        this.insQuals = quals;
        this.delQuals = quals;
        this.counts = counts;
        this.expectedCounts = expectedCounts;
    }

    public BaseIndex[] getBases() {
        return bases;
    }

    public Byte[] getQuals() {
        return quals;
    }

    public Byte[] getInsQuals() {
        return insQuals;
    }

    public Byte[] getDelQuals() {
        return delQuals;
    }

    public Byte[] getCounts() {
        return counts;
    }

    public byte[] getExpectedCounts() {
        return expectedCounts;
    }
}

}
