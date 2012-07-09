package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String DELETION_BAM = validationDataLocation + "filtered_deletion_for_reduce_reads.bam";
    final String STASH_BAM = validationDataLocation + "ReduceReadsStashBug.bam";
    final String STASH_L = " -L 14:73718184-73718284 -L 14:73718294-73718330 -L 14:73718360-73718556";
    final String L = " -L 20:10,100,000-10,120,000 ";

    private void RRTest(String testName, String args, String md5) {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, BAM) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("testDefaultCompression ", L, "887d05e3bdbe831a2584305088806a39");
    }

    @Test(enabled = true)
    public void testMultipleIntervals() {
        String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";
        RRTest("testMultipleIntervals ", intervals, "04208b8e2c2e65bcb9d0eb038702b006");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -cs 10 -minvar 0.3 -mindel 0.3 " + L, "3a607bc3ebaf84e9dc44e005c5f8a047");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -cs 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "afd39459c841b68a442abdd5ef5f8f27");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        RRTest("testIndelCompression ", " -cs 50 -L 20:10,100,500-10,100,600 ", "f7b9fa44c10bc4b2247813d2b8dc1973");
    }

    @Test(enabled = true)
    public void testFilteredDeletionCompression() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, DELETION_BAM) + " -o %s ";
        executeTest("testFilteredDeletionCompression", new WalkerTestSpec(base, Arrays.asList("891bd6dcda66611f343e8ff25f34aaeb")));
    }

    /**
     * Bug reported by Adam where a read that got clipped before actually belongs 2 intervals ahead
     * and a subsequent tail leaves only this read in the stash. The next read to come in is in fact 
     * before (alignment start) than this read, so the TreeSet breaks with a Key out of Range error
     * that was freaking hard to catch. 
     * 
     * This bam is simplified to replicate the exact bug with the three provided intervals.
     */
    @Test(enabled = true)
    public void testAddingReadAfterTailingTheStash() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s", STASH_L, REF, STASH_BAM) + " -o %s ";
        executeTest("testAddingReadAfterTailingTheStash", new WalkerTestSpec(base, Arrays.asList("886b43e1f26ff18425814dc7563931c6")));
    }
}

