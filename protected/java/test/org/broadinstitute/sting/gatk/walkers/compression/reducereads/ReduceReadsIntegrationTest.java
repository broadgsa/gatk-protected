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
    final String DIVIDEBYZERO_BAM = validationDataLocation + "ReduceReadsDivideByZeroBug.bam";
    final String DIVIDEBYZERO_L = " -L " + validationDataLocation + "ReduceReadsDivideByZeroBug.intervals";
    final String L = " -L 20:10,100,000-10,120,000 ";
    final String COREDUCTION_BAM_A = validationDataLocation + "coreduction.test.A.bam";
    final String COREDUCTION_BAM_B = validationDataLocation + "coreduction.test.B.bam";
    final String COREDUCTION_L = " -L 1:1,853,860-1,854,354 -L 1:1,884,131-1,892,057";
    final String OFFCONTIG_BAM = privateTestDir + "readOffb37contigMT.bam";

    private void RRTest(String testName, String args, String md5) {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, BAM) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("testDefaultCompression ", L, "98080d3c53f441564796fc143cf510da");
    }

    @Test(enabled = true)
    public void testMultipleIntervals() {
        String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";
        RRTest("testMultipleIntervals ", intervals, "c5dcdf4edf368b5b897d66f76034d9f0");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -cs 10 -minvar 0.3 -mindel 0.3 " + L, "27cb99e87eda5e46187e56f50dd37f26");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -cs 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "4e7f111688d49973c35669855b7a2eaf");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        RRTest("testIndelCompression ", " -cs 50 -L 20:10,100,500-10,100,600 ", "f6c9ea83608f35f113cf1f62a77ee6d0");
    }

    @Test(enabled = true)
    public void testFilteredDeletionCompression() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, DELETION_BAM) + " -o %s ";
        executeTest("testFilteredDeletionCompression", new WalkerTestSpec(base, Arrays.asList("122e4e60c4412a31d0aeb3cce879e841")));
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
        executeTest("testAddingReadAfterTailingTheStash", new WalkerTestSpec(base, Arrays.asList("647b0f0f95730de8e6bc4f74186ad4df")));
    }

    /**
     * Divide by zero bug reported by GdA and users in the forum. Happens when the downsampler goes over a region where all reads get
     * filtered out.
     */
    @Test(enabled = true)
    public void testDivideByZero() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s", DIVIDEBYZERO_L, REF, DIVIDEBYZERO_BAM) + " -o %s ";
        executeTest("testDivideByZero", new WalkerTestSpec(base, Arrays.asList("2c87985972dd43ee9dd50b463d93a511")));
    }

    @Test(enabled = true)
    public void testCoReduction() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s -I %s", COREDUCTION_L, REF, COREDUCTION_BAM_A, COREDUCTION_BAM_B) + " -o %s ";
        executeTest("testCoReduction", new WalkerTestSpec(base, Arrays.asList("5c30fde961a1357bf72c15144c01981b")));
    }

    /**
     * Bug happens when reads are soft-clipped off the  contig (usually in the MT). This test guarantees no changes to the upstream code will
     * break the current hard-clipping routine that protects reduce reads from such reads.
     */
    @Test(enabled = true)
    public void testReadOffContig() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, OFFCONTIG_BAM) + " -o %s ";
        executeTest("testReadOffContig", new WalkerTestSpec(base, Arrays.asList("2f17c1a78e9d0138217fdb83cede8f68")));
    }

}

