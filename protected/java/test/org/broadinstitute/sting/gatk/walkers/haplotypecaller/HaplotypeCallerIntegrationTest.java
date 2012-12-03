package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class HaplotypeCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String NA12878_BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String NA12878_CHR20_BAM = validationDataLocation + "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam";
    final String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final String NA12878_RECALIBRATED_BAM = privateTestDir + "NA12878.100kb.BQSRv2.example.bam";
    final String INTERVALS_FILE = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals";

    private void HCTest(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s -L %s", REF, bam, INTERVALS_FILE) + " --no_cmdline_in_header -o %s -minPruning 3";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCaller: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSample() {
        HCTest(CEUTRIO_BAM, "", "2b39732ff8e0de5bc2ae949aaf7a6f21");
    }

    @Test
    public void testHaplotypeCallerSingleSample() {
        HCTest(NA12878_BAM, "", "8b217638ff585effb9cc70e9a9aa544f");
    }

    // TODO -- add more tests for GGA mode, especially with input alleles that are complex variants and/or not trimmed
    @Test
    public void testHaplotypeCallerMultiSampleGGA() {
        HCTest(CEUTRIO_BAM, "--max_alternate_alleles 3 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf",
                "541aa8291f03ba33bd1ad3d731fd5657");
    }

    private void HCTestComplexVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10028767-10028967 -L 20:10431524-10431924 -L 20:10723661-10724061 -L 20:10903555-10903955 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerComplexVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSampleComplex() {
        HCTestComplexVariants(privateTestDir + "AFR.complex.variants.bam", "", "fd7170cbde7df04d4fbe1da7903c31c6");
    }

    private void HCTestSymbolicVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:5947969-5948369 -L 20:61091236-61091636 --no_cmdline_in_header -o %s -minPruning 1";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerSymbolicVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleSymbolic() {
        HCTestSymbolicVariants(NA12878_CHR20_BAM, "", "99456fc7207c1fe9f367a0d0afae87cd");
    }

    private void HCTestIndelQualityScores(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10,005,000-10,025,000 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerIndelQualityScores: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleIndelQualityScores() {
        HCTestIndelQualityScores(NA12878_RECALIBRATED_BAM, "", "6c1631785b3f832aecab1a99f0454762");
    }

    @Test
    public void HCTestProblematicReadsModifiedInActiveRegions() {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "haplotype-problem-4.bam") + " --no_cmdline_in_header -o %s -minPruning 3 -L 4:49139026-49139965";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("ec437d2d9f3ae07d155983be0155c8ed"));
        executeTest("HCTestProblematicReadsModifiedInActiveRegions: ", spec);
    }

    @Test
    public void HCTestStructuralIndels() {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "AFR.structural.indels.bam") + " --no_cmdline_in_header -o %s -minPruning 6 -L 20:8187565-8187800 -L 20:18670537-18670730";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("237601bbc39694c7413a332cbb656c8e"));
        executeTest("HCTestStructuralIndels: ", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing reduced reads
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void HCTestReducedBam() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam -o %s -L 1:67,225,396-67,288,518", 1,
                Arrays.asList("40bf739fb2b1743642498efe79ea6342"));
        executeTest("HC calling on a ReducedRead BAM", spec);
    }
}
