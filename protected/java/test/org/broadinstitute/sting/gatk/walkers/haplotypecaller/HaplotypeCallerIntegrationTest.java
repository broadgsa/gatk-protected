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
        HCTest(CEUTRIO_BAM, "", "7122d4f0ef94c5274aa3047cfebe08ed");
    }

    @Test
    public void testHaplotypeCallerSingleSample() {
        HCTest(NA12878_BAM, "", "6cd6e6787521c07a7bae98766fd628ab");
    }

    // TODO -- add more tests for GGA mode, especially with input alleles that are complex variants and/or not trimmed
    @Test
    public void testHaplotypeCallerMultiSampleGGA() {
        HCTest(CEUTRIO_BAM, "--max_alternate_alleles 3 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf",
                "44df2a9da4fbd2162ae44c3f2a6ef01f");
    }

    private void HCTestComplexVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10028767-10028967 -L 20:10431524-10431924 -L 20:10723661-10724061 -L 20:10903555-10903955 --no_cmdline_in_header -o %s -minPruning 4";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerComplexVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSampleComplex() {
        HCTestComplexVariants(privateTestDir + "AFR.complex.variants.bam", "", "4a413eeb7a75cab0ab5370b4c08dcf8e");
    }

    private void HCTestSymbolicVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:5947969-5948369 -L 20:61091236-61091636 --no_cmdline_in_header -o %s -minPruning 1";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerSymbolicVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleSymbolic() {
        HCTestSymbolicVariants(NA12878_CHR20_BAM, "", "77cf5b5273828dd1605bb23a5aeafcaa");
    }

    private void HCTestIndelQualityScores(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10,005,000-10,025,000 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerIndelQualityScores: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleIndelQualityScores() {
        HCTestIndelQualityScores(NA12878_RECALIBRATED_BAM, "", "87ca97f90e74caee35c35616c065821c");
    }

    @Test
    public void HCTestProblematicReadsModifiedInActiveRegions() {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "haplotype-problem-4.bam") + " --no_cmdline_in_header -o %s -minPruning 3 -L 4:49139026-49139965";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("3df42d0550b51eb9b55aac61e8b3c452"));
        executeTest("HCTestProblematicReadsModifiedInActiveRegions: ", spec);
    }

    @Test
    public void HCTestStructuralIndels() {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "AFR.structural.indels.bam") + " --no_cmdline_in_header -o %s -minPruning 6 -L 20:8187565-8187800 -L 20:18670537-18670730";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("4dbc72b72e3e2d9d812d5a398490e213"));
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
                Arrays.asList("f8c2745bf71f2659a57494fcaa2c103b"));
        executeTest("HC calling on a ReducedRead BAM", spec);
    }
}
