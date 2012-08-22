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
        HCTest(CEUTRIO_BAM, "", "e5b4a0627a1d69b9356f8a7cd2260e89");
    }

    @Test
    public void testHaplotypeCallerSingleSample() {
        HCTest(NA12878_BAM, "", "202d5b6edaf74f411c170099749f202f");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGA() {
        HCTest(CEUTRIO_BAM, "-gt_mode GENOTYPE_GIVEN_ALLELES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf", "561931ba3919808ec471e745cb3148c7");
    }

    private void HCTestComplexVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10028767-10028967 -L 20:10431524-10431924 -L 20:10723661-10724061 -L 20:10903555-10903955 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerComplexVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSampleComplex() {
        HCTestComplexVariants(CEUTRIO_BAM, "", "3424b398a9f47c8ac606a5c56eb7d8a7");
    }

    private void HCTestSymbolicVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:5947969-5948369 -L 20:61091236-61091636 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerSymbolicVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleSymbolic() {
        HCTestSymbolicVariants(NA12878_CHR20_BAM, "", "b71cfaea9390136c584c9671b149d573");
    }

    private void HCTestIndelQualityScores(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10,005,000-10,025,000 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerIndelQualityScores: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleIndelQualityScores() {
        HCTestIndelQualityScores(NA12878_RECALIBRATED_BAM, "", "e1f88fac91424740c0eaac1de48b3970");
    }
}
