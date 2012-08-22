package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class HaplotypeCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String NA12878_BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final String INTERVALS_FILE = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals";
    //final String RECAL_FILE = validationDataLocation + "NA12878.kmer.8.subset.recal_data.bqsr";

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
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10431524-10431924 -L 20:10723661-10724061 -L 20:10903555-10903955 --no_cmdline_in_header -o %s -minPruning 3";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerComplexVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSampleComplex() {
        HCTestComplexVariants(CEUTRIO_BAM, "", "316a0fb9c66c0a6aa40a170d5d8c0021");
    }
}

