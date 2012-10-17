package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 4/5/12
 * Time: 11:28 AM
 * To change this template use File | Settings | File Templates.
 */
public class UnifiedGenotyperGeneralPloidyIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String CEUTRIO_BAM = "/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37.list";
    final String LSV_BAM = validationDataLocation +"93pools_NA12878_ref_chr20_40m_41m.bam";
    final String REFSAMPLE_MT_CALLS = comparisonDataLocation + "Unvalidated/mtDNA/NA12878.snp.vcf";
    final String REFSAMPLE_NAME = "NA12878";
    final String MTINTERVALS = "MT:1-1000";
    final String LSVINTERVALS = "20:40,500,000-41,000,000";
    final String LSVINTERVALS_SHORT = "20:40,500,000-40,501,000";
    final String NA12891_CALLS = comparisonDataLocation + "Unvalidated/mtDNA/NA12891.snp.vcf";
    final String NA12878_WG_CALLS = comparisonDataLocation + "Unvalidated/NA12878/CEUTrio.HiSeq.WGS.b37_decoy.recal.ts_95.snp_indel_combined.vcf";
    final String LSV_ALLELES = validationDataLocation + "ALL.chr20_40m_41m.largeScaleValidationSites.vcf";

    private void PC_MT_Test(String bam, String args, String name, String md5) {
        final String base = String.format("-T UnifiedGenotyper -dcov 10000 -R %s -I %s -L %s --reference_sample_calls %s -refsample %s -ignoreLane ",
                REF, bam, MTINTERVALS, REFSAMPLE_MT_CALLS, REFSAMPLE_NAME) + " --no_cmdline_in_header -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testPoolCaller:"+name+" args=" + args, spec);
    }

    private void PC_LSV_Test(String args, String name, String model, String md5) {
        final String base = String.format("-T UnifiedGenotyper -dcov 10000 -R %s -I %s -L %s --reference_sample_calls %s -refsample %s -glm %s -ignoreLane ",
                REF, LSV_BAM, LSVINTERVALS, NA12878_WG_CALLS, REFSAMPLE_NAME, model) + " --no_cmdline_in_header -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testPoolCaller:"+name+" args=" + args, spec);
    }

    private void PC_LSV_Test_short(String args, String name, String model, String md5) {
        final String base = String.format("-T UnifiedGenotyper -dcov 10000 -R %s -I %s -L %s --reference_sample_calls %s -refsample %s -glm %s -ignoreLane ",
                REF, LSV_BAM, LSVINTERVALS_SHORT, NA12878_WG_CALLS, REFSAMPLE_NAME, model) + " --no_cmdline_in_header -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testPoolCaller:"+name+" args=" + args, spec);
    }

    private void PC_LSV_Test_NoRef(String args, String name, String model, String md5) {
        final String base = String.format("-T UnifiedGenotyper -dcov 10000 -R %s -I %s -L %s -glm %s -ignoreLane",
                REF, LSV_BAM, LSVINTERVALS, model) + " --no_cmdline_in_header -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testPoolCaller:"+name+" args=" + args, spec);
    }

    @Test(enabled = true)
    public void testSNP_ACS_Pools() {
        PC_LSV_Test_short(" -maxAltAlleles 1 -ploidy 6 -out_mode EMIT_ALL_CONFIDENT_SITES","LSV_SNP_ACS","SNP","ec19f0b7c7d57493cecfff988a4815c8");
    }

    @Test(enabled = true)
    public void testBOTH_GGA_Pools() {
        PC_LSV_Test(String.format(" -maxAltAlleles 2 -ploidy 24 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles %s",LSV_ALLELES),"LSV_BOTH_GGA","BOTH","67dabdbf1e6ed8a83d2e85766558a20a");
    }

    @Test(enabled = true)
    public void testINDEL_GGA_Pools() {
        PC_LSV_Test(String.format(" -maxAltAlleles 1 -ploidy 24 -gt_mode GENOTYPE_GIVEN_ALLELES  -out_mode EMIT_ALL_SITES -alleles %s",LSV_ALLELES),"LSV_INDEL_GGA","INDEL","d4bfae27f1b07923f381d708d8a34cf4");
    }

    @Test(enabled = true)
    public void testINDEL_maxAltAlleles2_ploidy3_Pools_noRef() {
        PC_LSV_Test_NoRef(" -maxAltAlleles 2 -ploidy 3","LSV_INDEL_DISC_NOREF_p3","INDEL","7d6f319b9edcb1ff8c290fef150a2df8");
    }

    @Test(enabled = true)
    public void testINDEL_maxAltAlleles2_ploidy1_Pools_noRef() {
        PC_LSV_Test_NoRef(" -maxAltAlleles 2 -ploidy 1","LSV_INDEL_DISC_NOREF_p1","INDEL","dd02890123e07e7412a49475cb6280f1");
    }

    @Test(enabled = true)
    public void testMT_SNP_DISCOVERY_sp4() {
         PC_MT_Test(CEUTRIO_BAM, " -maxAltAlleles 1 -ploidy 8", "MT_SNP_DISCOVERY_sp4","da84bf45f7080a46a7a78542b3a0629d");
    }

    @Test(enabled = true)
    public void testMT_SNP_GGA_sp10() {
        PC_MT_Test(CEUTRIO_BAM, String.format(" -maxAltAlleles 1 -ploidy 20 -gt_mode GENOTYPE_GIVEN_ALLELES  -out_mode EMIT_ALL_SITES -alleles %s",NA12891_CALLS), "MT_SNP_GGA_sp10", "f8ea18ec6a717a77fdf8c5f2482d8d8d");
    }
}
