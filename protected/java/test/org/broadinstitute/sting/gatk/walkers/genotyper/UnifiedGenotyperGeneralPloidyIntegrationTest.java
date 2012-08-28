package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;

import java.util.Arrays;
import org.testng.annotations.Test;

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

    private void PC_LSV_Test_NoRef(String args, String name, String model, String md5) {
        final String base = String.format("-T UnifiedGenotyper -dcov 10000 -R %s -I %s -L %s -glm %s -ignoreLane",
                REF, LSV_BAM, LSVINTERVALS, model) + " --no_cmdline_in_header -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testPoolCaller:"+name+" args=" + args, spec);
    }

    @Test(enabled = true)
    public void testBOTH_GGA_Pools() {
        PC_LSV_Test(String.format(" -maxAltAlleles 2 -ploidy 24 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles %s",LSV_ALLELES),"LSV_BOTH_GGA","BOTH","0ff90fa3882a3fb5089a7bba50dd8ae3");
    }

    @Test(enabled = true)
    public void testINDEL_GGA_Pools() {
        PC_LSV_Test(String.format(" -maxAltAlleles 1 -ploidy 24 -gt_mode GENOTYPE_GIVEN_ALLELES  -out_mode EMIT_ALL_SITES -alleles %s",LSV_ALLELES),"LSV_INDEL_GGA","INDEL","90af837f372e3d5143af30bf5c8c2b75");
    }

    @Test(enabled = true)
    public void testINDEL_maxAltAlleles2_ploidy3_Pools_noRef() {
        PC_LSV_Test_NoRef(" -maxAltAlleles 2 -ploidy 3","LSV_INDEL_DISC_NOREF_p3","INDEL","9514ed15c7030b6d47e04e6a3a2b0a3e");
    }

    @Test(enabled = true)
    public void testINDEL_maxAltAlleles2_ploidy1_Pools_noRef() {
        PC_LSV_Test_NoRef(" -maxAltAlleles 2 -ploidy 1","LSV_INDEL_DISC_NOREF_p1","INDEL","26598044436c8044f22ffa767b06a0f0");
    }

    @Test(enabled = true)
    public void testMT_SNP_DISCOVERY_sp4() {
         PC_MT_Test(CEUTRIO_BAM, " -maxAltAlleles 1 -ploidy 8", "MT_SNP_DISCOVERY_sp4","c32e10070e10d30d33e5b882c1f89413");
    }

    @Test(enabled = true)
    public void testMT_SNP_GGA_sp10() {
        PC_MT_Test(CEUTRIO_BAM, String.format(" -maxAltAlleles 1 -ploidy 20 -gt_mode GENOTYPE_GIVEN_ALLELES  -out_mode EMIT_ALL_SITES -alleles %s",NA12891_CALLS), "MT_SNP_GGA_sp10", "4d16d3c9475637bad70e9dc2eafe2da2");
    }
}
