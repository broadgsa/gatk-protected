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
    final String MTINTERVALS = "MT:1-3000";
    final String LSVINTERVALS = "20:40,000,000-41,000,000";
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

    @Test
    public void testBOTH_GGA_Pools() {
        PC_LSV_Test(String.format(" -maxAltAlleles 2 -ploidy 24 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles %s",LSV_ALLELES),"LSV_BOTH_GGA","BOTH","d76e3b910259da819f1e1b2adc68ba8d");
    }

    @Test
    public void testINDEL_GGA_Pools() {
        PC_LSV_Test(String.format(" -maxAltAlleles 1 -ploidy 24 -gt_mode GENOTYPE_GIVEN_ALLELES  -out_mode EMIT_ALL_SITES -alleles %s",LSV_ALLELES),"LSV_INDEL_GGA","INDEL","ffadcdaee613dab975197bed0fc78da3");
    }

    @Test
    public void testINDEL_maxAltAlleles2_ploidy3_Pools_noRef() {
        PC_LSV_Test_NoRef(" -maxAltAlleles 2 -ploidy 3","LSV_INDEL_DISC_NOREF_p3","INDEL","96087fe9240e3656cc2a4e0ff0174d5b");
    }

    @Test
    public void testINDEL_maxAltAlleles2_ploidy1_Pools_noRef() {
        PC_LSV_Test_NoRef(" -maxAltAlleles 2 -ploidy 1","LSV_INDEL_DISC_NOREF_p1","INDEL","6fdae7093831ecfc82a06dd707d62fe9");
    }

    @Test
    public void testMT_SNP_DISCOVERY_sp4() {
         PC_MT_Test(CEUTRIO_BAM, " -maxAltAlleles 1 -ploidy 8", "MT_SNP_DISCOVERY_sp4","6b27634214530d379db70391a9cfc2d7");
    }

    @Test
    public void testMT_SNP_GGA_sp10() {

        PC_MT_Test(CEUTRIO_BAM, String.format(" -maxAltAlleles 1 -ploidy 20 -gt_mode GENOTYPE_GIVEN_ALLELES  -out_mode EMIT_ALL_SITES -alleles %s",NA12891_CALLS), "MT_SNP_GGA_sp10", "e74d4c73ece45d7fb676b99364df4f1a");
    }

}
