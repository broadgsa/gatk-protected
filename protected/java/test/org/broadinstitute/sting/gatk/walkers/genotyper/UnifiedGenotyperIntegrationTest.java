package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

// ********************************************************************************** //
// Note that this class also serves as an integration test for the VariantAnnotator!  //
// ********************************************************************************** //

public class UnifiedGenotyperIntegrationTest extends WalkerTest {

    private final static String baseCommand = "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommandIndels = "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm INDEL -mbq 20 -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommandIndelsb37 = "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -mbq 20 --dbsnp " + b37dbSNP132;
    private final static String baseCommandNoCmdLineHeaderStdout = "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing normal calling
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000", 1,
                Arrays.asList("847605f4efafef89529fe0e496315edd"));
        executeTest("test MultiSample Pilot1", spec);
    }

    @Test
    public void testWithAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("5b31b811072a4df04524e13604015f9b"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec1);
    }

    @Test
    public void testWithAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("d9992e55381afb43742cc9b30fcd7538"));
        executeTest("test MultiSample Pilot2 with alleles passed in and emitting all sites", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("fea530fdc8677e10be4cc11625fa5376"));
        executeTest("test SingleSample Pilot2", spec);
    }

    @Test
    public void testMultipleSNPAlleles() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm BOTH --dbsnp " + b37dbSNP129 + " -I " + privateTestDir + "multiallelic.snps.bam -o %s -L " + privateTestDir + "multiallelic.snps.intervals", 1,
                Arrays.asList("d20c7a143b899f0239bf64b652ad3edb"));
        executeTest("test Multiple SNP alleles", spec);
    }

    @Test
    public void testBadRead() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm BOTH -I " + privateTestDir + "badRead.test.bam -o %s -L 1:22753424-22753464", 1,
                Arrays.asList("d915535c1458733f09f82670092fcab6"));
        executeTest("test bad read", spec);
    }

    @Test
    public void testReverseTrim() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -I " + validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam -o %s -L 20:10289124 -L 20:10090289", 1,
                Arrays.asList("e14c9b1f9f34d6c16de445bfa385be89"));
        executeTest("test reverse trim", spec);
    }

    @Test
    public void testMismatchedPLs() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -I " + privateTestDir + "mismatchedPLs.bam -o %s -L 1:24020341", 1,
                Arrays.asList("fb204e821a24d03bd3a671b6e01c449a"));
        executeTest("test mismatched PLs", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "5b8f477c287770b5769b05591e35bc2d";

    @Test
    public void testCompressedOutput() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("gz"), Arrays.asList(COMPRESSED_OUTPUT_MD5));
        executeTest("test compressed output", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parallelization
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testParallelization() {

        // Note that we need to turn off any randomization for this to work, so no downsampling and no annotations

        String md5 = "d408b4661b820ed86272415b8ea08780";

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -dt NONE -G none --contamination_fraction_to_filter 0.0 -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,075,000", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (single thread)", spec1);

        GenomeAnalysisEngine.resetRandomGenerator();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -dt NONE -G none --contamination_fraction_to_filter 0.0 -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,075,000 -nt 2", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (2 threads)", spec2);

        GenomeAnalysisEngine.resetRandomGenerator();

        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -dt NONE -G none --contamination_fraction_to_filter 0.0 -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,075,000 -nt 4", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (4 threads)", spec3);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parameters
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testMinBaseQualityScore() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 --min_base_quality_score 26", 1,
                Arrays.asList("6ee6537e9ebc1bfc7c6cf8f04b1582ff"));
        executeTest("test min_base_quality_score 26", spec);
    }

    @Test
    public void testSLOD() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --computeSLOD --no_cmdline_in_header -glm BOTH --dbsnp " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("55760482335497086458b09e415ecf54"));
        executeTest("test SLOD", spec);
    }

    @Test
    public void testNDA() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " --annotateNDA -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("938e888a40182878be4c3cc4859adb69"));
        executeTest("test NDA", spec);
    }

    @Test
    public void testCompTrack() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -comp:FOO " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("7dc186d420487e4e156a24ec8dea0951"));
        executeTest("test using comp track", spec);
    }

    @Test(enabled = false) // EB: for some reason this test crashes whenever I run it on my local machine
    public void testNoCmdLineHeaderStdout() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandNoCmdLineHeaderStdout + " -glm INDEL -L 1:67,225,396-67,288,518", 0,
                Collections.<String>emptyList());
        executeTest("testNoCmdLineHeaderStdout", spec);
    }

    @Test
    public void testOutputParameterSitesOnly() {
        testOutputParameters("-sites_only", "f99c7471127a6fb6f72e136bc873b2c9");
    }

    @Test
    public void testOutputParameterAllConfident() {
        testOutputParameters("--output_mode EMIT_ALL_CONFIDENT_SITES", "9dbc9389db39cf9697e93e0bf529314f");
    }

    @Test
    public void testOutputParameterAllSites() {
        testOutputParameters("--output_mode EMIT_ALL_SITES", "8b26088a035e579c4afd3b46737291e4");
    }

    private void testOutputParameters(final String args, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 " + args, 1,
                Arrays.asList(md5));
        executeTest(String.format("testParameter[%s]", args), spec);
    }

    @Test
    public void testConfidence() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_call_conf 10 ", 1,
                Arrays.asList("4af83a883ecc03a23b0aa6dd4b8f1ceb"));
        executeTest("test confidence 1", spec1);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing heterozygosity
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity1() {
        testHeterozosity( 0.01, "8dd37249e0a80afa86594c3f1e720760" );
    }

    @Test
    public void testHeterozyosity2() {
        testHeterozosity( 1.0 / 1850, "040d169e20fda56f8de009a6015eb384" );
    }

    private void testHeterozosity(final double arg, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000 --heterozygosity " + arg, 1,
                Arrays.asList(md5));
        executeTest(String.format("test heterozyosity[%s]", arg), spec);
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // testing calls with SLX, 454, and SOLID data
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiTechnologies() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,100,000",
                1,
                Arrays.asList("0e4713e4aa44f4f8fcfea7138295a627"));

        executeTest(String.format("test multiple technologies"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing calls with BAQ
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testCallingWithBAQ() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,100,000" +
                        " -baq CALCULATE_AS_NECESSARY",
                1,
                Arrays.asList("46ea5d1ceb8eed1d0db63c3577915d6c"));

        executeTest(String.format("test calling with BAQ"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing indel caller
    //
    // --------------------------------------------------------------------------------------------------------------
    // Basic indel testing with SLX data
    @Test
    public void testSimpleIndels() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,500,000",
                1,
                Arrays.asList("50329e15e5139be9e3b643f0b3ba8a53"));

        executeTest(String.format("test indel caller in SLX"), spec);
    }

    // Basic indel testing with SLX data
    @Test
    public void testIndelsWithLowMinAlleleCnt() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam" +
                        " -o %s" +
                        " -minIndelCnt 1" +
                        " -L 1:10,000,000-10,100,000",
                1,
                Arrays.asList("2b85e3bd6bf981afaf7324666740d74b"));

        executeTest(String.format("test indel caller in SLX with low min allele count"), spec);
    }

    @Test
    public void testMultiTechnologyIndels() {
         WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                 baseCommandIndels +
                         " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                         " -o %s" +
                         " -L 1:10,000,000-10,500,000",
                 1,
                 Arrays.asList("a6fd46eff78827060451a62cffd698a7"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    @Test
    public void testWithIndelAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("b8129bf754490cc3c76191d8cc4ec93f"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec);
    }

    @Test
    public void testWithIndelAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles "
                        + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("591332fa0b5b22778cf820ee257049d2"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in and emitting all sites", spec);
    }

    @Test
    public void testMultiSampleIndels1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("69df7a00f800204564ca3726e1871132"));
        List<File> result = executeTest("test MultiSample Pilot1 CEU indels", spec1).getFirst();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + result.get(0).getAbsolutePath() + " -I " + validationDataLocation +
                        "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("1256a7eceff2c2374c231ff981df486d"));
        executeTest("test MultiSample Pilot1 CEU indels using GENOTYPE_GIVEN_ALLELES", spec2);
    }

    @Test
    public void testGGAwithNoEvidenceInReads() {
        final String vcf = "small.indel.test.vcf";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndelsb37 + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + privateTestDir + vcf + " -I " + validationDataLocation +
                        "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam -o %s -L " + validationDataLocation + vcf, 1,
                Arrays.asList("d76eacc4021b78ccc0a9026162e814a7"));
        executeTest("test GENOTYPE_GIVEN_ALLELES with no evidence in reads", spec);
    }

    @Test
    public void testBaseIndelQualityScores() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndelsb37 +
                        " -I " + privateTestDir + "NA12878.100kb.BQSRv2.example.bam" +
                        " -o %s" +
                        " -L 20:10,000,000-10,100,000",
                1,
                Arrays.asList("1e0d2c15546c3b0959b00ffb75488b56"));

        executeTest(String.format("test UG with base indel quality scores"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing SnpEff
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testSnpEffAnnotationRequestedWithoutRodBinding() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000 " +
                "-A SnpEff",
                1,
                UserException.class);
        executeTest("testSnpEffAnnotationRequestedWithoutRodBinding", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing MinIndelFraction
    //
    // --------------------------------------------------------------------------------------------------------------

    final static String assessMinIndelFraction = baseCommandIndelsb37 + " -I " + validationDataLocation
            + "978604.bam -L 1:978,586-978,626 -o %s --sites_only -rf Sample -goodSM 7377 -goodSM 22-0022 -goodSM 134 -goodSM 344029-53 -goodSM 14030";

    @Test
    public void testMinIndelFraction0() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.0", 1,
                Arrays.asList("90adefd39ed67865b0cb275ad0f07383"));
        executeTest("test minIndelFraction 0.0", spec);
    }

    @Test
    public void testMinIndelFraction25() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.25", 1,
                Arrays.asList("2fded43949e258f8e9f68893c61c1bdd"));
        executeTest("test minIndelFraction 0.25", spec);
    }

    @Test
    public void testMinIndelFraction100() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 1", 1,
                Arrays.asList("3f07efb768e08650a7ce333edd4f9a52"));
        executeTest("test minIndelFraction 1.0", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing Ns in CIGAR
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testNsInCigar() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + validationDataLocation + "testWithNs.bam -o %s -L 8:141813600-141813700 -out_mode EMIT_ALL_SITES", 1,
                Arrays.asList("32f18ba50406cd8c8069ba07f2f89558"));
        executeTest("test calling on reads with Ns in CIGAR", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing reduced reads
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testReducedBam() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam -o %s -L 1:67,225,396-67,288,518", 1,
                Arrays.asList("c1077662411164182c5f75478344f83d"));
        executeTest("test calling on a ReducedRead BAM", spec);
    }

    @Test
    public void testReducedBamSNPs() {
        testReducedCalling("SNP", "dee6590e3b7079890bc3a9cb372c297e");
    }

    @Test
    public void testReducedBamINDELs() {
        testReducedCalling("INDEL", "3c02ee5187933bed44dc416a2e28511f");
    }


    private void testReducedCalling(final String model, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.reduced.bam -o %s -L 20:10,000,000-11,000,000 -glm " + model, 1,
                Arrays.asList(md5));
        executeTest("test calling on a ReducedRead BAM with " + model, spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing contamination down-sampling
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testContaminationDownsampling() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 --contamination_fraction_to_filter 0.20", 1,
                Arrays.asList("1f9071466fc40f4c6a0f58ac8e9135fb"));
        executeTest("test contamination_percentage_to_filter 0.20", spec);
    }


}
