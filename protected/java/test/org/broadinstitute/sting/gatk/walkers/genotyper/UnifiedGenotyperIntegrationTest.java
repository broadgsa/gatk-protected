/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

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
                Arrays.asList("dff4412a074940d26994f9552476b209"));
        executeTest("test SingleSample Pilot2", spec);
    }

    @Test
    public void testMultipleSNPAlleles() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm BOTH --dbsnp " + b37dbSNP129 + " -I " + privateTestDir + "multiallelic.snps.bam -o %s -L " + privateTestDir + "multiallelic.snps.intervals", 1,
                Arrays.asList("b41b95aaa2c453c9b75b3b29a9c2718e"));
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
                Arrays.asList("1e61de694b51d7c0f26da5179ee6bb0c"));
        executeTest("test reverse trim", spec);
    }

    @Test
    public void testMismatchedPLs() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -I " + privateTestDir + "mismatchedPLs.bam -o %s -L 1:24020341", 1,
                Arrays.asList("935ee705ffe8cc6bf1d9efcceea271c8"));
        executeTest("test mismatched PLs", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "e6e33f0ebabab027eabed51fe9a08da9";

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
        testHeterozosity( 0.01, "bdc8760d7ae1e01c0510b12c1e6fcfa3" );
    }

    @Test
    public void testHeterozyosity2() {
        testHeterozosity( 1.0 / 1850, "f508f06a47305e11e62776615cb14fe3" );
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
                Arrays.asList("13d91059f58fb50a07a6a34b9438a45b"));

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
                Arrays.asList("07d8b77a5f6697f3a47a4f1efb0dcf50"));

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
                Arrays.asList("0f026d2e568172cf32813cc54ea7ba23"));

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
                Arrays.asList("e7ad858e9d6617534761918561f3ed4c"));

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
                 Arrays.asList("39c7a813fd6ee82d3604f2a868b35b2a"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    @Test
    public void testWithIndelAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("9430fe36789a791fcff6162f768ae563"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec);
    }

    @Test
    public void testWithIndelAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles "
                        + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("8d8dbf483526b0b309f5728619a74a86"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in and emitting all sites", spec);
    }

    @Test
    public void testMultiSampleIndels1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("3d3c5691973a223209a1341272d881be"));
        List<File> result = executeTest("test MultiSample Pilot1 CEU indels", spec1).getFirst();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + result.get(0).getAbsolutePath() + " -I " + validationDataLocation +
                        "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("b6c1d5cd28ff584c5f5037afef4e883a"));
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
                Arrays.asList("db3026c49a3de7a5cb9a3d77635d0706"));
        executeTest("test minIndelFraction 0.0", spec);
    }

    @Test
    public void testMinIndelFraction25() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.25", 1,
                Arrays.asList("7ab8e5ee15ab98d6756b0eea0f4d3798"));
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
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "testWithNs.bam -o %s -L 8:141813600-141813700 -out_mode EMIT_ALL_SITES", 1,
                Arrays.asList("4d36969d4f8f1094f1fb6e7e085c19f6"));
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
                Arrays.asList("092e42a712afb660ec79ff11c55933e2"));
        executeTest("test calling on a ReducedRead BAM", spec);
    }

    @Test
    public void testReducedBamSNPs() {
        testReducedCalling("SNP", "c0de74ab8f4f14eb3a2c5d55c200ac5f");
    }

    @Test
    public void testReducedBamINDELs() {
        testReducedCalling("INDEL", "1c9aaf65ffaa12bb766855265a1c3f8e");
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
