/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class SelectVariantsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T SelectVariants -R " + b36KGReference + " -L 1 -o %s --no_cmdline_in_header" + args;
    }

    private static final String SAMPLE_EXCLUSION_MD5 = "2e52f21e7dcc67151a51630807a4eef2";
    private static final String INVERT_SELECTION_MD5 = "26d192b868746ab14133f145ae812e7c";
    private static final String MAX_FILTERED_GT_SELECTION_MD5 = "f83ac0deb7a8b022d6d40a85627a71ec";
    private static final String MIN_FILTERED_GT_SELECTION_MD5 = "346620b7a5d66dabf89d3f42d6e27db7";
    private static final String NO_CALL_FILTERING_KEEP_ONE = "6e2401190c5ada6a3bed2640c068f43b";
    private static final String NO_CALL_FILTERING_KEEP_TWO =  "6bced1ab6a3d58f1fd905b7f601987a3";

    @Test
    public void testDiscordanceNoSampleSpecified() {
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("9e08f761d2ba9a2bae9c279701aabc70")
        );
        spec.disableShadowBCF();

        executeTest("testDiscordanceNoSampleSpecified--" + testFile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = privateTestDir + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C --variant " + testfile),
                1,
                Arrays.asList("792962a5cc830e86dfc89caffbda1707")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("c9aa80cabf036a268a032a61d398cdd5")
        );
        spec.disableShadowBCF();

        executeTest("testDiscordance--" + testFile, spec);
    }

    @Test
    public void testComplexSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile),
                1,
                Arrays.asList("8090c349d12549b437a80e29c285fdd5")
        );
        spec.disableShadowBCF();
        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testComplexSelectionWithNonExistingSamples() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES -sn A -se '[CDH]' -sn Z -sn T -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile),
                1,
                Arrays.asList("8090c349d12549b437a80e29c285fdd5")
        );
        spec.disableShadowBCF();
        executeTest("testComplexSelectionWithNonExistingSamples--" + testfile, spec);
    }

    @Test
    public void testNonExistingFieldSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -env -ef -select 'foo!=0||DP>0' --variant " + testfile),
                1,
                Arrays.asList("e7ec1f2c8077d07b54721e68b603d42c")    // should yield empty vcf because the foo!=0 will yield complete expression false
        );
        spec.disableShadowBCF();
        executeTest("testNonExistingSelection--" + testfile, spec);
    }

    /**
     * Test excluding samples from file and sample name
     */
    @Test
    public void testSampleExclusionFromFileAndSeparateSample() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -xl_sn A -xl_sf " + samplesFile + " --variant " + testfile,
                1,
                Arrays.asList("30aabc865634bf887cad0c02cdcde042")
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusionFromFileAndSeparateSample--" + testfile, spec);
    }

    /**
     * Test excluding samples from file
     */
    @Test
    public void testSampleExclusionJustFromFile() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -xl_sf " + samplesFile + " --variant " + testfile,
                1,
                Arrays.asList("1afba8d53094bdef63db1e39d52be5aa")
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusionJustFromFile--" + testfile, spec);
    }

    /**
     * Test excluding samples from expression
     */
    @Test
    public void testSampleExclusionJustFromExpression() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -xl_se '[CDH]' --variant " + testfile,
                1,
                Arrays.asList(SAMPLE_EXCLUSION_MD5)
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusionJustFromExpression--" + testfile, spec);
    }

    /**
     * Test excluding samples from negation expression
     */
    @Test
    public void testSampleExclusionJustFromNegationExpression() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -se '[^CDH]' --variant " + testfile,
                1,
                Arrays.asList(SAMPLE_EXCLUSION_MD5)
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusionJustFromRegexExpression--" + testfile, spec);
    }

    /**
     * Test including samples that are not in the VCF
     */
    @Test
    public void testSampleInclusionWithNonexistingSamples() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -sn A -sn Z -sn Q -sf " + samplesFile + " --variant " + testfile,
                1,
                UserException.BadInput.class
        );
        spec.disableShadowBCF();

        executeTest("testSampleInclusionWithNonexistingSamples--" + testfile, spec);
    }

    @Test
    public void testConcordance() {
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc "
                        + b37hapmapGenotypes + " --variant " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("24114c01b81fc0052ee36523ccd1d338")
        );
        spec.disableShadowBCF();

        executeTest("testConcordance--" + testFile, spec);
    }

    /**
     * Test including variant types.
     */
    @Test
    public void testVariantTypeSelection() {
        String testFile = privateTestDir + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -restrictAllelesTo MULTIALLELIC -selectType MIXED --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("41dda9f4b9ec9f9b0f3593b2cbd82cd5")
        );

        executeTest("testVariantTypeSelection--" + testFile, spec);
    }

    /**
     * Test excluding indels that are larger than the specified size
     */
    @Test
    public void testMaxIndelLengthSelection() {
        String testFile = privateTestDir + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -selectType INDEL --variant " + testFile + " -o %s --no_cmdline_in_header --maxIndelSize 2",
                1,
                Arrays.asList("41dda9f4b9ec9f9b0f3593b2cbd82cd5")
        );

        executeTest("testMaxIndelLengthSelection--" + testFile, spec);
    }

    /**
     * Test excluding indels that are smaller than the specified size
     */
    @Test
    public void testMinIndelLengthSelection() {
        String testFile = privateTestDir + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -selectType INDEL --variant " + testFile + " -o %s --no_cmdline_in_header --minIndelSize 2",
                1,
                Arrays.asList("ed9dc00d0551630a2eed9e81a2a357d3")
        );

        executeTest("testMinIndelLengthSelection--" + testFile, spec);
    }

    @Test
    public void testUsingDbsnpName() {
        String testFile = privateTestDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("86d97e682b2dccff75d079f3b5d17f4b")
        );

        executeTest("testUsingDbsnpName--" + testFile, spec);
    }

    @Test
    public void testRemoveMLE() {
        String testFile = privateTestDir + "vcfexample.withMLE.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("86d97e682b2dccff75d079f3b5d17f4b")
        );

        executeTest("testRemoveMLE--" + testFile, spec);
    }

    @Test
    public void testKeepOriginalAC() {
        String testFile = privateTestDir + "vcfexample.loseAlleleInSelection.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --keepOriginalAC -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("6f26cf5a7fd20682e1de193e5bb5f61f")
        );

        executeTest("testKeepOriginalAC--" + testFile, spec);
    }

    @Test
    public void testKeepOriginalACAndENV() {
        String testFile = privateTestDir + "vcfexample.loseAlleleInSelection.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --keepOriginalAC -env -trimAlternates -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("e0ac2b37387048bf51ac2914bdd2e178")
        );

        executeTest("testKeepOriginalACAndENV--" + testFile, spec);
    }

    @Test
    public void testKeepOriginalDP() {
        String testFile = privateTestDir + "CEUtrioTest.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --keepOriginalDP -R " + b37KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("9ad02f0df308eecb0634b3cd386956e9")
        );

        executeTest("testKeepOriginalDP--" + testFile, spec);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() {
        String testFile = privateTestDir + "selectVariants.onePosition.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -select 'KG_FREQ < 0.5' --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("bfbfefbd4a84b093ee0b63eab8cc1be9")
        );

        executeTest("testMultipleRecordsAtOnePosition--" + testFile, spec);
    }

    @Test
    public void testNoGTs() {
        String testFile = privateTestDir + "vcf4.1.example.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("c78a65b41edbdd386211042e8f65220b")
        );

        executeTest("testNoGTs--" + testFile, spec);
    }

    @Test
    public void testSelectFromMultiAllelic() {
        String testfile = privateTestDir + "multi-allelic.bi-allelicInGIH.vcf";
        String samplesFile = privateTestDir + "GIH.samples.list";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " -o %s --no_cmdline_in_header -sf " + samplesFile + " --excludeNonVariants -trimAlternates --variant " + testfile,
                1,
                Arrays.asList("c963ca96d543ecccab8055295d2a4dab")
        );
        executeTest("test select from multi allelic with excludeNonVariants --" + testfile, spec);
    }

    @Test
    public void testMultiAllelicAnnotationOrdering() {
        String testfile = privateTestDir + "multi-allelic-ordering.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " -o %s --no_cmdline_in_header " +
                        "-sn SAMPLE-CC -sn SAMPLE-CT -sn SAMPLE-CA --excludeNonVariants --variant " + testfile,
                1,
                Arrays.asList("7f5484a74ab648608228eafea96f8ad3")
        );
        executeTest("test multi allelic annotation ordering --" + testfile, spec);
    }

    @Test()
    public void testFileWithoutInfoLineInHeader() {
        testFileWithoutInfoLineInHeader("testFileWithoutInfoLineInHeader", IllegalStateException.class);
    }

    @Test()
    public void testFileWithoutInfoLineInHeaderWithOverride() {
        testFileWithoutInfoLineInHeader("testFileWithoutInfoLineInHeaderWithOverride", null);
    }

    private void testFileWithoutInfoLineInHeader(final String name, final Class expectedException) {
        final String testFile = privateTestDir + "missingHeaderLine.vcf";
        final String cmd = "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp "
                + testFile + " -o %s --no_cmdline_in_header"
                + (expectedException == null ? " -U LENIENT_VCF_PROCESSING" : "");
        WalkerTestSpec spec =
                expectedException != null
                        ? new WalkerTestSpec(cmd, 1, expectedException)
                        : new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();

        executeTest(name, spec);
    }

    @Test
    public void testInvalidJexl() {
        // NOTE: JexlEngine singleton construction in VariantContextUtils sets silent to false.
        // However VariantFiltration.initialize() sets setSilent(true) on the shared instance.
        // Just in case this test runs after a VariantFiltration in the same VM, always set silent back to false.
        htsjdk.variant.variantcontext.VariantContextUtils.engine.get().setSilent(false);
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants" +
                        " -R "+b37KGReference +
                        " -V "+privateTestDir+"ac0.vcf" +
                        " -select 'vc.getGenotype(\"FAKE_SAMPLE\").isHomRef()'" +
                        " -o %s",
                1,
                UserException.class);
        executeTest("InvalidJexl", spec);
    }

    @Test
    public void testAlleleTrimming() {
        final String testFile = privateTestDir + "forHardLeftAlignVariantsTest.vcf";
        final String cmd = "-T SelectVariants -R " + b37KGReference + " -sn NA12878 -env -trimAlternates "
                + "-V " + testFile + " -o %s --no_cmdline_in_header";
        WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("21d83006b012eeea84c6612976348d3c"));
        executeTest("testAlleleTrimming", spec);
    }

    @DataProvider(name="unusedAlleleTrimmingProvider")
    public Object[][] unusedAlleleTrimmingProvider() {
        return new Object[][] {
                { privateTestDir+"forHardLeftAlignVariantsTest.vcf", "-trimAlternates", "21d83006b012eeea84c6612976348d3c"},
                { privateTestDir+"forHardLeftAlignVariantsTest.vcf", "", "8fc0c8a7de6bb579e1534b936f844090"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT", "595392b623b0869f1d87e46edf3de122"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT -env", "bba873b8eeeb4c01199140c37deb6f6b"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT -trimAlternates", "93858f706dac876a8581f6b89bb85cc5"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT -env -trimAlternates", "5d831401367eb8b0ab49ffa34e0dd278"}
        };
    }

    @Test(dataProvider = "unusedAlleleTrimmingProvider")
    public void testUnusedAlleleTrimming(final String vcf, final String extraArgs, final String md5) {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants" +
                        " -R "+b37KGReference +
                        " -V "+vcf +
                        " -o %s --no_cmdline_in_header" +
                        " "+extraArgs,
                1,
                Arrays.asList(md5)
        );
        executeTest(String.format("testUnusedAlleleTrimming: (%s,%s)", new File(vcf).getName(), extraArgs), spec);
    }

    /**
     *  Test with an empty VCF file
     */
    @Test
    public void testEmptyVcfException(){
        String testfile = privateTestDir + "reallyEmpty.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants" +
                        " -R " + b36KGReference +
                        " -V " + testfile +
                        " -o %s --no_cmdline_in_header",
                1,
                UserException.CommandLineException.class
        );
        spec.disableShadowBCF();

        executeTest("testEmptyVcfException--" + testfile, spec);
    }

    /**
     * Test with a VCF file that is not a file
     */
    @Test
    public void testNotFileVcfException(){
        String testfile = privateTestDir;

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants" +
                        " -R " + b36KGReference +
                        " -V " + testfile +
                        " -o %s --no_cmdline_in_header",
                1,
                UserException.CouldNotReadInputFile.class
        );
        spec.disableShadowBCF();

        executeTest("testNotFileVcfException--" + testfile, spec);
    }

    /**
     * Test with a VCF file that does not exist
     */
    @Test
    public void testMissingVcfException(){
        String testfile = "test.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants" +
                        " -R " + b36KGReference +
                        " -V " + testfile +
                        " -o %s --no_cmdline_in_header",
                1,
                UserException.CouldNotReadInputFile.class
        );
        spec.disableShadowBCF();

        executeTest("testMissingVcfException--" + testfile, spec);
    }

    /**
     * Test inverting the variant selection criteria by the -invertSelect argument
     */
    @Test
    public void testInvertSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 20000' -invertSelect --variant " + testfile),
                1,
                Arrays.asList(INVERT_SELECTION_MD5)
        );
        spec.disableShadowBCF();
        executeTest("testInvertSelection--" + testfile, spec);
    }

    /**
     * Test inverting the variant selection criteria by inverting the JEXL expression logic following -select
     */
    @Test
    public void testInvertJexlSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP >= 20000'--variant " + testfile),
                1,
                Arrays.asList(INVERT_SELECTION_MD5)
        );
        spec.disableShadowBCF();
        executeTest("testInvertJexlSelection--" + testfile, spec);
    }

    /**
     * Test selecting variants with IDs
     */
    @Test
    public void testKeepSelectionID() {
        String testFile = privateTestDir + "complexExample1.vcf";
        String idFile = privateTestDir + "complexExample1.vcf.id";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -IDs " + idFile + " --variant " + testFile),
                1,
                Arrays.asList("c6632b63617162455f02670174a2322a")
        );
        spec.disableShadowBCF();
        executeTest("testKeepSelectionID--" + testFile, spec);
    }

    /**
     * Test excluding variants with IDs
     * Also tests --forceValidOutput flag, which changes the GQ from floats to ints to match
     * header spec.
     */
    @Test
    public void testExcludeSelectionID() {
        String testFile = privateTestDir + "complexExample1.vcf";
        String idFile = privateTestDir + "complexExample1.vcf.id";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -xlIDs " + idFile + " --variant " + testFile + " --forceValidOutput"),
                1,
                Arrays.asList("613826acb415f05bf288536701a87855")
        );
        spec.disableShadowBCF();
        executeTest("testExcludeSelectionID--" + testFile, spec);
    }

    /**
     * Test excluding variant types
     */
    @Test
    public void testExcludeSelectionType() {
        String testFile = privateTestDir + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -xlSelectType SNP --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("ed9dc00d0551630a2eed9e81a2a357d3")
        );

        executeTest("testExcludeSelectionType--" + testFile, spec);
    }

    @Test
    public void testMendelianViolationSelection() {
        String testFile = privateTestDir + "CEUtrioTest.vcf";
        String pedFile = privateTestDir + "CEUtrio.ped";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R "+b37KGReference + " -mv -mvq 0 --variant  " + testFile + " -ped " + pedFile + " -o %s  --no_cmdline_in_header",
                1,
                Arrays.asList("c68779547b28dfef39792598df8a93e9"));

        executeTest("testMendelianViolationSelection--" + testFile, spec);
    }

    @Test
    public void testInvertMendelianViolationSelection() {
        String testFile = privateTestDir + "CEUtrioTest.vcf";
        String pedFile = privateTestDir + "CEUtrio.ped";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R "+b37KGReference + " -mv -mvq 0 -invMv --variant  " + testFile + " -ped " + pedFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("0ac6fda76228080bdb39c0e698440718"));

        executeTest("testInvertMendelianViolationSelection--" + testFile, spec);
    }

    @Test
    public void testMaxFilteredGenotypesSelection() {
        String testfile = privateTestDir + "filteredSamples.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --maxFilteredGenotypes 1 -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList(MAX_FILTERED_GT_SELECTION_MD5)
        );
        spec.disableShadowBCF();
        executeTest("testMaxFilteredGenotypesSelection--" + testfile, spec);
    }

    @Test
    public void testMinFilteredGenotypesSelection() {
        String testfile = privateTestDir + "filteredSamples.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --minFilteredGenotypes 2 -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList(MIN_FILTERED_GT_SELECTION_MD5)
        );
        spec.disableShadowBCF();
        executeTest("testMinFilteredGenotypesSelection--" + testfile, spec);
    }

    @Test
    public void testMaxFractionFilteredGenotypesSelection() {
        String testfile = privateTestDir + "filteredSamples.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --maxFractionFilteredGenotypes 0.4 -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList(MAX_FILTERED_GT_SELECTION_MD5)
        );
        spec.disableShadowBCF();
        executeTest("testMaxFractionFilteredGenotypesSelection--" + testfile, spec);
    }

    @Test
    public void testMinFractionFilteredGenotypesSelection() {
        String testfile = privateTestDir + "filteredSamples.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --minFractionFilteredGenotypes 0.6 -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList(MIN_FILTERED_GT_SELECTION_MD5)
        );
        spec.disableShadowBCF();
        executeTest("testMinFractionFilteredGenotypesSelection--" + testfile, spec);
    }

    @Test
    public void testSetFilteredGtoNocall() {
        String testfile = privateTestDir + "filteredSamples.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --setFilteredGtToNocall -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("410c6b7bb62fc43bb41eee627670f757")
        );

        spec.disableShadowBCF();
        executeTest("testSetFilteredGtoNocall--" + testfile, spec);
    }

    @Test
    public void testSetFilteredGtoNocallUpdateInfo() {
        String testfile = privateTestDir + "selectVariantsInfoField.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --setFilteredGtToNocall --removeUnusedAlternates --excludeNonVariants -R " + b37KGReference + " --variant " +
                        testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("349136d92f915f8c7ba8a2f92e51d6b7"));
        executeTest("testSetFilteredGtoNocallUpdateInfo", spec);
    }

    @Test
    public void testSACSimpleDiploid() {
        String testfile = privateTestDir + "261_S01_raw_variants_gvcf.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header -trimAlternates",
                1,
                Arrays.asList("d3bb7ea37a7c9dce8b34bf2020961619"));
        spec.disableShadowBCF();
        executeTest("testSACSimpleDiploid", spec);
    }

    @Test
    public void testSACDiploid() {
        String testfile = privateTestDir + "diploid-multisample-sac.g.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header -sn NA12891 -trimAlternates",
                1,
                Arrays.asList("67a92b4d4174ff41f6f88ddf5ab6e422"));
        spec.disableShadowBCF();
        executeTest("testSACDiploid", spec);
    }

    @Test
    public void testSACNonDiploid() {
        String testfile = privateTestDir + "tetraploid-multisample-sac.g.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header -sn NA12891 -trimAlternates",
                1,
                ReviewedGATKException.class);
        spec.disableShadowBCF();
        executeTest("testSACNonDiploid", spec);
    }

    @Test
    public void testMaxNoCall1() {
        final String testfile = privateTestDir + "vcfexample.forNoCallFiltering.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " --variant " + testfile + " -o %s --no_cmdline_in_header --maxNOCALLnumber 1",
                1,
                Arrays.asList(NO_CALL_FILTERING_KEEP_ONE));
        spec.disableShadowBCF();

        executeTest("testMaxNoCall1", spec);
    }

    @Test
    public void testMaxNoCall0_25()  {
        final String testfile = privateTestDir + "vcfexample.forNoCallFiltering.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " --variant " + testfile + " -o %s --no_cmdline_in_header --maxNOCALLfraction 0.25",
                1,
                Arrays.asList(NO_CALL_FILTERING_KEEP_ONE));
        spec.disableShadowBCF();

        executeTest("testMaxNoCall0_25", spec);
    }

    @Test
    public void testMaxNoCall2() {
        final String testfile = privateTestDir + "vcfexample.forNoCallFiltering.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " --variant " + testfile + " -o %s --no_cmdline_in_header --maxNOCALLnumber 2",
                1,
                Arrays.asList(NO_CALL_FILTERING_KEEP_TWO));
        spec.disableShadowBCF();

        executeTest("testMaxNoCall2", spec);
    }

    @Test
    public void testMaxNoCall0_5() {
        final String testfile = privateTestDir + "vcfexample.forNoCallFiltering.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " --variant " + testfile + " -o %s --no_cmdline_in_header --maxNOCALLfraction 0.5",
                1,
                Arrays.asList(NO_CALL_FILTERING_KEEP_TWO));
        spec.disableShadowBCF();

        executeTest("testMaxNoCall0_5", spec);
    }
}
