/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
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
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class SelectVariantsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T SelectVariants -R " + b36KGReference + " -L 1 -o %s --no_cmdline_in_header" + args;
    }

    private static final String sampleExclusionMD5 = "eea22fbf1e490e59389a663c3d6a6537";
    private static final String invertSelectionMD5 = "831bc0a5a723b0681a910d668ff3757b";

    @Test
    public void testDiscordanceNoSampleSpecified() {
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("954415f84996d27b07d00855e96d33a2")
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
                Arrays.asList("125d1c9fa111cd38dfa2ff3900f16b57")
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
                Arrays.asList("ca1b5226eaeaffb78d4abd9d2ee10c43")
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
                Arrays.asList("4386fbb258dcef4437495a37f5a83c53")
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
                Arrays.asList("4386fbb258dcef4437495a37f5a83c53")
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
                Arrays.asList("44e77cea624cfff2b8acc3a4b30485cb")    // should yield empty vcf because the foo!=0 will yield complete expression false
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
                Arrays.asList("1f5c72951a35667c4bdf1be153787e27")
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
                Arrays.asList("875d7e00ac8081e87ab9fb1b20c83677")
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
                Arrays.asList(sampleExclusionMD5)
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
                Arrays.asList(sampleExclusionMD5)
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
                Arrays.asList("946e7f2e0ae08dc0e931c1634360fc46")
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
                Arrays.asList("2c50ab2ae96fae40bfc2b8398fc5e54e")
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
                Arrays.asList("2c50ab2ae96fae40bfc2b8398fc5e54e")
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
                Arrays.asList("fa5f3eb4f0fc5cedc93e6c519c0c8bcb")
        );

        executeTest("testMinIndelLengthSelection--" + testFile, spec);
    }

    @Test
    public void testUsingDbsnpName() {
        String testFile = privateTestDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("a554459c9ccafb9812ff6d8c06c11726")
        );

        executeTest("testUsingDbsnpName--" + testFile, spec);
    }

    @Test
    public void testRemoveMLE() {
        String testFile = privateTestDir + "vcfexample.withMLE.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("a554459c9ccafb9812ff6d8c06c11726")
        );

        executeTest("testRemoveMLE--" + testFile, spec);
    }

    @Test
    public void testKeepOriginalAC() {
        String testFile = privateTestDir + "vcfexample.loseAlleleInSelection.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --keepOriginalAC -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("ad7e8b25e431a3229a78cec063876559")
        );

        executeTest("testKeepOriginalAC--" + testFile, spec);
    }

    @Test
    public void testKeepOriginalACAndENV() {
        String testFile = privateTestDir + "vcfexample.loseAlleleInSelection.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --keepOriginalAC -env -trimAlternates -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("4695c99d96490ed4e5b1568c5b52dea6")
        );

        executeTest("testKeepOriginalACAndENV--" + testFile, spec);
    }

    @Test
    public void testKeepOriginalDP() {
        String testFile = privateTestDir + "CEUtrioTest.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --keepOriginalDP -R " + b37KGReference + " -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("e897097a47aee5516dc4f1c0b9d69037")
        );

        executeTest("testKeepOriginalDP--" + testFile, spec);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() {
        String testFile = privateTestDir + "selectVariants.onePosition.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -select 'KG_FREQ < 0.5' --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("44f7c47395ca5b2afef5313f592c8cea")
        );

        executeTest("testMultipleRecordsAtOnePosition--" + testFile, spec);
    }

    @Test
    public void testNoGTs() {
        String testFile = privateTestDir + "vcf4.1.example.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("ef3c5f75074a5dd2b2cd2715856a2542")
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
                Arrays.asList("69862fb97e8e895fe65c7abb14b03cee")
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
                Arrays.asList("8fe7cdca8638461909262cb0769b2527")
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
        WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("9df942000eb18b12d9008c7d9b5c4178"));
        executeTest("testAlleleTrimming", spec);
    }

    @DataProvider(name="unusedAlleleTrimmingProvider")
    public Object[][] unusedAlleleTrimmingProvider() {
        return new Object[][] {
                { privateTestDir+"forHardLeftAlignVariantsTest.vcf", "-trimAlternates", "9df942000eb18b12d9008c7d9b5c4178"},
                { privateTestDir+"forHardLeftAlignVariantsTest.vcf", "", "981b757e3dc6bf3864ac7e493cf9d30d"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT", "8ded359dd87fd498ff38736ea0fa4c28"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT -env", "a7e7288dcd779cfac6983069de45b79c"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT -trimAlternates", "2e726d06a8d317199e8dda74691948a3"},
                { privateTestDir+"multi-allelic-ordering.vcf", "-sn SAMPLE-CC -sn SAMPLE-CT -env -trimAlternates", "1e5585f86c347da271a79fbfc61ac849"}
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
     * Test inverting the variant selection criteria
     */
    @Test
    public void testInvertSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 20000' --invert_selection --variant " + testfile),
                1,
                Arrays.asList(invertSelectionMD5)
        );
        spec.disableShadowBCF();
        executeTest("testInvertSelection--" + testfile, spec);
    }

    /**
     * Test inverting the variant JEXL selection criteria
     */
    @Test
    public void testInvertJexlSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP >= 20000'--variant " + testfile),
                1,
                Arrays.asList(invertSelectionMD5)
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
                Arrays.asList("2c50ab2ae96fae40bfc2b8398fc5e54e")
        );
        spec.disableShadowBCF();
        executeTest("testKeepSelectionID--" + testFile, spec);
    }

    /**
     * Test excluding variants with IDs
     */
    @Test
    public void testExcludeSelectionID() {
        String testFile = privateTestDir + "complexExample1.vcf";
        String idFile = privateTestDir + "complexExample1.vcf.id";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -xlIDs " + idFile + " --variant " + testFile),
                1,
                Arrays.asList("77514a81233e1bbc0f5e47b0fb76a89a")
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
                Arrays.asList("fa5f3eb4f0fc5cedc93e6c519c0c8bcb")
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
                Arrays.asList("406243096074a417d2aa103bd3d13e01"));

        executeTest("testMendelianViolationSelection--" + testFile, spec);
    }

    @Test
    public void testInvertMendelianViolationSelection() {
        String testFile = privateTestDir + "CEUtrioTest.vcf";
        String pedFile = privateTestDir + "CEUtrio.ped";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R "+b37KGReference + " -mv -mvq 0 -inv_mv --variant  " + testFile + " -ped " + pedFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("35921fb2dedca0ead83027a66b725794"));

        executeTest("testMendelianViolationSelection--" + testFile, spec);
    }
}
