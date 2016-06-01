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

import htsjdk.samtools.util.TestUtil;
import org.apache.commons.io.FileUtils;
import org.apache.log4j.Level;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

public class ValidateVariantsIntegrationTest extends WalkerTest {

    protected static final String EMPTY_MD5 = "d41d8cd98f00b204e9800998ecf8427e";
    protected static final String DEFAULT_REGION = "1:10001292-10001303";


    public static String baseTestString(final String file, String type) {
        return baseTestString(file,type,DEFAULT_REGION,b36KGReference);
    }

    public static String baseTestString(String file, String type, String region, String reference) {
        final String typeArgString = type.startsWith("-") ? " --validationTypeToExclude " + type.substring(1) : excludeValidationTypesButString(type);
        return "-T ValidateVariants -R " + reference + " -L " + region + " --variant:vcf " + privateTestDir + file + typeArgString;
    }

    private static String excludeValidationTypesButString(String type) {
        if (type == "ALL")
            return "";
        final ValidateVariants.ValidationType vtype = ValidateVariants.ValidationType.valueOf(type);
        final StringBuilder sbuilder = new StringBuilder();
        for (final ValidateVariants.ValidationType t : ValidateVariants.ValidationType.CONCRETE_TYPES)
            if (t != vtype)
                sbuilder.append(" --validationTypeToExclude " + t.toString());
        return sbuilder.toString();
    }

    @Test
    public void testGoodFile() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleGood.vcf", "ALL"),
                0,
                Arrays.asList(EMPTY_MD5)
        );

        executeTest("test good file", spec);
    }
    
    @Test
    public void testBadRefBase1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "REF"),
                0,
                UserException.FailsStrictValidation.class
        );

        executeTest("test bad ref base #1", spec);
    }

    @Test
    public void testBadRefBase2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad2.vcf", "REF"),
                0,
                UserException.FailsStrictValidation.class
        );

        executeTest("test bad ref base #2", spec);
    }

    @Test
    public void testBadChrCount1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "CHR_COUNTS"),
                0,
                UserException.FailsStrictValidation.class
        );

        executeTest("test bad chr counts #1", spec);
    }

    @Test
    public void testBadChrCount2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad2.vcf", "CHR_COUNTS"),
                0,
                UserException.FailsStrictValidation.class
        );

        executeTest("test bad chr counts #2", spec);
    }

    @Test
    public void testBadID() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "IDS") + " --dbsnp " + b36dbSNP129,
                0,
                UserException.FailsStrictValidation.class
        );
        executeTest("test bad RS ID", spec);
    }

    @Test
    public void testBadAllele() {
        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString("validationExampleBad.vcf", "ALLELES"),
            0,
            UserException.FailsStrictValidation.class
        );

        executeTest("test bad alt allele", spec);
    }

    @Test
    public void testBadAllele2() {
        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString("validationExampleBad3.vcf", "REF"),
            0,
            UserException.FailsStrictValidation.class
        );

        executeTest("test bad ref allele in deletion", spec);
    }

    @Test
    public void testNoValidation() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationExampleBad.vcf", "-ALL"),
                0,
                Arrays.asList(EMPTY_MD5)
        );

        executeTest("test no validation", spec);
    }

    @Test
    public void testComplexEvents() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("complexEvents.vcf", "ALL", DEFAULT_REGION, b37KGReference),
                0,
                Arrays.asList(EMPTY_MD5)
        );

        executeTest("test validating complex events", spec);
    }

    @Test(description = "Checks out of order header contigs")
    public void testOutOfOrderHeaderContigs() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("complexEvents-outOfOrder.vcf", "ALL", DEFAULT_REGION, b37KGReference),
                0,
                Arrays.asList(EMPTY_MD5));
        executeTest("test out of order header contigs", spec);
    }

    @Test(description = "Fixes '''bug''' reported in story https://www.pivotaltracker.com/story/show/68725164")
    public void testUnusedAlleleFix() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationUnusedAllelesBugFix.vcf","-ALLELES","1:1-739000",b37KGReference),0,Arrays.asList(EMPTY_MD5));
        executeTest("test unused allele bug fix", spec);
    }

    @Test(description = "Checks '''bug''' reported in story https://www.pivotaltracker.com/story/show/68725164")
    public void testUnusedAlleleError() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("validationUnusedAllelesBugFix.vcf","ALLELES","1:1-739000",b37KGReference),0, UserException.FailsStrictValidation.class);
        executeTest("test unused allele bug fix", spec);
    }

    @Test(description = "Checks '''bug''' reported in issue https://github.com/broadinstitute/gsa-unstable/issues/963")
    public void testLargeReferenceAlleleError() throws IOException {
        // Need to see log INFO messages
        Level level = logger.getLevel();
        logger.setLevel(Level.INFO);

        File logFile = createTempFile("testLargeReferenceAlleleError.log", ".tmp");
        String logFileName = logFile.getAbsolutePath();

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("longAlleles.vcf", "ALL", "1", b37KGReference) + " -log " + logFileName,
                0, Arrays.asList(EMPTY_MD5));
        executeTest("test long reference allele bug error", spec);

        // Make sure the "reference allele too long" message is in the log
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains(ValidateVariants.REFERENCE_ALLELE_TOO_LONG_MSG));
        
        // Set the log level back
        logger.setLevel(level);
    }

    @Test(description = "Checks '''bug''' is fixed, reported in issue https://github.com/broadinstitute/gsa-unstable/issues/963")
    public void testLargeReferenceAlleleFix() throws IOException {
        // Need to see log INFO messages
        Level level = logger.getLevel();
        logger.setLevel(Level.INFO);

        File logFile = createTempFile("testLargeReferenceAllele.log", ".tmp");
        String logFileName = logFile.getAbsolutePath();

        // expand window for the large reference allele
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("longAlleles.vcf","ALL","1",b37KGReference) + " --reference_window_stop 208 -log " + logFileName,
                0, Arrays.asList(EMPTY_MD5));
        executeTest("test long reference allele bug fix", spec);

        // Make sure the "reference allele too long" message is not in the log
        Assert.assertFalse(FileUtils.readFileToString(logFile).contains(ValidateVariants.REFERENCE_ALLELE_TOO_LONG_MSG));

        // All of the validation tests have passed since UserException.FailsStrictValidation is not thrown.

        // Set the log level back
        logger.setLevel(level);
    }

    @Test(description = "Checks '''issue''' reported in issue https://github.com/broadinstitute/gsa-unstable/issues/964")
    public void testWrongContigHeaderLengthError()  {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("longAlleles-wrongLength.vcf", "ALL", "1", b37KGReference),
                0, UserException.IncompatibleSequenceDictionaries.class);
        executeTest("test wrong header contig length error", spec);
    }

    @Test
    public void testAllowWrongContigHeaderLengthDictIncompat()  {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("longAlleles-wrongLength.vcf", "ALL", "1", b37KGReference) + "  --reference_window_stop 208 -U ALLOW_SEQ_DICT_INCOMPATIBILITY ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("test to allow wrong header contig length, not checking dictionary incompatibility", spec);
    }

    @Test
    public void testAllowWrongContigHeaderLength()  {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("longAlleles-wrongLength.vcf", "ALL", "1", b37KGReference) + "  --reference_window_stop 208 -U ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("test to allow wrong header contig length, no compatibility checks", spec);
    }

    @Test
    public void testGoodGvcf()  {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12891.AS.chr20snippet.g.vcf", "ALL", "20:10433000-10437000", b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList("d41d8cd98f00b204e9800998ecf8427e"));
        executeTest("tests correct gvcf", spec);
    }

    @Test
    public void testGoodGvcfExcludingAlleles()  {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12891.AS.chr20snippet.g.vcf", "-ALLELES", "20:10433000-10437000", b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList("d41d8cd98f00b204e9800998ecf8427e"));
        executeTest("tests correct gvcf", spec);
    }


    @Test(expectedExceptions = RuntimeException.class )
    public void testBadGvcfMissingNON_REF()  {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12891.AS.chr20snippet.BAD_MISSING_NON_REF.g.vcf", "-ALLELES", "20:10433000-10437000", b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("tests capture of missing NON_REF allele", spec);
    }

    @Test(expectedExceptions = RuntimeException.class )
    public void testBadGvcfRegions() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("diploid-gvcf.bad-IncompleteRegion.vcf", "-ALLELES", "20:10433000-10437000", b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("tests capture of non-complete region", spec);
    }

   @Test(expectedExceptions = RuntimeException.class )
    public void testNonOverlappingRegions()  {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12891.AS.chr20snippet_BAD_INCOMPLETE_REGION.g.vcf", "-ALLELES", "Y:4966254-4967190", b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("tests capture of non-complete region", spec);
    }

    @Test
    public void testNonOverlappingRegionsBP_RESOLUTION()  {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("gvcf.basepairResolution.vcf", "-ALLELES", "20:10000000-10010000", b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("tests capture of non-complete region, on BP_RESOLUTION gvcf", spec);
    }
    @Test
    public void testCorrectCreationOfBlocks() throws IOException {
        final File tempDir = IOUtils.tempDir("RefBlocks", "test", new File(privateTestDir));
        tempDir.mkdir();
        tempDir.deleteOnExit();
        final File output = File.createTempFile("RefBlocks", ".g.vcf", tempDir);
        String baseIntervals = " 1:1-100 -L 5:1-200 ";
        String intervalString = " -L " + baseIntervals;
        final WalkerTestSpec hc = new WalkerTestSpec("-T HaplotypeCaller " + intervalString + " -I " + privateTestDir + "NA12878.4.snippet.bam " +
                " -R /humgen/1kg/reference/human_g1k_v37_decoy.fasta -ERC GVCF -o " + output, Collections.singletonList(EMPTY_MD5));
        executeTest("running hc", hc);

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(tempDir.getName() + "/" + output.getName(), "-ALLELES", baseIntervals, b37KGReference) + " -gvcf  --reference_window_stop 208 -U ",
                0, Collections.singletonList(EMPTY_MD5));
        executeTest("testing the correct creation of reference blocks", spec);

        TestUtil.recursiveDelete(tempDir);
    }

}
