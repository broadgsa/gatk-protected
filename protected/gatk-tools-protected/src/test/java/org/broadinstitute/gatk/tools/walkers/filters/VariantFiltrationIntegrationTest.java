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

package org.broadinstitute.gatk.tools.walkers.filters;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b36KGReference;
    }


    @Test
    public void testNoAction() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("cc2398d5a69e4d3ead25e70e8c786dfa"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("92d9d04451f8fa510dd34db837a8314b"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask " + privateTestDir + "vcfexample2.vcf --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("e1e65e8d66e02715a73137a2cbb2d2d1"));
        executeTest("test mask all", spec1);
    }

    @Test
    public void testMask2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF " + privateTestDir + "vcfMask.vcf --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("270e06d7e83e40e6a273668be890fe56"));
        executeTest("test mask some", spec2);
    }

    @Test
    public void testMask3() {
        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName foo -maskExtend 10 --mask:VCF " + privateTestDir + "vcfMask.vcf --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("e4d1c66da07d1d0ccfe674ed2b013849"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testMaskReversed() {
        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName outsideGoodSites -filterNotInMask --mask:BED " + privateTestDir + "goodMask.bed --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("42e6dd3f30618753dc7411509c83c14f"));
        executeTest("test filter sites not in mask", spec3);
    }

    @Test
    public void testIllegalFilterName() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName 'foo < foo' --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                UserException.class);
        executeTest("test illegal filter name", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("f635fa722b962a6cdf429f84b0aa8b12"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("0a7cdf3117b4e8bcfcc41fde34ef0a9e"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7f798982ec0d90aa32484e35a89c7b4f"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testInvertFilter() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000 --invertFilterExpression", 1,
                Arrays.asList("01dc0130f85ef609c6f4f6e7235ef644"));
        executeTest("test inversion of selection of filter with separate names #2", spec);
    }

    @Test
    public void testInvertJexlFilter() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance >= 0.7' --filterName FSF -filter 'FisherStrand != 1.4' --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("23b6c0f2f682b744501ffdd1141d62ca")); // Differs from testInvertFilter() because their VCF header FILTER description uses the -filter argument. Their filter statuses are identical.
        executeTest("test inversion of selection of filter via JEXL with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("fbb22221a531d68a6394042ebf247aa0"));
        executeTest("test genotype filter #1", spec1);
    }

    @Test
    public void testGenotypeFilters2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'isHomVar == 1' -G_filterName foo --variant " + privateTestDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("99eacba26598ab31501e09604ff88be2"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo --variant:VCF " + privateTestDir + "twoDeletions.vcf", 1,
                Arrays.asList("ff84050b969d1068eb9ef39058d0e007"));
        executeTest("test deletions", spec);
    }

    @Test
    public void testUnfilteredBecomesFilteredAndPass() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                    + " --filterExpression 'FS > 60.0' --filterName SNP_FS -V " + privateTestDir + "unfilteredForFiltering.vcf", 1,
                Arrays.asList("3f828fd8768ec0e1708ae509fdb5e998"));
        executeTest("testUnfilteredBecomesFilteredAndPass", spec);
    }

    @Test
    public void testFilteringDPfromINFO() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --filterExpression 'DP < 8' --filterName lowDP -V " + privateTestDir + "filteringDepthInFormat.vcf", 1,
                Arrays.asList("397aa0a403d8c2aa47beb14788d81517"));
        executeTest("testFilteringDPfromINFO", spec);
    }

    @Test
    public void testFilteringDPfromFORMAT() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --genotypeFilterExpression 'DP < 8' --genotypeFilterName lowDP -V " + privateTestDir + "filteringDepthInFormat.vcf", 1,
                Arrays.asList("8c82cc0ec68b2397512d7bf255c17730"));
        executeTest("testFilteringDPfromFORMAT", spec);
    }

    // The current htsjdk implementation of JEXL matching on genotype fields is buggy. When the filter uses an
    // annotation that is present in both FORMAT and INFO, and the FORMAT value is missing, the current code (Dec 2016)
    // will look up the INFO value. Here we use a made-up annotation Z instead of DP to avoid having to rig the test
    // so that the INFO value will give the same matching results as the FORMAT value.
    @Test
    public void testFilteringZfromFORMATWithMissing() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --genotypeFilterExpression 'Z < 10' --genotypeFilterName lowZ -V " + privateTestDir + "filteringDepthInFormatWithMissing.vcf", 1,
                Arrays.asList("fa1bf49c94a36fa79c99f741a9fad14e"));
        executeTest("testFilteringDPfromFORMATWithMissing", spec);
    }

    // Same comment as above.
    @Test
    public void testFilteringZfromFORMATAndFailMissing() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --missingValuesInExpressionsShouldEvaluateAsFailing --genotypeFilterExpression 'Z < 10' --genotypeFilterName lowZ -V " + privateTestDir + "filteringDepthInFormatWithMissing.vcf", 1,
                Arrays.asList("0bf39f791e1a4d9a935fca86407b1c0f"));
        executeTest("testFilteringDPfromFORMATAndFailMissing", spec);
    }

    @Test
    public void testInvertGenotypeFilterExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --genotypeFilterExpression 'DP < 8' --genotypeFilterName highDP -V " + privateTestDir + "filteringDepthInFormat.vcf --invertGenotypeFilterExpression", 1,
                Arrays.asList("cc8eab3f7cacc6d7effc9df9ba98210a"));
        executeTest("testInvertGenotypeFilterExpression", spec);
    }

    @Test
    public void testInvertJexlGenotypeFilterExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --genotypeFilterExpression 'DP >= 8' --genotypeFilterName highDP -V " + privateTestDir + "filteringDepthInFormat.vcf", 1,
                Arrays.asList("bb22c3ca06bf22ca244504b2b3e27259")); // Differs from testInvertFilter because FILTER description uses the -genotypeFilterExpression argument
        executeTest("testInvertJexlGenotypeFilterExpression", spec);
    }

    @Test
    public void testSetFilteredGtoNocall() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " --genotypeFilterExpression 'DP < 8' --genotypeFilterName lowDP -V " + privateTestDir + "filteringDepthInFormat.vcf --setFilteredGtToNocall", 1,
                Arrays.asList("7a9d29b6a002aa4fba11bcaf24401a70"));
        executeTest("testSetFilteredGtoNocall", spec);
    }


    @Test
    public void testSetFilteredGtoNocallUpdateInfo() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b37KGReference
                        + " -G_filter 'GQ < 20' -G_filterName lowDP -G_filter 'DP<10' -G_filterName lowGQ -V " + privateTestDir + "variantFiltrationInfoField.vcf --setFilteredGtToNocall",
                1,
                Arrays.asList("fd39d39082d8a1535003a9d216905e73"));
        executeTest("testSetFilteredGtoNocallUpdateInfo", spec);
    }

    @Test
    public void testSetVcfFilteredGtoNocall() {
        String testfile = privateTestDir + "filteredSamples.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants --setFilteredGtToNocall -R " + b37KGReference + " --variant " + testfile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("fb3030518f7c4120989d37b2cca9abe6")
        );

        spec.disableShadowBCF();
        executeTest("testSetVcfFilteredGtoNocall--" + testfile, spec);
    }
}
