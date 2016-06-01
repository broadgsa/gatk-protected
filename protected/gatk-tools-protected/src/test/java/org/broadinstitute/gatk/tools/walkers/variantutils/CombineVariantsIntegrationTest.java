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
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Tests CombineVariants
 */
public class CombineVariantsIntegrationTest extends WalkerTest {
    //
    // TODO TODO TODO TODO TODO TODO TODO TODO
    // TODO TODO TODO TODO TODO TODO TODO TODO
    //
    // TODO WHEN THE HC EMITS VALID VCF HEADERS ENABLE BCF AND REMOVE lenientVCFProcessing ARGUMENTS
    //
    // TODO TODO TODO TODO TODO TODO TODO TODO
    // TODO TODO TODO TODO TODO TODO TODO TODO
    // TODO TODO TODO TODO TODO TODO TODO TODO
    //
    private static String baseTestString(String args) {
        return baseTestString(args, b36KGReference);
    }

    private static String baseTestString(String args, String ref) {
        return "-T CombineVariants --no_cmdline_in_header -L 1:1-50,000,000 -o %s -R " + ref + args;
        //return "-T CombineVariants --no_cmdline_in_header -L 1:1-50,000,000 -o %s -U LENIENT_VCF_PROCESSING -R " + b36KGReference + args;
    }

    private void cvExecuteTest(final String name, final WalkerTestSpec spec, final boolean parallel) {
        spec.disableShadowBCF();
        if ( parallel )
            executeTestParallel(name, spec);
        else
            executeTest(name, spec);
    }

    public void test1InOut(String file, String md5) {
        test1InOut(file, md5, "");
    }

    public void test1InOut(String file, String md5, String args) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1 -V:v1 " + validationDataLocation + file + args),
                 1,
                 Arrays.asList(md5));
         cvExecuteTest("testInOut1--" + file, spec, true);
    }

    public void combine2(String file1, String file2, String args, String md5) {
        combine2(file1, file2, args, md5, true);
    }

    public void combine2(String file1, String file2, String args, String md5, final boolean parallel) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 baseTestString(" -priority v1,v2 -V:v1 " + validationDataLocation + file1 + " -V:v2 "+ validationDataLocation + file2 + args),
                 1,
                 Arrays.asList(md5));
         cvExecuteTest("combine2 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec, parallel);
    }

    public void combineSites(String args, String md5) {
        String file1 = "1000G_omni2.5.b37.sites.vcf";
        String file2 = "hapmap_3.3.b37.sites.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants --no_cmdline_in_header -o %s -R " + b37KGReference
                        + " -L 1:1-10,000,000 -V:omni " + validationDataLocation + file1
                        + " -V:hm3 " + validationDataLocation + file2 + args,
                1,
                Arrays.asList(md5));
        cvExecuteTest("combineSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec, true);
    }

    public void combinePLs(String file1, String file2, String md5) {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T CombineVariants --no_cmdline_in_header -o %s -R " + b36KGReference + " -priority v1,v2 -V:v1 " + privateTestDir + file1 + " -V:v2 " + privateTestDir + file2,
                 1,
                 Arrays.asList(md5));
         cvExecuteTest("combine PLs 1:" + new File(file1).getName() + " 2:" + new File(file2).getName(), spec, true);
    }

    @Test public void test1SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "e3dbdfa14aefb2f6bd1213287d34a2e5", " -U LENIENT_VCF_PROCESSING"); }
    @Test public void test2SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "d727fab83b4265859c4a902f6e66ac3d", " -setKey foo -U LENIENT_VCF_PROCESSING"); }
    @Test public void test3SNP() { test1InOut("pilot2.snps.vcf4.genotypes.vcf", "42fc3d2c68415a61ff15e594a63d9349", " -setKey null -U LENIENT_VCF_PROCESSING"); }
    @Test public void testOfficialCEUPilotCalls() { test1InOut("CEU.trio.2010_03.genotypes.vcf.gz", "a3994d6145bb3813950939238db4c592"); } // official project VCF files in tabix format

    @Test public void test1Indel1() { test1InOut("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "e7fd959312e2aff0b4231963ee690aec"); }
    @Test public void test1Indel2() { test1InOut("CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "23439a1f0108b57a14e18efe9482cc88"); }

    @Test public void combineWithPLs() { combinePLs("combine.3.vcf", "combine.4.vcf", "27aa46cdb022be3959e7240a0d7ac794"); }

    @Test public void combineTrioCalls() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", "", "9bdda937754e1407183406808f560723"); } // official project VCF files in tabix format
    @Test public void combineTrioCallsMin() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "YRI.trio.2010_03.genotypes.vcf.gz", " -minimalVCF", "6344953a82a422115bd647ec1d696b94"); } // official project VCF files in tabix format
    @Test public void combine2Indels() { combine2("CEU.dindel.vcf4.trio.2010_06.indel.genotypes.vcf", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "51cf4543e46c8434e32c426f59507d1e"); }

    @Test public void combineSNPsAndIndels() { combine2("CEU.trio.2010_03.genotypes.vcf.gz", "CEU.dindel.vcf4.low_coverage.2010_06.indel.genotypes.vcf", "", "f9d1d7e6246f0ce9e493357d5b320323"); }

    @Test public void uniqueSNPs() {
        // parallelism must be disabled because the input VCF is malformed (DB=0) and parallelism actually fixes this which breaks the md5s
        //both of these files have the YRI trio and merging of duplicate samples without priority must be specified with UNSORTED merge type
        combine2("pilot2.snps.vcf4.genotypes.vcf", "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", " -genotypeMergeOptions UNSORTED", "5aece78046bfb7d6ee8dc4d551542e3a", true);
    }

    @Test public void omniHM3Union() { combineSites(" -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED", "0897efcc0046bd94760315838d4d0fa5"); }
    @Test public void omniHM3Intersect() { combineSites(" -filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED", "8b12b09a6ec4e3fde2352bbf82637f1e"); }

    @Test public void threeWayWithRefs() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:NA19240_BGI "+validationDataLocation+"NA19240.BGI.RG.vcf" +
                        " -V:NA19240_ILLUMINA "+validationDataLocation+"NA19240.ILLUMINA.RG.vcf" +
                        " -V:NA19240_WUGSC "+validationDataLocation+"NA19240.WUGSC.RG.vcf" +
                        " -V:denovoInfo "+validationDataLocation+"yri_merged_validation_data_240610.annotated.b36.vcf" +
                        " -setKey centerSet" +
                        " -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED" +
                        " -U LENIENT_VCF_PROCESSING" +
                        " -priority NA19240_BGI,NA19240_ILLUMINA,NA19240_WUGSC,denovoInfo" +
                        " -genotypeMergeOptions UNIQUIFY -L 1"),
                1,
                Arrays.asList("8b75e835ed19c06c358a2185cd0e14db"));
        cvExecuteTest("threeWayWithRefs", spec, true);
    }

    // complex examples with filtering, indels, and multiple alleles
    public void combineComplexSites(String args, String md5) {
        String file1 = "combine.1.vcf";
        String file2 = "combine.2.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants --no_cmdline_in_header -o %s -R " + b37KGReference
                        + " -V:one " + privateTestDir + file1
                        + " -V:two " + privateTestDir + file2 + args,
                1,
                Arrays.asList(md5));
        cvExecuteTest("combineComplexSites 1:" + new File(file1).getName() + " 2:" + new File(file2).getName() + " args = " + args, spec, true);
    }

    @Test public void complexTestFull() { combineComplexSites("", "80c3b7ba39c8a3f3511fc1ea61ecd4da"); }
    @Test public void complexTestMinimal() { combineComplexSites(" -minimalVCF", "e7da95bbcf3890a4debbfa07cbd646e5"); }
    @Test public void complexTestSitesOnly() { combineComplexSites(" -sites_only", "d7ee8da6ceee4dd1212122c2e9cab2a6"); }
    @Test public void complexTestSitesOnlyMinimal() { combineComplexSites(" -sites_only -minimalVCF", "d7ee8da6ceee4dd1212122c2e9cab2a6"); }

    @Test
    public void combineDBSNPDuplicateSites() {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T CombineVariants --no_cmdline_in_header -L 1:902000-903000 -o %s -R " + b37KGReference + " -V:v1 " + b37dbSNP132,
                 1,
                 Arrays.asList("b0d4b86702b44fc4faa527c34adf6239"));
         cvExecuteTest("combineDBSNPDuplicateSites:", spec, true);
    }

    @Test
    public void combineLeavesUnfilteredRecordsUnfiltered() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants --no_cmdline_in_header -o %s "
                        + " -R " + b37KGReference
                        + " -V " + privateTestDir + "combineVariantsLeavesRecordsUnfiltered.vcf",
                1,
                Arrays.asList("0f221847e76521250de1abcba535e49c"));
        cvExecuteTest("combineLeavesUnfilteredRecordsUnfiltered: ", spec, false);
    }

    @Test
    public void combiningGVCFsFails() {
        try {
            WalkerTestSpec spec = new WalkerTestSpec(
                    "-T CombineVariants --no_cmdline_in_header -o %s "
                            + " -R " + b37KGReference
                            + " -V " + privateTestDir + "gvcfExample1.vcf",
                    1,
                    Arrays.asList("FAILFAILFAILFAILFAILFAILFAILFAIL"));
            executeTest("combiningGVCFsFails", spec);
        } catch (Exception e) { } // do nothing
    }

    @Test
    public void combineSymbolicVariants() {
        // Just checking that this does not fail, hence no output files and MD5
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants --no_cmdline_in_header -o %s "
                        + " -R " + hg19ReferenceWithChrPrefixInChromosomeNames
                        + " -V " + privateTestDir + "WES-chr1.DEL.vcf"
                        + " -V " + privateTestDir + "WGS-chr1.DEL.vcf"
                        + " -genotypeMergeOptions UNIQUIFY",
                0,
                Arrays.asList(""));
        executeTest("combineSymbolicVariants: ", spec);
    }

    @Test
    public void combineSpanningDels() {
        // Just checking that this does not fail, hence no output files and MD5
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineVariants --no_cmdline_in_header -o %s "
                        + " -R " + b37KGReference
                        + " -V " + privateTestDir + "test.spanningdel.combine.1.vcf "
                        + " -V " + privateTestDir + "test.spanningdel.combine.2.vcf "
                        + " -genotypeMergeOptions UNIQUIFY",
                0,
                Arrays.asList(""));
        executeTest("combineSpanningDels: ", spec);
    }
}
