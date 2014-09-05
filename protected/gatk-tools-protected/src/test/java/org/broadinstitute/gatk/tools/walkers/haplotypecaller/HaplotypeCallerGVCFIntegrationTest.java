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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HaplotypeCallerGVCFIntegrationTest extends WalkerTest {

    @DataProvider(name = "MyDataProviderHaploid")
    public Object[][] makeMyDataProviderHaploid() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "646ec07bd026da1e72b5e789f5aa3a3d"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "9acbee336e91cbfc1abeebd41bbcc9dd"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "e696ffc927af7f7a36dc7d49dad2c4f8"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "c223d6fe112d2bb698811600c3b7f6af"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "dacb94af2632e4dc4a1948306dd1661c"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, "ef9093880efedac09b78c8fb26420e84"});

        return tests.toArray(new Object[][]{});
    }


    @DataProvider(name = "MyDataProvider")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "f8bc02b296d89369fe7ddf3923743c26"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "ad0b6274991d451a6491f3118f0e7b2b"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "956938fd1b5c7bf1b717326a7dd46d91"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "1c3570461e96ad6d66c6abb0fd6ee865"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "66019a0914f905522da6bd3b557a57d1"});

        final String NA12878bandedResolutionMD5 = "de763f6c9a5aec4586a2671941e4c96d";
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, NA12878bandedResolutionMD5});
        tests.add(new Object[]{NA12878_WEx + " -I " + privateTestDir + "NA20313.highCoverageRegion.bam -sn NA12878",
                ReferenceConfidenceMode.GVCF, WExIntervals, NA12878bandedResolutionMD5});

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "MyDataProviderTetraploid")
    public Object[][] makeMyDataProviderTetraploid() {
        List<Object[]> tests = new ArrayList<>();

        final String PCRFreeIntervals = "-L 20:10,000,000-10,010,000";
        final String WExIntervals = "-L 20:10,000,000-10,100,000 -isr INTERSECTION -L " + hg19Chr20Intervals;

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.NONE, PCRFreeIntervals, "df5b803750e604bedd6d713c59f46c33"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.BP_RESOLUTION, PCRFreeIntervals, "8e06e6502d958cab79b9275460c06389"});
        tests.add(new Object[]{NA12878_PCRFREE, ReferenceConfidenceMode.GVCF, PCRFreeIntervals, "f32fb9a2485c583bc22ff4f09217318c"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.NONE, WExIntervals, "eac793500fbc7de46259000dbbcdd27d"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.BP_RESOLUTION, WExIntervals, "4a431e0e387de3f791318f67d8855b0b"});
        tests.add(new Object[]{NA12878_WEx, ReferenceConfidenceMode.GVCF, WExIntervals, "1058d3fe6553e07f002f994759c9647d"});

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MyDataProvider")
    public void testHCWithGVCF(String bam, ReferenceConfidenceMode mode, String intervals, String md5) {
        final String commandLine = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCF bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MyDataProviderHaploid", enabled=true)
    public void testHCWithGVCFHaploid(final String bam, final ReferenceConfidenceMode mode, final String intervals, final String md5) {
        final String commandLine = String.format("-T HaplotypeCaller -ploidy 1 --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCFHaploid bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MyDataProviderTetraploid", enabled=true)
    public void testHCWithGVCFTetraploid(final String bam, final ReferenceConfidenceMode mode, final String intervals, final String md5) {
        final String commandLine = String.format("-T HaplotypeCaller -ploidy 4 --disableDithering --pcr_indel_model NONE -R %s -I %s %s -ERC %s --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, bam, intervals, mode, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final String name = "testHCWithGVCFTetraploid bam=" + bam + " intervals= " + intervals + " gvcf= " + mode;
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(md5));
        executeTest(name, spec);
    }

    @Test
    public void testERCRegionWithNoCalledHaplotypes() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001", HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest("testERCRegionWithNoCalledHaplotypes", spec);
    }

    @Test()
    public void testMissingGVCFIndexException() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001");
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.GVCFIndexException.class);
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    @Test()
    public void testWrongParameterGVCFIndexException() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001", HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER + 1);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.GVCFIndexException.class);
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    @Test()
    public void testWrongTypeGVCFIndexException() {
        // ensure non-optimal, if optimal changes
        GATKVCFIndexType type = GATKVCFIndexType.DYNAMIC_SEEK;
        if (HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE == GATKVCFIndexType.DYNAMIC_SEEK)
            type = GATKVCFIndexType.DYNAMIC_SIZE;

        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, privateTestDir + "noCallRefModel.bam", "20:17000001-18000001", type, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", 1, UserException.GVCFIndexException.class);
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    private final static String WRONG_GVCF_RECORD_ORDER_BUGFIX_INTERVALS = privateTestDir + "gvcf_unsorted_records_bug.interval_list";
    private final static String WRONG_GVCF_RECORD_ORDER_BUGFIX_BAM = privateTestDir + "gvcf_unsorted_records_bug.bam";

    @Test()
    public void testWrongGVCFNonVariantRecordOrderBugFix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, WRONG_GVCF_RECORD_ORDER_BUGFIX_BAM, WRONG_GVCF_RECORD_ORDER_BUGFIX_INTERVALS, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("e6a4e571abb59b925d59a38d244f0abe"));
        spec.disableShadowBCF();
        executeTest("testMissingGVCFIndexingStrategyException", spec);
    }

    private static final String NOCALL_GVCF_BUGFIX_INTERVALS = privateTestDir + "gvcf_nocall_bug.interval_list";
    private static final String NOCALL_GVCF_BUGFIX_BAM = privateTestDir + "gvcf_nocall_bug.bam";
    private static final String GENERAL_PLOIDY_BUGFIX1_BAM = privateTestDir + "general-ploidy-arrayindex-bug-1.bam";
    private static final String GENERAL_PLOIDY_BUGFIX1_INTERVALS = privateTestDir + "general-ploidy-arrayindex-bug-1.intervals";
    private static final String GENERAL_PLOIDY_BUGFIX2_BAM = privateTestDir + "general-ploidy-arrayindex-bug-2.bam";
    private static final String GENERAL_PLOIDY_BUGFIX2_INTERVALS = privateTestDir + "general-ploidy-arrayindex-bug-2.intervals";


    @Test
    public void testNoCallGVCFMissingPLsBugFix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d",
                b37KGReference, NOCALL_GVCF_BUGFIX_BAM, NOCALL_GVCF_BUGFIX_INTERVALS, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("187e43b1034478edced1eb5c664bac34"));
        spec.disableShadowBCF();
        executeTest("testNoCallGVCFMissingPLsBugFix", spec);
    }

    /**
     * Checking that the bug's Exception is no longer been thrown (thus no md5).
     */
    @Test(enabled=true)
    public void testGeneralPloidyArrayIndexBug1Fix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 1 -maxAltAlleles 2 -isr INTERSECTION -L 1:23696115-23696189",
                b37KGReference, GENERAL_PLOIDY_BUGFIX1_BAM, GENERAL_PLOIDY_BUGFIX1_INTERVALS, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest(" testGeneralPloidyArrayIndexBug1Fix", spec);
    }

    /**
     * Checking that the bug's Exception is no longer been thrown (thus no md5).
     */
    @Test(enabled=true)
    public void testGeneralPloidyArrayIndexBug2Fix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -R %s -I %s -L %s -ERC GVCF --no_cmdline_in_header -variant_index_type %s -variant_index_parameter %d -ploidy 2 -maxAltAlleles 2 -A DepthPerSampleHC -A StrandBiasBySample -L 1:38052860-38052937",
                b37KGReference, GENERAL_PLOIDY_BUGFIX2_BAM, GENERAL_PLOIDY_BUGFIX2_INTERVALS, HaplotypeCaller.OPTIMAL_GVCF_INDEX_TYPE, HaplotypeCaller.OPTIMAL_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest(" testGeneralPloidyArrayIndexBug2Fix", spec);
    }

}
