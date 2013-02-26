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

import java.util.Arrays;
import java.util.Random;

public class BiasedDownsamplingIntegrationTest extends WalkerTest {

    private final static String baseCommand1 = "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommand2 = "-T UnifiedGenotyper -R " + hg19Reference + " --no_cmdline_in_header -glm BOTH -L 20:1,000,000-5,000,000";
    private final static String baseCommand3 = "-T UnifiedGenotyper -R " + hg19Reference + " --no_cmdline_in_header -glm BOTH -L 20:4,000,000-5,000,000";
    private final String ArtificalBAMLocation = privateTestDir + "ArtificallyContaminatedBams/";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing UnifiedGenotyper contamination down-sampling
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testContaminationDownsamplingFlat() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand1 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -contamination 0.20", 1,
                Arrays.asList("1f9071466fc40f4c6a0f58ac8e9135fb"));
        executeTest("test contamination_percentage_to_filter 0.20", spec);
    }

    @Test
    public void testContaminationDownsamplingFlatAndPerSample() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand1 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 --contamination_fraction_per_sample_file " + ArtificalBAMLocation + "NA12878.NA19240.contam.txt --contamination_fraction_to_filter 0.10", 1,
                Arrays.asList("53395814dd6990448a01a294ccd69bd2"));
        executeTest("test contamination_percentage_to_filter per-sample and .20 overall", spec);
    }

    @Test
    public void testContaminationDownsamplingPerSampleOnly() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand1 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -contaminationFile " + ArtificalBAMLocation + "NA19240.contam.txt", 1,
                Arrays.asList("4af83a883ecc03a23b0aa6dd4b8f1ceb"));
        executeTest("test contamination_percentage_to_filter per-sample", spec);
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // testing UnifiedGenotyper contamination down-sampling on BAMs with artificially created contaminated.
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    private void testDefaultContamination() {
        final String bam1 = "NA11918.with.1.NA12842.reduced.bam";
        final String bam2 = "NA12842.with.1.NA11918.reduced.bam";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand2 + " -I " + ArtificalBAMLocation + bam1 + " -I " + ArtificalBAMLocation + bam2 + " -o %s ", 1,
                Arrays.asList("e2e5a8dd313f8d7e382e7d49dfac59a2"));
        executeTest("test contamination on Artificial Contamination (flat) on " + bam1 + " and " + bam2 + " with default downsampling.", spec);
    }

    private void testFlatContamination(final String bam1, final String bam2, final Double downsampling, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand2 + " -I " + ArtificalBAMLocation + bam1 + " -I " + ArtificalBAMLocation + bam2 + " -o %s -contamination " + downsampling.toString(), 1,
                Arrays.asList(md5));
        executeTest("test contamination on Artificial Contamination (flat) on " + bam1 + " and " + bam2 + " downsampling " + downsampling.toString(), spec);
    }

    @Test
    public void testFlatContaminationCase1() {
        testFlatContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.05, "e2e5a8dd313f8d7e382e7d49dfac59a2");
    }

    @Test
    public void testFlatContaminationCase2() {
        testFlatContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.1, "549737002f98775fea8f46e7ea174dde");
    }

    @Test
    public void testFlatContaminationCase3() {
        testFlatContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.2, "529d82c2a33fcc303a5dc55de2d56979");
    }

    @Test
    public void testFlatContaminationCase4() {
        testFlatContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.2.NA11918.reduced.bam", 0.1, "b5689972fbb7d230a372ee5f0da1c6d7");
    }

    @Test
    public void testFlatContaminationCase5() {
        testFlatContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.2.NA11918.reduced.bam", 0.2, "9dceee2e921b53fbc1ce137a7e0b7b74");
    }

    @Test
    public void testFlatContaminationCase6() {
        testFlatContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.2.NA11918.reduced.bam", 0.3, "d6a74061033503af80dcaea065bfa075");
    }

    @Test
    public void testFlatContaminationCase7() {
        testFlatContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.1, "7d1b5efab58a1b8f9d99fcf5af82f15a");
    }

    @Test
    public void testFlatContaminationCase8() {
        testFlatContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.2, "a7f8d5c79626aff59d7f426f79d8816e");
    }

    @Test
    public void testFlatContaminationCase9() {
        testFlatContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.3, "fcf482398b7c908e3e2d1e4d5da6377b");
    }

    private void testPerSampleContamination(String bam1, String bam2, String persampleFile, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand2 + " -I " + ArtificalBAMLocation + bam1 + " -I " + ArtificalBAMLocation + bam2 + " -o %s  -contaminationFile " + persampleFile, 1,
                Arrays.asList(md5));
        executeTest("test contamination on Artificial Contamination (per-sample) on " + bam1 + " and " + bam2 + " with " + persampleFile, spec);
    }

    @Test
    public void testPerSampleContaminationCase1() {
        testPerSampleContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.1.txt", "e00278527a294833259e9e411728e395");
    }

    @Test
    public void testPerSampleContaminationCase2() {
        testPerSampleContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.2.txt", "a443e793f0b0e2ffce1b751634d706e2");
    }

    @Test
    public void testPerSampleContaminationCase3() {
        testPerSampleContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.3.txt", "e11d83a7815ce757afbcf7689568cb25");
    }

    @Test
    public void testPerSampleContaminationCase4() {
        testPerSampleContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.4.txt", "615042eeeffe042bd1c86279d34f80b6");
    }

    @Test
    public void testPerSampleContaminationCase5() {
        testPerSampleContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.1.txt", "9bc99fc79ca34744bf26cb19ee4ef44d");
    }

    @Test
    public void testPerSampleContaminationCase6() {
        testPerSampleContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.2.txt", "143626fe5fce765d6c997a64f058a813");
    }

    @Test
    public void testPerSampleContaminationCase7() {
        testPerSampleContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.3.txt", "f2593674cef894eda4e0be9cf3158f57");
    }

    @Test
    public void testPerSampleContaminationCase8() {
        testPerSampleContamination("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.4.txt", "fb7ce0740767ae3896b3e552026da1e4");
    }


    private void testPerSampleEqualsFlat(final String bam1, final String bam2, final String persampleFile, final Double downsampling, final String md5) {
        final String command =  baseCommand3 + " -I " + ArtificalBAMLocation + bam1 + " -I " + ArtificalBAMLocation + bam2 + " -o %s  ";

        WalkerTestSpec spec = new WalkerTestSpec( command +" -contaminationFile " + persampleFile, 1, Arrays.asList(md5));
        final Random rnd = GenomeAnalysisEngine.getRandomGenerator();

        rnd.setSeed(123451); // so that the two test cases have a hope of giving the same result
        executeTest("test contamination on Artificial Contamination, with per-sample file on " + bam1 + " and " + bam2 + " with " + persampleFile, spec);

        spec = new WalkerTestSpec(command + "-contamination " + downsampling.toString(), 1, Arrays.asList(md5));

        rnd.setSeed(123451); // so that the two test cases have a hope of giving the same result
        executeTest("test contamination on Artificial Contamination, with flat contamination on " + bam1 + " and " + bam2 + " with " + downsampling.toString(), spec);

    }

    // verify that inputing a file with an effectively flat contamination level is equivalent to handing in a flat contamination level

    @Test
    public void testPerSampleEqualsFlatContaminationCase1() {
        testPerSampleEqualsFlat("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.6.txt", 0.0, "");
    }

    @Test
    public void testPerSampleEqualsFlatContaminationCase2() {
        testPerSampleEqualsFlat("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.7.txt", 0.15, "");
    }

    @Test
    public void testPerSampleEqualsFlatContaminationCase3() {
        testPerSampleEqualsFlat("NA11918.with.2.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", ArtificalBAMLocation + "contamination.case.8.txt", 0.3, "");
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // testing HaplotypeCaller Contamination Removal
    //
    // --------------------------------------------------------------------------------------------------------------


    @Test
    public void testHCContaminationDownsamplingFlat() {
        final String baseCommand = "-T HaplotypeCaller -R " + b36KGReference + " --no_cmdline_in_header --dbsnp " + b36dbSNP129;
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -contamination 0.20", 1,
                Arrays.asList("c3a253467ead7b1cfe9fd9dd310828b1"));
        executeTest("HC calling with contamination_percentage_to_filter 0.20", spec);
    }

    //  HaplotypeCaller can only (currently) use flat contamination reduction, not per-sample. Until that is implemented, this test
    @Test
    public void testHCCannotProcessPerSampleContamination() {
        final String baseCommand = "-T HaplotypeCaller -R " + hg19Reference + " --no_cmdline_in_header  -L 20:3,000,000-5,000,000";
        final String bam1 = "NA11918.with.1.NA12842.reduced.bam";
        final String perSampleFile = ArtificalBAMLocation + "contamination.case.1.txt";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand + " -I " + ArtificalBAMLocation + bam1 + " -o %s  -contaminationFile " + perSampleFile, 1,
                UserException.class);
        executeTest("HC should fail on per-Sample contamination removal.", spec);
    }


    private void testHCFlatContamination(final String bam1, final String bam2, final Double downsampling, final String md5) {
        final String baseCommand = "-T HaplotypeCaller -R " + hg19Reference + " --no_cmdline_in_header -L 20:3,000,000-5,000,000";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseCommand + " -I " + ArtificalBAMLocation + bam1 + " -I " + ArtificalBAMLocation + bam2 + " -o %s -contamination " + downsampling.toString(), 1,
                Arrays.asList(md5));
        executeTest("HC test contamination on Artificial Contamination (flat) on " + bam1 + " and " + bam2 + " downsampling " + downsampling.toString(), spec);
    }

    @Test
    public void testHCFlatContaminationCase1() {
        testHCFlatContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.05, "c3e695381d8627e3922d8c642b66c3ce");
    }

    @Test
    public void testHCFlatContaminationCase2() {
        testHCFlatContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.1, "002d2b45336d88d7c04e19f9f26e29d9");
    }

    @Test
    public void testHCFlatContaminationCase3() {
        testHCFlatContamination("NA11918.with.1.NA12842.reduced.bam", "NA12842.with.1.NA11918.reduced.bam", 0.2, "1809a33ac112d1a3bd7a071c566794dd");
    }

}
