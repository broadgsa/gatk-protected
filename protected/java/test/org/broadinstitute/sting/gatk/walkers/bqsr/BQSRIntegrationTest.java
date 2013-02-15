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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author ebanks
 * @since 7/16/12
 */
public class BQSRIntegrationTest extends WalkerTest {

    private static class BQSRTest {
        final String reference;
        final String interval;
        final String bam;
        final String args;
        final String md5;

        private BQSRTest(String reference, String bam, String interval, String args, String md5) {
            this.reference = reference;
            this.bam = bam;
            this.interval = interval;
            this.args = args;
            this.md5 = md5;
        }

        public String getCommandLine() {
            return " -T BaseRecalibrator" +
                    " -R " + reference +
                    " -I " + bam +
                    " -L " + interval +
                    args +
                    " -knownSites " + (reference.equals(b36KGReference) ? b36dbSNP129 : hg18dbSNP132) +
                    " --allow_potentially_misencoded_quality_scores" +  // TODO -- remove me when we get new SOLiD bams
                    " -o %s" +
                    " -sortAllCols";
        }

        @Override
        public String toString() {
            return String.format("BQSR(bam='%s', args='%s')", bam, args);
        }
    }

    private static final String HiSeqBam = privateTestDir + "HiSeq.1mb.1RG.bam";
    private static final String HiSeqInterval = "chr1:10,000,000-10,100,000";

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, "", "61fd466b5e94d2d67e116f6f67c9f939")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov ContextCovariate", "e08b5bcdb64f4beea03730e5631a14ca")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov CycleCovariate", "448a45dc154c95d1387cb5cdddb67071")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --indels_context_size 4", "c1e7999e445d51bbe2e775dac5325643")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --low_quality_tail 5", "a57c16918cdfe12d55a89c21bf195279")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --quantizing_levels 6", "836dccacf48ccda6b2843d07e8f1ef4d")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --mismatches_context_size 4", "0fb2aedc2f8d66b5821cb570f15a8c4d")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", "", "c9953f020a65c1603a6d71aeeb1b95f3")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-10,200,000", "", "85a120b7d86b61597b86b9e93decbdfc")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.1RG.bam", "1:10,000,000-10,200,000", "", "5248dc49aec0323c74b496bb4928c73c")},
                {new BQSRTest(b36KGReference, validationDataLocation + "originalQuals.1kg.chr1.1-1K.1RG.bam", "1:1-1,000", " -OQ", "cb52f267e0010f849f50b0bf1de474a1")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-20,000,000", " --solid_recal_mode REMOVE_REF_BIAS", "1425a5063ee757dbfc013df24e65a67a")},
                {new BQSRTest(b36KGReference, privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam", "1:50,000-80,000", " --solid_nocall_strategy LEAVE_READ_UNRECALIBRATED", "c1c3cda8caceed619d3d439c3990cd26")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:anyNameABCD,VCF " + privateTestDir + "vcfexample3.vcf", "c9953f020a65c1603a6d71aeeb1b95f3")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:bed " + validationDataLocation + "bqsrKnownTest.bed", "5bfff0c699345cca12a9b33acf95588f")},
        };
    }

    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.md5));
        executeTest("testBQSR-"+params.args, spec).getFirst();
    }

    @Test
    public void testBQSRFailWithoutDBSNP() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam" +
                        " -L 1:10,000,000-10,200,000" +
                        " -o %s",
                1, // just one output file
                UserException.CommandLineException.class);
        executeTest("testBQSRFailWithoutDBSNP", spec);
    }

    @Test
    public void testBQSRCSV() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam" +
                        " -knownSites " + b36dbSNP129 +
                        " -L 1:10,000,000-10,200,000" +
                        " -o /dev/null" +
                        " -sortAllCols" +
                        " --plot_pdf_file /dev/null" +
                        " --intermediate_csv_file %s",
                Arrays.asList("dd6e0e1e3f53f8ae0c8f5de21ded6ee9"));
        executeTest("testBQSR-CSVfile", spec);
    }

    @Test
    public void testBQSRFailWithSolidNoCall() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam" +
                        " -L 1:50,000-80,000" +
                        " --allow_potentially_misencoded_quality_scores" +  // TODO -- remove me when we get new SOLiD bams
                        " -o %s" +
                        " -sortAllCols",
                1, // just one output file
                UserException.class);
        executeTest("testBQSRFailWithSolidNoCall", spec);
    }

    @Test
    public void testBQSRFailWithReducedBam() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b37KGReference +
                        " -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam" +
                        " -L 1:67,225,396-67,288,518" +
                        " -o /dev/null",
                0,
                UserException.class);
        executeTest("testBQSRFailWithReducedBam", spec);
    }

    private static class PRTest {
        final String args;
        final String md5;

        private PRTest(String args, String md5) {
            this.args = args;
            this.md5 = md5;
        }

        @Override
        public String toString() {
            return String.format("PrintReads(args='%s')", args);
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createPRTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{1, new PRTest(" -qq -1", "fcc136b877fbde38791533b0f1ae39e4")});
        tests.add(new Object[]{1, new PRTest(" -qq 6", "f21b537c1689b8051b878ea5cc9b61a0")});
        tests.add(new Object[]{1, new PRTest(" -DIQ", "1d04a242bf825177d6a45eff9fbed647")});

        for ( final int nct : Arrays.asList(1, 2, 4) ) {
            tests.add(new Object[]{nct, new PRTest("", "b6f343ac69c63cdb49205c13e67297fc")});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PRTest")
    public void testPR(final int nct, PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + privateTestDir + "HiSeq.1mb.1RG.bam" +
                        " -nct " + nct +
                        " --no_pg_tag" +
                        " -BQSR " + privateTestDir + "HiSeq.20mb.1RG.table" +
                        params.args +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }

    @Test
    public void testPRNoFailWithHighMaxCycle() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + HiSeqBam +
                        " -L " + HiSeqInterval +
                        " --no_pg_tag" +
                        " -BQSR " + privateTestDir + "HiSeq.1mb.1RG.highMaxCycle.table" +
                        " -o /dev/null",
                0,
                Arrays.<String>asList());
        executeTest("testPRNoFailWithHighMaxCycle", spec);
    }

    @Test
    public void testPRFailWithLowMaxCycle() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + HiSeqBam +
                        " -L " + HiSeqInterval +
                        " --no_pg_tag" +
                        " -BQSR " + privateTestDir + "HiSeq.1mb.1RG.lowMaxCycle.table" +
                        " -o /dev/null",
                0,
                UserException.class);
        executeTest("testPRFailWithLowMaxCycle", spec);
    }
}
