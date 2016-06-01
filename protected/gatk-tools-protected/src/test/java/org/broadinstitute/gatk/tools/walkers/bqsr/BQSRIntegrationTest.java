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

package org.broadinstitute.gatk.tools.walkers.bqsr;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.commandline.ArgumentException;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
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
                    " -knownSites " + (reference.equals(b36KGReference) ? b36dbSNP129 : (reference.equals(b37KGReference) ? b37dbSNP129 : hg18dbSNP132)) +
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
    private static final String SimpleCigarMatchMismatchBam = privateTestDir + "simpleCigarMatchMismatch.bam";
    private static final String SimpleCigarMatchMismatchInterval = "1:1-60";

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, "", "fc9df1faf67bab70d32f89bcf4fa39db")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov ContextCovariate", "73ec38eb23b1739ecef8194cbb1132a3")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov CycleCovariate", "2d5721193ed4410d1a7d8db467a1fa05")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --indels_context_size 4", "16df7f1745f17f190c9fc33c475b91d8")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --low_quality_tail 5", "01811003ae811ee74c4b8d3eb5e992fe")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --quantizing_levels 6", "5e0eea6b0b300fbd2edabc3506ad3a60")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --mismatches_context_size 4", "8500b9747c16cb8eb17082163bdb8069")},
                {new BQSRTest(b36KGReference, SimpleCigarMatchMismatchBam, SimpleCigarMatchMismatchInterval, "", "56dfb2918a4cdae3ef9d705a43e85194")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", "", "0b5a8e259e997e4c7b5836d4c28e6f4d")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-10,200,000", "", "281682124584ab384f23359934df0c3b")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.1RG.bam", "1:10,000,000-10,200,000", "", "0a92fdff5fd26227c29d34eda5a32f49")},
                {new BQSRTest(hg18Reference, privateTestDir + "originalQuals.1kg.chr1.1-1K.1RG.bam", "chr1:1-1,000", " -OQ", "90d8c24077e8ae9a0037a9aad5f09e31")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-20,000,000", " --solid_recal_mode REMOVE_REF_BIAS", "c41ef02c640ef1fed4bfc03b9b33b616")},
                {new BQSRTest(b36KGReference, privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam", "1:50,000-80,000", " --solid_nocall_strategy LEAVE_READ_UNRECALIBRATED", "b577cd1d529425f66db49620db09fdca")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:anyNameABCD,VCF " + privateTestDir + "vcfexample3.vcf", "0b5a8e259e997e4c7b5836d4c28e6f4d")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:bed " + validationDataLocation + "bqsrKnownTest.bed", "9ad49269c0156f8ab1173261bf23e600")},
                // make sure we work with ION torrent bam
                {new BQSRTest(b37KGReference, privateTestDir + "iontorrent.bam", "20:10,000,000-10,200,000", "", "04bfa4760767022e7f5252e6e4432cc1")},
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

        tests.add(new Object[]{1, new PRTest(" -qq -1", "8a38828e3b14ce067614d4248e3ea95a")});
        tests.add(new Object[]{1, new PRTest(" -qq 6", "e4f23250b2c87f0d68d042cc3d2ec1d3")});
        tests.add(new Object[]{1, new PRTest(" -DIQ", "2dfa45f004d3a371fd290ed67fbdf573")});
        tests.add(new Object[]{1, new PRTest(" --useOriginalQualities -SQQ 10 -SQQ 20 -SQQ 30", "4882354d9e603f9bbe7c9591bba0a573")});
        tests.add(new Object[]{1, new PRTest(" --useOriginalQualities -SQQ 10 -SQQ 20 -SQQ 30 -RDQ", "6ffdfc4593e83f7c234b6249412433af")});
        for ( final int nct : Arrays.asList(1, 2, 4) ) {
            tests.add(new Object[]{nct, new PRTest("", "6451093cadfc14d7359617b2a7ea6db8")});
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

    @Test
    public void testPRFailWithBadPL() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b37KGReference +
                        " -I " + privateTestDir + "badPLForBQSR.bam" +
                        " -L 1:10,000,000-10,200,000" +
                        " -o %s",
                1,
                UserException.class);
        executeTest("testPRFailWithBadPL", spec);
    }

    @Test
    public void testPRWithConflictingArguments_qqAndSQQ() {
        // -qq and -SQQ shouldn't be able to be run in the same command
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + HiSeqBam +
                        " -L " + HiSeqInterval +
                        " -qq 4 -SQQ 9" +
                        " -BQSR " + privateTestDir + "HiSeq.1mb.1RG.highMaxCycle.table" +
                        " -o /dev/null",
                0,
                ArgumentException.class);
        executeTest("testPRWithConflictingArguments_qqAndSQQ", spec);
    }

    @Test
    public void testPRWithConflictingArguments_qlAndSQQ() {
        // Arguments -SQQ and -ql conflict and should throw an exception
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + HiSeqBam +
                        " -L " + HiSeqInterval +
                        " -SQQ 4 -ql 4 " +
                        " -BQSR " + privateTestDir + "HiSeq.1mb.1RG.lowMaxCycle.table" +
                        " -o /dev/null",
                0,
                ArgumentException.class);
        executeTest("testPRWithConflictingArguments_qlAndSQQ", spec);
    }
}
