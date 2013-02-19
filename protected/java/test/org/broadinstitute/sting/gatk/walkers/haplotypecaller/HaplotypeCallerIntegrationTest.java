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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

public class HaplotypeCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String NA12878_BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String NA12878_CHR20_BAM = validationDataLocation + "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam";
    final String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final String NA12878_RECALIBRATED_BAM = privateTestDir + "NA12878.100kb.BQSRv2.example.bam";
    final String INTERVALS_FILE = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals";

    private void HCTest(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s -L %s", REF, bam, INTERVALS_FILE) + " --no_cmdline_in_header -o %s -minPruning 3";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCaller: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSample() {
        HCTest(CEUTRIO_BAM, "", "1e49fd927d79594a993ea6c4a1d10004");
    }

    @Test
    public void testHaplotypeCallerSingleSample() {
        HCTest(NA12878_BAM, "", "1b39ac32c9cbba26ed60c6b06be81359");
    }

    @Test(enabled = false)
    public void testHaplotypeCallerSingleSampleWithDbsnp() {
        HCTest(NA12878_BAM, "-D " + b37dbSNP132, "");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGA() {
        HCTest(CEUTRIO_BAM, "--max_alternate_alleles 3 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf",
                "f751363288740c6fd9179a487be61fb4");
    }

    private void HCTestComplexGGA(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " --no_cmdline_in_header -o %s -minPruning 3 -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerComplexGGA: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGAComplex() {
        HCTestComplexGGA(NA12878_CHR20_BAM, "-L 20:119673-119823 -L 20:121408-121538",
                "719402122fe92cfe7a3fa6b7cdb66f26");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGAMultiAllelic() {
        HCTestComplexGGA(NA12878_CHR20_BAM, "-L 20:133041-133161 -L 20:300207-300337",
                "71ef8d0217c1a73dd360413dccd05f4d");
    }

    private void HCTestComplexVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10028767-10028967 -L 20:10431524-10431924 -L 20:10723661-10724061 -L 20:10903555-10903955 --no_cmdline_in_header -o %s -minPruning 4";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerComplexVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerMultiSampleComplex() {
        HCTestComplexVariants(privateTestDir + "AFR.complex.variants.bam", "", "88ae6e7b34514043bfc78b1ecf29a341");
    }

    private void HCTestSymbolicVariants(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:5947969-5948369 -L 20:61091236-61091636 --no_cmdline_in_header -o %s -minPruning 1";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerSymbolicVariants: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleSymbolic() {
        HCTestSymbolicVariants(NA12878_CHR20_BAM, "", "855827f901b63b41dcd37dd49dd3a1ac");
    }

    private void HCTestIndelQualityScores(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, bam) + " -L 20:10,005,000-10,025,000 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerIndelQualityScores: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleIndelQualityScores() {
        HCTestIndelQualityScores(NA12878_RECALIBRATED_BAM, "", "c0ac5a1f75c66052b19684eb37c088cb");
    }

    // That problem bam came from a user on the forum and it spotted a problem where the ReadClipper
    // was modifying the GATKSamRecord and that was screwing up the traversal engine from map call to
    // map call. So the test is there for consistency but not for correctness. I'm not sure we can trust
    // any of the calls in that region because it is so messy. The only thing I would maybe be worried about is
    // that the three calls that are missing happen to all be the left most calls in the region
    @Test
    public void HCTestProblematicReadsModifiedInActiveRegions() {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "haplotype-problem-4.bam") + " --no_cmdline_in_header -o %s -minPruning 3 -L 4:49139026-49139965";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("598d245498c0d0b55e263f0a061a77e3"));
        executeTest("HCTestProblematicReadsModifiedInActiveRegions: ", spec);
    }

    @Test
    public void HCTestStructuralIndels() {
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "AFR.structural.indels.bam") + " --no_cmdline_in_header -o %s -minPruning 6 -L 20:8187565-8187800 -L 20:18670537-18670730";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("a4e74226b16a7d8c5999620c2f6be1ba"));
        executeTest("HCTestStructuralIndels: ", spec);
    }

    @Test
    public void HCTestDoesNotFailOnBadRefBase() {
        // don't care about the output - just want to make sure it doesn't fail
        final String base = String.format("-T HaplotypeCaller -R %s -I %s", REF, privateTestDir + "NA12878.readsOverBadBase.chr3.bam") + " --no_cmdline_in_header -o /dev/null -L 3:60830000-60840000 --minPruning 3 -stand_call_conf 2 -stand_emit_conf 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Collections.<String>emptyList());
        executeTest("HCTestDoesNotFailOnBadRefBase: ", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing reduced reads
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void HCTestReducedBam() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam -o %s -L 1:67,225,396-67,288,518", 1,
                Arrays.asList("9f0bb0b97857c66937de39670e195d00"));
        executeTest("HC calling on a ReducedRead BAM", spec);
    }

    @Test
    public void testReducedBamWithReadsNotFullySpanningDeletion() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "reduced.readNotFullySpanningDeletion.bam -o %s -L 1:167871297", 1,
                Arrays.asList("87bd7ac2f7d65580838c7c956ccf52b7"));
        executeTest("test calling on a ReducedRead BAM where the reads do not fully span a deletion", spec);
    }
}
