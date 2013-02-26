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

package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String DELETION_BAM = validationDataLocation + "filtered_deletion_for_reduce_reads.bam";
    final String STASH_BAM = validationDataLocation + "ReduceReadsStashBug.bam";
    final String STASH_L = " -L 14:73718184-73718284 -L 14:73718294-73718330 -L 14:73718360-73718556";
    final String DIVIDEBYZERO_BAM = validationDataLocation + "ReduceReadsDivideByZeroBug.bam";
    final String DIVIDEBYZERO_L = " -L " + validationDataLocation + "ReduceReadsDivideByZeroBug.intervals";
    final String L = " -L 20:10,100,000-10,120,000 ";
    final String COREDUCTION_BAM_A = validationDataLocation + "coreduction.test.A.bam";
    final String COREDUCTION_BAM_B = validationDataLocation + "coreduction.test.B.bam";
    final String COREDUCTION_L = " -L 1:1,853,860-1,854,354 -L 1:1,884,131-1,892,057";
    final String OFFCONTIG_BAM = privateTestDir + "readOffb37contigMT.bam";
    final String INSERTIONS_AT_EDGE_OF_CONSENSUS_BAM = privateTestDir + "rr-too-many-insertions.bam";

    private void RRTest(String testName, String args, String md5) {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, BAM) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("testDefaultCompression ", L, "17908e8515217c4693d303ed68108ccc");
    }

    @Test(enabled = true)
    public void testInsertionsAtEdgeOfConsensus() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, INSERTIONS_AT_EDGE_OF_CONSENSUS_BAM) + " -o %s ";
        executeTest("testInsertionsAtEdgeOfConsensus", new WalkerTestSpec(base, Arrays.asList("3103667fc68c3136a8cfa8e22429f94e")));
    }

    @Test(enabled = true)
    public void testMultipleIntervals() {
        String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";
        RRTest("testMultipleIntervals ", intervals, "497c5e36c2beaad2fcdbd02a0b9c121b");
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -cs 10 -minvar 0.3 -mindel 0.3 " + L, "0ff4142e4d7b6a9a9c76012246ad9e2d");
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -cs 30 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "7890a37444a0e05b902f63a83238ce37");
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        RRTest("testIndelCompression ", " -cs 50 -L 20:10,100,500-10,100,600 ", "f58ae2154e0e5716be0e850b7605856e");
    }

    @Test(enabled = true)
    public void testFilteredDeletionCompression() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, DELETION_BAM) + " -o %s ";
        executeTest("testFilteredDeletionCompression", new WalkerTestSpec(base, Arrays.asList("bfe0693aea74634f1035a9bd11302517")));
    }

    /**
     * Bug reported by Adam where a read that got clipped before actually belongs 2 intervals ahead
     * and a subsequent tail leaves only this read in the stash. The next read to come in is in fact 
     * before (alignment start) than this read, so the TreeSet breaks with a Key out of Range error
     * that was freaking hard to catch. 
     * 
     * This bam is simplified to replicate the exact bug with the three provided intervals.
     */
    @Test(enabled = true)
    public void testAddingReadAfterTailingTheStash() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s", STASH_L, REF, STASH_BAM) + " -o %s ";
        executeTest("testAddingReadAfterTailingTheStash", new WalkerTestSpec(base, Arrays.asList("f118e83c394d21d901a24230379864fc")));
    }

    /**
     * Divide by zero bug reported by GdA and users in the forum. Happens when the downsampler goes over a region where all reads get
     * filtered out.
     */
    @Test(enabled = true)
    public void testDivideByZero() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s", DIVIDEBYZERO_L, REF, DIVIDEBYZERO_BAM) + " -o %s ";
        executeTest("testDivideByZero", new WalkerTestSpec(base, Arrays.asList("bd5198a3e21034887b741faaaa3964bf")));
    }

    @Test(enabled = true)
    public void testCoReduction() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s -I %s", COREDUCTION_L, REF, COREDUCTION_BAM_A, COREDUCTION_BAM_B) + " -o %s ";
        executeTest("testCoReduction", new WalkerTestSpec(base, Arrays.asList("13c44a9afa92ae728bf55b7075cc5de3")));
    }

    /**
     * Bug happens when reads are soft-clipped off the  contig (usually in the MT). This test guarantees no changes to the upstream code will
     * break the current hard-clipping routine that protects reduce reads from such reads.
     */
    @Test(enabled = true)
    public void testReadOffContig() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, OFFCONTIG_BAM) + " -o %s ";
        executeTest("testReadOffContig", new WalkerTestSpec(base, Arrays.asList("922be8b1151dd0d92602af93b77f7a51")));
    }

}

