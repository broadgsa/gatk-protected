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
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final static String DBSNP = b37dbSNP132;
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
    final String BOTH_ENDS_OF_PAIR_IN_VARIANT_REGION_BAM = privateTestDir + "bothEndsOfPairInVariantRegion.bam";
    final String INSERTIONS_AT_EDGE_OF_CONSENSUS_BAM = privateTestDir + "rr-too-many-insertions.bam";

    final static String emptyFileMd5 = "d41d8cd98f00b204e9800998ecf8427e";

    protected Pair<List<File>, List<String>> executeTest(final String name, final WalkerTestSpec spec) {
        return executeTest(name, spec, emptyFileMd5);
    }

    protected Pair<List<File>, List<String>> executeTest(final String name, final WalkerTestSpec spec, final String qualsTestMD5) {
        final Pair<List<File>, List<String>> result = super.executeTest(name, spec);

        // perform some Reduce Reads specific testing now
        if ( result != null ) {

            // generate a new command-line based on the old one
            spec.disableImplicitArgs();
            final String[] originalArgs = spec.getArgsWithImplicitArgs().split(" ");

            final StringBuilder reducedInputs = new StringBuilder();
            for ( final File file : result.getFirst() ) {
                reducedInputs.append(" -I:reduced ");
                reducedInputs.append(file.getAbsolutePath());
            }

            // the coverage test is a less stricter version of the quals test so we can safely ignore it for now
            //final String coverageCommand = createCommandLine("AssessReducedCoverage", originalArgs);
            //super.executeTest(name + " : COVERAGE_TEST", new WalkerTestSpec(coverageCommand + reducedInputs.toString(), Arrays.asList(emptyFileMd5)));

            // run the quals test
            final String qualsCommand = createCommandLine("AssessReducedQuals", originalArgs);
            super.executeTest(name + " : QUALS_TEST", new WalkerTestSpec(qualsCommand + reducedInputs.toString(), Arrays.asList(qualsTestMD5)));
        }

        return result;
    }

    /*
     * Generate a new command-line based on the old one
     *
     * @param walkerName    the new walker name to use
     * @param originalArgs  the original arguments used for the test
     * @return the new command line
     */
    private String createCommandLine(final String walkerName, final String[] originalArgs) {

        final StringBuilder newArgs = new StringBuilder();

        for ( int i = 0; i < originalArgs.length; i++ ) {
            final String arg = originalArgs[i];

            if ( arg.equals("-T") ) {
                newArgs.append("-T ");
                newArgs.append(walkerName);
            } else if ( arg.startsWith("-I") ) {
                newArgs.append("-I:original ");
                newArgs.append(originalArgs[++i]);
            } else if ( arg.equals("-R") || arg.equals("-L") ) {
                newArgs.append(arg);
                newArgs.append(" ");
                newArgs.append(originalArgs[++i]);
            }

            // always add a trailing space
            newArgs.append(" ");
        }

        newArgs.append("-o %s");

        return newArgs.toString();
    }

    protected Pair<List<File>, List<String>> executeTestWithoutAdditionalRRTests(final String name, final WalkerTestSpec spec) {
        return super.executeTest(name, spec);
    }

    private void RRTest(final String testName, final String args, final String md5, final boolean useKnowns) {
        this.RRTest(testName, args, md5, useKnowns, emptyFileMd5);
    }

    private void RRTest(final String testName, final String args, final String md5, final boolean useKnowns, final String qualsTestMD5) {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, BAM) + " -o %s" + (useKnowns ? " -known " + DBSNP : "") + " ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList("bam"), Arrays.asList(md5));
        executeTest(testName, spec, qualsTestMD5);
    }

        @Test(enabled = true)
    public void testDefaultCompression() {
        RRTest("testDefaultCompression ", L, "fa1cffc4539e0c20b818a11da5dba5b9", false);
    }

    @Test(enabled = true)
    public void testDefaultCompressionWithKnowns() {
        RRTest("testDefaultCompressionWithKnowns ", L, "d1b5fbc402810d9cdc020bb3503f1325", true);
    }

    private final String intervals = "-L 20:10,100,000-10,100,500 -L 20:10,200,000-10,200,500 -L 20:10,300,000-10,300,500 -L 20:10,400,000-10,500,000 -L 20:10,500,050-10,500,060 -L 20:10,600,000-10,600,015 -L 20:10,700,000-10,700,110";

    @Test(enabled = true)
    public void testMultipleIntervals() {
        RRTest("testMultipleIntervals ", intervals, "7e9dcd157ad742d4ebae7e56bc4af663", false);
    }

    @Test(enabled = true)
    public void testMultipleIntervalsWithKnowns() {
        RRTest("testMultipleIntervalsWithKnowns ", intervals, "dbb1e95e1bcad956701142afac763717", true);
    }

    @Test(enabled = true)
    public void testHighCompression() {
        RRTest("testHighCompression ", " -cs 10 -min_pvalue 0.3 -minvar 0.3 -mindel 0.3 " + L, "8f8fd1a53fa0789116f45e4cf2625906", false);
    }

    @Test(enabled = true)
    public void testHighCompressionWithKnowns() {
        RRTest("testHighCompressionWithKnowns ", " -cs 10 -min_pvalue 0.3 -minvar 0.3 -mindel 0.3 " + L, "52fd2a77802a4677b604abb18e15d96a", true);
    }

    @Test(enabled = true)
    public void testLowCompression() {
        RRTest("testLowCompression ", " -cs 30 -min_pvalue 0.001 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "79c6543d5ce84ebc2ca74404498edbd1", false);
    }

    @Test(enabled = true)
    public void testLowCompressionWithKnowns() {
        RRTest("testLowCompressionWithKnowns ", " -cs 30 -min_pvalue 0.001 -minvar 0.01 -mindel 0.01 -minmap 5 -minqual 5 " + L, "271aec358b309603291a974b5ba3bd60", true);
    }

    @Test(enabled = true)
    public void testBadPvalueInput() {
        final String cmd = String.format("-T ReduceReads -npt -R %s -I %s ", REF, BAM) + "-o %s -min_pvalue -0.01";
        WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, UserException.BadArgumentValue.class);
        executeTest("testBadPvalueInput", spec);
    }

    @Test(enabled = true)
    public void testIndelCompression() {
        final String md5 = "d20e6012300898a0315c795cab7583d8";
        RRTest("testIndelCompression ", " -cs 50 -L 20:10,100,500-10,100,600 ", md5, false);
        RRTest("testIndelCompressionWithKnowns ", " -cs 50 -L 20:10,100,500-10,100,600 ", md5, true);
    }

    @Test(enabled = true)
    public void testFilteredDeletionCompression() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, DELETION_BAM) + " -o %s ";
        executeTest("testFilteredDeletionCompression", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("e5da09662708f562c0c617ba73cf4763")), "4f916da29d91852077f0a2fdbdd2c7f6");
    }

    private static final String COREDUCTION_QUALS_TEST_MD5 = "26d84a2bd549a01a63fcebf8847a1b7d";

    @Test(enabled = true)
    public void testCoReduction() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s -I %s", COREDUCTION_L, REF, COREDUCTION_BAM_A, COREDUCTION_BAM_B) + " -o %s ";
        executeTest("testCoReduction", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("5f4d2c1d9c010dfd6865aeba7d0336fe")), COREDUCTION_QUALS_TEST_MD5);
    }

    @Test(enabled = true)
    public void testCoReductionWithKnowns() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s -I %s -known %s", COREDUCTION_L, REF, COREDUCTION_BAM_A, COREDUCTION_BAM_B, DBSNP) + " -o %s ";
        executeTest("testCoReductionWithKnowns", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("ca48dd972bf57595c691972c0f887cb4")), COREDUCTION_QUALS_TEST_MD5);
    }

    @Test(enabled = true)
    public void testInsertionsAtEdgeOfConsensus() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, INSERTIONS_AT_EDGE_OF_CONSENSUS_BAM) + " -o %s ";
        executeTest("testInsertionsAtEdgeOfConsensus", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("760500a5b036b987f84099f45f26a804")));
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
        executeTest("testAddingReadAfterTailingTheStash", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("67f8a3a647f8ec5212104bdaafd8c862")), "3eab32c215ba68e75efd5ab7e9f7a2e7");
    }

    /**
     * Divide by zero bug reported by GdA and users in the forum. Happens when the downsampler goes over a region where all reads get
     * filtered out.
     */
    @Test(enabled = true)
    public void testDivideByZero() {
        String base = String.format("-T ReduceReads %s -npt -R %s -I %s", DIVIDEBYZERO_L, REF, DIVIDEBYZERO_BAM) + " -o %s ";
        // we expect to lose coverage due to the downsampling so don't run the systematic tests
        executeTestWithoutAdditionalRRTests("testDivideByZero", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("c459a6153a17c2cbf8441e1918fda9c8")));
    }

    /**
     * Bug happens when reads are soft-clipped off the contig (usually in the MT). This test guarantees no changes to the upstream code will
     * break the current hard-clipping routine that protects reduce reads from such reads.
     */
    @Test(enabled = true)
    public void testReadOffContig() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", REF, OFFCONTIG_BAM) + " -o %s ";
        executeTest("testReadOffContig", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("0ce693b4ff925998867664e4099f3248")));
    }

    /**
     * Confirm that if both ends of pair are in same variant region, compressed names of both ends of pair are the same.
     */
    @Test(enabled = true)
    public void testPairedReadsInVariantRegion() {
        String base = String.format("-T ReduceReads -npt -R %s -I %s ", hg19Reference, BOTH_ENDS_OF_PAIR_IN_VARIANT_REGION_BAM) +
                " -o %s  --downsample_coverage 250 -dcov 50  ";
        executeTest("testPairedReadsInVariantRegion", new WalkerTestSpec(base, Arrays.asList("bam"), Arrays.asList("7e7b358443827ca239db3b98f299aec6")), "2af063d1bd3c322b03405dbb3ecf59a9");
    }
}

