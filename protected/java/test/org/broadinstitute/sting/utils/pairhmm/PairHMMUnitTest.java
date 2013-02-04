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

package org.broadinstitute.sting.utils.pairhmm;


// the imports for unit testing.

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class PairHMMUnitTest extends BaseTest {
    final static boolean EXTENSIVE_TESTING = true;
    PairHMM exactHMM = new ExactPairHMM(); // the log truth implementation
    PairHMM originalHMM = new OriginalPairHMM(); // the reference implementation
    PairHMM cachingHMM = new CachingPairHMM();
    PairHMM loglessHMM = new LoglessCachingPairHMM();

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class BasicLikelihoodTestProvider extends TestDataProvider {
        final String ref, read;
        final byte[] refBasesWithContext, readBasesWithContext;
        final int baseQual, insQual, delQual, gcp;
        final int expectedQual;
        final static String CONTEXT = "ACGTAATGACGATTGCA";
        final static String LEFT_FLANK = "GATTTATCATCGAGTCTGC";
        final static String RIGHT_FLANK = "CATGGATCGTTATCAGCTATCTCGAGGGATTCACTTAACAGTTTTA";

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp ) {
            this(ref, read, baseQual, insQual, delQual, expectedQual, gcp, false, false);
        }

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp, final boolean left, final boolean right) {
            super(BasicLikelihoodTestProvider.class, String.format("ref=%s read=%s b/i/d/c quals = %d/%d/%d/%d l/r flank = %b/%b e[qual]=%d", ref, read, baseQual, insQual, delQual, gcp, left, right, expectedQual));
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
            this.gcp = gcp;
            this.read = read;
            this.ref = ref;
            this.expectedQual = expectedQual;

            refBasesWithContext = asBytes(ref, left, right);
            readBasesWithContext = asBytes(read, false, false);
        }

        public double expectedLogL() {
            return (expectedQual / -10.0) + 0.03 ;
        }

        public double toleranceFromTheoretical() {
            return 0.2;
        }

        public double toleranceFromReference() {
            return 1E-4;
        }

        public double toleranceFromExact() {
            return 1E-9;
        }

        public double calcLogL( final PairHMM pairHMM, boolean anchorIndel ) {
            pairHMM.initialize(readBasesWithContext.length, refBasesWithContext.length);
            return pairHMM.computeReadLikelihoodGivenHaplotypeLog10(
                    refBasesWithContext, readBasesWithContext,
                    qualAsBytes(baseQual, false, anchorIndel), qualAsBytes(insQual, true, anchorIndel), qualAsBytes(delQual, true, anchorIndel),
                    qualAsBytes(gcp, false, anchorIndel), 0, true);
        }

        private final byte[] asBytes(final String bases, final boolean left, final boolean right) {
            return ( (left ? LEFT_FLANK : "") + CONTEXT + bases + CONTEXT + (right ? RIGHT_FLANK : "")).getBytes();
        }

        private byte[] qualAsBytes(final int phredQual, final boolean doGOP, final boolean anchorIndel) {
            final byte phredQuals[] = new byte[readBasesWithContext.length];

            if( anchorIndel ) {
                // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
                Arrays.fill(phredQuals, (byte)100);

                // update just the bases corresponding to the provided micro read with the quality scores
                if( doGOP ) {
                    phredQuals[0 + CONTEXT.length()] = (byte)phredQual;
                } else {
                    for ( int i = 0; i < read.length(); i++)
                        phredQuals[i + CONTEXT.length()] = (byte)phredQual;
                }
            } else {
                Arrays.fill(phredQuals, (byte)phredQual);
            }

            return phredQuals;
        }
    }

    @DataProvider(name = "BasicLikelihoodTestProvider")
    public Object[][] makeBasicLikelihoodTests() {
        // context on either side is ACGTTGCA REF ACGTTGCA
        // test all combinations
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30, 40, 50) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(20, 30, 40, 50) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(8, 10, 20) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(2,3,4,5,7,8,9,10,20,30,35) : Arrays.asList(2);

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {

                    // test substitutions
                    for ( final byte refBase : BaseUtils.BASES ) {
                        for ( final byte readBase : BaseUtils.BASES ) {
                            final String ref  = new String(new byte[]{refBase});
                            final String read = new String(new byte[]{readBase});
                            final int expected = refBase == readBase ? 0 : baseQual;
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                        }
                    }

                    // test insertions and deletions
                    for ( final int size : sizes ) {
                        for ( final byte base : BaseUtils.BASES ) {
                            final int expected = indelQual + (size - 2) * gcp;

                            for ( boolean insertionP : Arrays.asList(true, false)) {
                                final String small = Utils.dupString((char)base, 1);
                                final String big = Utils.dupString((char) base, size);

                                final String ref = insertionP ? small : big;
                                final String read = insertionP ? big : small;

                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, true, false);
                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, false, true);
                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, true, true);
                            }
                        }
                    }
                }
            }
        }

        return BasicLikelihoodTestProvider.getTests(BasicLikelihoodTestProvider.class);
    }

    final Random random = new Random(87860573);
    @DataProvider(name = "OptimizedLikelihoodTestProvider")
    public Object[][] makeOptimizedLikelihoodTests() {
        // context on either side is ACGTTGCA REF ACGTTGCA
        // test all combinations
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 30, 40, 60) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(20, 40, 60) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(3, 20, 50, 90, 160) : Arrays.asList(2);

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {
                    for ( final int refSize : sizes ) {
                        for ( final int readSize : sizes ) {
                            String ref = "";
                            String read = "";
                            for( int iii = 0; iii < refSize; iii++) {
                                ref += (char) BaseUtils.BASES[random.nextInt(4)];
                            }
                            for( int iii = 0; iii < readSize; iii++) {
                                read += (char) BaseUtils.BASES[random.nextInt(4)];
                            }
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, -0, gcp);
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, -0, gcp, true, false);
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, -0, gcp, false, true);
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, -0, gcp, true, true);
                        }
                    }
                }
            }
        }

        return BasicLikelihoodTestProvider.getTests(BasicLikelihoodTestProvider.class);
    }

    @Test(dataProvider = "BasicLikelihoodTestProvider", enabled = true)
    public void testBasicLikelihoods(BasicLikelihoodTestProvider cfg) {
        double exactLogL = cfg.calcLogL( exactHMM, true );
        double calculatedLogL = cfg.calcLogL( originalHMM, true );
        double optimizedLogL = cfg.calcLogL( cachingHMM, true );
        double loglessLogL = cfg.calcLogL( loglessHMM, true );
        double expectedLogL = cfg.expectedLogL();
        //logger.warn(String.format("Test: logL calc=%.2f optimized=%.2f logless=%.2f expected=%.2f for %s", calculatedLogL, optimizedLogL, loglessLogL, expectedLogL, cfg.toString()));
        Assert.assertEquals(exactLogL, expectedLogL, cfg.toleranceFromTheoretical());
        Assert.assertEquals(calculatedLogL, expectedLogL, cfg.toleranceFromTheoretical());
        Assert.assertEquals(optimizedLogL, calculatedLogL, cfg.toleranceFromReference());
        Assert.assertEquals(loglessLogL, exactLogL, cfg.toleranceFromExact());
    }

    @Test(dataProvider = "OptimizedLikelihoodTestProvider", enabled = true)
    public void testOptimizedLikelihoods(BasicLikelihoodTestProvider cfg) {
        double exactLogL = cfg.calcLogL( exactHMM, false );
        double calculatedLogL = cfg.calcLogL( originalHMM, false );
        double optimizedLogL = cfg.calcLogL( cachingHMM, false );
        double loglessLogL = cfg.calcLogL( loglessHMM, false );
        //logger.warn(String.format("Test: logL calc=%.2f optimized=%.2f logless=%.2f expected=%.2f for %s", calculatedLogL, optimizedLogL, loglessLogL, expectedLogL, cfg.toString()));
        Assert.assertEquals(optimizedLogL, calculatedLogL, cfg.toleranceFromReference());
        Assert.assertEquals(loglessLogL, exactLogL, cfg.toleranceFromExact());
    }

    @Test
    public void testMismatchInEveryPositionInTheReadWithCenteredHaplotype() {
        byte[] haplotype1 = "TTCTCTTCTGTTGTGGCTGGTT".getBytes();

        final int offset = 2;
        byte[] gop = new byte[haplotype1.length - 2 * offset];
        Arrays.fill(gop, (byte) 80);
        byte[] gcp = new byte[haplotype1.length - 2 * offset];
        Arrays.fill(gcp, (byte) 80);

        for( int k = 0; k < haplotype1.length - 2 * offset; k++ ) {
            byte[] quals = new byte[haplotype1.length - 2 * offset];
            Arrays.fill(quals, (byte) 90);
            // one read mismatches the haplotype
            quals[k] = 20;

            byte[] mread = Arrays.copyOfRange(haplotype1,offset,haplotype1.length-offset);
            // change single base at position k to C. If it's a C, change to T
            mread[k] = ( mread[k] == (byte)'C' ? (byte)'T' : (byte)'C');
            originalHMM.initialize(mread.length, haplotype1.length);
            double res1 = originalHMM.computeReadLikelihoodGivenHaplotypeLog10(
                    haplotype1, mread,
                    quals, gop, gop,
                    gcp, 0, false);

            System.out.format("H:%s\nR:  %s\n Pos:%d Result:%4.2f\n",new String(haplotype1), new String(mread), k,res1);

            Assert.assertEquals(res1, -2.0, 1e-2);
        }
    }

    @Test
    public void testMismatchInEveryPositionInTheRead() {
        byte[] haplotype1 = "TTCTCTTCTGTTGTGGCTGGTT".getBytes();

        final int offset = 2;
        byte[] gop = new byte[haplotype1.length - offset];
        Arrays.fill(gop, (byte) 80);
        byte[] gcp = new byte[haplotype1.length - offset];
        Arrays.fill(gcp, (byte) 80);

        for( int k = 0; k < haplotype1.length - offset; k++ ) {
            byte[] quals = new byte[haplotype1.length - offset];
            Arrays.fill(quals, (byte) 90);
            // one read mismatches the haplotype
            quals[k] = 20;

            byte[] mread = Arrays.copyOfRange(haplotype1,offset,haplotype1.length);
            // change single base at position k to C. If it's a C, change to T
            mread[k] = ( mread[k] == (byte)'C' ? (byte)'T' : (byte)'C');
            originalHMM.initialize(mread.length, haplotype1.length);
            double res1 = originalHMM.computeReadLikelihoodGivenHaplotypeLog10(
                    haplotype1, mread,
                    quals, gop, gop,
                    gcp, 0, false);

            System.out.format("H:%s\nR:  %s\n Pos:%d Result:%4.2f\n",new String(haplotype1), new String(mread), k,res1);

            Assert.assertEquals(res1, -2.0, 1e-2);
        }
    }
}