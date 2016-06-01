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

package org.broadinstitute.gatk.utils.pairhmm;


// the imports for unit testing.

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class PairHMMUnitTest extends BaseTest {
    private final static boolean ALLOW_READS_LONGER_THAN_HAPLOTYPE = true;
    private final static boolean DEBUG = false;
    final static boolean EXTENSIVE_TESTING = true;
    final N2MemoryPairHMM exactHMM = new Log10PairHMM(true); // the log truth implementation
    final N2MemoryPairHMM originalHMM = new Log10PairHMM(false); // the reference implementation
    final N2MemoryPairHMM loglessHMM = new LoglessPairHMM();
    final PairHMM arrayHMM = new ArrayLoglessPairHMM();
    final N2MemoryPairHMM fastloglessHMM = new FastLoglessPairHMM((byte) 10);

    @BeforeClass
    public void initialize() {
        exactHMM.doNotUseTristateCorrection();
        originalHMM.doNotUseTristateCorrection();
        loglessHMM.doNotUseTristateCorrection();
        arrayHMM.doNotUseTristateCorrection();
        fastloglessHMM.doNotUseTristateCorrection();
    }

    private List<N2MemoryPairHMM> getHMMs() {
        return Arrays.asList(exactHMM, originalHMM, loglessHMM, fastloglessHMM);
    }

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class BasicLikelihoodTestProvider {
        final String ref, nextRef, read;
        final byte[] refBasesWithContext, nextRefBasesWithContext, readBasesWithContext;
        final int baseQual, insQual, delQual, gcp;
        final int expectedQual;
        final boolean left, right;
        final static String CONTEXT = "ACGTAATGACGATTGCA";
        final static String LEFT_FLANK = "GATTTATCATCGAGTCTGC";
        final static String RIGHT_FLANK = "CATGGATCGTTATCAGCTATCTCGAGGGATTCACTTAACAGTTTTA";

        public BasicLikelihoodTestProvider(final String ref, final String nextRef, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp ) {
            this(ref, nextRef, read, baseQual, insQual, delQual, expectedQual, gcp, false, false);
        }

        public BasicLikelihoodTestProvider(final String ref, final String nextRef, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp, final boolean left, final boolean right) {
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
            this.gcp = gcp;
            this.read = read;
            this.ref = ref;
            this.nextRef = nextRef;
            this.expectedQual = expectedQual;
            this.left = left;
            this.right = right;

            refBasesWithContext = asBytes(ref, left, right);
            nextRefBasesWithContext = asBytes(nextRef, left, right);
            readBasesWithContext = asBytes(read, false, false);
        }

        @Override
        public String toString() {
            return String.format("ref=%s nextRef=%s read=%s b/i/d/c quals = %d/%d/%d/%d l/r flank = %b/%b e[qual]=%d", ref, nextRef, read, baseQual, insQual, delQual, gcp, left, right, expectedQual);
        }

        public double expectedLogL() {
            return (expectedQual / -10.0) + 0.03 + Math.log10(1.0/refBasesWithContext.length);
        }

        public double getTolerance(final PairHMM hmm) {
            if ( hmm instanceof LoglessPairHMM || hmm instanceof ArrayLoglessPairHMM)
                return toleranceFromExact();
            if ( hmm instanceof Log10PairHMM ) {
                return ((Log10PairHMM)hmm).isDoingExactLog10Calculations() ? toleranceFromExact() : toleranceFromReference();
            } else
                return toleranceFromTheoretical();
        }

        public double toleranceFromTheoretical() {
            return 0.2;
        }

        public double toleranceFromReference() {
            return 1E-3; // has to be very tolerant -- this approximation is quite approximate
        }

        public double toleranceFromExact() {
            return 1E-9;
        }

        public double calcLogL( final PairHMM pairHMM, boolean anchorIndel ) {
            pairHMM.initialize(readBasesWithContext.length, refBasesWithContext.length);
            return pairHMM.computeReadLikelihoodGivenHaplotypeLog10(
                    refBasesWithContext, readBasesWithContext,
                    qualAsBytes(baseQual, false, anchorIndel), qualAsBytes(insQual, true, anchorIndel), qualAsBytes(delQual, true, anchorIndel),
                    qualAsBytes(gcp, false, anchorIndel), true, nextRefBasesWithContext);
        }

        private byte[] asBytes(final String bases, final boolean left, final boolean right) {
            if(bases == null)
                return null;
            else
                return ( (left ? LEFT_FLANK : "") + CONTEXT + bases + CONTEXT + (right ? RIGHT_FLANK : "")).getBytes();
        }

        private byte[] qualAsBytes(final int phredQual, final boolean doGOP, final boolean anchorIndel) {
            final byte phredQuals[] = new byte[readBasesWithContext.length];

            if( anchorIndel ) {
                // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
                Arrays.fill(phredQuals, (byte)100);

                // update just the bases corresponding to the provided micro read with the quality scores
                if( doGOP ) {
                    phredQuals[CONTEXT.length()] = (byte)phredQual;
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

        final List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {

                    // test substitutions
                    for ( final byte refBase : BaseUtils.BASES ) {
                        for ( final byte readBase : BaseUtils.BASES ) {
                            final String ref  = new String(new byte[]{refBase});
                            final String read = new String(new byte[]{readBase});
                            final int expected = refBase == readBase ? 0 : baseQual;
                            // runBasicLikelihoodTests uses calcLogL(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                            tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp)});
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

                                // runBasicLikelihoodTests uses calcLogL(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp)});
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp, true, false)});
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp, false, true)});
                                tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, expected, gcp, true, true)});
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "OptimizedLikelihoodTestProvider")
    public Object[][] makeOptimizedLikelihoodTests() {
        Utils.resetRandomGenerator();
        final Random random = Utils.getRandomGenerator();
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 30, 40, 60) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(20, 40, 60) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(3, 20, 50, 90, 160) : Arrays.asList(2);

        final List<Object[]> tests = new ArrayList<Object[]>();

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

                            for ( final boolean leftFlank : Arrays.asList(true, false) )
                                for ( final boolean rightFlank : Arrays.asList(true, false) )
                                    // runOptimizedLikelihoodTests uses calcLogL(), which runs HMM with recacheReads=true. Since we will not cache, should pass null in place of a nextRef
                                    tests.add(new Object[]{new BasicLikelihoodTestProvider(ref, null, read, baseQual, indelQual, indelQual, -0, gcp, leftFlank, rightFlank)});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "BasicLikelihoodTestProvider")
    public void testBasicLikelihoods(BasicLikelihoodTestProvider cfg) {
        if ( ALLOW_READS_LONGER_THAN_HAPLOTYPE || cfg.read.length() <= cfg.ref.length() ) {
            final double exactLogL = cfg.calcLogL( exactHMM, true );
            for ( final PairHMM hmm : getHMMs() ) {
                final double actualLogL = cfg.calcLogL( hmm, true );
                final double expectedLogL = cfg.expectedLogL();

                // compare to our theoretical expectation with appropriate tolerance
                Assert.assertEquals(actualLogL, expectedLogL, cfg.toleranceFromTheoretical(), "Failed with hmm " + hmm);
                // compare to the exact reference implementation with appropriate tolerance
                Assert.assertEquals(actualLogL, exactLogL, cfg.getTolerance(hmm), "Failed with hmm " + hmm);
                Assert.assertTrue(MathUtils.goodLog10Probability(actualLogL), "Bad log10 likelihood " + actualLogL);
            }
        }
    }

    @Test(enabled = !DEBUG, dataProvider = "OptimizedLikelihoodTestProvider")
    public void testOptimizedLikelihoods(BasicLikelihoodTestProvider cfg) {
        if ( ALLOW_READS_LONGER_THAN_HAPLOTYPE || cfg.read.length() <= cfg.ref.length() ) {
            final double exactLogL = cfg.calcLogL( exactHMM, false );

            for ( final PairHMM hmm : getHMMs() ) {
                final double calculatedLogL = cfg.calcLogL( hmm, false );
                // compare to the exact reference implementation with appropriate tolerance
                Assert.assertEquals(calculatedLogL, exactLogL, cfg.getTolerance(hmm), String.format("Test: logL calc=%.2f expected=%.2f for %s with hmm %s", calculatedLogL, exactLogL, cfg.toString(), hmm));
                Assert.assertTrue(MathUtils.goodLog10Probability(calculatedLogL), "Bad log10 likelihood " + calculatedLogL);
            }
        }
    }


    @Test(enabled = !DEBUG)
    public void testMismatchInEveryPositionInTheReadWithCenteredHaplotype() {
        final byte[] haplotype1 = "TTCTCTTCTGTTGTGGCTGGTT".getBytes();
        final byte matchQual = 90;
        final byte mismatchQual = 20;
        final byte indelQual = 80;
        final int offset = 2;
        final byte[] gop = new byte[haplotype1.length - 2 * offset];
        Arrays.fill(gop, indelQual);
        final byte[] gcp = new byte[haplotype1.length - 2 * offset];
        Arrays.fill(gcp, indelQual);
        loglessHMM.initialize(gop.length, haplotype1.length);

        for( int k = 0; k < haplotype1.length - 2 * offset; k++ ) {
            final byte[] quals = new byte[haplotype1.length - 2 * offset];
            Arrays.fill(quals, matchQual);
            // one base mismatches the haplotype
            quals[k] = mismatchQual;
            final byte[] mread = Arrays.copyOfRange(haplotype1,offset,haplotype1.length-offset);
            // change single base at position k to C. If it's a C, change to T
            mread[k] = ( mread[k] == (byte)'C' ? (byte)'T' : (byte)'C');
            final double res1 = loglessHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotype1, mread, quals, gop, gop, gcp, true, null);
            final double expected = Math.log10(1.0/haplotype1.length * Math.pow(QualityUtils.qualToProb(matchQual), mread.length-1) * QualityUtils.qualToErrorProb(mismatchQual));
            Assert.assertEquals(res1, expected, 1e-2);
        }
    }

    @Test(enabled = ! DEBUG)
    public void testMismatchInEveryPositionInTheRead() {
        final byte[] haplotype1 = "TTCTCTTCTGTTGTGGCTGGTT".getBytes();
        final byte matchQual = 90;
        final byte mismatchQual = 20;
        final byte indelQual = 80;

        final int offset = 2;
        final byte[] gop = new byte[haplotype1.length - offset];
        Arrays.fill(gop, indelQual);
        final byte[] gcp = new byte[haplotype1.length - offset];
        Arrays.fill(gcp, indelQual);
        loglessHMM.initialize(gop.length, haplotype1.length);

        for( int k = 0; k < haplotype1.length - offset; k++ ) {
            final byte[] quals = new byte[haplotype1.length - offset];
            Arrays.fill(quals, matchQual);
            // one base mismatches the haplotype with low qual
            quals[k] = mismatchQual;
            final byte[] mread = Arrays.copyOfRange(haplotype1,offset,haplotype1.length);
            // change single base at position k to C. If it's a C, change to T
            mread[k] = ( mread[k] == (byte)'C' ? (byte)'T' : (byte)'C');
            final double res1 = loglessHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotype1, mread, quals, gop, gop, gcp, true , null);
            final double expected = Math.log10(1.0/haplotype1.length * Math.pow(QualityUtils.qualToProb(matchQual), mread.length-1) * QualityUtils.qualToErrorProb(mismatchQual));
            Assert.assertEquals(res1, expected, 1e-2);
        }
    }

    @DataProvider(name = "HMMProvider")
    public Object[][] makeHMMProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int readSize : Arrays.asList(1, 2, 5, 10) ) {
            for ( final int refSize : Arrays.asList(1, 2, 5, 10) ) {
                if ( refSize > readSize ) {
                    for ( final PairHMM hmm : getHMMs() )
                        tests.add(new Object[]{hmm, readSize, refSize});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "HMMProvider")
    void testMultipleReadMatchesInHaplotype(final PairHMM hmm, final int readSize, final int refSize) {
        final byte[] readBases =  Utils.dupBytes((byte)'A', readSize);
        final byte[] refBases = ("CC" + new String(Utils.dupBytes((byte)'A', refSize)) + "GGA").getBytes();
        final byte baseQual = 20;
        final byte insQual = 37;
        final byte delQual = 37;
        final byte gcp = 10;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        final double d = hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
        Assert.assertTrue(d <= 0.0, "Likelihoods should be <= 0 but got "+ d);
    }

    @Test(enabled = !DEBUG, dataProvider = "HMMProvider")
    void testAllMatchingRead(final PairHMM hmm, final int readSize, final int refSize) {
        final byte[] readBases =  Utils.dupBytes((byte)'A', readSize);
        final byte[] refBases = Utils.dupBytes((byte)'A', refSize);
        final byte baseQual = 20;
        final byte insQual = 100;
        final byte delQual = 100;
        final byte gcp = 100;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
        double expected =  0;
        final double initialCondition = ((double) Math.abs(refBases.length-readBases.length+1))/refBases.length;
        if (readBases.length < refBases.length) {
            expected = Math.log10(initialCondition * Math.pow(QualityUtils.qualToProb(baseQual), readBases.length));
        } else if (readBases.length > refBases.length) {
            expected = Math.log10(initialCondition * Math.pow(QualityUtils.qualToProb(baseQual), refBases.length) * Math.pow(QualityUtils.qualToErrorProb(insQual), readBases.length - refBases.length));
        }
        Assert.assertEquals(d, expected, 1e-3, "Likelihoods should sum to just the error prob of the read " + String.format("readSize=%d refSize=%d", readSize, refSize));
    }

    @DataProvider(name = "HMMProviderWithBigReads")
    public Object[][] makeBigReadHMMProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String read1 = "ACCAAGTAGTCACCGT";
        final String ref1  = "ACCAAGTAGTCACCGTAACG";

        for ( final int nReadCopies : Arrays.asList(1, 2, 10, 20, 50) ) {
            for ( final int nRefCopies : Arrays.asList(1, 2, 10, 20, 100) ) {
                if ( nRefCopies > nReadCopies ) {
                    for ( final PairHMM hmm : getHMMs() ) {
                        final String read = Utils.dupString(read1, nReadCopies);
                        final String ref  = Utils.dupString(ref1, nRefCopies);
                        tests.add(new Object[]{hmm, read, ref});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "HMMProviderWithBigReads")
    void testReallyBigReads(final PairHMM hmm, final String read, final String ref) {
        final byte[] readBases =  read.getBytes();
        final byte[] refBases = ref.getBytes();
        final byte baseQual = 30;
        final byte insQual = 40;
        final byte delQual = 40;
        final byte gcp = 10;
        hmm.initialize(readBases.length, refBases.length);
        // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
        hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
    }

    @Test(enabled = !DEBUG)
    void testPreviousBadValue() {
        final byte[] readBases = "A".getBytes();
        final byte[] refBases =  "AT".getBytes();
        final byte baseQual = 30;
        final byte insQual = 40;
        final byte delQual = 40;
        final byte gcp = 10;
        // running HMMs with no haplotype caching. Should therefore pass null in place of nextRef bases
        exactHMM.initialize(readBases.length, refBases.length);
        exactHMM.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);

        loglessHMM.initialize(readBases.length, refBases.length);
        loglessHMM.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);

        arrayHMM.initialize(readBases.length, refBases.length);
        arrayHMM.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                Utils.dupBytes(baseQual, readBases.length),
                Utils.dupBytes(insQual, readBases.length),
                Utils.dupBytes(delQual, readBases.length),
                Utils.dupBytes(gcp, readBases.length), true, null);
    }

    @DataProvider(name = "JustHMMProvider")
    public Object[][] makeJustHMMProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final PairHMM hmm : getHMMs() ) {
            tests.add(new Object[]{hmm});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "JustHMMProvider")
    void testMaxLengthsBiggerThanProvidedRead(final PairHMM hmm) {
        final byte[] readBases = "CTATCTTAGTAAGCCCCCATACCTGCAAATTTCAGGATGTCTCCTCCAAAAATCAACA".getBytes();
        final byte[] refBases =  "CTATCTTAGTAAGCCCCCATACCTGCAAATTTCAGGATGTCTCCTCCAAAAATCAAAACTTCTGAGAAAAAAAAAAAAAATTAAATCAAACCCTGATTCCTTAAAGGTAGTAAAAAAACATCATTCTTTCTTAGTGGAATAGAAACTAGGTCAAAAGAACAGTGATTC".getBytes();

        final byte[] quals = new byte[]{35,34,31,32,35,34,32,31,36,30,31,32,36,34,33,32,32,32,33,32,30,35,33,35,36,36,33,33,33,32,32,32,37,33,36,35,33,32,34,31,36,35,35,35,35,33,34,31,31,30,28,27,26,29,26,25,29,29};
        final byte[] insQual = new byte[]{46,46,46,46,46,47,45,46,45,48,47,44,45,48,46,43,43,42,48,48,45,47,47,48,48,47,48,45,38,47,45,39,47,48,47,47,48,46,49,48,49,48,46,47,48,44,44,43,39,32,34,36,46,48,46,44,45,45};
        final byte[] delQual = new byte[]{44,44,44,43,45,44,43,42,45,46,45,43,44,47,45,40,40,40,45,46,43,45,45,44,46,46,46,43,35,44,43,36,44,45,46,46,44,44,47,43,47,45,45,45,46,45,45,46,44,35,35,35,45,47,45,44,44,43};
        final byte[] gcp = Utils.dupBytes((byte) 10, delQual.length);
        hmm.initialize(readBases.length + 100, refBases.length + 100);
        for ( int nExtraMaxSize = 0; nExtraMaxSize < 100; nExtraMaxSize++ ) {
            // running HMM with no haplotype caching. Should therefore pass null in place of nextRef bases
            hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases, quals, insQual, delQual, gcp, true, null);
        }
    }

    @DataProvider(name = "HaplotypeIndexingProvider")
    public Object[][] makeHaplotypeIndexingProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // First difference (root2, root3) is the base position immediately following first difference (root1, root2)
        final String root1    = "ACGTGTCAAACCGGGTT";
        final String root2    = "ACGTGTCACACTGGGTT"; // differs in two locations from root1
        final String root3    = "ACGTGTCACTCCGCGTT"; // differs in two locations from root2

        final String read1    = "ACGTGTCACACTGGATT"; // 1 diff from 2, 2 diff from root1, 2 diff from root3
        final String read2    = root1; // same as root1
        final String read3    = root2; // same as root2
        final String read4    = "ACGTGTCACACTGGATTCGAT";
        final String read5    = "CCAGTAACGTGTCACACTGGATTCGAT";

//        for ( final String read : Arrays.asList(read2) ) {
        for ( final String read : Arrays.asList(read1, read2, read3, read4, read5) ) {
            for ( final PairHMM hmm : getHMMs() ) {
//                int readLength = read.length(); {
                for ( int readLength = 10; readLength < read.length(); readLength++ ) {
                    final String myRead = read.substring(0, readLength);
                    tests.add(new Object[]{hmm, root1, root2, root3, myRead});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "HaplotypeIndexingProvider")
    void testHaplotypeIndexing(final PairHMM hmm, final String root1, final String root2, final String root3, final String read) {
        final double TOLERANCE = 1e-9;
        final String prefix   = "AACCGGTTTTTGGGCCCAAACGTACGTACAGTTGGTCAACATCGATCAGGTTCCGGAGTAC";

        final int maxReadLength = read.length();
        final int maxHaplotypeLength = prefix.length() + root1.length();

        // the initialization occurs once, at the start of the evalution of reads
        hmm.initialize(maxReadLength, maxHaplotypeLength);

        for ( int prefixStart = prefix.length(); prefixStart >= 0; prefixStart-- ) {
            final String myPrefix = prefix.substring(prefixStart, prefix.length());
            final String hap1 = myPrefix + root1;
            final String hap2 = myPrefix + root2;
            final String hap3 = myPrefix + root3;

            final int hapStart = PairHMM.findFirstPositionWhereHaplotypesDiffer(hap1.getBytes(), hap2.getBytes());

            // Run the HMM on the first haplotype, peaking ahead the second, to set up caching
            // Then run on the second haplotype in both cached and uncached mode, and verify that results are the same
            // When evaluating actual2, it is important that we both apply old caching from hap1 and set up new caching for hap3, to ensure read/write operations do not cause conflicts
            final double actual1 = testHaplotypeIndexingCalc(hmm, hap1, hap2, read, 0, true);
            final double actual2 = testHaplotypeIndexingCalc(hmm, hap2, hap3, read, hapStart, false);
            final double expected2 = testHaplotypeIndexingCalc(hmm, hap2, null, read, 0, true);
            Assert.assertEquals(actual2, expected2, TOLERANCE, "HMM " + hmm.getClass() + " Caching calculation failed for read " + read + " against haplotype with prefix '" + myPrefix
                    + "' expected " + expected2 + " but got " + actual2 + " with hapStart of " + hapStart);
        }
    }

    private double testHaplotypeIndexingCalc(final PairHMM hmm, final String hap, final String nextHap, final String read, final int hapStart, final boolean recache) {
        final byte[] readBases = read.getBytes();
        // if not peaking ahead to capture info for a future cache run, the next haplotype will be null, and this should be passed to HMM
        final byte[] nextHapBases = nextHap == null ? null : nextHap.getBytes();
        final byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);
        final byte[] insQuals = Utils.dupBytes((byte)45, readBases.length);
        final byte[] delQuals = Utils.dupBytes((byte)40, readBases.length);
        final byte[] gcp = Utils.dupBytes((byte)10, readBases.length);
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10(hap.getBytes(), readBases, baseQuals, insQuals, delQuals, gcp, recache, nextHapBases);
        Assert.assertTrue(MathUtils.goodLog10Probability(d), "Likelihoods = " + d + " was bad for read " + read + " and ref " + hap + " with hapStart " + hapStart);
        return d;
    }

    @Test(enabled = !DEBUG)
    public void testFindFirstPositionWhereHaplotypesDiffer() {
        for ( int haplotypeSize1 = 10; haplotypeSize1 < 30; haplotypeSize1++ ) {
            for ( int haplotypeSize2 = 10; haplotypeSize2 < 50; haplotypeSize2++ ) {
                final int maxLength = Math.max(haplotypeSize1, haplotypeSize2);
                final int minLength = Math.min(haplotypeSize1, haplotypeSize2);
                for ( int differingSite = 0; differingSite < maxLength + 1; differingSite++) {
                    for ( final boolean oneIsDiff : Arrays.asList(true, false) ) {
                        final byte[] hap1 = Utils.dupBytes((byte)'A', haplotypeSize1);
                        final byte[] hap2 = Utils.dupBytes((byte)'A', haplotypeSize2);
                        final int expected = oneIsDiff
                                ? makeDiff(hap1, differingSite, minLength)
                                : makeDiff(hap2, differingSite, minLength);
                        final int actual = PairHMM.findFirstPositionWhereHaplotypesDiffer(hap1, hap2);
                        Assert.assertEquals(actual, expected, "Bad differing site for " + new String(hap1) + " vs. " + new String(hap2));
                    }
                }
            }
        }
    }

    private int makeDiff(final byte[] bytes, final int site, final int minSize) {
        if ( site < bytes.length ) {
            bytes[site] = 'C';
            return Math.min(site, minSize);
        } else
            return minSize;
    }

    @DataProvider(name = "UninitializedHMMs")
    public Object[][] makeUninitializedHMMs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final LoglessPairHMM myLoglessPairHMM = new LoglessPairHMM();
        myLoglessPairHMM.doNotUseTristateCorrection();
        tests.add(new Object[]{myLoglessPairHMM});

        final ArrayLoglessPairHMM myArrayLoglessPairHMM = new ArrayLoglessPairHMM();
        myArrayLoglessPairHMM.doNotUseTristateCorrection();
        tests.add(new Object[]{myArrayLoglessPairHMM});

        final Log10PairHMM myLog10PairHMM = new Log10PairHMM(true);
        myLog10PairHMM.doNotUseTristateCorrection();
        tests.add(new Object[]{myLog10PairHMM});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, expectedExceptions = IllegalStateException.class, dataProvider = "UninitializedHMMs")
    public void testNoInitializeCall(final PairHMM hmm) {
        byte[] readBases = "A".getBytes();
        byte[] refBases =  "AT".getBytes();
        byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);

        // didn't call initialize => should exception out
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                baseQuals, baseQuals, baseQuals, baseQuals, true, null);
    }

    @Test(enabled = true, expectedExceptions = IllegalArgumentException.class, dataProvider = "JustHMMProvider")
    public void testHapTooLong(final PairHMM hmm) {
        byte[] readBases = "AAA".getBytes();
        byte[] refBases =  "AAAT".getBytes();
        byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);

        hmm.initialize(3, 3);
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                baseQuals, baseQuals, baseQuals, baseQuals, true, null);
    }

    @Test(enabled = true, expectedExceptions = IllegalArgumentException.class, dataProvider = "JustHMMProvider")
    public void testReadTooLong(final PairHMM hmm) {
        byte[] readBases = "AAA".getBytes();
        byte[] refBases =  "AAAT".getBytes();
        byte[] baseQuals = Utils.dupBytes((byte)30, readBases.length);

        hmm.initialize(2, 3);
        double d = hmm.computeReadLikelihoodGivenHaplotypeLog10( refBases, readBases,
                baseQuals, baseQuals, baseQuals, baseQuals, true, null);
    }
}