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

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;


/**
 * Unit tests for {@link PairHMMModel}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class PairHMMModelUnitTest extends BaseTest {

    final double TOLERANCE = 1E-9;

    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbs(final int insQual, final int delQual, final int gcp, final double[] expected) {
        final double[] actual = PairHMMModel.qualToTransProbs((byte)insQual,(byte)delQual,(byte)gcp);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, PairHMMModel.TRANS_PROB_ARRAY_LENGTH);
        assertEqualsDoubleArray(actual,expected,TOLERANCE);
        Assert.assertEquals(actual.length, PairHMMModel.TRANS_PROB_ARRAY_LENGTH);
    }

    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbsLog10(final int insQuals, final int delQual, final int gcp, final double[] expected) {
        final double[] logExpected = new double[expected.length];
        for (int i = 0; i < logExpected.length; i++)
            logExpected[i] = Math.log10(expected[i]);
        final double[] actual = PairHMMModel.qualToTransProbsLog10((byte)insQuals,(byte)delQual,(byte)gcp);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length, PairHMMModel.TRANS_PROB_ARRAY_LENGTH);
        assertEqualsDoubleArray(actual,logExpected,TOLERANCE);
    }

    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbsFill(final int insQual, final int delQual, final int gcp, final double[] expected) {
        final double[] actual = new double[PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
        PairHMMModel.qualToTransProbs(actual, (byte) insQual, (byte) delQual, (byte) gcp);
        assertEqualsDoubleArray(actual,expected,TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbs(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.qualToTransProbs(insQuals,delQuals,gapQuals);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length,expected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length,expected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],expected[i],TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbsLog10(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.qualToTransProbsLog10(insQuals,delQuals,gapQuals);
        final double[][] logExpected = new double[expected.length][expected[0].length];
        for (int i = 1; i < expected.length; i++)
            for (int j = 0; j < expected[0].length; j++)
                logExpected[i][j] = Math.log10(expected[i][j]);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length,logExpected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length,logExpected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],logExpected[i],TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbsLog10Fill(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.createTransitionMatrix(insQuals.length);
        PairHMMModel.qualToTransProbsLog10(actual,insQuals,delQuals,gapQuals);
        final double[][] logExpected = new double[expected.length][expected[0].length];
        for (int i = 1; i < expected.length; i++)
            for (int j = 0; j < expected[0].length; j++)
                logExpected[i][j] = Math.log10(expected[i][j]);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length,logExpected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length,logExpected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],logExpected[i],TOLERANCE);
    }

    @Test(dataProvider="qualToTransDataProvider")
    public void testQualsToTransProbsFill(final byte[] insQuals, final byte[] delQuals, final byte[] gapQuals, final double[][] expected) {
        final double[][] actual = PairHMMModel.createTransitionMatrix(insQuals.length);
        PairHMMModel.qualToTransProbs(actual,insQuals,delQuals,gapQuals);
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.length,expected.length);
        Assert.assertNotNull(actual[0]);
        Assert.assertEquals(actual[0].length,expected[0].length);
        for (int i = 0; i < actual.length ; i++)
            assertEqualsDoubleArray(actual[i],expected[i],TOLERANCE);
    }
    @Test(dataProvider="qualToProbsDataProvider")
    public void testQualToProbsLog10Fill(final int insQuals, final int delQual, final int gcp, final double[] expected) {
        final double[] logExpected = new double[expected.length];
        for (int i = 0; i < logExpected.length; i++)
            logExpected[i] = Math.log10(expected[i]);
        final double[] actual = new double[PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
        PairHMMModel.qualToTransProbsLog10(actual, (byte) insQuals, (byte) delQual, (byte) gcp);
        assertEqualsDoubleArray(actual,logExpected,TOLERANCE);
    }


    @DataProvider(name="qualToTransDataProvider")
    public Iterator<Object[]> qualToTransDataProvider() {
        return new Iterator<Object[]>() {

            private final Iterator<Integer> readLengthIterator = readLengthIterator();
            private Iterator<int[]> qualsIterator = qualIterator();

            @Override
            public boolean hasNext() {
                return readLengthIterator.hasNext();
            }

            @Override
            public Object[] next() {
                final int readLength = readLengthIterator.next();
                double[][] matrix = new double[readLength+1][PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
                final byte[] insQuals = new byte[readLength];
                final byte[] delQuals = new byte[readLength];
                final byte[] gapQuals = new byte[readLength];
                for (int i = 0; i < readLength; i++) {
                    if (!qualsIterator.hasNext())
                        qualsIterator = qualIterator();
                    final int[] quals = qualsIterator.next();
                    final int insQual = quals[0];
                    final int delQual = quals[1];
                    final int gapQual = quals[2];
                    final double[] trans = qualsToProbs(insQual, delQual, gapQual);
                    matrix[i+1] = trans;
                    insQuals[i] = (byte)insQual;
                    delQuals[i] = (byte)delQual;
                    gapQuals[i] = (byte)gapQual;
                }

                return new Object[] { insQuals, delQuals, gapQuals, matrix };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }


    @DataProvider(name="qualToProbsDataProvider")
    public Iterator<Object[]> qualToProbsDataProvider() {
        return new Iterator<Object[]>() {
            private final Iterator<int[]> qualsIterator = qualIterator();

            @Override
            public boolean hasNext() {
                return qualsIterator.hasNext();
            }

            @Override
            public Object[] next() {
                final int[] quals = qualsIterator.next();
                final int insQual = quals[0];
                final int delQual = quals[1];
                final int gapQual = quals[2];

                final double[] trans = qualsToProbs(insQual, delQual, gapQual);


                return new Object[] { insQual, delQual, gapQual, trans };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private double[] qualsToProbs(final int insQual, final int delQual, final int gapQual) {
        final double[] trans = new double[PairHMMModel.TRANS_PROB_ARRAY_LENGTH];
        final double matchToMatch = PairHMMModel.matchToMatchProb(insQual, delQual);
        final double matchToInsert = QualityUtils.qualToErrorProb(insQual);
        final double matchToDeletion = QualityUtils.qualToErrorProb(delQual);
        final double indelToMatch = QualityUtils.qualToProb(gapQual);
        final double indelToIndel = QualityUtils.qualToErrorProb(gapQual);

        trans[PairHMMModel.matchToMatch] = matchToMatch;
        trans[PairHMMModel.matchToInsertion] = matchToInsert;
        trans[PairHMMModel.matchToDeletion] = matchToDeletion;
        trans[PairHMMModel.indelToMatch] = indelToMatch;
        trans[PairHMMModel.deletionToDeletion] = trans[PairHMMModel.insertionToInsertion] = indelToIndel;
        return trans;
    }

    private Iterator<Integer> readLengthIterator() {
        return Arrays.asList(READ_LENGTHS).iterator();
    }

    private Iterator<int[]> qualIterator() {
        final int totalCount = INS_QUALS.length * DEL_QUALS.length * GAP_QUALS.length;

        return new Iterator<int[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < totalCount;
            }

            @Override
            public int[] next() {
                final int gap = i % GAP_QUALS.length;
                final int indelGroup = i / GAP_QUALS.length;
                final int del = indelGroup % DEL_QUALS.length;
                final int ins = indelGroup % DEL_QUALS.length;
                i++;
                return new int[] { INS_QUALS[ins], DEL_QUALS[del], GAP_QUALS[gap]};
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }



    @Test(dataProvider = "dualTestDataProvider")
    public void testDoubleQualToProb(final int insQual, final int delQual, final double log10Expected, final double expected) {
        Assert.assertEquals(PairHMMModel.matchToMatchProb(insQual, delQual),expected,TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProbLog10(insQual, delQual),log10Expected,TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProb((byte) insQual, (byte) delQual),expected,TOLERANCE);
        Assert.assertEquals(PairHMMModel.matchToMatchProbLog10((byte) insQual, (byte) delQual),log10Expected,TOLERANCE);
    }

    @DataProvider(name = "dualTestDataProvider")
    private Iterator<Object[]> dualTestDataProvider() {
        final int[] testQuals = new int[] { 0, 1, 2, 5, 10, 13, 17, 20, 23, 27, 30, 43, 57, 70, 100, 200, 254};

        return new Iterator<Object[]>() {
            private int i = 0;
            private int j = 0;

            @Override
            public Object[] next() {

                final int qual1 =  testQuals[i];
                final int qual2 =  testQuals[j];

                final double errorProb1 = Math.pow(10,- 0.1 * qual1);
                final double errorProb2 = Math.pow(10,- 0.1 * qual2);
                final double expected = Math.max(0, (1 - (errorProb1 + errorProb2)));
                final Object[] result = new Object[] { qual1, qual2,Math.log10(Math.min(1,expected)),Math.min(1, expected)};

                if (++j >= testQuals.length) {
                    i++;
                    j = i;
                }
                return result;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }

            @Override
            public boolean hasNext() {
                return i < testQuals.length;
            }
        };
    }


    private static int[] INS_QUALS = {30, 45, 20, 10, 5, 60, 123 };

    private static int[] DEL_QUALS = {30, 45, 20, 10, 5, 60, 123 };

    private static int[] GAP_QUALS = {10, 20, 5};

    private static Integer[] READ_LENGTHS = { 0, 1, 5, 20, 100, 250};
}
