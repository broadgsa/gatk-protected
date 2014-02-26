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

package org.broadinstitute.sting.utils.smithwaterman;

import net.sf.samtools.TextCigarCodec;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GlobalEdgeGreedySWPairwiseAlignmentUnitTest extends BaseTest {

    private final static boolean DEBUG = false;

    @Test(enabled = !DEBUG)
    public void testReadAlignedToRefComplexAlignment() {
        final String reference = "AAAGGACTGACTG";
        final String read      = "ACTGACTGACTG";
        final GlobalEdgeGreedySWPairwiseAlignment sw = new GlobalEdgeGreedySWPairwiseAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getCigar().toString(), "1M1D11M");
    }

    @Test(enabled = !DEBUG)
    public void testIndelsAtStartAndEnd() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match + "GGG";
        final int expectedStart = 0;
        final String expectedCigar = "3D5M3I";
        final GlobalEdgeGreedySWPairwiseAlignment sw = new GlobalEdgeGreedySWPairwiseAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = !DEBUG)
    public void testDegenerateAlignmentWithIndelsAtBothEnds() {
        logger.warn("testDegenerateAlignmentWithIndelsAtBothEnds");
        final String ref = "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final String alt =               "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final int expectedStart = 0;
        final String expectedCigar = "6I45M";
        final GlobalEdgeGreedySWPairwiseAlignment sw = new GlobalEdgeGreedySWPairwiseAlignment(ref.getBytes(), alt.getBytes(), SWParameterSet.STANDARD_NGS);
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = !DEBUG)
    public void testAlignReallyLongDeletion() {
        final String ref = "CGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAGGTGATCCACTCGCCTCGGTCTCCCAAAGTGTTGGGATTACAGGCATGAACCACTGCACCTGGCCTAGTGTTTGGGAAAACTATACTAGGAAAAGAATAGTTGCTTTAAGTCATTCTTTGATTATTCTGAGAATTGGCATATAGCTGCCATTATAACCTACTTTTGCTAAATATAATAATAATAATCATTATTTTTATTTTTTGAGACAGGGTCTTGTTTTGTCACCCCGGCTGGAGTGAAGTGGCGCAATCTCGGCTCACTGCAACCTCCACCTCCGGGTGCAAGCAATTCTCCTGCCTCAGCCTCTTGAGTAGCTAGGATTACAGGCACAAGCCATCATGCCCAGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGTCAGGCTGGTCTTGAACTCCTGACCTCAGGT";
        final String alt = "CGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCATGTTGGTCAGGCTGGTCTTGAACTCCTGACCTCAGGT";

        final GlobalEdgeGreedySWPairwiseAlignment sw = new GlobalEdgeGreedySWPairwiseAlignment(ref.getBytes(), alt.getBytes(), SWParameterSet.STANDARD_NGS);
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), 0);
        Assert.assertEquals(sw.getCigar().toString(), "47M419D31M");
    }

    public static final Parameters params = new Parameters(20.0, -10.0, -26.0, -1.1);
    @DataProvider(name = "SWData")
    public Object[][] makeSWData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // simple cases
        tests.add(new Object[]{"A", "C", "1M"});
        tests.add(new Object[]{"AAA", "AAA", "3M"});
        tests.add(new Object[]{"AAA", "AGA", "3M"});
        tests.add(new Object[]{"AAA", "GAA", "3M"});
        tests.add(new Object[]{"AAA", "AAG", "3M"});

        // small single indels
        tests.add(new Object[]{"ACACACAC", "ACACAC", "6M2D"});
        tests.add(new Object[]{"ACACAC", "ACACACAC", "6M2I"});
        tests.add(new Object[]{"XXACACACXX", "XXACACACACXX", "8M2I2M"});
        tests.add(new Object[]{"XXACACACXX", "XXACACXX", "6M2D2M"});
        tests.add(new Object[]{"ACGT", "AACGT", "1I4M"});
        tests.add(new Object[]{"ACGT", "ACCGT", "2M1I2M"});
        tests.add(new Object[]{"ACGT", "ACGGT", "3M1I1M"});
        tests.add(new Object[]{"ACGT", "ACGTT", "4M1I"});
        tests.add(new Object[]{"ACGT", "CGT", "1D3M"});
        tests.add(new Object[]{"ACGT", "AGT", "1M1D2M"});
        tests.add(new Object[]{"ACGT", "ACT", "2M1D1M"});
        tests.add(new Object[]{"ACGT", "ACG", "3M1D"});

        // mismatches through out the sequences
        final String ref = "ACGTAACCGGTT";
        for ( int diff = 0; diff < ref.length(); diff++ ) {
            final byte[] altBases = ref.getBytes();
            altBases[diff] = 'N';
            tests.add(new Object[]{ref, new String(altBases), ref.length() + "M"});
        }
        for ( int diff1 = 0; diff1 < ref.length(); diff1++ ) {
            for ( int diff2 = 0; diff2 < ref.length(); diff2++ ) {
                final byte[] altBases = ref.getBytes();
                altBases[diff1] = 'N';
                altBases[diff2] = 'N';
                tests.add(new Object[]{ref, new String(altBases), ref.length() + "M"});
            }
        }

        // prefixes and suffixes matching
        final String totalPrefix = "ACG";
        final String totalSuffix = "GCT";
        for ( int prefixSize = 0; prefixSize < totalPrefix.length(); prefixSize++) {
            for ( int suffixSize = 0; suffixSize < totalPrefix.length(); suffixSize++) {
                if ( prefixSize + suffixSize == 0 )
                    continue;
                for ( int indelSize = 1; indelSize < 50; indelSize++ ) {
                    final String prefix = totalPrefix.substring(0, prefixSize);
                    final String suffix = totalSuffix.substring(0, suffixSize);
                    final String insert = Utils.dupString("N", indelSize);
                    tests.add(new Object[]{prefix + suffix, prefix + insert + suffix, prefix.length() + "M" + indelSize + "I" + suffix.length() + "M"});
                    tests.add(new Object[]{prefix + insert + suffix, prefix + suffix, prefix.length() + "M" + indelSize + "D" + suffix.length() + "M"});
                }
            }
        }

        // larger indels with prefixes/suffixes
        tests.add(new Object[]{"ACTGTTTTGAACATCAGTTATTTTAAACTTTTAAGTTGTTAGCACAGCAAAAGCAACAAAATTCTAAGTGCAGTAATCACTTTACTGCGTGGTCATATGAAATCAAGGCAATGTTATGAGTATTACTGGAAAGCTGGACAGAGTAACGGGAAAAGTGACTAAAACTATGC", "CCTGTTTTGAACATCAGTTATTTTAAACTTTTAAGTTGTTAGCACAGCAAAAGCAACAAAATTCTAAGTGCAGTAATCACTTTACTGCGTGGTCATATGAAATCAAGGCAATGTTATGAGTATTACTGGAAAGCTGGACAGAGTAACGGGAAAAGTGACT", "160M10D"});
        tests.add(new Object[]{"LLLLLTATTAAGTAGTGCTCTATGTTGTCAACTAATTTATTTCCCATTTCAAACATTAGTTGACATGTTTTCATTTCTCTTTTGGAAGGAAACAACTAAATATGTTATCAATCCATCATTTACTTGTACAATAAATAAAGTTCTAAATCACTGCACAGTGTAAAATGGCAAATAGACTTCCCCATAACACAAAGCCATCCTGAAAAGTTTTGTTCATTTTAGAAGRRRRR", "LLLLLARRRRR", "5M219D6M"});
        tests.add(new Object[]{"LLLLLTATTTTTTRRRRR", "LLLLLARRRRR", "5M7D6M"});

        // systematic testing
        for ( final int forwardMatches : Arrays.asList(0, 1, 5, 10)) {
            for ( final int forwardMismatches : Arrays.asList(0, 1, 2)) {
                for ( final int middleMatches : Arrays.asList(0, 1, 5, 10)) {
                    for ( final int delSize : Arrays.asList(0, 1, 2, 3 )) {
                        for ( final int insSize : Arrays.asList(0, 1, 2, 3 )) {
                            for ( final int reverseMismatches : Arrays.asList(0, 1, 2)) {
                                for ( final int reverseMatches : Arrays.asList(0, 1, 5, 10)) {
                                    // if there is an insertion and deletion, they should cancel each other out (at least partially)
                                    final int overlap = Math.min(delSize, insSize);
                                    final int myDelSize = delSize - overlap;
                                    final int myInsSize = insSize - overlap;

                                    // this case is too difficult to create a CIGAR for because SW will (legitimately) prefer to switch the indel and mismatches
                                    final int totalMismatches = forwardMismatches + reverseMismatches;
                                    if ( (myDelSize > 0 || myInsSize > 0 ) && (totalMismatches >= myDelSize || totalMismatches >= myInsSize) )
                                        continue;

                                    final StringBuilder refBuilder = new StringBuilder();
                                    final StringBuilder altBuilder = new StringBuilder();
                                    final StringBuilder cigarBuilder = new StringBuilder();

                                    refBuilder.append(Utils.dupString('A', forwardMatches + forwardMismatches + middleMatches));
                                    altBuilder.append(Utils.dupString('A', forwardMatches));
                                    altBuilder.append(Utils.dupString('C', forwardMismatches));
                                    altBuilder.append(Utils.dupString('A', middleMatches));
                                    cigarBuilder.append(forwardMatches + forwardMismatches + middleMatches);
                                    cigarBuilder.append("M");

                                    if ( myDelSize > 0 ) {
                                        refBuilder.append(Utils.dupString('G', myDelSize));
                                        cigarBuilder.append(myDelSize);
                                        cigarBuilder.append("D");
                                    }
                                    if ( myInsSize > 0 ) {
                                        altBuilder.append(Utils.dupString('T', myInsSize));
                                        cigarBuilder.append(myInsSize);
                                        cigarBuilder.append("I");
                                    }
                                    if ( overlap > 0 ) {
                                        refBuilder.append(Utils.dupString('G', overlap));
                                        altBuilder.append(Utils.dupString('T', overlap));
                                        cigarBuilder.append(overlap);
                                        cigarBuilder.append("M");
                                    }
                                    if ( delSize > 0 || insSize > 0 ) {
                                        refBuilder.append(Utils.dupString('A', middleMatches));
                                        altBuilder.append(Utils.dupString('A', middleMatches));
                                        cigarBuilder.append(middleMatches);
                                        cigarBuilder.append("M");
                                    }

                                    refBuilder.append(Utils.dupString('A', reverseMismatches + reverseMatches));
                                    altBuilder.append(Utils.dupString('C', reverseMismatches));
                                    altBuilder.append(Utils.dupString('A', reverseMatches));
                                    cigarBuilder.append(reverseMismatches + reverseMatches);
                                    cigarBuilder.append("M");

                                    if ( refBuilder.length() > 0 && altBuilder.length() > 0 )
                                        tests.add(new Object[]{refBuilder.toString(), altBuilder.toString(), cigarBuilder.toString()});
                                }
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SWData", enabled = !DEBUG)
    public void testSW(final String seq1, final String seq2, final String expectedCigar) {
        final GlobalEdgeGreedySWPairwiseAlignment alignment = new GlobalEdgeGreedySWPairwiseAlignment(seq1.getBytes(), seq2.getBytes(), new Parameters(5.0, -5.0, -25.0, -1.0));
        Assert.assertEquals(alignment.getCigar(), AlignmentUtils.consolidateCigar(TextCigarCodec.getSingleton().decode(expectedCigar)));
    }

    /**
     * For debugging purposes only
     */
    @Test(enabled = DEBUG)
    public void testDebugging() {
        final String ref = "A";
        final String alt = "C";

        final GlobalEdgeGreedySWPairwiseAlignment sw = new GlobalEdgeGreedySWPairwiseAlignment(ref.getBytes(), alt.getBytes(), new Parameters(5.0, -5.0, -25.0, -1.0));
        Assert.assertEquals(sw.getCigar().toString(), "1M");
    }
}
