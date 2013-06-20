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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller.readthreading;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.Kmer;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class ReadThreadingGraphUnitTest extends BaseTest {
    private final static boolean DEBUG = false;

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    private void assertNonUniques(final ReadThreadingGraph assembler, String ... nonUniques) {
        final Set<String> actual = new HashSet<>();
        assembler.buildGraphIfNecessary();
        for ( final Kmer kmer : assembler.getNonUniqueKmers() ) actual.add(kmer.baseString());
        final Set<String> expected = new HashSet<>(Arrays.asList(nonUniques));
        Assert.assertEquals(actual, expected);
    }

    @Test(enabled = ! DEBUG)
    public void testNonUniqueMiddle() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "GACACACAGTCA";
        final String read1 = "GACAC---GTCA";
        final String read2 =   "CAC---GTCA";
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
        assertNonUniques(assembler, "ACA", "CAC");
    }

    @Test(enabled = ! DEBUG)
    public void testReadsCreateNonUnique() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "GCAC--GTCA"; // CAC is unique
        final String read1 = "GCACACGTCA"; // makes CAC non unique because it has a duplication
        final String read2 =    "CACGTCA"; // shouldn't be allowed to match CAC as start
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
//        assembler.convertToSequenceGraph().printGraph(new File("test.dot"), 0);

        assertNonUniques(assembler, "CAC");
        //assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }

    @Test(enabled = ! DEBUG)
         public void testCountingOfStartEdges() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "NNNGTCAAA"; // ref has some bases before start
        final String read1 =    "GTCAAA"; // starts at first non N base

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.buildGraphIfNecessary();
//        assembler.printGraph(new File("test.dot"), 0);

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final boolean headerVertex = source.getSuffix() == 'N' || target.getSuffix() == 'N';
            if ( headerVertex ) {
                Assert.assertEquals(edge.getMultiplicity(), 1, "Bases in the unique reference header should have multiplicity of 1");
            } else {
                Assert.assertEquals(edge.getMultiplicity(), 2, "Should have multiplicity of 2 for any edge outside the ref header but got " + edge + " " + source + " -> " + target);
            }
        }
    }

    @Test(enabled = !DEBUG)
    public void testCountingOfStartEdgesWithMultiplePrefixes() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        assembler.increaseCountsThroughBranches = true;
        final String ref   = "NNNGTCAXX"; // ref has some bases before start
        final String alt1  = "NNNCTCAXX"; // alt1 has SNP right after N
        final String read  =     "TCAXX"; // starts right after SNP, but merges right before branch

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(alt1), false);
        assembler.addSequence(getBytes(read), false);
        assembler.buildGraphIfNecessary();
        assembler.printGraph(new File("test.dot"), 0);

        final List<String> oneCountVertices = Arrays.asList("NNN", "NNG", "NNC", "NGT", "NCT");
        final List<String> threeCountVertices = Arrays.asList("CAX", "AXX");

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final int expected = oneCountVertices.contains(target.getSequenceString()) ? 1 : (threeCountVertices.contains(target.getSequenceString()) ? 3 : 2);
            Assert.assertEquals(edge.getMultiplicity(), expected, "Bases at edge " + edge + " from " + source + " to " + target + " has bad multiplicity");
        }
    }

    @Test(enabled = !DEBUG)
    public void testCyclesInGraph() {

        // b37 20:12655200-12655850
        final String ref = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        // SNP at 20:12655528 creates a cycle for small kmers
        final String alt = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        final List<GATKSAMRecord> reads = new ArrayList<>();
        for ( int index = 0; index < alt.length() - 100; index += 20 )
            reads.add(ArtificialSAMUtils.createArtificialRead(Arrays.copyOfRange(alt.getBytes(), index, index + 100), Utils.dupBytes((byte) 30, 100), 100 + "M"));

        // test that there are cycles detected for small kmer
        final ReadThreadingGraph rtgraph25 = new ReadThreadingGraph(25);
        rtgraph25.addSequence("ref", ref.getBytes(), null, true);
        for ( final GATKSAMRecord read : reads )
            rtgraph25.addRead(read);
        rtgraph25.buildGraphIfNecessary();
        Assert.assertTrue(rtgraph25.hasCycles());

        // test that there are no cycles detected for large kmer
        final ReadThreadingGraph rtgraph75 = new ReadThreadingGraph(75);
        rtgraph75.addSequence("ref", ref.getBytes(), null, true);
        for ( final GATKSAMRecord read : reads )
            rtgraph75.addRead(read);
        rtgraph75.buildGraphIfNecessary();
        Assert.assertFalse(rtgraph75.hasCycles());
    }

    @Test(enabled = !DEBUG)
    public void testNsInReadsAreNotUsedForGraph() {

        final int length = 100;
        final byte[] ref = Utils.dupBytes((byte)'A', length);

        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(25);
        rtgraph.addSequence("ref", ref, null, true);

        // add reads with Ns at any position
        for ( int i = 0; i < length; i++ ) {
            final byte[] bases = ref.clone();
            bases[i] = 'N';
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, Utils.dupBytes((byte) 30, length), length + "M");
            rtgraph.addRead(read);
        }
        rtgraph.buildGraphIfNecessary();

        final SeqGraph graph = rtgraph.convertToSequenceGraph();
        final KBestPaths<SeqVertex,BaseEdge> pathFinder = new KBestPaths<>(false);
        Assert.assertEquals(pathFinder.getKBestPaths(graph, length, graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex()).size(), 1);
    }

    @DataProvider(name = "DanglingTails")
    public Object[][] makeDanglingTailsData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // add 1M to the expected CIGAR because it includes the previous (common) base too
        tests.add(new Object[]{"AAAAAAAAAA", "CAAA", "5M", true, 3});                  // incomplete haplotype
        tests.add(new Object[]{"AAAAAAAAAA", "CAAAAAAAAAA", "1M1I10M", true, 10});     // insertion
        tests.add(new Object[]{"CCAAAAAAAAAA", "AAAAAAAAAA", "1M2D10M", true, 10});    // deletion
        tests.add(new Object[]{"AAAAAAAA", "CAAAAAAA", "9M", true, 7});                // 1 snp
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2});                // several snps
        tests.add(new Object[]{"AAAAA", "C", "1M4D1M", true, -1});                     // funky SW alignment
        tests.add(new Object[]{"AAAAA", "CA", "1M3D2M", true, 1});                     // very little data
        tests.add(new Object[]{"AAAAAAA", "CAAAAAC", "8M", true, -1});                 // ends in mismatch
        tests.add(new Object[]{"AAAAAA", "CGAAAACGAA", "1M2I4M2I2M", false, 0});       // alignment is too complex

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingTails", enabled = !DEBUG)
    public void testDanglingTails(final String refEnd,
                                  final String altEnd,
                                  final String cigar,
                                  final boolean cigarIsGood,
                                  final int mergePointDistanceFromSink) {

        final int kmerSize = 15;

        // construct the haplotypes
        final String commonPrefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        final String ref = commonPrefix + refEnd;
        final String alt = commonPrefix + altEnd;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), null, true);
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        rtgraph.addRead(read);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling tail
        MultiDeBruijnVertex altSink = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSink(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSink == null, "We found more than one non-reference sink");
                altSink = v;
            }
        }

        Assert.assertTrue(altSink != null, "We did not find a non-reference sink");

        // confirm that the SW alignment agrees with our expectations
        final ReadThreadingGraph.DanglingTailMergeResult result = rtgraph.generateCigarAgainstReferencePath(altSink);
        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the goodness of the cigar agrees with our expectations
        Assert.assertEquals(rtgraph.cigarIsOkayToMerge(result.cigar), cigarIsGood);

        // confirm that the tail merging works as expected
        if ( cigarIsGood ) {
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 1 || mergePointDistanceFromSink == -1);

            // confirm that we created the appropriate edge
            if ( mergePointDistanceFromSink >= 0 ) {
                MultiDeBruijnVertex v = altSink;
                for ( int i = 0; i < mergePointDistanceFromSink; i++ ) {
                    if ( rtgraph.inDegreeOf(v) != 1 )
                        Assert.fail("Encountered vertex with multiple sources");
                    v = rtgraph.getEdgeSource(rtgraph.incomingEdgeOf(v));
                }
                Assert.assertTrue(rtgraph.outDegreeOf(v) > 1);
            }
        }
    }


// TODO -- update to use determineKmerSizeAndNonUniques directly
//    @DataProvider(name = "KmerSizeData")
//    public Object[][] makeKmerSizeDataProvider() {
//        List<Object[]> tests = new ArrayList<Object[]>();
//
//        // this functionality can be adapted to provide input data for whatever you might want in your data
//        tests.add(new Object[]{3, 3, 3, Arrays.asList("ACG"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 3, Arrays.asList("CAGACG"), Arrays.asList()});
//
//        tests.add(new Object[]{3, 3, 3, Arrays.asList("AAAAC"), Arrays.asList("AAA")});
//        tests.add(new Object[]{3, 4, 4, Arrays.asList("AAAAC"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 4, Arrays.asList("AAAAC"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 3, Arrays.asList("CAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 4, Arrays.asList("CAAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 4, Arrays.asList("CAAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 5, Arrays.asList("ACGAAAAACG"), Arrays.asList()});
//
//        for ( int maxSize = 3; maxSize < 20; maxSize++ ) {
//            for ( int dupSize = 3; dupSize < 20; dupSize++ ) {
//                final int expectedSize = Math.min(maxSize, dupSize);
//                final String dup = Utils.dupString("C", dupSize);
//                final List<String> nonUnique = dupSize > maxSize ? Arrays.asList(Utils.dupString("C", maxSize)) : Collections.<String>emptyList();
//                tests.add(new Object[]{3, maxSize, expectedSize, Arrays.asList("ACGT", "A" + dup + "GT"), nonUnique});
//                tests.add(new Object[]{3, maxSize, expectedSize, Arrays.asList("A" + dup + "GT", "ACGT"), nonUnique});
//            }
//        }
//
//        return tests.toArray(new Object[][]{});
//    }
//
//    /**
//     * Example testng test using MyDataProvider
//     */
//    @Test(dataProvider = "KmerSizeData")
//    public void testDynamicKmerSizing(final int min, final int max, final int expectKmer, final List<String> seqs, final List<String> expectedNonUniques) {
//        final ReadThreadingGraph assembler = new ReadThreadingGraph(min, max);
//        for ( String seq : seqs ) assembler.addSequence(seq.getBytes(), false);
//        assembler.buildGraphIfNecessary();
//        Assert.assertEquals(assembler.getKmerSize(), expectKmer);
//        assertNonUniques(assembler, expectedNonUniques.toArray(new String[]{}));
//    }


}
