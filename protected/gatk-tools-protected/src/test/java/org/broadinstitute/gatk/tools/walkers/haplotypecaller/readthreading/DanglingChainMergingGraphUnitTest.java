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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.KBestHaplotype;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.KBestHaplotypeFinder;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class DanglingChainMergingGraphUnitTest extends BaseTest {

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    @DataProvider(name = "DanglingTails")
    public Object[][] makeDanglingTailsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the previous (common) base too
        tests.add(new Object[]{"AAAAAAAAAA", "CAAA", "5M", true, 3});                  // incomplete haplotype
        tests.add(new Object[]{"AAAAAAAAAA", "CAAAAAAAAAA", "1M1I10M", true, 10});     // insertion
        tests.add(new Object[]{"CCAAAAAAAAAA", "AAAAAAAAAA", "1M2D10M", true, 10});    // deletion
        tests.add(new Object[]{"AAAAAAAA", "CAAAAAAA", "9M", true, 7});                // 1 snp
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2});                // several snps
        tests.add(new Object[]{"AAAAA", "C", "1M4D1M", false, -1});                    // funky SW alignment
        tests.add(new Object[]{"AAAAA", "CA", "1M3D2M", false, 1});                    // very little data
        tests.add(new Object[]{"AAAAAAA", "CAAAAAC", "8M", true, -1});                 // ends in mismatch
        tests.add(new Object[]{"AAAAAA", "CGAAAACGAA", "1M2I4M2I2M", false, 0});       // alignment is too complex
        tests.add(new Object[]{"AAAAA", "XXXXX", "1M5I", false, -1});                  // insertion

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingTails")
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
        rtgraph.addSequence("ref", ref.getBytes(), true);
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
        final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink, 0, 4);

        if ( result == null ) {
            Assert.assertFalse(cigarIsGood);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the goodness of the cigar agrees with our expectations
        Assert.assertEquals(rtgraph.cigarIsOkayToMerge(result.cigar, false, true), cigarIsGood);

        // confirm that the tail merging works as expected
        if ( cigarIsGood ) {
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 1 || mergePointDistanceFromSink == -1);

            // confirm that we created the appropriate edge
            if ( mergePointDistanceFromSink >= 0 ) {
                MultiDeBruijnVertex v = altSink;
                for ( int i = 0; i < mergePointDistanceFromSink; i++ ) {
                    if ( rtgraph.inDegreeOf(v) != 1 )
                        Assert.fail("Encountered vertex with multiple edges");
                    v = rtgraph.getEdgeSource(rtgraph.incomingEdgeOf(v));
                }
                Assert.assertTrue(rtgraph.outDegreeOf(v) > 1);
            }
        }
    }

    @Test
    public void testWholeTailIsInsertion() {
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(10);
        final ReadThreadingGraph.DanglingChainMergeHelper result = new ReadThreadingGraph.DanglingChainMergeHelper(null, null, "AXXXXX".getBytes(), "AAAAAA".getBytes(), new TextCigarCodec().decode("5I1M"));
        final int mergeResult = rtgraph.mergeDanglingTail(result);
        Assert.assertEquals(mergeResult, 0);
    }

    @Test
    public void testGetBasesForPath() {

        final int kmerSize = 4;
        final String testString = "AATGGGGCAATACTA";

        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
        graph.addSequence(testString.getBytes(), true);
        graph.buildGraphIfNecessary();

        final List<MultiDeBruijnVertex> vertexes = new ArrayList<>();
        MultiDeBruijnVertex v = graph.getReferenceSourceVertex();
        while ( v != null ) {
            vertexes.add(v);
            v = graph.getNextReferenceVertex(v);
        }

        final String resultForTails = new String(graph.getBasesForPath(vertexes, false));
        Assert.assertEquals(resultForTails, testString.substring(kmerSize-1));
        final String resultForHeads = new String(graph.getBasesForPath(vertexes, true));
        Assert.assertEquals(resultForHeads, "GTAAGGGCAATACTA");  // because the source node will be reversed
    }

    @DataProvider(name = "DanglingHeads")
    public Object[][] makeDanglingHeadsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the last (common) base too
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true});        // 1 snp
        tests.add(new Object[]{"XXXAACCGGTTACGT", "XAAACCGGTTACGT", "7M", false});         // 1 snp
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "XAACGGTTACGT", "4M1D4M", false});   // deletion
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYYCGGTTACGT", "8M", true});        // 2 snps
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", true});    // 2 snps
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYCGGTTACGT", "7M", true});         // very little data
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "YCCGGTTACGT", "6M", true});         // begins in mismatch

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingHeads")
    public void testDanglingHeads(final String ref,
                                  final String alt,
                                  final String cigar,
                                  final boolean shouldBeMerged) {

        final int kmerSize = 5;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        rtgraph.addRead(read);
        rtgraph.setMaxMismatchesInDanglingHead(10);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling head
        MultiDeBruijnVertex altSource = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSource(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSource == null, "We found more than one non-reference source");
                altSource = v;
            }
        }

        Assert.assertTrue(altSource != null, "We did not find a non-reference source");

        // confirm that the SW alignment agrees with our expectations
        final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstUpwardsReferencePath(altSource, 0, 1);

        if ( result == null ) {
            Assert.assertFalse(shouldBeMerged);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the tail merging works as expected
        final int mergeResult = rtgraph.mergeDanglingHead(result);
        Assert.assertTrue(mergeResult > 0 || !shouldBeMerged);

        // confirm that we created the appropriate bubble in the graph only if expected
        rtgraph.cleanNonRefPaths();
        final SeqGraph seqGraph = rtgraph.convertToSequenceGraph();
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(seqGraph, seqGraph.getReferenceSourceVertex(), seqGraph.getReferenceSinkVertex());
        Assert.assertEquals(paths.size(), shouldBeMerged ? 2 : 1);
    }
}
