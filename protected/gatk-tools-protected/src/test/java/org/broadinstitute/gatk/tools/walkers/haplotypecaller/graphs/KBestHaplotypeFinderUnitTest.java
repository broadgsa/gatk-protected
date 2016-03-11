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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 1/31/13
 */
public class KBestHaplotypeFinderUnitTest extends BaseTest {

    @DataProvider(name = "BasicPathFindingData")
    public Object[][] makeBasicPathFindingData() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int nStartNodes : Arrays.asList(1, 2, 3) ) {
            for ( final int nBranchesPerBubble : Arrays.asList(2, 3) ) {
                for ( final int nEndNodes : Arrays.asList(1, 2, 3) ) {
                    tests.add(new Object[]{nStartNodes, nBranchesPerBubble, nEndNodes});
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    private static int weight = 1;
    final Set<SeqVertex> createVertices(final SeqGraph graph, final int n, final SeqVertex source, final SeqVertex target) {
        final List<String> seqs = Arrays.asList("A", "C", "G", "T");
        final Set<SeqVertex> vertices = new LinkedHashSet<>();
        for ( int i = 0; i < n; i++ ) {
            final SeqVertex v = new SeqVertex(seqs.get(i));
            graph.addVertex(v);
            vertices.add(v);
            if ( source != null ) graph.addEdge(source, v, new BaseEdge(false, weight++));
            if ( target != null ) graph.addEdge(v, target, new BaseEdge(false, weight++));
        }
        return vertices;
    }

    @Test(dataProvider = "BasicPathFindingData")
    public void testBasicPathFinding(final int nStartNodes, final int nBranchesPerBubble, final int nEndNodes) {
        final SeqGraph graph = new SeqGraph(11);

        final SeqVertex middleTop = new SeqVertex("GTAC");
        final SeqVertex middleBottom = new SeqVertex("ACTG");
        graph.addVertices(middleTop, middleBottom);
        final Set<SeqVertex> starts = createVertices(graph, nStartNodes, null, middleTop);
        @SuppressWarnings("unused")
        final Set<SeqVertex> bubbles = createVertices(graph, nBranchesPerBubble, middleTop, middleBottom);
        final Set<SeqVertex> ends = createVertices(graph, nEndNodes, middleBottom, null);

        // enumerate all possible paths
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph, starts, ends);

        final int expectedNumOfPaths = nStartNodes * nBranchesPerBubble * nEndNodes;
        Assert.assertEquals(paths.size(), expectedNumOfPaths, "Didn't find the expected number of paths");

        double lastScore = 0;
        for ( final KBestHaplotype kbh : paths ) {
            final Path<SeqVertex,BaseEdge> path = kbh.path();
            Assert.assertTrue(kbh.score() <= lastScore, "Paths out of order.   Path " + path + " has score " + path.getScore() + " above previous " + lastScore);
            lastScore = kbh.score();
        }

        // get the best path, and make sure it's the same as our optimal path overall
        final Path<SeqVertex,BaseEdge> best = paths.get(0).path();
        final List<KBestHaplotype> justOne = new KBestHaplotypeFinder(graph,starts, ends).subList(0,1);
        Assert.assertEquals(justOne.size(), 1);

        Assert.assertTrue(justOne.get(0).path().pathsAreTheSame(best), "Best path from complete enumerate " + best + " not the same as from k = 1 search " + justOne.get(0));
    }

    @DataProvider(name = "BasicBubbleDataProvider")
    public Object[][] makeBasicBubbleDataProvider() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int refBubbleLength : Arrays.asList(1, 5, 10) ) {
            for ( final int altBubbleLength : Arrays.asList(1, 5, 10) ) {
                tests.add(new Object[]{refBubbleLength, altBubbleLength});
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BasicBubbleDataProvider")
    public void testBasicBubbleData(final int refBubbleLength, final int altBubbleLength) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(3);
        final String preRef = "ATGG";
        final String postRef = "GGGGC";

        SeqVertex v = new SeqVertex(preRef);
        SeqVertex v2Ref = new SeqVertex(Utils.dupString('A', refBubbleLength));
        SeqVertex v2Alt = new SeqVertex(Utils.dupString('A', altBubbleLength-1) + "T");
        SeqVertex v3 = new SeqVertex(postRef);

        graph.addVertex(v);
        graph.addVertex(v2Ref);
        graph.addVertex(v2Alt);
        graph.addVertex(v3);
        graph.addEdge(v, v2Ref, new BaseEdge(true, 10));
        graph.addEdge(v2Ref, v3, new BaseEdge(true, 10));
        graph.addEdge(v, v2Alt, new BaseEdge(false, 5));
        graph.addEdge(v2Alt, v3, new BaseEdge(false, 5));

        // Construct the test path
        Path<SeqVertex,BaseEdge> path = new Path<>(v, graph);
        path = new Path<>(path, graph.getEdge(v, v2Alt));
        path = new Path<>(path, graph.getEdge(v2Alt, v3));

        // Construct the actual cigar string implied by the test path
        Cigar expectedCigar = new Cigar();
        expectedCigar.add(new CigarElement(preRef.length(), CigarOperator.M));
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength, CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength,CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(postRef.length(), CigarOperator.M));

        final String ref = preRef + v2Ref.getSequenceString() + postRef;
        Assert.assertEquals(path.calculateCigar(ref.getBytes()).toString(), AlignmentUtils.consolidateCigar(expectedCigar).toString(), "Cigar string mismatch");
    }

    @DataProvider(name = "GetBasesData")
    public Object[][] makeGetBasesData() {
        List<Object[]> tests = new ArrayList<>();

        final List<String> frags = Arrays.asList("ACT", "GAC", "CAT");

        for ( int n = 1; n <= frags.size(); n++ ) {
            for ( final List<String> comb : Utils.makePermutations(frags, n, false) ) {
                tests.add(new Object[]{comb});
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetBasesData")
    public void testGetBases(final List<String> frags) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(3);

        SeqVertex prev = null;
        for (final String s : frags) {
            SeqVertex v = new SeqVertex(s);
            graph.addVertex(v);
            if ( prev != null )
                graph.addEdge(prev, v);
            prev = v;
        }

        // enumerate all possible paths
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph,graph.getSources(),graph.getSinks());
        Assert.assertEquals(paths.size(), 1);
        final Path<SeqVertex,BaseEdge> path = paths.get(0).path();
        Assert.assertEquals(new String(path.getBases()), Utils.join("", frags), "Path doesn't have the expected sequence");
    }

    @DataProvider(name = "TripleBubbleDataProvider")
    public Object[][] makeTripleBubbleDataProvider() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int refBubbleLength : Arrays.asList(1, 5, 10) ) {
            for ( final int altBubbleLength : Arrays.asList(1, 5, 10) ) {
                for ( final boolean offRefEnding : Arrays.asList(true, false) ) {
                    for ( final boolean offRefBeginning : Arrays.asList(false) ) {
                        tests.add(new Object[]{refBubbleLength, altBubbleLength, offRefBeginning, offRefEnding});
                    }
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TripleBubbleDataProvider")
    public void testTripleBubbleData(final int refBubbleLength, final int altBubbleLength, final boolean offRefBeginning, final boolean offRefEnding) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(11);
        final String preAltOption = "ATCGATCGATCGATCGATCG";
        final String postAltOption = "CCCC";
        final String preRef = "ATGG";
        final String postRef = "GGCCG";
        final String midRef1 = "TTCCT";
        final String midRef2 = "CCCAAAAAAAAAAAA";

        SeqVertex preV = new SeqVertex(preAltOption);
        SeqVertex v = new SeqVertex(preRef);
        SeqVertex v2Ref = new SeqVertex(Utils.dupString('A', refBubbleLength));
        SeqVertex v2Alt = new SeqVertex(Utils.dupString('A', altBubbleLength-1) + "T");
        SeqVertex v4Ref = new SeqVertex(Utils.dupString('C', refBubbleLength));
        SeqVertex v4Alt = new SeqVertex(Utils.dupString('C', altBubbleLength-1) + "T");
        SeqVertex v6Ref = new SeqVertex(Utils.dupString('G', refBubbleLength));
        SeqVertex v6Alt = new SeqVertex(Utils.dupString('G', altBubbleLength-1) + "T");
        SeqVertex v3 = new SeqVertex(midRef1);
        SeqVertex v5 = new SeqVertex(midRef2);
        SeqVertex v7 = new SeqVertex(postRef);
        SeqVertex postV = new SeqVertex(postAltOption);

        final String ref = preRef + v2Ref.getSequenceString() + midRef1 + v4Ref.getSequenceString() + midRef2 + v6Ref.getSequenceString() + postRef;

        graph.addVertex(preV);
        graph.addVertex(v);
        graph.addVertex(v2Ref);
        graph.addVertex(v2Alt);
        graph.addVertex(v3);
        graph.addVertex(v4Ref);
        graph.addVertex(v4Alt);
        graph.addVertex(v5);
        graph.addVertex(v6Ref);
        graph.addVertex(v6Alt);
        graph.addVertex(v7);
        graph.addVertex(postV);
        graph.addEdge(preV, v, new BaseEdge(false, 1));
        graph.addEdge(v, v2Ref, new BaseEdge(true, 10));
        graph.addEdge(v2Ref, v3, new BaseEdge(true, 10));
        graph.addEdge(v, v2Alt, new BaseEdge(false, 5));
        graph.addEdge(v2Alt, v3, new BaseEdge(false, 5));
        graph.addEdge(v3, v4Ref, new BaseEdge(true, 10));
        graph.addEdge(v4Ref, v5, new BaseEdge(true, 10));
        graph.addEdge(v3, v4Alt, new BaseEdge(false, 5));
        graph.addEdge(v4Alt, v5, new BaseEdge(false, 5));
        graph.addEdge(v5, v6Ref, new BaseEdge(true, 11));
        graph.addEdge(v6Ref, v7, new BaseEdge(true, 11));
        graph.addEdge(v5, v6Alt, new BaseEdge(false, 55));
        graph.addEdge(v6Alt, v7, new BaseEdge(false, 55));
        graph.addEdge(v7, postV, new BaseEdge(false, 1));

        // Construct the test path
        Path<SeqVertex,BaseEdge> path = new Path<>( (offRefBeginning ? preV : v), graph);
        if( offRefBeginning )
            path = new Path<>(path, graph.getEdge(preV, v));
        path = new Path<>(path, graph.getEdge(v, v2Alt));
        path = new Path<>(path, graph.getEdge(v2Alt, v3));
        path = new Path<>(path, graph.getEdge(v3, v4Ref));
        path = new Path<>(path, graph.getEdge(v4Ref, v5));
        path = new Path<>(path, graph.getEdge(v5, v6Alt));
        path = new Path<>(path, graph.getEdge(v6Alt, v7));
        if( offRefEnding )
            path = new Path<>(path, graph.getEdge(v7,postV));

        // Construct the actual cigar string implied by the test path
        Cigar expectedCigar = new Cigar();
        if( offRefBeginning ) {
            expectedCigar.add(new CigarElement(preAltOption.length(), CigarOperator.I));
        }
        expectedCigar.add(new CigarElement(preRef.length(), CigarOperator.M));
        // first bubble
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength,CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength,CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength,CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(midRef1.length(), CigarOperator.M));
        // second bubble is ref path
        expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        expectedCigar.add(new CigarElement(midRef2.length(), CigarOperator.M));
        // third bubble
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength,CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength,CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength,CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(postRef.length(), CigarOperator.M));
        if( offRefEnding ) {
            expectedCigar.add(new CigarElement(postAltOption.length(), CigarOperator.I));
        }

        Assert.assertEquals(path.calculateCigar(ref.getBytes()).toString(),
                AlignmentUtils.consolidateCigar(expectedCigar).toString(),
                "Cigar string mismatch: ref = " + ref + " alt " + new String(path.getBases()));
    }

    @Test
    public void testIntraNodeInsertionDeletion() {
        // Construct the assembly graph
        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex top = new SeqVertex("T");
        final SeqVertex bot = new SeqVertex("T");
        final SeqVertex alt = new SeqVertex("AAACCCCC");
        final SeqVertex ref = new SeqVertex("CCCCCGGG");

        graph.addVertices(top, bot, alt, ref);
        graph.addEdges(new BaseEdge(true, 1), top, ref, bot);
        graph.addEdges(new BaseEdge(false, 1), top, alt, bot);

        @SuppressWarnings("all")
        final KBestHaplotypeFinder bestPathFinder = new KBestHaplotypeFinder(graph,top,bot);

        Assert.assertEquals(bestPathFinder.size(), 2);

        final Path<SeqVertex,BaseEdge> refPath = bestPathFinder.get(0).path();
        final Path<SeqVertex,BaseEdge> altPath = bestPathFinder.get(1).path();

        final String refString = top.getSequenceString() + ref.getSequenceString() + bot.getSequenceString();
        Assert.assertEquals(refPath.calculateCigar(refString.getBytes()).toString(), "10M");
        Assert.assertEquals(altPath.calculateCigar(refString.getBytes()).toString(), "1M3I5M3D1M");
    }

    @Test
    public void testHardSWPath() {
        // Construct the assembly graph
        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex top = new SeqVertex( "NNN" );
        final SeqVertex bot = new SeqVertex( "NNN" );
        final SeqVertex alt = new SeqVertex( "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA" );
        final SeqVertex ref = new SeqVertex( "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA" );
        graph.addVertices(top, bot, alt, ref);
        graph.addEdges(new BaseEdge(true, 1), top, ref, bot);
        graph.addEdges(new BaseEdge(false, 1), top, alt, bot);

        @SuppressWarnings("all")
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph, top, bot);

        Assert.assertEquals(paths.size(), 2);

        final Path<SeqVertex,BaseEdge> refPath = paths.get(0).path();
        final Path<SeqVertex,BaseEdge> altPath = paths.get(1).path();

        final String refString = top.getSequenceString() + ref.getSequenceString() + bot.getSequenceString();

        logger.warn("RefPath : " + refPath + " cigar " + refPath.calculateCigar(refString.getBytes()));
        logger.warn("AltPath : " + altPath + " cigar " + altPath.calculateCigar(refString.getBytes()));

        Assert.assertEquals(refPath.calculateCigar(refString.getBytes()).toString(), "51M");
        Assert.assertEquals(altPath.calculateCigar(refString.getBytes()).toString(), "3M6I48M");
    }

    // -----------------------------------------------------------------
    //
    // Systematic tests to ensure that we get the correct SW result for
    // a variety of variants in the ref vs alt bubble
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "SystematicRefAltSWTestData")
    public Object[][] makeSystematicRefAltSWTestData() {
        final List<Object[]> tests = new ArrayList<>();

        final List<List<String>> allDiffs = Arrays.asList(
                Arrays.asList("G", "C", "1M"),
                Arrays.asList("G", "", "1D"),
                Arrays.asList("", "C", "1I"),
                Arrays.asList("AAA", "CGT", "3M"),
                Arrays.asList("TAT", "CAC", "3M"),
                Arrays.asList("GCTG", "GTCG", "4M"),
                Arrays.asList("AAAAA", "", "5D"),
                Arrays.asList("", "AAAAA", "5I"),
                Arrays.asList("AAAAACC", "CCGGGGGG", "5D2M6I")
        );

        for ( final String prefix : Arrays.asList("", "X", "XXXXXXXXXXXXX")) {
            for ( final String end : Arrays.asList("", "X", "XXXXXXXXXXXXX")) {
                for ( final List<String> diffs : allDiffs )
                    tests.add(new Object[]{prefix, end, diffs.get(0), diffs.get(1), diffs.get(2)});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SystematicRefAltSWTestData")
    public void testRefAltSW(final String prefix, final String end, final String refMid, final String altMid, final String midCigar) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(11);

        final int padSize = 0;
        SeqVertex top = new SeqVertex(Utils.dupString("N", padSize));
        SeqVertex ref = new SeqVertex(prefix + refMid + end);
        SeqVertex alt = new SeqVertex(prefix + altMid + end);
        SeqVertex bot = new SeqVertex(Utils.dupString("N", padSize));

        graph.addVertices(top, ref, alt, bot);
        graph.addEdges(new BaseEdge(true, 1), top, ref, bot);
        graph.addEdges(new BaseEdge(false, 1), top, alt, bot);

        // Construct the test path
        Path<SeqVertex,BaseEdge> path = Path.makePath(Arrays.asList(top, alt, bot), graph);

        Cigar expected = new Cigar();
        expected.add(new CigarElement(padSize, CigarOperator.M));
        if ( ! prefix.equals("") ) expected.add(new CigarElement(prefix.length(), CigarOperator.M));
        for ( final CigarElement elt : TextCigarCodec.decode(midCigar).getCigarElements() ) expected.add(elt);
        if ( ! end.equals("") ) expected.add(new CigarElement(end.length(), CigarOperator.M));
        expected.add(new CigarElement(padSize, CigarOperator.M));
        expected = AlignmentUtils.consolidateCigar(expected);

        final String refString = top.getSequenceString() + ref.getSequenceString() + bot.getSequenceString();
        final Cigar pathCigar = path.calculateCigar(refString.getBytes());

        logger.warn("diffs: " + ref + " vs. " + alt + " cigar " + midCigar);
        logger.warn("Path " + path + " with cigar " + pathCigar);
        logger.warn("Expected cigar " + expected);

        Assert.assertEquals(pathCigar, expected, "Cigar mismatch: ref = " + refString + " vs alt = " + new String(path.getBases()));
    }

    @Test
    public void testLeftAlignCigarSequentially() {
        String preRefString = "GATCGATCGATC";
        String postRefString = "TTT";
        String refString = "ATCGAGGAGAGCGCCCCG";
        String indelString1 = "X";
        String indelString2 = "YZ";
        int refIndel1 = 10;
        int refIndel2 = 12;

        for ( final int indelSize1 : Arrays.asList(1, 2, 3, 4) ) {
            for ( final int indelOp1 : Arrays.asList(1, -1) ) {
                for ( final int indelSize2 : Arrays.asList(1, 2, 3, 4) ) {
                    for ( final int indelOp2 : Arrays.asList(1, -1) ) {

                        Cigar expectedCigar = new Cigar();
                        expectedCigar.add(new CigarElement(refString.length(), CigarOperator.M));
                        expectedCigar.add(new CigarElement(indelSize1, (indelOp1 > 0 ? CigarOperator.I : CigarOperator.D)));
                        expectedCigar.add(new CigarElement((indelOp1 < 0 ? refIndel1 - indelSize1 : refIndel1), CigarOperator.M));
                        expectedCigar.add(new CigarElement(refString.length(), CigarOperator.M));
                        expectedCigar.add(new CigarElement(indelSize2 * 2, (indelOp2 > 0 ? CigarOperator.I : CigarOperator.D)));
                        expectedCigar.add(new CigarElement((indelOp2 < 0 ? (refIndel2 - indelSize2) * 2 : refIndel2 * 2), CigarOperator.M));
                        expectedCigar.add(new CigarElement(refString.length(), CigarOperator.M));

                        Cigar givenCigar = new Cigar();
                        givenCigar.add(new CigarElement(refString.length() + refIndel1/2, CigarOperator.M));
                        givenCigar.add(new CigarElement(indelSize1, (indelOp1 > 0 ? CigarOperator.I : CigarOperator.D)));
                        givenCigar.add(new CigarElement((indelOp1 < 0 ? (refIndel1/2 - indelSize1) : refIndel1/2) + refString.length() + refIndel2/2 * 2, CigarOperator.M));
                        givenCigar.add(new CigarElement(indelSize2 * 2, (indelOp2 > 0 ? CigarOperator.I : CigarOperator.D)));
                        givenCigar.add(new CigarElement((indelOp2 < 0 ? (refIndel2/2 - indelSize2) * 2 : refIndel2/2 * 2) + refString.length(), CigarOperator.M));

                        String theRef = preRefString + refString + Utils.dupString(indelString1, refIndel1) + refString + Utils.dupString(indelString2, refIndel2) + refString + postRefString;
                        String theRead = refString + Utils.dupString(indelString1, refIndel1 + indelOp1 * indelSize1) + refString + Utils.dupString(indelString2, refIndel2 + indelOp2 * indelSize2) + refString;

                        Cigar calculatedCigar = CigarUtils.leftAlignCigarSequentially(AlignmentUtils.consolidateCigar(givenCigar), theRef.getBytes(), theRead.getBytes(), preRefString.length(), 0);
                        Assert.assertEquals(AlignmentUtils.consolidateCigar(calculatedCigar).toString(), AlignmentUtils.consolidateCigar(expectedCigar).toString(), "Cigar strings do not match!");
                    }
                }
            }
        }
    }

    @Test(enabled = true)
    public void testLeftAlignCigarSequentiallyAdjacentID() {
        final String ref = "GTCTCTCTCTCTCTCTCTATATATATATATATATTT";
        final String hap = "GTCTCTCTCTCTCTCTCTCTCTATATATATATATTT";
        final Cigar originalCigar = TextCigarCodec.decode("18M4I12M4D2M");

        final Cigar result = CigarUtils.leftAlignCigarSequentially(originalCigar, ref.getBytes(), hap.getBytes(), 0, 0);
        logger.warn("Result is " + result);
        Assert.assertEquals(originalCigar.getReferenceLength(), result.getReferenceLength(), "Reference lengths are different");
    }
}
