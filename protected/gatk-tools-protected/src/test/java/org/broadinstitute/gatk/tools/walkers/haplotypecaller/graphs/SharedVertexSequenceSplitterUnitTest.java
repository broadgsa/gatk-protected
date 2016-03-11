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

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class SharedVertexSequenceSplitterUnitTest extends BaseTest {
    private final static boolean PRINT_GRAPHS = false;

    @DataProvider(name = "PrefixSuffixData")
    public Object[][] makePrefixSuffixData() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{Arrays.asList("A", "C"), 0, 0});
        tests.add(new Object[]{Arrays.asList("C", "C"), 1, 0});
        tests.add(new Object[]{Arrays.asList("ACT", "AGT"), 1, 1});
        tests.add(new Object[]{Arrays.asList("ACCT", "AGT"), 1, 1});
        tests.add(new Object[]{Arrays.asList("ACT", "ACT"), 3, 0});
        tests.add(new Object[]{Arrays.asList("ACTA", "ACT"), 3, 0});
        tests.add(new Object[]{Arrays.asList("ACTA", "ACTG"), 3, 0});
        tests.add(new Object[]{Arrays.asList("ACTA", "ACTGA"), 3, 1});
        tests.add(new Object[]{Arrays.asList("GCTGA", "ACTGA"), 0, 4});

        tests.add(new Object[]{Arrays.asList("A", "C", "A"), 0, 0});
        tests.add(new Object[]{Arrays.asList("A", "A", "A"), 1, 0});
        tests.add(new Object[]{Arrays.asList("A", "AA", "A"), 1, 0});
        tests.add(new Object[]{Arrays.asList("A", "ACA", "A"), 1, 0});
        tests.add(new Object[]{Arrays.asList("ACT", "ACAT", "ACT"), 2, 1});
        tests.add(new Object[]{Arrays.asList("ACT", "ACAT", "ACGT"), 2, 1});
        tests.add(new Object[]{Arrays.asList("AAAT", "AAA", "CAAA"), 0, 0});
        tests.add(new Object[]{Arrays.asList("AACTTT", "AAGTTT", "AAGCTTT"), 2, 3});
        tests.add(new Object[]{Arrays.asList("AAA", "AAA", "CAAA"), 0, 3});
        tests.add(new Object[]{Arrays.asList("AAA", "AAA", "AAA"), 3, 0});

        tests.add(new Object[]{Arrays.asList("AC", "ACA", "AC"), 2, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PrefixSuffixData")
    public void testPrefixSuffix(final List<String> strings, int expectedPrefixLen, int expectedSuffixLen) {
        final List<byte[]> bytes = new ArrayList<>();
        int min = Integer.MAX_VALUE;
        for ( final String s : strings ) {
            bytes.add(s.getBytes());
            min = Math.min(min, s.length());
        }

        final int actualPrefixLen = GraphUtils.compPrefixLen(bytes, min);
        Assert.assertEquals(actualPrefixLen, expectedPrefixLen, "Failed prefix test");

        final int actualSuffixLen = GraphUtils.compSuffixLen(bytes, min - actualPrefixLen);
        Assert.assertEquals(actualSuffixLen, expectedSuffixLen, "Failed suffix test");
    }

    @Test(dataProvider = "PrefixSuffixData")
    public void testPrefixSuffixVertices(final List<String> strings, int expectedPrefixLen, int expectedSuffixLen) {
        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : strings ) {
            v.add(new SeqVertex(s));
        }

        final String expectedPrefix = strings.get(0).substring(0, expectedPrefixLen);
        final String expectedSuffix = strings.get(0).substring(strings.get(0).length() - expectedSuffixLen);

        final Pair<SeqVertex, SeqVertex> result = SharedVertexSequenceSplitter.commonPrefixAndSuffixOfVertices(v);
        Assert.assertEquals(result.getFirst().getSequenceString(), expectedPrefix, "Failed suffix test");
        Assert.assertEquals(result.getSecond().getSequenceString(), expectedSuffix, "Failed suffix test");

        Assert.assertEquals(result.getFirst().isEmpty(), expectedPrefix.isEmpty());
        Assert.assertEquals(result.getSecond().isEmpty(), expectedSuffix.isEmpty());
    }

    @Test(dataProvider = "PrefixSuffixData")
    public void testSplitter(final List<String> strings, int expectedPrefixLen, int expectedSuffixLen) {
        final SeqGraph graph = new SeqGraph(11);

        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : strings ) {
            v.add(new SeqVertex(s));
        }

        graph.addVertices(v.toArray(new SeqVertex[v.size()]));

        final String expectedPrefix = strings.get(0).substring(0, expectedPrefixLen);
        final String expectedSuffix = strings.get(0).substring(strings.get(0).length() - expectedSuffixLen);

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(graph, v);
        splitter.split();

        Assert.assertEquals(splitter.prefixV.getSequenceString(), expectedPrefix);
        Assert.assertEquals(splitter.suffixV.getSequenceString(), expectedSuffix);

        Assert.assertTrue(splitter.splitGraph.outDegreeOf(splitter.prefixV) <= strings.size());
        Assert.assertEquals(splitter.splitGraph.inDegreeOf(splitter.prefixV), 0);

        Assert.assertTrue(splitter.splitGraph.inDegreeOf(splitter.suffixV) <= strings.size());
        Assert.assertEquals(splitter.splitGraph.outDegreeOf(splitter.suffixV), 0);

        for ( final SeqVertex mid : splitter.newMiddles ) {
            Assert.assertNotNull(splitter.splitGraph.getEdge(splitter.prefixV, mid));
            Assert.assertNotNull(splitter.splitGraph.getEdge(mid, splitter.suffixV));
        }
    }

    @DataProvider(name = "CompleteCycleData")
    public Object[][] makeCompleteCycleData() {
        List<Object[]> tests = new ArrayList<>();

        for ( final boolean hasTop : Arrays.asList(true, false) ) {
            for ( final boolean hasBot : Arrays.asList(true, false) ) {
                if ( ! hasTop && ! hasBot ) continue;
                tests.add(new Object[]{Arrays.asList("A", "A"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "AC"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "CA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("AC", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("AT", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("ATA", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("ATAA", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("ATAACA", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CCCAAA", "AAA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CCCAAAAAA", "AAA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CCCAAAAAA", "CCCAAA"), hasTop, hasBot});

                tests.add(new Object[]{Arrays.asList("A", "A", "A"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "A", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "C", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("AC", "C", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CA", "C", "C"), hasTop, hasBot});
                // all merged
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "AGA"), hasTop, hasBot});
                // prefix and suffix
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "ACA"), hasTop, hasBot});
                // 2 -> prefix, leave C
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "AGAC"), hasTop, hasBot});
                // 2 -> prefix, leave CCC
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "AGACCC"), hasTop, hasBot});
                // 2 -> suffix, leave A/T
                tests.add(new Object[]{Arrays.asList("TAGA", "TAGA", "AAGA"), hasTop, hasBot});
                // 2 -> suffix, leave T, delete 1
                tests.add(new Object[]{Arrays.asList("TAGA", "TAGA", "AGA"), hasTop, hasBot});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CompleteCycleData")
    public void testSplitterCompleteCycle(final List<String> strings, final boolean hasTop, final boolean hasBot) {
        final SeqGraph graph = new SeqGraph(11);

        int edgeWeight = 1;
        final SeqVertex top = hasTop ? new SeqVertex("AAAAAAAA") : null;
        final SeqVertex bot = hasBot ? new SeqVertex("GGGGGGGG") : null;
        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : strings ) {
            v.add(new SeqVertex(s));
        }
        graph.addVertices(v.toArray(new SeqVertex[v.size()]));
        final SeqVertex first = v.get(0);

        if ( hasTop ) {
            graph.addVertex(top);
            for ( final SeqVertex vi : v )
                graph.addEdge(top, vi, new BaseEdge(vi == first, edgeWeight++));
        }

        if ( hasBot ) {
            graph.addVertex(bot);
            for ( final SeqVertex vi : v )
                graph.addEdge(vi, bot, new BaseEdge(vi == first, edgeWeight++));
        }

        final Set<String> haplotypes = new HashSet<>();
        final KBestHaplotypeFinder originalPaths = new KBestHaplotypeFinder((SeqGraph) graph.clone(),graph.getSources(),graph.getSinks());
        for ( final KBestHaplotype path : originalPaths )
            haplotypes.add(new String(path.bases()));

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(graph, v);
        splitter.split();
        if ( PRINT_GRAPHS ) graph.printGraph(new File(Utils.join("_", strings) + "_" + hasTop + "_" + hasBot + ".original.dot"), 0);
        if ( PRINT_GRAPHS ) splitter.splitGraph.printGraph(new File(Utils.join("_", strings) + "_" + hasTop + "_" + hasBot + ".split.dot"), 0);
        splitter.updateGraph(top, bot);
        if ( PRINT_GRAPHS ) graph.printGraph(new File(Utils.join("_", strings) + "_" + hasTop + "_" + hasBot + ".updated.dot"), 0);

        final KBestHaplotypeFinder splitPaths = new KBestHaplotypeFinder(graph,graph.getSources(),graph.getSinks());
        for ( final KBestHaplotype path : splitPaths ) {
            final String h = new String(path.bases());
            Assert.assertTrue(haplotypes.contains(h), "Failed to find haplotype " + h);
        }


        final List<byte[]> sortedOriginalPaths = new ArrayList<>(originalPaths.size());
        for (final KBestHaplotype kbh : originalPaths.unique())
            sortedOriginalPaths.add(kbh.bases());
        Collections.sort(sortedOriginalPaths, BaseUtils.BASES_COMPARATOR);
        final List<byte[]> sortedSplitPaths = new ArrayList<>(splitPaths.size());
        for (final KBestHaplotype kbh : splitPaths.unique())
            sortedSplitPaths.add(kbh.bases());
        Collections.sort(sortedSplitPaths, BaseUtils.BASES_COMPARATOR);

        Assert.assertEquals(sortedSplitPaths,sortedOriginalPaths,Utils.join("_", strings) + "_" + hasTop + "_" + hasBot);
    }

    @DataProvider(name = "MeetsMinSequenceData")
    public Object[][] makeMeetsMinSequenceData() {
        final List<Object[]> tests = new ArrayList<>();

        final boolean prefixBiased = SharedVertexSequenceSplitter.prefersPrefixMerging();
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 0, true, true});
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 1, prefixBiased, ! prefixBiased});
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 2, prefixBiased, ! prefixBiased});
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 3, false, false});
        tests.add(new Object[]{Arrays.asList("A", "AC"), 1, true, false});
        tests.add(new Object[]{Arrays.asList("A", "AC"), 2, false, false});
        tests.add(new Object[]{Arrays.asList("AT", "AC"), 1, true, false});
        tests.add(new Object[]{Arrays.asList("AAT", "AAC"), 1, true, false});
        tests.add(new Object[]{Arrays.asList("AAT", "AAC"), 2, true, false});
        tests.add(new Object[]{Arrays.asList("AAT", "AAC"), 3, false, false});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 1, true, true});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 2, true, true});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 3, false, true});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 4, false, false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MeetsMinSequenceData")
    public void testSplitterCompleteCycle(final List<String> mids, final int minSeqLength, final boolean prefixMeets, final boolean suffixMeets) {
        final SeqGraph graph = new SeqGraph(11);

        final SeqVertex top = new SeqVertex("AAAAAAAA");
        final SeqVertex bot = new SeqVertex("GGGGGGGG");
        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : mids ) { v.add(new SeqVertex(s)); }
        graph.addVertices(v.toArray(new SeqVertex[v.size()]));
        graph.addVertices(top, bot);
        for ( final SeqVertex vi : v ) { graph.addEdge(top, vi); graph.addEdge(vi, bot); }

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(graph, v);
        Assert.assertEquals(splitter.meetsMinMergableSequenceForPrefix(minSeqLength), prefixMeets, "Prefix failed");
        Assert.assertEquals(splitter.meetsMinMergableSequenceForSuffix(minSeqLength), suffixMeets, "Suffix failed");
        Assert.assertEquals(splitter.meetsMinMergableSequenceForEitherPrefixOrSuffix(minSeqLength), suffixMeets || prefixMeets, "Either prefix or suffix failed");
    }
}
