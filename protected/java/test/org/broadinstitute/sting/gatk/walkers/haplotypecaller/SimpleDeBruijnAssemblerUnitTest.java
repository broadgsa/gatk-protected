package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/27/12
 */

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class SimpleDeBruijnAssemblerUnitTest extends BaseTest {


    private class MergeNodesWithNoVariationTestProvider extends TestDataProvider {
        public byte[] sequence;
        public int KMER_LENGTH;

        public MergeNodesWithNoVariationTestProvider(String seq, int kmer) {
            super(MergeNodesWithNoVariationTestProvider.class, String.format("Merge nodes with no variation test. kmer = %d, seq = %s", kmer, seq));
            sequence = seq.getBytes();
            KMER_LENGTH = kmer;
        }

        public DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> expectedGraph() {
            DeBruijnVertex v = new DeBruijnVertex(sequence, 0);
            DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            graph.addVertex(v);
            return graph;
        }

        public DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> calcGraph() {

            DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
            for (int i = 0; i < kmersInSequence - 1; i++) {
                // get the kmers
                final byte[] kmer1 = new byte[KMER_LENGTH];
                System.arraycopy(sequence, i, kmer1, 0, KMER_LENGTH);
                final byte[] kmer2 = new byte[KMER_LENGTH];
                System.arraycopy(sequence, i+1, kmer2, 0, KMER_LENGTH);

                SimpleDeBruijnAssembler.addKmersToGraph(graph, kmer1, kmer2, false);
            }
            SimpleDeBruijnAssembler.mergeNodes(graph);
            return graph;
        }
    }

    @DataProvider(name = "MergeNodesWithNoVariationTestProvider")
    public Object[][] makeMergeNodesWithNoVariationTests() {
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 3);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 4);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 5);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 6);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 7);
        new MergeNodesWithNoVariationTestProvider("GGTTAACCATGCAGACGGGAGGCTGAGCGAGAGTTTT", 6);
        new MergeNodesWithNoVariationTestProvider("AATACCATTGGAGTTTTTTTCCAGGTTAAGATGGTGCATTGAATCCACCCATCTACTTTTGCTCCTCCCAAAACTCACTAAAACTATTATAAAGGGATTTTGTTTAAAGACACAAACTCATGAGGACAGAGAGAACAGAGTAGACAATAGTGGGGGAAAAATAAGTTGGAAGATAGAAAACAGATGGGTGAGTGGTAATCGACTCAGCAGCCCCAAGAAAGCTGAAACCCAGGGAAAGTTAAGAGTAGCCCTATTTTCATGGCAAAATCCAAGGGGGGGTGGGGAAAGAAAGAAAAACAGAAAAAAAAATGGGAATTGGCAGTCCTAGATATCTCTGGTACTGGGCAAGCCAAAGAATCAGGATAACTGGGTGAAAGGTGATTGGGAAGCAGTTAAAATCTTAGTTCCCCTCTTCCACTCTCCGAGCAGCAGGTTTCTCTCTCTCATCAGGCAGAGGGCTGGAGAT", 66);
        new MergeNodesWithNoVariationTestProvider("AATACCATTGGAGTTTTTTTCCAGGTTAAGATGGTGCATTGAATCCACCCATCTACTTTTGCTCCTCCCAAAACTCACTAAAACTATTATAAAGGGATTTTGTTTAAAGACACAAACTCATGAGGACAGAGAGAACAGAGTAGACAATAGTGGGGGAAAAATAAGTTGGAAGATAGAAAACAGATGGGTGAGTGGTAATCGACTCAGCAGCCCCAAGAAAGCTGAAACCCAGGGAAAGTTAAGAGTAGCCCTATTTTCATGGCAAAATCCAAGGGGGGGTGGGGAAAGAAAGAAAAACAGAAAAAAAAATGGGAATTGGCAGTCCTAGATATCTCTGGTACTGGGCAAGCCAAAGAATCAGGATAACTGGGTGAAAGGTGATTGGGAAGCAGTTAAAATCTTAGTTCCCCTCTTCCACTCTCCGAGCAGCAGGTTTCTCTCTCTCATCAGGCAGAGGGCTGGAGAT", 76);

        return MergeNodesWithNoVariationTestProvider.getTests(MergeNodesWithNoVariationTestProvider.class);
    }

    @Test(dataProvider = "MergeNodesWithNoVariationTestProvider", enabled = true)
    public void testMergeNodesWithNoVariation(MergeNodesWithNoVariationTestProvider cfg) {
        logger.warn(String.format("Test: %s", cfg.toString()));
        Assert.assertTrue(graphEquals(cfg.calcGraph(), cfg.expectedGraph()));
    }

    @Test(enabled = true)
    public void testPruneGraph() {
        DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
        DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> expectedGraph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        DeBruijnVertex v = new DeBruijnVertex("ATGG".getBytes(), 0);
        DeBruijnVertex v2 = new DeBruijnVertex("ATGGA".getBytes(), 0);
        DeBruijnVertex v3 = new DeBruijnVertex("ATGGT".getBytes(), 0);
        DeBruijnVertex v4 = new DeBruijnVertex("ATGGG".getBytes(), 0);
        DeBruijnVertex v5 = new DeBruijnVertex("ATGGC".getBytes(), 0);
        DeBruijnVertex v6 = new DeBruijnVertex("ATGGCCCCCC".getBytes(), 0);

        graph.addVertex(v);
        graph.addVertex(v2);
        graph.addVertex(v3);
        graph.addVertex(v4);
        graph.addVertex(v5);
        graph.addVertex(v6);
        graph.addEdge(v, v2, new DeBruijnEdge(false, 1));
        graph.addEdge(v2, v3, new DeBruijnEdge(false, 3));
        graph.addEdge(v3, v4, new DeBruijnEdge(false, 5));
        graph.addEdge(v4, v5, new DeBruijnEdge(false, 3));
        graph.addEdge(v5, v6, new DeBruijnEdge(false, 2));

        expectedGraph.addVertex(v2);
        expectedGraph.addVertex(v3);
        expectedGraph.addVertex(v4);
        expectedGraph.addVertex(v5);
        expectedGraph.addEdge(v2, v3, new DeBruijnEdge(false, 3));
        expectedGraph.addEdge(v3, v4, new DeBruijnEdge(false, 5));
        expectedGraph.addEdge(v4, v5, new DeBruijnEdge(false, 3));

        SimpleDeBruijnAssembler.pruneGraph(graph, 2);

        Assert.assertTrue(graphEquals(graph, expectedGraph));

        graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
        expectedGraph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        graph.addVertex(v);
        graph.addVertex(v2);
        graph.addVertex(v3);
        graph.addVertex(v4);
        graph.addVertex(v5);
        graph.addVertex(v6);
        graph.addEdge(v, v2, new DeBruijnEdge(true, 1));
        graph.addEdge(v2, v3, new DeBruijnEdge(false, 3));
        graph.addEdge(v3, v4, new DeBruijnEdge(false, 5));
        graph.addEdge(v4, v5, new DeBruijnEdge(false, 3));

        expectedGraph.addVertex(v);
        expectedGraph.addVertex(v2);
        expectedGraph.addVertex(v3);
        expectedGraph.addVertex(v4);
        expectedGraph.addVertex(v5);
        expectedGraph.addEdge(v, v2, new DeBruijnEdge(true, 1));
        expectedGraph.addEdge(v2, v3, new DeBruijnEdge(false, 3));
        expectedGraph.addEdge(v3, v4, new DeBruijnEdge(false, 5));
        expectedGraph.addEdge(v4, v5, new DeBruijnEdge(false, 3));

        SimpleDeBruijnAssembler.pruneGraph(graph, 2);

        Assert.assertTrue(graphEquals(graph, expectedGraph));
    }

    @Test(enabled = true)
    public void testEliminateNonRefPaths() {
        DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
        DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> expectedGraph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        DeBruijnVertex v = new DeBruijnVertex("ATGG".getBytes(), 0);
        DeBruijnVertex v2 = new DeBruijnVertex("ATGGA".getBytes(), 0);
        DeBruijnVertex v3 = new DeBruijnVertex("ATGGT".getBytes(), 0);
        DeBruijnVertex v4 = new DeBruijnVertex("ATGGG".getBytes(), 0);
        DeBruijnVertex v5 = new DeBruijnVertex("ATGGC".getBytes(), 0);
        DeBruijnVertex v6 = new DeBruijnVertex("ATGGCCCCCC".getBytes(), 0);
        
        graph.addVertex(v);
        graph.addVertex(v2);
        graph.addVertex(v3);
        graph.addVertex(v4);
        graph.addVertex(v5);
        graph.addVertex(v6);
        graph.addEdge(v, v2, new DeBruijnEdge(false));
        graph.addEdge(v2, v3, new DeBruijnEdge(true));
        graph.addEdge(v3, v4, new DeBruijnEdge(true));
        graph.addEdge(v4, v5, new DeBruijnEdge(true));
        graph.addEdge(v5, v6, new DeBruijnEdge(false));

        expectedGraph.addVertex(v2);
        expectedGraph.addVertex(v3);
        expectedGraph.addVertex(v4);
        expectedGraph.addVertex(v5);
        expectedGraph.addEdge(v2, v3, new DeBruijnEdge());
        expectedGraph.addEdge(v3, v4, new DeBruijnEdge());
        expectedGraph.addEdge(v4, v5, new DeBruijnEdge());

        SimpleDeBruijnAssembler.eliminateNonRefPaths(graph);

        Assert.assertTrue(graphEquals(graph, expectedGraph));

        
        
        
        graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
        expectedGraph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        graph.addVertex(v);
        graph.addVertex(v2);
        graph.addVertex(v3);
        graph.addVertex(v4);
        graph.addVertex(v5);
        graph.addVertex(v6);
        graph.addEdge(v, v2, new DeBruijnEdge(true));
        graph.addEdge(v2, v3, new DeBruijnEdge(true));
        graph.addEdge(v4, v5, new DeBruijnEdge(false));
        graph.addEdge(v5, v6, new DeBruijnEdge(false));

        expectedGraph.addVertex(v);
        expectedGraph.addVertex(v2);
        expectedGraph.addVertex(v3);
        expectedGraph.addEdge(v, v2, new DeBruijnEdge());
        expectedGraph.addEdge(v2, v3, new DeBruijnEdge());

        SimpleDeBruijnAssembler.eliminateNonRefPaths(graph);

        Assert.assertTrue(graphEquals(graph, expectedGraph));



        graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
        expectedGraph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        graph.addVertex(v);
        graph.addVertex(v2);
        graph.addVertex(v3);
        graph.addVertex(v4);
        graph.addVertex(v5);
        graph.addVertex(v6);
        graph.addEdge(v, v2, new DeBruijnEdge(true));
        graph.addEdge(v2, v3, new DeBruijnEdge(true));
        graph.addEdge(v4, v5, new DeBruijnEdge(false));
        graph.addEdge(v5, v6, new DeBruijnEdge(false));
        graph.addEdge(v4, v2, new DeBruijnEdge(false));

        expectedGraph.addVertex(v);
        expectedGraph.addVertex(v2);
        expectedGraph.addVertex(v3);
        expectedGraph.addEdge(v, v2, new DeBruijnEdge());
        expectedGraph.addEdge(v2, v3, new DeBruijnEdge());

        SimpleDeBruijnAssembler.eliminateNonRefPaths(graph);

        Assert.assertTrue(graphEquals(graph, expectedGraph));
    }

    private boolean graphEquals(DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> g1, DefaultDirectedGraph<DeBruijnVertex,DeBruijnEdge> g2) {
        if( !(g1.vertexSet().containsAll(g2.vertexSet()) && g2.vertexSet().containsAll(g1.vertexSet())) ) {
            return false;
        }
        for( DeBruijnEdge e1 : g1.edgeSet() ) {
            boolean found = false;
            for( DeBruijnEdge e2 : g2.edgeSet() ) {
                if( e1.equals(g1, e2, g2) ) { found = true; break; }
            }
            if( !found ) { return false; }
        }
        for( DeBruijnEdge e2 : g2.edgeSet() ) {
            boolean found = false;
            for( DeBruijnEdge e1 : g1.edgeSet() ) {
                if( e2.equals(g2, e1, g1) ) { found = true; break; }
            }
            if( !found ) { return false; }
        }
        return true;
    }
}
