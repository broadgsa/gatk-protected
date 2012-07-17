package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */
// Class for finding the K best paths (as determined by the sum of multiplicities of the edges) in a graph.
// This is different from most graph traversals because we want to test paths from any source node to any sink node.
public class KBestPaths {

    // static access only
    protected KBestPaths() { }
    private static int MAX_PATHS_TO_HOLD = 100;

    protected static class MyInt { public int val = 0; }

    // class to keep track of paths
    protected static class Path {

        // the last vertex seen in the path
        private DeBruijnVertex lastVertex;

        // the list of edges comprising the path
        private ArrayList<DeBruijnEdge> edges;

        // the scores for the path
        private int totalScore = 0, lowestEdge = -1;

        public Path( final DeBruijnVertex initialVertex ) {
            lastVertex = initialVertex;
            edges = new ArrayList<DeBruijnEdge>(0);
        }

        public Path( final Path p, final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final DeBruijnEdge edge ) {
            lastVertex = graph.getEdgeTarget(edge);
            edges = new ArrayList<DeBruijnEdge>(p.edges);
            edges.add(edge);
            totalScore = p.totalScore + edge.getMultiplicity();
            lowestEdge = ( p.lowestEdge == -1 ) ? edge.getMultiplicity() : Math.min(p.lowestEdge, edge.getMultiplicity());
        }

        public boolean containsEdge( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final DeBruijnEdge edge ) {
            final DeBruijnVertex targetVertex = graph.getEdgeTarget(edge);
            for( final DeBruijnEdge e : edges ) {
                if( e.equals(graph, edge) || graph.getEdgeTarget(e).equals(targetVertex) ) {
                    return true;
                }
            }
            
            return false;
        }

        public ArrayList<DeBruijnEdge> getEdges() { return edges; }

        public int getScore() { return totalScore; }

        public int getLowestEdge() { return lowestEdge; }

        public DeBruijnVertex getLastVertexInPath() { return lastVertex; }

        public byte[] getBases( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
            if( edges.size() == 0 ) { return lastVertex.getSequence(); }
            
            byte[] bases = graph.getEdgeSource( edges.get(0) ).getSequence();
            for( final DeBruijnEdge e : edges ) {
                bases = ArrayUtils.addAll(bases, graph.getEdgeTarget( e ).getSuffix());
            }
            return bases;
        }
    }

    protected static class PathComparatorTotalScore implements Comparator<Path> {
        public int compare(final Path path1, final Path path2) {
            return path1.totalScore - path2.totalScore;
        }
    }

    protected static class PathComparatorLowestEdge implements Comparator<Path> {
        public int compare(final Path path1, final Path path2) {
            return path2.lowestEdge - path1.lowestEdge;
        }
    }

    public static List<Path> getKBestPaths( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final int k ) {
        if( k > MAX_PATHS_TO_HOLD/2 ) { throw new ReviewedStingException("Asked for more paths than MAX_PATHS_TO_HOLD!"); }
        final ArrayList<Path> bestPaths = new ArrayList<Path>();
        
        // run a DFS for best paths
        for( final DeBruijnVertex v : graph.vertexSet() ) {
            if( graph.inDegreeOf(v) == 0 ) {
                findBestPaths(graph, new Path(v), bestPaths);
            }
        }

        Collections.sort(bestPaths, new PathComparatorLowestEdge() );
        Collections.reverse(bestPaths);
        return bestPaths.subList(0, Math.min(k, bestPaths.size()));
    }

    private static void findBestPaths( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final Path path, final List<Path> bestPaths ) {
        findBestPaths(graph, path, bestPaths, new MyInt());
    }

    private static void findBestPaths( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final Path path, final List<Path> bestPaths, MyInt n ) {

        // did we hit the end of a path?
        if ( allOutgoingEdgesHaveBeenVisited(graph, path) ) {
            if ( bestPaths.size() >= MAX_PATHS_TO_HOLD ) {
                // clean out some low scoring paths
                Collections.sort(bestPaths, new PathComparatorLowestEdge() );
                for(int iii = 0; iii < 20; iii++) { bestPaths.remove(0); }
            }
            bestPaths.add(path);
        } else if( n.val > 10000) {
            // do nothing, just return
        } else {
            // recursively run DFS
            final ArrayList<DeBruijnEdge> edgeArrayList = new ArrayList<DeBruijnEdge>();
            edgeArrayList.addAll(graph.outgoingEdgesOf(path.lastVertex));
            Collections.sort(edgeArrayList);
            Collections.reverse(edgeArrayList);
            for ( final DeBruijnEdge edge : edgeArrayList ) {
                // make sure the edge is not already in the path
                if ( path.containsEdge(graph, edge) )
                    continue;

                final Path newPath = new Path(path, graph, edge);
                n.val++;
                findBestPaths(graph, newPath, bestPaths, n);
            }
        }
    }

    private static boolean allOutgoingEdgesHaveBeenVisited( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final Path path ) {
        for( final DeBruijnEdge edge : graph.outgoingEdgesOf(path.lastVertex) ) {
            if( !path.containsEdge(graph, edge) ) {
                return false;
            }
        }
        return true;
    }
}
