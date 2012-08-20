package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.jgrapht.graph.DefaultDirectedGraph;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */

// simple edge class for connecting nodes in the graph
public class DeBruijnEdge {

    private int multiplicity;
    private boolean isRef;

    public DeBruijnEdge() {
        multiplicity = 1;
        isRef = false;
    }

    public DeBruijnEdge( final boolean isRef ) {
        multiplicity = 1;
        this.isRef = isRef;
    }

    public DeBruijnEdge( final boolean isRef, final int multiplicity ) {
        this.multiplicity = multiplicity;
        this.isRef = isRef;
    }

    public int getMultiplicity() {
        return multiplicity;
    }

    public void setMultiplicity( final int value ) {
        multiplicity = value;
    }

    public boolean getIsRef() {
        return isRef;
    }

    public void setIsRef( final boolean isRef ) {
        this.isRef = isRef;
    }

    public boolean equals( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final DeBruijnEdge edge ) {
        return (graph.getEdgeSource(this).equals(graph.getEdgeSource(edge))) && (graph.getEdgeTarget(this).equals(graph.getEdgeTarget(edge)));
    }

    public boolean equals( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final DeBruijnEdge edge, final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph2 ) {
        return (graph.getEdgeSource(this).equals(graph2.getEdgeSource(edge))) && (graph.getEdgeTarget(this).equals(graph2.getEdgeTarget(edge)));
    }

    public static class EdgeWeightComparator implements Comparator<DeBruijnEdge> {
        @Override
        public int compare(final DeBruijnEdge edge1, final DeBruijnEdge edge2) {
            return edge1.multiplicity - edge2.multiplicity;
        }
    }
}
