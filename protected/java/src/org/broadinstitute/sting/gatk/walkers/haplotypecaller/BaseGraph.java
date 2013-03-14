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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.jgrapht.EdgeFactory;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.DepthFirstIterator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 2/6/13
 */

public class BaseGraph<T extends BaseVertex> extends DefaultDirectedGraph<T, BaseEdge> {
    protected final static Logger logger = Logger.getLogger(BaseGraph.class);
    private final int kmerSize;

    /**
     * Construct an empty BaseGraph
     */
    public BaseGraph() {
        this(11);
    }

    /**
     * Edge factory that creates non-reference multiplicity 1 edges
     * @param <T> the new of our vertices
     */
    private static class MyEdgeFactory<T extends BaseVertex> implements EdgeFactory<T, BaseEdge> {
        @Override
        public BaseEdge createEdge(T sourceVertex, T targetVertex) {
            return new BaseEdge(false, 1);
        }
    }

    /**
     * Construct a DeBruijnGraph with kmerSize
     * @param kmerSize
     */
    public BaseGraph(final int kmerSize) {
        super(new MyEdgeFactory<T>());

        if ( kmerSize < 1 ) throw new IllegalArgumentException("kmerSize must be >= 1 but got " + kmerSize);
        this.kmerSize = kmerSize;
    }

    /**
     * How big of a kmer did we use to create this graph?
     * @return
     */
    public int getKmerSize() {
        return kmerSize;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    public boolean isReferenceNode( final T v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        for( final BaseEdge e : edgesOf(v) ) {
            if( e.isRef() ) { return true; }
        }
        return false;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node (in degree == 0)
     */
    public boolean isSource( final T v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        return inDegreeOf(v) == 0;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a sink node (out degree == 0)
     */
    public boolean isSink( final T v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        return outDegreeOf(v) == 0;
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    @Ensures({"result != null"})
    public byte[] getAdditionalSequence( final T v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to pull sequence from a null vertex."); }
        return v.getAdditionalSequence(isSource(v));
    }

    /**
     * @param e the edge to test
     * @return  true if this edge is a reference source edge
     */
    public boolean isRefSource( final BaseEdge e ) {
        if( e == null ) { throw new IllegalArgumentException("Attempting to test a null edge."); }
        for( final BaseEdge edgeToTest : incomingEdgesOf(getEdgeSource(e)) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    public boolean isRefSource( final T v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        for( final BaseEdge edgeToTest : incomingEdgesOf(v) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @param e the edge to test
     * @return  true if this edge is a reference sink edge
     */
    public boolean isRefSink( final BaseEdge e ) {
        if( e == null ) { throw new IllegalArgumentException("Attempting to test a null edge."); }
        for( final BaseEdge edgeToTest : outgoingEdgesOf(getEdgeTarget(e)) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference sink
     */
    public boolean isRefSink( final T v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        for( final BaseEdge edgeToTest : outgoingEdgesOf(v) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public T getReferenceSourceVertex( ) {
        for( final T v : vertexSet() ) {
            if( isReferenceNode(v) && isRefSource(v) ) {
                return v;
            }
        }
        return null;
    }

    /**
     * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public T getReferenceSinkVertex( ) {
        for( final T v : vertexSet() ) {
            if( isReferenceNode(v) && isRefSink(v) ) {
                return v;
            }
        }
        return null;
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the next reference vertex if it exists
     */
    public T getNextReferenceVertex( final T v ) {
        if( v == null ) { return null; }
        for( final BaseEdge edgeToTest : outgoingEdgesOf(v) ) {
            if( edgeToTest.isRef() ) {
                return getEdgeTarget(edgeToTest);
            }
        }
        return null;
    }

    /**
     * Traverse the graph and get the previous reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the previous reference vertex if it exists
     */
    public T getPrevReferenceVertex( final T v ) {
        if( v == null ) { return null; }
        for( final BaseEdge edgeToTest : incomingEdgesOf(v) ) {
            if( isReferenceNode(getEdgeSource(edgeToTest)) ) {
                return getEdgeSource(edgeToTest);
            }
        }
        return null;
    }

    /**
     * Does a reference path exist between the two vertices?
     * @param fromVertex    from this vertex, can be null
     * @param toVertex      to this vertex, can be null
     * @return              true if a reference path exists in the graph between the two vertices
     */
    public boolean referencePathExists(final T fromVertex, final T toVertex) {
        T v = fromVertex;
        if( v == null ) {
            return false;
        }
        v = getNextReferenceVertex(v);
        if( v == null ) {
            return false;
        }
        while( !v.equals(toVertex) ) {
            v = getNextReferenceVertex(v);
            if( v == null ) {
                return false;
            }
        }
        return true;
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    public byte[] getReferenceBytes( final T fromVertex, final T toVertex, final boolean includeStart, final boolean includeStop ) {
        if( fromVertex == null ) { throw new IllegalArgumentException("Starting vertex in requested path cannot be null."); }
        if( toVertex == null ) { throw  new IllegalArgumentException("From vertex in requested path cannot be null."); }

        byte[] bytes = null;
        T v = fromVertex;
        if( includeStart ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        v = getNextReferenceVertex(v); // advance along the reference path
        while( v != null && !v.equals(toVertex) ) {
            bytes = ArrayUtils.addAll( bytes, getAdditionalSequence(v) );
            v = getNextReferenceVertex(v); // advance along the reference path
        }
        if( includeStop && v != null && v.equals(toVertex)) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        return bytes;
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    public void addVertices(final T ... vertices) {
        for ( final T v : vertices )
            addVertex(v);
    }

    /**
     * Get the set of vertices connected by outgoing edges of V
     * @param v a non-null vertex
     * @return a set of vertices connected by outgoing edges from v
     */
    public Set<T> outgoingVerticesOf(final T v) {
        final Set<T> s = new HashSet<T>();
        for ( final BaseEdge e : outgoingEdgesOf(v) ) {
            s.add(getEdgeTarget(e));
        }
        return s;
    }

    /**
     * Get the set of vertices connected to v by incoming edges
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    public Set<T> incomingVerticesOf(final T v) {
        final Set<T> s = new HashSet<T>();
        for ( final BaseEdge e : incomingEdgesOf(v) ) {
            s.add(getEdgeSource(e));
        }
        return s;
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param destination File to write to
     */
    public void printGraph(final File destination, final int pruneFactor) {
        PrintStream stream = null;

        try {
            stream = new PrintStream(new FileOutputStream(destination));
            printGraph(stream, true, pruneFactor);
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException(e);
        } finally {
            if ( stream != null ) stream.close();
        }
    }

    // TODO -- generalize to support both types of graphs.  Need some kind of display string function
    public void printGraph(final PrintStream graphWriter, final boolean writeHeader, final int pruneFactor) {
        if ( writeHeader )
            graphWriter.println("digraph assemblyGraphs {");

        for( final BaseEdge edge : edgeSet() ) {
//            if( edge.getMultiplicity() > PRUNE_FACTOR ) {
            graphWriter.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [" + (edge.getMultiplicity() <= pruneFactor ? "style=dotted,color=grey," : "") + "label=\"" + edge.getMultiplicity() + "\"];");
//            }
            if( edge.isRef() ) {
                graphWriter.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [color=red];");
            }
            //if( !edge.isRef() && edge.getMultiplicity() <= PRUNE_FACTOR ) { System.out.println("Graph pruning warning!"); }
        }

        for( final T v : vertexSet() ) {
            graphWriter.println("\t" + v.toString() + " [label=\"" + new String(getAdditionalSequence(v)) + "\",shape=box]");
        }

        if ( writeHeader )
            graphWriter.println("}");
    }

    protected void cleanNonRefPaths() {
        if( getReferenceSourceVertex() == null || getReferenceSinkVertex() == null ) {
            return;
        }

        // Remove non-ref edges connected before and after the reference path
        final Set<BaseEdge> edgesToCheck = new HashSet<BaseEdge>();
        edgesToCheck.addAll(incomingEdgesOf(getReferenceSourceVertex()));
        while( !edgesToCheck.isEmpty() ) {
            final BaseEdge e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( incomingEdgesOf(getEdgeSource(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        edgesToCheck.addAll(outgoingEdgesOf(getReferenceSinkVertex()));
        while( !edgesToCheck.isEmpty() ) {
            final BaseEdge e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( outgoingEdgesOf(getEdgeTarget(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        // Run through the graph and clean up singular orphaned nodes
        final List<T> verticesToRemove = new LinkedList<T>();
        for( final T v : vertexSet() ) {
            if( inDegreeOf(v) == 0 && outDegreeOf(v) == 0 ) {
                verticesToRemove.add(v);
            }
        }
        removeAllVertices(verticesToRemove);
    }

    protected void pruneGraph( final int pruneFactor ) {
        final List<BaseEdge> edgesToRemove = new ArrayList<BaseEdge>();
        for( final BaseEdge e : edgeSet() ) {
            if( e.getMultiplicity() <= pruneFactor && !e.isRef() ) { // remove non-reference edges with weight less than or equal to the pruning factor
                edgesToRemove.add(e);
            }
        }
        removeAllEdges(edgesToRemove);

        // Run through the graph and clean up singular orphaned nodes
        final List<T> verticesToRemove = new ArrayList<T>();
        for( final T v : vertexSet() ) {
            if( inDegreeOf(v) == 0 && outDegreeOf(v) == 0 ) {
                verticesToRemove.add(v);
            }
        }

        removeAllVertices(verticesToRemove);
    }

    public void removeVerticesNotConnectedToRef() {
        final HashSet<T> toRemove = new HashSet<T>(vertexSet());
        final HashSet<T> visited = new HashSet<T>();

        final LinkedList<T> toVisit = new LinkedList<T>();
        final T refV = getReferenceSourceVertex();
        if ( refV != null ) {
            toVisit.add(refV);
            while ( ! toVisit.isEmpty() ) {
                final T v = toVisit.pop();
                if ( ! visited.contains(v) ) {
                    toRemove.remove(v);
                    visited.add(v);
                    for ( final T prev : incomingVerticesOf(v) ) toVisit.add(prev);
                    for ( final T next : outgoingVerticesOf(v) ) toVisit.add(next);
                }
            }

//            for ( final T remove : toRemove )
//                logger.info("Cleaning up nodes not attached to any reference node: " + remove.toString());

            removeAllVertices(toRemove);
        }
    }

    public static <T extends BaseVertex> boolean graphEquals(final BaseGraph<T> g1, BaseGraph<T> g2) {
        if( !(g1.vertexSet().containsAll(g2.vertexSet()) && g2.vertexSet().containsAll(g1.vertexSet())) ) {
            return false;
        }
        for( BaseEdge e1 : g1.edgeSet() ) {
            boolean found = false;
            for( BaseEdge e2 : g2.edgeSet() ) {
                if( e1.equals(g1, e2, g2) ) { found = true; break; }
            }
            if( !found ) { return false; }
        }
        for( BaseEdge e2 : g2.edgeSet() ) {
            boolean found = false;
            for( BaseEdge e1 : g1.edgeSet() ) {
                if( e2.equals(g2, e1, g1) ) { found = true; break; }
            }
            if( !found ) { return false; }
        }
        return true;
    }
}
