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

import com.google.java.contract.Ensures;
import htsjdk.samtools.Cigar;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.gatk.utils.sam.CigarUtils;

import java.util.*;

/**
 * A path thought a BaseGraph
 *
 * class to keep track of paths
 *
 * User: depristo
 * Date: 3/19/13
 * Time: 2:34 PM
 *
 */
public class Path<T extends BaseVertex, E extends BaseEdge> {

    // the last vertex seen in the path
    protected final T lastVertex;
    protected final List<E> edgesInOrder;

    // the scores for the path
    protected final int totalScore;

    // the graph from which this path originated
    protected final BaseGraph<T, E> graph;

    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    public Path(final T initialVertex, final BaseGraph<T, E> graph) {
        if ( initialVertex == null ) throw new IllegalArgumentException("initialVertex cannot be null");
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( ! graph.containsVertex(initialVertex) ) throw new IllegalArgumentException("Vertex " + initialVertex + " must be part of graph " + graph);

        lastVertex = initialVertex;
        edgesInOrder = Collections.emptyList();
        totalScore = 0;
        this.graph = graph;
    }

    protected Path(final BaseGraph<T, E> graph, final List<E> edges, final T lastVertex, final int totalScore) {
        this.graph = graph;
        edgesInOrder = edges;
        this.lastVertex = lastVertex;
        this.totalScore = totalScore;
    }

    /**
     * Convenience constructor for testing that creates a path through vertices in graph
     *
     * @deprecated use {@link PathBuilder}.
     */
    protected static <T extends BaseVertex, E extends BaseEdge> Path<T,E> makePath(final List<T> vertices, final BaseGraph<T, E> graph) {
        Path<T,E> path = new Path<T,E>(vertices.get(0), graph);
        for ( int i = 1; i < vertices.size(); i++ )
            path = new Path<T,E>(path, graph.getEdge(path.lastVertex, vertices.get(i)));
        return path;
    }

    /**
     * Create a new path with the same field values.
     *
     * @param p the template path.
     *
     * @throws NullPointerException if {@code p} is {@code null}.
     */
    protected Path(final Path<T,E> p) {
        this.edgesInOrder = p.edgesInOrder;
        this.lastVertex = p.lastVertex;
        this.totalScore = p.totalScore;
        this.graph = p.graph;
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edge the edge to extend path with.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
     */
    public Path(final Path<T,E> p, final E edge) {
        if ( p == null ) throw new IllegalArgumentException("Path cannot be null");
        if ( edge == null ) throw new IllegalArgumentException("Edge cannot be null");
        if ( ! p.graph.containsEdge(edge) ) throw new IllegalArgumentException("Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeSource(edge).equals(p.lastVertex) ) { throw new IllegalStateException("Edges added to path must be contiguous."); }

        graph = p.graph;
        lastVertex = p.graph.getEdgeTarget(edge);
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.addAll(p.edgesInOrder);
        edgesInOrder.add(edge);
        totalScore = p.totalScore + edge.getMultiplicity();
    }

    /**
     * Length of the path in edges.
     *
     * @return {@code 0} or greater.
     */
    public int length() {
        return edgesInOrder.size();
    }

    /**
     * Prepend a path with an edge.
     *
     * @param edge the extending edge.
     * @param p the original path.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a target the first vertex in {@code p}.
     */
    public Path(final E edge, final Path<T,E> p) {
        if ( p == null ) throw new IllegalArgumentException("Path cannot be null");
        if ( edge == null ) throw new IllegalArgumentException("Edge cannot be null");
        if ( ! p.graph.containsEdge(edge) ) throw new IllegalArgumentException("Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeTarget(edge).equals(p.getFirstVertex())) { throw new IllegalStateException("Edges added to path must be contiguous."); }
        graph = p.graph;
        lastVertex = p.lastVertex;
        edgesInOrder = new ArrayList<>(p.length() + 1);
        edgesInOrder.add(edge);
        edgesInOrder.addAll(p.getEdges());
        totalScore = p.totalScore + edge.getMultiplicity();
    }

   /**
     * Checks whether a given path is a suffix of this path.
     *
     * @param other the path to compare against.
     * @throws IllegalArgumentException if <code>other</code> is <code>null</code>, or the come from
     *   different graphs.
     * @return true if <code>other</code> is a suffix of this path.
     */
    public boolean isSuffix(final Path<T, E> other) {
        if ( other == null ) throw new IllegalArgumentException("path cannot be null");
        if (other.getGraph() != this.getGraph()) throw new IllegalArgumentException("the other path most belong to the same path");
        if (!lastVertex.equals(other.lastVertex))
          return false;
        final ListIterator<E> myIt = edgesInOrder.listIterator(edgesInOrder.size());
        final ListIterator<E> otherIt = other.edgesInOrder.listIterator(other.edgesInOrder.size());
        while (myIt.hasPrevious() && otherIt.hasPrevious())
            if (otherIt.previous() != myIt.previous())
                return false;
        return !otherIt.hasPrevious();
    }

    /**
     * Check that two paths have the same edges and total score
     * @param path the other path we might be the same as
     * @return true if this and path are the same
     */
    protected boolean pathsAreTheSame(Path<T,E> path) {
        return totalScore == path.totalScore && edgesInOrder.equals(path.edgesInOrder);
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("Path{score=" + totalScore + ", path=");
        boolean first = true;
        for ( final T v : getVertices() ) {
            if ( first )
                first = false;
            else
                b.append(" -> ");
            b.append(v.getSequenceString());
        }
        b.append('}');
        return b.toString();
    }

    /**
     * Get the graph of this path
     * @return a non-null graph
     */
    @Ensures("result != null")
    public BaseGraph<T, E> getGraph() {
        return graph;
    }

    /**
     * Get the edges of this path in order
     * @return a non-null list of edges
     */
    @Ensures("result != null")
    public List<E> getEdges() { return edgesInOrder; }

    /**
     * Get the list of vertices in this path in order defined by the edges of the path
     * @return a non-null, non-empty list of vertices
     */
    @Ensures({"result != null", "!result.isEmpty()"})
    public List<T> getVertices() {
        if ( getEdges().isEmpty() )
            return Collections.singletonList(lastVertex);
        else {
            final LinkedList<T> vertices = new LinkedList<T>();
            boolean first = true;
            for ( final E e : getEdges() ) {
                if ( first ) {
                    vertices.add(graph.getEdgeSource(e));
                    first = false;
                }
                vertices.add(graph.getEdgeTarget(e));
            }
            return vertices;
        }
    }

    /**
     * Get the total score of this path (bigger is better)
     * @return a positive integer
     */
    @Ensures("result >= 0")
    public int getScore() { return totalScore; }

    /**
     * Get the final vertex of the path
     * @return a non-null vertex
     */
    @Ensures("result != null")
    public T getLastVertex() { return lastVertex; }

    /**
     * Get the first vertex in this path
     * @return a non-null vertex
     */
    public T getFirstVertex() {
        if (edgesInOrder.size() == 0) {
            return lastVertex;
        } else {
            return getGraph().getEdgeSource(edgesInOrder.get(0));
        }
    }

    /**
     * The base sequence for this path. Pull the full sequence for source nodes and then the suffix for all subsequent nodes
     * @return  non-null sequence of bases corresponding to this path
     */
    @Ensures({"result != null"})
    public byte[] getBases() {
        if( getEdges().isEmpty() ) { return graph.getAdditionalSequence(lastVertex); }

        byte[] bases = graph.getAdditionalSequence(graph.getEdgeSource(edgesInOrder.get(0)));
        for( final E e : edgesInOrder ) {
            bases = ArrayUtils.addAll(bases, graph.getAdditionalSequence(graph.getEdgeTarget(e)));
        }
        return bases;
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public  Cigar calculateCigar(final byte[] refSeq) {
        return CigarUtils.calculateCigar(refSeq,getBases());
    }

    /**
     * Tests that this and other have the same score and vertices in the same order with the same seq
     * @param other the other path to consider.  Cannot be null
     * @return true if this and path are equal, false otherwise
     */
    public boolean equalScoreAndSequence(final Path<T,E> other) {
        if ( other == null ) throw new IllegalArgumentException("other cannot be null");
        return getScore() == other.getScore() && equalSequence(other);
    }

    /**
     * Tests that this and other have the same vertices in the same order with the same seq
     * @param other the other path to consider.  Cannot be null
     * @return true if this and path are equal, false otherwise
     */
    public boolean equalSequence(final Path<T,E> other) {
        final List<T> mine = getVertices();
        final List<T> yours = other.getVertices();
        if ( mine.size() == yours.size() ) { // hehehe
            for ( int i = 0; i < mine.size(); i++ )
                if ( ! mine.get(i).seqEquals(yours.get(i)) )
                    return false;
        }
        return true;
    }

}
