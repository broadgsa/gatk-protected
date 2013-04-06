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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

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
public class Path<T extends BaseVertex> {
    private final static int MAX_CIGAR_ELEMENTS_BEFORE_FAILING_SW = 20;

    // the last vertex seen in the path
    private final T lastVertex;

    // the list of edges comprising the path
    private Set<BaseEdge> edgesAsSet = null;
    private final LinkedList<BaseEdge> edgesInOrder;

    // the scores for the path
    private final int totalScore;

    // the graph from which this path originated
    private final BaseGraph<T> graph;

    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    private static final double SW_MATCH = 20.0;
    private static final double SW_MISMATCH = -15.0;
    private static final double SW_GAP = -26.0;
    private static final double SW_GAP_EXTEND = -1.1;
    private static final byte[] STARTING_SW_ANCHOR_BYTES = "XXXXXXXXX".getBytes();

    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path with follow through
     */
    public Path(final T initialVertex, final BaseGraph<T> graph) {
        if ( initialVertex == null ) throw new IllegalArgumentException("initialVertex cannot be null");
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( ! graph.containsVertex(initialVertex) ) throw new IllegalArgumentException("Vertex " + initialVertex + " must be part of graph " + graph);

        lastVertex = initialVertex;
        edgesInOrder = new LinkedList<BaseEdge>();
        totalScore = 0;
        this.graph = graph;
    }

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend
     * @param edge the edge to extend path by
     */
    public Path(final Path<T> p, final BaseEdge edge) {
        if ( p == null ) throw new IllegalArgumentException("Path cannot be null");
        if ( edge == null ) throw new IllegalArgumentException("Edge cannot be null");
        if ( ! p.graph.containsEdge(edge) ) throw new IllegalArgumentException("Graph must contain edge " + edge + " but it doesn't");
        if ( ! p.graph.getEdgeSource(edge).equals(p.lastVertex) ) { throw new IllegalStateException("Edges added to path must be contiguous."); }

        graph = p.graph;
        lastVertex = p.graph.getEdgeTarget(edge);
        edgesInOrder = new LinkedList<BaseEdge>(p.getEdges());
        edgesInOrder.add(edge);
        totalScore = p.totalScore + edge.getMultiplicity();
    }

    /**
     * Get the collection of edges leaving the last vertex of this path
     * @return a non-null collection
     */
    public Collection<BaseEdge> getOutgoingEdgesOfLastVertex() {
        return getGraph().outgoingEdgesOf(getLastVertex());
    }

    /**
     * Does this path contain the given edge
     * @param edge  the given edge to test
     * @return      true if the edge is found in this path
     */
    public boolean containsEdge( final BaseEdge edge ) {
        if( edge == null ) { throw new IllegalArgumentException("Attempting to test null edge."); }
        if ( edgesInOrder.isEmpty() ) return false;

        // initialize contains cache if necessary
        if ( edgesAsSet == null ) edgesAsSet = new HashSet<BaseEdge>(edgesInOrder);
        return edgesAsSet.contains(edge);
    }

    /**
     * Does this path contain the given vertex?
     *
     * @param v a non-null vertex
     * @return true if v occurs within this path, false otherwise
     */
    public boolean containsVertex(final T v) {
        if ( v == null ) throw new IllegalArgumentException("Vertex cannot be null");

        // TODO -- warning this is expensive.  Need to do vertex caching
        return getVertices().contains(v);
    }

    /**
     * Check that two paths have the same edges and total score
     * @param path the other path we might be the same as
     * @return true if this and path are the same
     */
    protected boolean pathsAreTheSame(Path<T> path) {
        return totalScore == path.totalScore && edgesInOrder.equals(path.edgesInOrder);
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("Path{score=" + totalScore + ", path=");
        boolean first = true;
        for ( final T v : getVertices() ) {
            if ( first ) {
                first = false;
            } else {
                b.append(" -> ");
            }
            b.append(v.getSequenceString());
        }
        return b.toString();
    }

    /**
     * Get the graph of this path
     * @return a non-null graph
     */
    @Ensures("result != null")
    public BaseGraph<T> getGraph() {
        return graph;
    }

    /**
     * Get the edges of this path in order
     * @return a non-null list of edges
     */
    @Ensures("result != null")
    public List<BaseEdge> getEdges() { return edgesInOrder; }

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
            for ( final BaseEdge e : getEdges() ) {
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
     * The base sequence for this path. Pull the full sequence for source nodes and then the suffix for all subsequent nodes
     * @return  non-null sequence of bases corresponding to this path
     */
    @Ensures({"result != null"})
    public byte[] getBases() {
        if( getEdges().isEmpty() ) { return graph.getAdditionalSequence(lastVertex); }

        byte[] bases = graph.getAdditionalSequence(graph.getEdgeSource(edgesInOrder.getFirst()));
        for( final BaseEdge e : edgesInOrder ) {
            bases = ArrayUtils.addAll(bases, graph.getAdditionalSequence(graph.getEdgeTarget(e)));
        }
        return bases;
    }

    /**
     * Calculate the cigar string for this path using a bubble traversal of the assembly graph and running a Smith-Waterman alignment on each bubble
     * @return  non-null Cigar string with reference length equal to the refHaplotype's reference length
     */
    @Ensures("result != null")
    public Cigar calculateCigar() {
        final Cigar cigar = new Cigar();
        // special case for paths that start on reference but not at the reference source node
        if( edgesInOrder.getFirst().isRef() && !graph.isRefSource(edgesInOrder.getFirst()) ) {
            for( final CigarElement ce : calculateCigarForCompleteBubble(null, null, graph.getEdgeSource(edgesInOrder.getFirst())).getCigarElements() ) {
                cigar.add(ce);
            }
        }

        // reset the bubble state machine
        final BubbleStateMachine<T> bsm = new BubbleStateMachine<T>(cigar);

        for( final BaseEdge e : getEdges() ) {
            if ( e.hasSameSourceAndTarget(graph, edgesInOrder.getFirst()) ) {
                advanceBubbleStateMachine( bsm, graph.getEdgeSource(e), null );
            }
            advanceBubbleStateMachine( bsm, graph.getEdgeTarget(e), e );
        }

        // special case for paths that don't end on reference
        if( bsm.inBubble ) {
            for( final CigarElement ce : calculateCigarForCompleteBubble(bsm.bubbleBytes, bsm.lastSeenReferenceNode, null).getCigarElements() ) {
                bsm.cigar.add(ce);
            }
        } else if( edgesInOrder.getLast().isRef() && !graph.isRefSink(edgesInOrder.getLast()) ) { // special case for paths that end of the reference but haven't completed the entire reference circuit
            for( final CigarElement ce : calculateCigarForCompleteBubble(bsm.bubbleBytes, graph.getEdgeTarget(edgesInOrder.getLast()), null).getCigarElements() ) {
                bsm.cigar.add(ce);
            }
        }

        return AlignmentUtils.consolidateCigar(bsm.cigar);
    }

    /**
     * Advance the bubble state machine by incorporating the next node in the path.
     * @param bsm   the current bubble state machine
     * @param node  the node to be incorporated
     * @param e     the edge which generated this node in the path
     */
    @Requires({"bsm != null", "graph != null", "node != null"})
    private void advanceBubbleStateMachine( final BubbleStateMachine<T> bsm, final T node, final BaseEdge e ) {
        if( graph.isReferenceNode( node ) ) {
            if( !bsm.inBubble ) { // just add the ref bases as M's in the Cigar string, and don't do anything else
                if( e !=null && !e.isRef() ) {
                    if( graph.referencePathExists( graph.getEdgeSource(e), node) ) {
                        for( final CigarElement ce : calculateCigarForCompleteBubble(null, graph.getEdgeSource(e), node).getCigarElements() ) {
                            bsm.cigar.add(ce);
                        }
                        bsm.cigar.add( new CigarElement( graph.getAdditionalSequence(node).length, CigarOperator.M) );
                    } else if ( graph.getEdgeSource(e).equals(graph.getEdgeTarget(e)) ) { // alt edge at ref node points to itself
                        bsm.cigar.add( new CigarElement( graph.getAdditionalSequence(node).length, CigarOperator.I) );
                    } else {
                        bsm.inBubble = true;
                        bsm.bubbleBytes = null;
                        bsm.lastSeenReferenceNode = graph.getEdgeSource(e);
                        bsm.bubbleBytes = ArrayUtils.addAll( bsm.bubbleBytes, graph.getAdditionalSequence(node) );
                    }
                } else {
                    bsm.cigar.add( new CigarElement( graph.getAdditionalSequence(node).length, CigarOperator.M) );
                }
            } else if( bsm.lastSeenReferenceNode != null && !graph.referencePathExists( bsm.lastSeenReferenceNode, node ) ) { // add bases to the bubble string until we get back to the reference path
                bsm.bubbleBytes = ArrayUtils.addAll( bsm.bubbleBytes, graph.getAdditionalSequence(node) );
            } else { // close the bubble and use a local SW to determine the Cigar string
                for( final CigarElement ce : calculateCigarForCompleteBubble(bsm.bubbleBytes, bsm.lastSeenReferenceNode, node).getCigarElements() ) {
                    bsm.cigar.add(ce);
                }
                bsm.inBubble = false;
                bsm.bubbleBytes = null;
                bsm.lastSeenReferenceNode = null;
                bsm.cigar.add( new CigarElement( graph.getAdditionalSequence(node).length, CigarOperator.M) );
            }
        } else { // non-ref vertex
            if( bsm.inBubble ) { // just keep accumulating until we get back to the reference path
                bsm.bubbleBytes = ArrayUtils.addAll( bsm.bubbleBytes, graph.getAdditionalSequence(node) );
            } else { // open up a bubble
                bsm.inBubble = true;
                bsm.bubbleBytes = null;
                bsm.lastSeenReferenceNode = (e != null ? graph.getEdgeSource(e) : null );
                bsm.bubbleBytes = ArrayUtils.addAll( bsm.bubbleBytes, graph.getAdditionalSequence(node) );
            }
        }
    }

    /**
     * Now that we have a completed bubble run a Smith-Waterman alignment to determine the cigar string for this bubble
     * @param bubbleBytes   the bytes that comprise the alternate allele path in this bubble
     * @param fromVertex    the vertex that marks the beginning of the reference path in this bubble (null indicates ref source vertex)
     * @param toVertex      the vertex that marks the end of the reference path in this bubble (null indicates ref sink vertex)
     * @return              the cigar string generated by running a SW alignment between the reference and alternate paths in this bubble
     */
    @Requires({"graph != null"})
    @Ensures({"result != null"})
    private Cigar calculateCigarForCompleteBubble( final byte[] bubbleBytes, final T fromVertex, final T toVertex ) {
        final byte[] refBytes = graph.getReferenceBytes(fromVertex == null ? graph.getReferenceSourceVertex() : fromVertex, toVertex == null ? graph.getReferenceSinkVertex() : toVertex, fromVertex == null, toVertex == null);

        final Cigar returnCigar = new Cigar();

        // add padding to anchor ref/alt bases in the SW matrix
        byte[] padding = STARTING_SW_ANCHOR_BYTES;
        boolean goodAlignment = false;
        SWPairwiseAlignment swConsensus = null;
        while( !goodAlignment && padding.length < 1000 ) {
            padding = ArrayUtils.addAll(padding, padding); // double the size of the padding each time
            final byte[] reference = ArrayUtils.addAll( ArrayUtils.addAll(padding, refBytes), padding );
            final byte[] alternate = ArrayUtils.addAll( ArrayUtils.addAll(padding, bubbleBytes), padding );
            swConsensus = new SWPairwiseAlignment( reference, alternate, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
            if( swConsensus.getAlignmentStart2wrt1() == 0 && !swConsensus.getCigar().toString().contains("S") && swConsensus.getCigar().getReferenceLength() == reference.length ) {
                goodAlignment = true;
            }
        }
        if( !goodAlignment ) {
            returnCigar.add(new CigarElement(1, CigarOperator.N));
            return returnCigar;
        }

        final Cigar swCigar = swConsensus.getCigar();
        if( swCigar.numCigarElements() > MAX_CIGAR_ELEMENTS_BEFORE_FAILING_SW ) { // this bubble is too divergent from the reference
            returnCigar.add(new CigarElement(1, CigarOperator.N));
        } else {
            for( int iii = 0; iii < swCigar.numCigarElements(); iii++ ) {
                // now we need to remove the padding from the cigar string
                int length = swCigar.getCigarElement(iii).getLength();
                if( iii == 0 ) { length -= padding.length; }
                if( iii == swCigar.numCigarElements() - 1 ) { length -= padding.length; }
                if( length > 0 ) {
                    returnCigar.add(new CigarElement(length, swCigar.getCigarElement(iii).getOperator()));
                }
            }
            if( (refBytes == null && returnCigar.getReferenceLength() != 0) || ( refBytes != null && returnCigar.getReferenceLength() != refBytes.length ) ) {
                throw new IllegalStateException("SmithWaterman cigar failure: " + (refBytes == null ? "-" : new String(refBytes)) + " against " + new String(bubbleBytes) + " = " + swConsensus.getCigar());
            }
        }

        return returnCigar;
    }

    // class to keep track of the bubble state machine
    private static class BubbleStateMachine<T extends BaseVertex> {
        public boolean inBubble = false;
        public byte[] bubbleBytes = null;
        public T lastSeenReferenceNode = null;
        public Cigar cigar = null;

        public BubbleStateMachine( final Cigar initialCigar ) {
            inBubble = false;
            bubbleBytes = null;
            lastSeenReferenceNode = null;
            cigar = initialCigar;
        }
    }

    /**
     * Tests that this and other have the same score and vertices in the same order with the same seq
     * @param other the other path to consider.  Cannot be null
     * @return true if this and path are equal, false otherwise
     */
    public boolean equalScoreAndSequence(final Path<T> other) {
        if ( other == null ) throw new IllegalArgumentException("other cannot be null");
        return getScore() == other.getScore() && equalSequence(other);
    }

    /**
     * Tests that this and other have the same vertices in the same order with the same seq
     * @param other the other path to consider.  Cannot be null
     * @return true if this and path are equal, false otherwise
     */
    public boolean equalSequence(final Path<T> other) {
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
