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
import com.google.java.contract.Requires;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;

import java.io.File;
import java.util.*;

/**
 * A graph that contains base sequence at each node
 *
 * @author: depristo
 * @since 03/2013
 */
public class SeqGraph extends BaseGraph<SeqVertex> {
    /**
     * Construct an empty SeqGraph
     */
    public SeqGraph() {
        super();
    }

    /**
     * Construct an empty SeqGraph where we'll add nodes based on a kmer size of kmer
     *
     * The kmer size is purely information.  It is useful when converting a Debruijn graph -> SeqGraph
     * for us to track the kmer used to make the transformation.
     *
     * @param kmer kmer
     */
    public SeqGraph(final int kmer) {
        super(kmer);
    }

    /**
     * Simplify this graph, merging vertices together and restructuring the graph in an
     * effort to minimize the number of overall vertices in the graph without changing
     * in any way the sequences implied by a complex enumeration of all paths through the graph.
     */
    public void simplifyGraph() {
        zipLinearChains();
        mergeBranchingNodes();
        zipLinearChains();
    }

    /**
     * Zip up all of the simple linear chains present in this graph.
     */
    protected void zipLinearChains() {
        while( zipOneLinearChain() ) {
            // just keep going until zipOneLinearChain says its done
        }
    }

    /**
     * Merge together two vertices in the graph v1 -> v2 into a single vertex v' containing v1 + v2 sequence
     *
     * Only works on vertices where v1's only outgoing edge is to v2 and v2's only incoming edge is from v1.
     *
     * If such a pair of vertices is found, they are merged and the graph is update.  Otherwise nothing is changed.
     *
     * @return true if any such pair of vertices could be found, false otherwise
     */
    protected boolean zipOneLinearChain() {
        for( final BaseEdge e : edgeSet() ) {
            final SeqVertex outgoingVertex = getEdgeTarget(e);
            final SeqVertex incomingVertex = getEdgeSource(e);
            if( !outgoingVertex.equals(incomingVertex)
                    && outDegreeOf(incomingVertex) == 1 && inDegreeOf(outgoingVertex) == 1
                    && isReferenceNode(incomingVertex) == isReferenceNode(outgoingVertex) ) {

                final Set<BaseEdge> outEdges = outgoingEdgesOf(outgoingVertex);
                final Set<BaseEdge> inEdges = incomingEdgesOf(incomingVertex);
                if( inEdges.size() == 1 && outEdges.size() == 1 ) {
                    inEdges.iterator().next().setMultiplicity( inEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() / 2 ) );
                    outEdges.iterator().next().setMultiplicity( outEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() / 2 ) );
                } else if( inEdges.size() == 1 ) {
                    inEdges.iterator().next().setMultiplicity( inEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() - 1 ) );
                } else if( outEdges.size() == 1 ) {
                    outEdges.iterator().next().setMultiplicity( outEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() - 1 ) );
                }

                final SeqVertex addedVertex = new SeqVertex( ArrayUtils.addAll(incomingVertex.getSequence(), outgoingVertex.getSequence()) );
                addVertex(addedVertex);
                for( final BaseEdge edge : outEdges ) {
                    addEdge(addedVertex, getEdgeTarget(edge), new BaseEdge(edge.isRef(), edge.getMultiplicity()));
                }
                for( final BaseEdge edge : inEdges ) {
                    addEdge(getEdgeSource(edge), addedVertex, new BaseEdge(edge.isRef(), edge.getMultiplicity()));
                }

                removeVertex(incomingVertex);
                removeVertex(outgoingVertex);
                return true;
            }
        }

        return false;
    }

    /**
     * Perform as many branch simplifications and merging operations as possible on this graph,
     * modifying it in place.
     */
    protected void mergeBranchingNodes() {
        boolean foundNodesToMerge = true;
        while( foundNodesToMerge ) {
            foundNodesToMerge = false;

            for( final SeqVertex v : vertexSet() ) {
                foundNodesToMerge = simplifyDiamondIfPossible(v);
                if ( foundNodesToMerge )
                    break;
            }
        }
    }

    /**
     * A simple structure that looks like:
     *
     *      v
     *    / |  \      \
     *   m1 m2 m3 ... mn
     *    \ |  /      /
     *      b
     *
     * Only returns true if all outgoing edges of v go to vertices that all only connect to
     * a single bottom node, and that all middle nodes have only the single edge
     *
     * @param v the vertex to test if its the top of a diamond pattern
     * @return true if v is the root of a diamond
     */
    protected boolean isRootOfDiamond(final SeqVertex v) {
        final Set<BaseEdge> ve = outgoingEdgesOf(v);
        if ( ve.size() <= 1 )
            return false;

        SeqVertex bottom = null;
        for ( final BaseEdge e : ve ) {
            final SeqVertex mi = getEdgeTarget(e);

            // all nodes must have at least 1 connection
            if ( outDegreeOf(mi) < 1 )
                return false;

            // can only have 1 incoming node, the root vertex
            if ( inDegreeOf(mi) != 1 )
                return false;

            // make sure that all outgoing vertices of mi go only to the bottom node
            for ( final SeqVertex mt : outgoingVerticesOf(mi) ) {
                if ( bottom == null )
                    bottom = mt;
                else if ( ! bottom.equals(mt) )
                    return false;
            }
        }

        // bottom has some connections coming in from other nodes, don't allow
        if ( inDegreeOf(bottom) != ve.size() )
            return false;

        return true;
    }

    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return
     */
    @Requires("!middleVertices.isEmpty()")
    private byte[] commonSuffixOfEdgeTargets(final Set<SeqVertex> middleVertices) {
        final String[] kmers = new String[middleVertices.size()];

        int i = 0;
        for ( final SeqVertex v : middleVertices ) {
            kmers[i++] = (StringUtils.reverse(v.getSequenceString()));
        }

        final String commonPrefix = StringUtils.getCommonPrefix(kmers);
        return commonPrefix.equals("") ? null : StringUtils.reverse(commonPrefix).getBytes();
    }

    /**
     * Get the node that is the bottom of a diamond configuration in the graph starting at top
     *
     * @param top
     * @return
     */
    @Requires("top != null")
    @Ensures({"result != null"})
    private SeqVertex getDiamondBottom(final SeqVertex top) {
        final BaseEdge topEdge = outgoingEdgesOf(top).iterator().next();
        final SeqVertex middle = getEdgeTarget(topEdge);
        final BaseEdge middleEdge = outgoingEdgesOf(middle).iterator().next();
        return getEdgeTarget(middleEdge);
    }

    /**
     * Get the set of vertices that are in the middle of a diamond starting at top
     * @param top
     * @return
     */
    @Requires("top != null")
    @Ensures({"result != null", "!result.isEmpty()"})
    final Set<SeqVertex> getMiddleVertices(final SeqVertex top) {
        final Set<SeqVertex> middles = new HashSet<SeqVertex>();
        for ( final BaseEdge topToMiddle : outgoingEdgesOf(top) ) {
            middles.add(getEdgeTarget(topToMiddle));
        }
        return middles;
    }

    /**
     * Simply a diamond configuration in the current graph starting at top, if possible
     *
     * If top is actually the top of a diamond that can be simplified (i.e., doesn't have any
     * random edges or other structure that would cause problems with the transformation), then this code
     * performs the following transformation on this graph (modifying it):
     *
     * A -> M1 -> B, A -> M2 -> B, A -> Mn -> B
     *
     * becomes
     *
     * A -> M1' -> B', A -> M2' -> B', A -> Mn' -> B'
     *
     * where B' is composed of the longest common suffix of all Mi nodes + B, and Mi' are each
     * middle vertex without their shared suffix.
     *
     * @param top a proposed vertex in this graph that might start a diamond (but doesn't have to)
     * @return true top actually starts a diamond and it could be simplified
     */
    private boolean simplifyDiamondIfPossible(final SeqVertex top) {
        if ( ! isRootOfDiamond(top) )
            return false;

        final SeqVertex diamondBottom = getDiamondBottom(top);
        final Set<SeqVertex> middleVertices = getMiddleVertices(top);
        final List<SeqVertex> verticesToRemove = new LinkedList<SeqVertex>();
        final List<BaseEdge> edgesToRemove = new LinkedList<BaseEdge>();

        // all of the edges point to the same sink, so it's time to merge
        final byte[] commonSuffix = commonSuffixOfEdgeTargets(middleVertices);
        if ( commonSuffix != null ) {
            final BaseEdge botToNewBottom = new BaseEdge(false, 0);
            final BaseEdge elimMiddleNodeEdge = new BaseEdge(false, 0);
            final SeqVertex newBottomV = new SeqVertex(commonSuffix);
            addVertex(newBottomV);

            for ( final SeqVertex middle : middleVertices ) {
                final SeqVertex withoutSuffix = middle.withoutSuffix(commonSuffix);
                final BaseEdge topToMiddleEdge = getEdge(top, middle);
                final BaseEdge middleToBottomE = getEdge(middle, diamondBottom);

                // clip out the two edges, since we'll be replacing them later
                edgesToRemove.add(topToMiddleEdge);
                edgesToRemove.add(middleToBottomE);

                if ( withoutSuffix != null ) { // this node is a deletion
                    addVertex(withoutSuffix);
                    // update edge from top -> middle to be top -> without suffix
                    addEdge(top, withoutSuffix, new BaseEdge(topToMiddleEdge));
                    addEdge(withoutSuffix, newBottomV, new BaseEdge(middleToBottomE));
                } else {
                    // this middle node is == the common suffix, wo we're removing the edge
                    elimMiddleNodeEdge.add(topToMiddleEdge);
                }
                // include the ref and multi of mid -> bot in our edge from new bot -> bot
                botToNewBottom.add(middleToBottomE);
                verticesToRemove.add(middle);
            }

            // add an edge from top to new bottom, because some middle nodes were removed
            if ( elimMiddleNodeEdge.getMultiplicity() > 0 )
                addEdge(top, newBottomV, elimMiddleNodeEdge);

            addEdge(newBottomV, diamondBottom, botToNewBottom);

            removeAllEdges(edgesToRemove);
            removeAllVertices(verticesToRemove);
            return true;
        } else {
            return false;
        }
    }
}
