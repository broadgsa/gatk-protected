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
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.collections.Pair;

import java.util.*;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V
 *
 * This algorithm creates a new SeqGraph with the following configuration
 *
 * prefix -> has outgoing edges to all seq_i
 * suffix -> has incoming edges for all seq_i
 *
 * There are a few special cases that must be handled.  First, Vi could be simply
 * == to the prefix or the suffix.  These generate direct connections between
 * the prefix and suffix nodes, and they are handled internally by the algorithm.
 *
 * Note that for convenience, we will always create newTop and newBottom nodes, but
 * these may be empty node (i.e., they contain no sequence).  That allows them to be
 * trivially merged, if desired, when the graph is incorporated into an overall
 * graph.
 *
 * The product of this operation is a SeqGraph that contains the split.  There's a
 * function to merge reconnect this graph into the graph that contains the middle nodes
 *
 * The process guarentees a few things about the output:
 *
 * -- Preserves the paths and weights among all vertices
 *
 * It produces a graph that has some unusual properties
 *
 * -- May add nodes with no sequence (isEmpty() == true) to preserve connectivity among the graph
 * -- May introduce edges with no multiplicity to preserve paths through the graph
 *
 * The overall workflow of using this class is simple:
 *
 * find vertices V in graph that you want to split out
 * s = new SharedVertexSequenceSplitter(graph, V)
 * s.updateGraph(graph)
 *
 * to update the graph with the modifications created by this splitter
 *
 * User: depristo
 * Date: 3/22/13
 * Time: 8:31 AM
 */
public class SharedVertexSequenceSplitter {
    final private SeqGraph outer;
    final protected SeqVertex prefixV, suffixV;
    final protected Collection<SeqVertex> toSplits;

    // updated in split routine
    protected SeqGraph splitGraph = null;
    protected Collection<SeqVertex> newMiddles = null;
    protected List<BaseEdge> edgesToRemove = null;

    /**
     * Create a new graph that contains the vertices in toSplitsArg with their shared suffix and prefix
     * sequences extracted out.
     *
     * @param graph the graph containing the vertices in toSplitsArg
     * @param toSplitsArg a collection of vertices to split.  Must be contained within graph, and have only connections
     *                    from a single shared top and/or bottom node
     */
    public SharedVertexSequenceSplitter(final SeqGraph graph, final Collection<SeqVertex> toSplitsArg) {
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( toSplitsArg == null ) throw new IllegalArgumentException("toSplitsArg cannot be null");
        if ( toSplitsArg.size() < 2 ) throw new IllegalArgumentException("Can only split at least 2 vertices but only got " + toSplitsArg);
        if ( ! graph.vertexSet().containsAll(toSplitsArg) ) throw new IllegalArgumentException("graph doesn't contain all of the vertices to split");

        this.outer = graph;
        this.toSplits = toSplitsArg;

        // all of the edges point to the same sink, so it's time to merge
        final Pair<SeqVertex, SeqVertex> prefixAndSuffix = commonPrefixAndSuffixOfVertices(toSplits);
        prefixV = prefixAndSuffix.getFirst();
        suffixV = prefixAndSuffix.getSecond();
    }

    /**
     * Given sequencing that are all equal, does this splitter make those into prefix or suffix nodes?
     * @return true if we merge equal nodes into prefix nodes or suffix nodes
     */
    protected static boolean prefersPrefixMerging() {
        return true;
    }

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @see #updateGraph(SeqVertex, SeqVertex) for a full description of top and bottom
     *
     * @param top the top vertex, may be null
     * @param bottom the bottom vertex, may be null
     * @return true if some useful splitting was done, false otherwise
     */
    public boolean splitAndUpdate(final SeqVertex top, final SeqVertex bottom) {
        split();
        updateGraph(top, bottom);
        return true;
    }

    /**
     * Does either the common suffix or prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if either suffix or prefix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForEitherPrefixOrSuffix(final int minCommonSequence) {
        return meetsMinMergableSequenceForPrefix(minCommonSequence) || meetsMinMergableSequenceForSuffix(minCommonSequence);
    }

    /**
     * Does the common prefix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if prefix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForPrefix(final int minCommonSequence) {
        return prefixV.length() >= minCommonSequence;
    }

    /**
     * Does the common suffix have at least minCommonSequence bases in it?
     * @param minCommonSequence a minimum length of the common sequence, must be >= 0
     * @return true if suffix length >= minCommonSequence
     */
    public boolean meetsMinMergableSequenceForSuffix(final int minCommonSequence) {
        return suffixV.length() >= minCommonSequence;
    }

    /**
     * Actually do the splitting up of the vertices
     *
     * Must be called before calling updateGraph
     */
    public void split() {
        splitGraph = new SeqGraph(outer.getKmerSize());
        newMiddles = new LinkedList<SeqVertex>();
        edgesToRemove = new LinkedList<BaseEdge>();

        splitGraph.addVertices(prefixV, suffixV);

        for ( final SeqVertex mid : toSplits ) {
            final BaseEdge toMid = processEdgeToRemove(mid, outer.incomingEdgeOf(mid));
            final BaseEdge fromMid = processEdgeToRemove(mid, outer.outgoingEdgeOf(mid));

            final SeqVertex remaining = mid.withoutPrefixAndSuffix(prefixV.getSequence(), suffixV.getSequence());
            if ( remaining != null ) {
                // there's some sequence prefix + seq + suffix, so add the node and make edges
                splitGraph.addVertex(remaining);
                newMiddles.add(remaining);
                // update edge from top -> middle to be top -> without suffix
                splitGraph.addEdge(prefixV, remaining, toMid);
                splitGraph.addEdge(remaining, suffixV, fromMid);
            } else {
                // prefix + suffix completely explain this node
                splitGraph.addOrUpdateEdge(prefixV, suffixV, toMid.copy().add(fromMid));
            }
        }
    }

    /**
     * Update graph outer, replacing the previous middle vertices that were split out with the new
     * graph structure of the split, linking this subgraph into the graph at top and bot (the
     * vertex connecting the middle nodes and the vertex outgoing of all middle node)
     *
     * @param top an optional top node that must have outgoing edges to all split vertices.  If null, this subgraph
     *            will be added without any incoming edges
     * @param bot an optional bottom node that must have incoming edges to all split vertices.  If null, this subgraph
     *            will be added without any outgoing edges to the rest of the graph
     */
    public void updateGraph(final SeqVertex top, final SeqVertex bot) {
        if ( ! outer.vertexSet().containsAll(toSplits) ) throw new IllegalArgumentException("graph doesn't contain all of the original vertices to split");
        if ( top == null && bot == null ) throw new IllegalArgumentException("Cannot update graph without at least one top or bot vertex, but both were null");
        if ( top != null && ! outer.containsVertex(top) ) throw new IllegalArgumentException("top " + top + " not found in graph " + outer);
        if ( bot != null && ! outer.containsVertex(bot) ) throw new IllegalArgumentException("bot " + bot + " not found in graph " + outer);
        if ( splitGraph == null ) throw new IllegalStateException("Cannot call updateGraph until split() has been called");

        outer.removeAllVertices(toSplits);
        outer.removeAllEdges(edgesToRemove);

        outer.addVertices(newMiddles);

        final boolean hasPrefixSuffixEdge = splitGraph.getEdge(prefixV, suffixV) != null;
        final boolean hasOnlyPrefixSuffixEdges = hasPrefixSuffixEdge && splitGraph.outDegreeOf(prefixV) == 1;
        final boolean needPrefixNode = ! prefixV.isEmpty() || (top == null && ! hasOnlyPrefixSuffixEdges);
        final boolean needSuffixNode = ! suffixV.isEmpty() || (bot == null && ! hasOnlyPrefixSuffixEdges);

        // if prefix / suffix are needed, keep them
        final SeqVertex topForConnect = needPrefixNode ? prefixV : top;
        final SeqVertex botForConnect = needSuffixNode ? suffixV : bot;

        if ( needPrefixNode ) {
            outer.addVertex(prefixV);
            if ( top != null ) outer.addEdge(top, prefixV, BaseEdge.orRef(splitGraph.outgoingEdgesOf(prefixV), 1));
        }

        if ( needSuffixNode ) {
            outer.addVertex(suffixV);
            if ( bot != null ) outer.addEdge(suffixV, bot, BaseEdge.orRef(splitGraph.incomingEdgesOf(suffixV), 1));
        }

        if ( topForConnect != null ) {
            for ( final BaseEdge e : splitGraph.outgoingEdgesOf(prefixV) ) {
                final SeqVertex target = splitGraph.getEdgeTarget(e);

                if ( target == suffixV ) { // going straight from prefix -> suffix
                    if ( botForConnect != null )
                        outer.addEdge(topForConnect, botForConnect, e);
                } else {
                    outer.addEdge(topForConnect, target, e);
                }
            }
        }

        if ( botForConnect != null ) {
            for ( final BaseEdge e : splitGraph.incomingEdgesOf(suffixV) ) {
                outer.addEdge(splitGraph.getEdgeSource(e), botForConnect, e);
            }
        }
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
    protected static Pair<SeqVertex, SeqVertex> commonPrefixAndSuffixOfVertices(final Collection<SeqVertex> middleVertices) {
        final List<byte[]> kmers = new ArrayList<byte[]>(middleVertices.size());

        int min = Integer.MAX_VALUE;
        for ( final SeqVertex v : middleVertices ) {
            kmers.add(v.getSequence());
            min = Math.min(min, v.getSequence().length);
        }

        final int prefixLen = GraphUtils.compPrefixLen(kmers, min);
        final int suffixLen = GraphUtils.compSuffixLen(kmers, min - prefixLen);

        final byte[] kmer = kmers.get(0);
        final byte[] prefix = Arrays.copyOfRange(kmer, 0, prefixLen);
        final byte[] suffix = Arrays.copyOfRange(kmer, kmer.length - suffixLen, kmer.length);
        return new Pair<SeqVertex, SeqVertex>(new SeqVertex(prefix), new SeqVertex(suffix));
    }

    /**
     * Helper function that returns an edge that we should use for splitting
     *
     * If e is null, creates a new 0 multiplicity edge, set to ref is any edges to V are ref
     * If e is not null, returns a new copy of e, and schedules e for removal
     *
     * @param e a non-null edge
     * @return a non-null edge
     */
    @Requires("v != null")
    @Ensures("result != null")
    private BaseEdge processEdgeToRemove(final SeqVertex v, final BaseEdge e) {
        if ( e == null ) {
            // there's no edge, so we return a newly allocated one and don't schedule e for removal
            // the weight must be 0 to preserve sum through the diamond
            return new BaseEdge(outer.isReferenceNode(v), 0);
        } else {
            // schedule edge for removal, and return a freshly allocated one for our graph to use
            edgesToRemove.add(e);
            return e.copy();
        }
    }
}