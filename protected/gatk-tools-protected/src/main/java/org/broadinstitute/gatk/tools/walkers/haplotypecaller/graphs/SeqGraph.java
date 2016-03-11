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
import org.jgrapht.EdgeFactory;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * A graph that contains base sequence at each node
 *
 * @author: depristo
 * @since 03/2013
 */
public class SeqGraph extends BaseGraph<SeqVertex, BaseEdge> {
    /**
     * Edge factory that creates non-reference multiplicity 1 edges
     */
    private static class MyEdgeFactory implements EdgeFactory<SeqVertex, BaseEdge> {
        @Override
        public BaseEdge createEdge(SeqVertex sourceVertex, SeqVertex targetVertex) {
            return new BaseEdge(false, 1);
        }
    }

    private final static boolean PRINT_SIMPLIFY_GRAPHS = false;

    /**
     * The minimum number of common bp from the prefix (head merging) or suffix (tail merging)
     * required before we'll merge in such configurations.  A large value here is critical to avoid
     * merging inappropriate head or tail nodes, which introduces large insertion / deletion events
     * as the merge operation creates a link among the non-linked sink / source vertices
     */
    protected final static int MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES = 10;

    /**
     * How many cycles of the graph simplifications algorithms will we run before
     * thinking something has gone wrong and throw an exception?
     */
    private final static int MAX_REASONABLE_SIMPLIFICATION_CYCLES = 100;

    /**
     * Construct an empty SeqGraph where we'll add nodes based on a kmer size of kmer
     *
     * The kmer size is purely information.  It is useful when converting a Debruijn graph -> SeqGraph
     * for us to track the kmer used to make the transformation.
     *
     * @param kmer kmer
     */
    public SeqGraph(final int kmer) {
        super(kmer, new MyEdgeFactory());
    }

    /**
     * Simplify this graph, merging vertices together and restructuring the graph in an
     * effort to minimize the number of overall vertices in the graph without changing
     * in any way the sequences implied by a complex enumeration of all paths through the graph.
     */
    public void simplifyGraph() {
        simplifyGraph(Integer.MAX_VALUE);
    }

    protected void simplifyGraph(final int maxCycles) {
        // start off with one round of zipping of chains for performance reasons
        zipLinearChains();

        SeqGraph prevGraph = null;
        for( int i = 0; i < maxCycles; i++ ) {
            if ( i > MAX_REASONABLE_SIMPLIFICATION_CYCLES ) {
                logger.warn("Infinite loop detected in simpliciation routines.  Writing current graph to debugMeMark.dot");
                printGraph(new File("debugMeMark.dot"), 0);
                throw new IllegalStateException("Infinite loop detected in simplification routines for kmer graph " + getKmerSize());
            }

            final boolean didSomeWork = simplifyGraphOnce(i);
            if ( ! didSomeWork )
                // no simplification algorithm could run, so stop
                break;

            // we get five cycles before we start looking for changes in the graph
            // by cloning ourselves and then checking for any changes
            if ( i > 5 ) {
                // the previous graph and this graph have the same structure, so the simplification
                // algorithms are looping endless between states.  Just break and consider ourselves done
                if ( prevGraph != null && graphEquals(prevGraph, this) )
                    break;

                prevGraph = (SeqGraph)clone();
            }
        }
    }

    /**
     * Run one full cycle of the graph simplification algorithms
     * @return true if any algorithms said they did some simplification
     */
    private boolean simplifyGraphOnce(final int iteration) {
        //logger.info("simplifyGraph iteration " + i);
        // iterate until we haven't don't anything useful
        boolean didSomeWork = false;
        printGraphSimplification(new File("simplifyGraph." + iteration + ".1.dot"));
        didSomeWork |= new MergeDiamonds().transformUntilComplete();
        didSomeWork |= new MergeTails().transformUntilComplete();
        printGraphSimplification(new File("simplifyGraph." + iteration + ".2.diamonds_and_tails.dot"));

        didSomeWork |= new SplitCommonSuffices().transformUntilComplete();
        printGraphSimplification(new File("simplifyGraph." + iteration + ".3.split_suffix.dot"));
        didSomeWork |= new MergeCommonSuffices().transformUntilComplete();
        printGraphSimplification(new File("simplifyGraph." + iteration + ".4.merge_suffix.dot"));

        didSomeWork |= zipLinearChains();
        return didSomeWork;
    }

    /**
     * Print simplication step of this graph, if PRINT_SIMPLIFY_GRAPHS is enabled
     * @param file the destination for the graph DOT file
     */
    private void printGraphSimplification(final File file) {
        if ( PRINT_SIMPLIFY_GRAPHS )
            subsetToNeighbors(getReferenceSourceVertex(), 5).printGraph(file, 0);
    }

    /**
     * Zip up all of the simple linear chains present in this graph.
     *
     * Merges together all pairs of vertices in the graph v1 -> v2 into a single vertex v' containing v1 + v2 sequence
     *
     * Only works on vertices where v1's only outgoing edge is to v2 and v2's only incoming edge is from v1.
     *
     * If such a pair of vertices is found, they are merged and the graph is update.  Otherwise nothing is changed.
     *
     * @return true if any such pair of vertices could be found, false otherwise
     */
    public boolean zipLinearChains() {
        // create the list of start sites [doesn't modify graph yet]
        final List<SeqVertex> zipStarts = new LinkedList<SeqVertex>();
        for ( final SeqVertex source : vertexSet() ) {
            if ( isLinearChainStart(source) )
                zipStarts.add(source);
        }

        if ( zipStarts.isEmpty() ) // nothing to do, as nothing could start a chain
            return false;

        // At this point, zipStarts contains all of the vertices in this graph that might start some linear
        // chain of vertices.  We walk through each start, building up the linear chain of vertices and then
        // zipping them up with mergeLinearChain, if possible
        boolean mergedOne = false;
        for ( final SeqVertex zipStart : zipStarts ) {
            final LinkedList<SeqVertex> linearChain = traceLinearChain(zipStart);

            // merge the linearized chain, recording if we actually did some useful work
            mergedOne |= mergeLinearChain(linearChain);
        }

        return mergedOne;
    }

    /**
     * Is source vertex potentially a start of a linear chain of vertices?
     *
     * We are a start of a zip chain if our out degree is 1 and either the
     * the vertex has no incoming connections or 2 or more (we must start a chain) or
     * we have exactly one incoming vertex and that one has out-degree > 1 (i.e., source's incoming
     * vertex couldn't be a start itself
     *
     * @param source a non-null vertex
     * @return true if source might start a linear chain
     */
    @Requires("source != null")
    private boolean isLinearChainStart(final SeqVertex source) {
        return outDegreeOf(source) == 1
                && ( inDegreeOf(source) != 1
                     || outDegreeOf(incomingVerticesOf(source).iterator().next()) > 1 );
    }

    /**
     * Get all of the vertices in a linear chain of vertices starting at zipStart
     *
     * Build a list of vertices (in order) starting from zipStart such that each sequential pair of vertices
     * in the chain A and B can be zipped together.
     *
     * @param zipStart a vertex that starts a linear chain
     * @return a list of vertices that comprise a linear chain starting with zipStart.  The resulting
     *         list will always contain at least zipStart as the first element.
     */
    @Requires("isLinearChainStart(zipStart)")
    @Ensures({"result != null", "result.size() >= 1"})
    private LinkedList<SeqVertex> traceLinearChain(final SeqVertex zipStart) {
        final LinkedList<SeqVertex> linearChain = new LinkedList<SeqVertex>();
        linearChain.add(zipStart);

        boolean lastIsRef = isReferenceNode(zipStart); // remember because this calculation is expensive
        SeqVertex last = zipStart;
        while (true) {
            if ( outDegreeOf(last) != 1 )
                // cannot extend a chain from last if last has multiple outgoing branches
                break;

            // there can only be one (outgoing edge of last) by contract
            final SeqVertex target = getEdgeTarget(outgoingEdgeOf(last));

            if ( inDegreeOf(target) != 1 || last.equals(target) )
                // cannot zip up a target that has multiple incoming nodes or that's a cycle to the last node
                break;

            final boolean targetIsRef = isReferenceNode(target);
            if ( lastIsRef != targetIsRef ) // both our isRef states must be equal
                break;

            linearChain.add(target); // extend our chain by one

            // update our last state to be the current state, and continue
            last = target;
            lastIsRef = targetIsRef;
        }

        return linearChain;
    }

    /**
     * Merge a linear chain of vertices into a single combined vertex, and update this graph to such that
     * the incoming edges into the first element of the linearChain and the outgoing edges from linearChain.getLast()
     * all point to this new combined vertex.
     *
     * @param linearChain a non-empty chain of vertices that can be zipped up into a single vertex
     * @return true if we actually merged at least two vertices together
     */
    protected boolean mergeLinearChain(final LinkedList<SeqVertex> linearChain) {
        if ( linearChain.isEmpty() ) throw new IllegalArgumentException("BUG: cannot have linear chain with 0 elements but got " + linearChain);

        final SeqVertex first = linearChain.getFirst();
        final SeqVertex last = linearChain.getLast();

        if ( first == last ) return false; // only one element in the chain, cannot be extended

        // create the combined vertex, and add it to the graph
        // TODO -- performance problem -- can be optimized if we want

        final SeqVertex addedVertex = mergeLinearChainVertices(linearChain);
        addVertex(addedVertex);

        // update the incoming and outgoing edges to point to the new vertex
        for( final BaseEdge edge : outgoingEdgesOf(last) ) { addEdge(addedVertex, getEdgeTarget(edge), edge.copy()); }
        for( final BaseEdge edge : incomingEdgesOf(first) )  { addEdge(getEdgeSource(edge), addedVertex, edge.copy()); }

        removeAllVertices(linearChain);
        return true;
    }

    protected SeqVertex mergeLinearChainVertices(final List<SeqVertex> vertices) {
        final List<byte[]> seqs = new LinkedList<byte[]>();
        for ( SeqVertex v : vertices ) seqs.add(v.getSequence());
        final byte[] seqsCat = org.broadinstitute.gatk.utils.Utils.concat(seqs.toArray(new byte[][]{}));
        return new SeqVertex( seqsCat );
    }

    /**
     * Base class for transformation operations that need to iterate over proposed vertices, where
     * each proposed vertex is a seed vertex for a potential transformation.
     *
     * transformUntilComplete will iteratively apply the tryToTransform function on each vertex in the graph
     * until no vertex can be found that can be transformed.
     *
     * Note that in order to eventually terminate tryToTransform must transform the graph such that eventually
     * no vertices are candidates for further transformations.
     */
    private abstract class VertexBasedTransformer {
        /**
         * For testing purposes we sometimes want to test that can be transformed capabilities are working
         * without actually modifying the graph */
        private boolean dontModifyGraphEvenIfPossible = false;

        public boolean dontModifyGraphEvenIfPossible() { return dontModifyGraphEvenIfPossible; }
        public void setDontModifyGraphEvenIfPossible() { this.dontModifyGraphEvenIfPossible = true; }

        /**
         * Merge until the graph has no vertices that are candidates for merging
         */
        public boolean transformUntilComplete() {
            boolean didAtLeastOneTransform = false;
            boolean foundNodesToMerge = true;
            while( foundNodesToMerge ) {
                foundNodesToMerge = false;

                for( final SeqVertex v : vertexSet() ) {
                    foundNodesToMerge = tryToTransform(v);
                    if ( foundNodesToMerge ) {
                        didAtLeastOneTransform = true;
                        break;
                    }
                }
            }

            return didAtLeastOneTransform;
        }

        /**
         * Merge, if possible, seeded on the vertex v
         * @param v the proposed seed vertex to merge
         * @return true if some useful merging happened, false otherwise
         */
        abstract boolean tryToTransform(final SeqVertex v);
    }

    /**
     * Merge diamond configurations:
     *
     * Performance the transformation:
     *
     * { A -> x + S_i + y -> Z }
     *
     * goes to:
     *
     * { A -> x -> S_i -> y -> Z }
     *
     * for all nodes that match this configuration.
     */
    protected class MergeDiamonds extends VertexBasedTransformer {
        @Override
        protected boolean tryToTransform(final SeqVertex top) {
            final Set<SeqVertex> middles = outgoingVerticesOf(top);
            if ( middles.size() <= 1 )
                // we can only merge if there's at least two middle nodes
                return false;

            SeqVertex bottom = null;
            for ( final SeqVertex mi : middles ) {
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
            if ( inDegreeOf(bottom) != middles.size() )
                return false;

            if ( dontModifyGraphEvenIfPossible() ) return true;

            // actually do the merging, returning true if at least 1 base was successfully split
            final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(SeqGraph.this, middles);
            if (splitter.meetsMinMergableSequenceForEitherPrefixOrSuffix(1))
                return splitter.splitAndUpdate(top, bottom);
            else
                return false;
        }
    }

    /**
     * Merge tail configurations:
     *
     * Performs the transformation:
     *
     * { A -> x + S_i + y }
     *
     * goes to:
     *
     * { A -> x -> S_i -> y }
     *
     * for all nodes that match this configuration.
     *
     * Differs from the diamond transform in that no bottom node is required
     */
    protected class MergeTails extends VertexBasedTransformer {
        @Override
        protected boolean tryToTransform(final SeqVertex top) {
            final Set<SeqVertex> tails = outgoingVerticesOf(top);
            if ( tails.size() <= 1 )
                return false;

            for ( final SeqVertex t : tails )
                if ( ! isSink(t) || inDegreeOf(t) > 1 )
                    return false;

            if ( dontModifyGraphEvenIfPossible() ) return true;

            final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(SeqGraph.this, tails);

            if (splitter.meetsMinMergableSequenceForSuffix(MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES))
                return splitter.splitAndUpdate(top, null);
            else
                return false;
        }
    }

    /**
     * Merge headless configurations:
     *
     * Performs the transformation:
     *
     * { x + S_i -> y -> Z }
     *
     * goes to:
     *
     * { x -> S_i -> y + Z }
     *
     * for all nodes that match this configuration.
     */
    protected class MergeCommonSuffices extends VertexBasedTransformer {
        @Override
        boolean tryToTransform(final SeqVertex bottom) {
            return new SharedSequenceMerger().merge(SeqGraph.this, bottom);
        }
    }

    /**
     * Performs the transformation:
     *
     * { x + S_i + y -> Z }
     *
     * goes to:
     *
     * { x -> S_i -> y -> Z }
     *
     * for all nodes that match this configuration.
     *
     * Differs from the diamond transform in that no top node is required
     */
    protected class SplitCommonSuffices extends VertexBasedTransformer {
        final Set<SeqVertex> alreadySplit = new HashSet<SeqVertex>();

        @Override
        boolean tryToTransform(final SeqVertex bottom) {
            if ( alreadySplit.contains(bottom) )
                return false;
            else {
                alreadySplit.add(bottom);
                return new CommonSuffixSplitter().split(SeqGraph.this, bottom);
            }
        }
    }
}
