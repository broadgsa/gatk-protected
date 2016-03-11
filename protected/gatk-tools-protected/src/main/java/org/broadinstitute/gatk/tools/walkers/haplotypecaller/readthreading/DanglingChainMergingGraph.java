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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading;

import com.google.java.contract.Ensures;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.smithwaterman.*;
import org.jgrapht.EdgeFactory;

import java.util.*;

public abstract class DanglingChainMergingGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> {

    private static final int MAX_CIGAR_COMPLEXITY = 3;
    private int maxMismatchesInDanglingHead = -1;

    protected boolean alreadyBuilt;

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    protected DanglingChainMergingGraph(final int kmerSize, final EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> edgeFactory) {
        super(kmerSize, edgeFactory);
    }

    protected void setMaxMismatchesInDanglingHead(final int maxMismatchesInDanglingHead) {
        this.maxMismatchesInDanglingHead = maxMismatchesInDanglingHead;
    }

    /**
     * Edge factory that encapsulates the numPruningSamples assembly parameter
     */
    protected static class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
        final int numPruningSamples;

        public MyEdgeFactory(int numPruningSamples) {
            this.numPruningSamples = numPruningSamples;
        }

        @Override
        public MultiSampleEdge createEdge(final MultiDeBruijnVertex sourceVertex, final MultiDeBruijnVertex targetVertex) {
            return new MultiSampleEdge(false, 1, numPruningSamples);
        }

        public MultiSampleEdge createEdge(final boolean isRef, final int multiplicity) {
            return new MultiSampleEdge(isRef, multiplicity, numPruningSamples);
        }

    }

    /**
     * Class to keep track of the important dangling chain merging data
     */
    protected static final class DanglingChainMergeHelper {
        final List<MultiDeBruijnVertex> danglingPath, referencePath;
        final byte[] danglingPathString, referencePathString;
        final Cigar cigar;

        public DanglingChainMergeHelper(final List<MultiDeBruijnVertex> danglingPath,
                                        final List<MultiDeBruijnVertex> referencePath,
                                        final byte[] danglingPathString,
                                        final byte[] referencePathString,
                                        final Cigar cigar) {
            this.danglingPath = danglingPath;
            this.referencePath = referencePath;
            this.danglingPathString = danglingPathString;
            this.referencePathString = referencePathString;
            this.cigar = cigar;
        }
    }

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingTails requires the graph be already built");

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( outDegreeOf(v) == 0 && ! isRefSink(v) ) {
                attempted++;
                nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength);
            }
        }

        logger.debug("Recovered " + nRecovered + " of " + attempted + " dangling tails");
    }

    /**
     * Try to recover dangling heads
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     */
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingHeads requires the graph be already built");

        // we need to build a list of dangling heads because that process can modify the graph (and otherwise generate
        // a ConcurrentModificationException if we do it while iterating over the vertexes)
        final List<MultiDeBruijnVertex> danglingHeads = new ArrayList<>();

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( inDegreeOf(v) == 0 && ! isRefSource(v) )
                danglingHeads.add(v);
        }

        // now we can try to recover the dangling heads
        for ( final MultiDeBruijnVertex v : danglingHeads ) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength);
        }

        logger.debug("Recovered " + nRecovered + " of " + attempted + " dangling heads");
    }

    /**
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    protected int recoverDanglingTail(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( outDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling tail for " + vertex + " but it has out-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingTailMergeResult == null || ! cigarIsOkayToMerge(danglingTailMergeResult.cigar, false, true) )
            return 0;

        // merge
        return mergeDanglingTail(danglingTailMergeResult);
    }

    /**
     * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
     * @return 1 if we successfully recovered a vertex and 0 otherwise
     */
    protected int recoverDanglingHead(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        if ( inDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling head for " + vertex + " but it has in-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingChainMergeHelper danglingHeadMergeResult = generateCigarAgainstUpwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingHeadMergeResult == null || ! cigarIsOkayToMerge(danglingHeadMergeResult.cigar, true, false) )
            return 0;

        // merge
        return mergeDanglingHead(danglingHeadMergeResult);
    }

    /**
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar    the cigar to analyze
     * @param requireFirstElementM if true, require that the first cigar element be an M operator in order for it to be okay
     * @param requireLastElementM  if true, require that the last cigar element be an M operator in order for it to be okay
     * @return true if it's okay to merge, false otherwise
     */
    protected boolean cigarIsOkayToMerge(final Cigar cigar, final boolean requireFirstElementM, final boolean requireLastElementM) {

        final List<CigarElement> elements = cigar.getCigarElements();
        final int numElements = elements.size();

        // don't allow more than a couple of different ops
        if ( numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY )
            return false;

        // the first element must be an M
        if ( requireFirstElementM && elements.get(0).getOperator() != CigarOperator.M )
            return false;

        // the last element must be an M
        if ( requireLastElementM && elements.get(numElements - 1).getOperator() != CigarOperator.M )
            return false;

        // note that there are checks for too many mismatches in the dangling branch later in the process

        return true;
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult   the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    protected int mergeDanglingTail(final DanglingChainMergeHelper danglingTailMergeResult) {

        final List<CigarElement> elements = danglingTailMergeResult.cigar.getCigarElements();
        final CigarElement lastElement = elements.get(elements.size() - 1);
        if ( lastElement.getOperator() != CigarOperator.M )
            throw new IllegalArgumentException("The last Cigar element must be an M");

        final int lastRefIndex = danglingTailMergeResult.cigar.getReferenceLength() - 1;
        final int matchingSuffix = Math.min(GraphUtils.longestSuffixMatch(danglingTailMergeResult.referencePathString, danglingTailMergeResult.danglingPathString, lastRefIndex), lastElement.getLength());
        if ( matchingSuffix == 0 )
            return 0;

        final int altIndexToMerge = Math.max(danglingTailMergeResult.cigar.getReadLength() - matchingSuffix - 1, 0);

        // there is an important edge condition that we need to handle here: Smith-Waterman correctly calculates that there is a
        // deletion, that deletion is left-aligned such that the LCA node is part of that deletion, and the rest of the dangling
        // tail is a perfect match to the suffix of the reference path.  In this case we need to push the reference index to merge
        // down one position so that we don't incorrectly cut a base off of the deletion.
        final boolean firstElementIsDeletion = elements.get(0).getOperator() == CigarOperator.D;
        final boolean mustHandleLeadingDeletionCase =  firstElementIsDeletion && (elements.get(0).getLength() + matchingSuffix == lastRefIndex + 1);
        final int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

        // another edge condition occurs here: if Smith-Waterman places the whole tail into an insertion then it will try to
        // merge back to the LCA, which results in a cycle in the graph.  So we do not want to merge in such a case.
        if ( refIndexToMerge == 0 )
            return 0;

        // it's safe to merge now
        addEdge(danglingTailMergeResult.danglingPath.get(altIndexToMerge), danglingTailMergeResult.referencePath.get(refIndexToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Actually merge the dangling head if possible
     *
     * @param danglingHeadMergeResult   the result from generating a Cigar for the dangling head against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    protected int mergeDanglingHead(final DanglingChainMergeHelper danglingHeadMergeResult) {

        final List<CigarElement> elements = danglingHeadMergeResult.cigar.getCigarElements();
        final CigarElement firstElement = elements.get(0);
        if ( firstElement.getOperator() != CigarOperator.M )
            throw new IllegalArgumentException("The first Cigar element must be an M");

        final int indexesToMerge = bestPrefixMatch(danglingHeadMergeResult.referencePathString, danglingHeadMergeResult.danglingPathString, firstElement.getLength());
        if ( indexesToMerge <= 0 )
            return 0;

        // we can't push back the reference path
        if ( indexesToMerge >= danglingHeadMergeResult.referencePath.size() - 1 )
            return 0;

        // but we can manipulate the dangling path if we need to
        if ( indexesToMerge >= danglingHeadMergeResult.danglingPath.size() &&
                ! extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge - danglingHeadMergeResult.danglingPath.size() + 2) )
            return 0;

        addEdge(danglingHeadMergeResult.referencePath.get(indexesToMerge+1), danglingHeadMergeResult.danglingPath.get(indexesToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));

        return 1;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param vertex   the sink of the dangling chain
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    protected DanglingChainMergeHelper generateCigarAgainstDownwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {
        final int minTailPathLength = Math.max(1, minDanglingBranchLength); // while heads can be 0, tails absolutely cannot

        // find the lowest common ancestor path between this vertex and the diverging master path if available
        final List<MultiDeBruijnVertex> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor);
        if ( altPath == null || isRefSource(altPath.get(0)) || altPath.size() < minTailPathLength + 1 ) // add 1 to include the LCA
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.downwards, Arrays.asList(incomingEdgeOf(altPath.get(1))));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, false);
        final byte[] altBases = getBasesForPath(altPath, false);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the source) and the reference path.
     *
     * @param vertex   the source of the dangling head
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    protected DanglingChainMergeHelper generateCigarAgainstUpwardsReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor, final int minDanglingBranchLength) {

        // find the highest common descendant path between vertex and the reference source if available
        final List<MultiDeBruijnVertex> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor);
        if ( altPath == null || isRefSink(altPath.get(0)) || altPath.size() < minDanglingBranchLength + 1 ) // add 1 to include the LCA
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0), TraversalDirection.upwards, Collections.<MultiSampleEdge>emptyList());

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath, true);
        final byte[] altBases = getBasesForPath(altPath, true);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingChainMergeHelper(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Finds the path upwards in the graph from this vertex to the first diverging node, including that (lowest common ancestor) vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto another path or
     *  has an ancestor with multiple incoming edges before hitting the reference path
     */
    protected List<MultiDeBruijnVertex> findPathUpwardsToLowestCommonAncestor(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( inDegreeOf(v) == 1 && outDegreeOf(v) < 2 ) {
            final MultiSampleEdge edge = incomingEdgeOf(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor )
                path.clear();
            // otherwise it is safe to use
            else
                path.addFirst(v);
            v = getEdgeSource(edge);
        }
        path.addFirst(v);

        return outDegreeOf(v) > 1 ? path : null;
    }

    /**
     * Finds the path downwards in the graph from this vertex to the reference sequence, including the highest common descendant vertex.
     * However note that the path is reversed so that this vertex ends up at the end of the path.
     * Also note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     *  has a descendant with multiple outgoing edges before hitting the reference path
     */
    protected List<MultiDeBruijnVertex> findPathDownwardsToHighestCommonDescendantOfReference(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( ! isReferenceNode(v) && outDegreeOf(v) == 1 ) {
            final MultiSampleEdge edge = outgoingEdgeOf(v);
            // if it has too low a weight, don't use it (or previous vertexes) for the path
            if ( edge.getPruningMultiplicity() < pruneFactor )
                path.clear();
                // otherwise it is safe to use
            else
                path.addFirst(v);
            v = getEdgeTarget(edge);
        }
        path.addFirst(v);

        return isReferenceNode(v) ? path : null;
    }

    private enum TraversalDirection {
        downwards,
        upwards
    }

    /**
     * Finds the path in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start   the reference vertex to start from
     * @param direction describes which direction to move in the graph (i.e. down to the reference sink or up to the source)
     * @param blacklistedEdges edges to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the path (non-null, non-empty)
     */
    protected List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start,
                                                         final TraversalDirection direction,
                                                         final Collection<MultiSampleEdge> blacklistedEdges) {

        final List<MultiDeBruijnVertex> path = new ArrayList<>();

        MultiDeBruijnVertex v = start;
        while ( v != null ) {
            path.add(v);
            v = (direction == TraversalDirection.downwards ? getNextReferenceVertex(v, true, blacklistedEdges) : getPrevReferenceVertex(v));
        }

        return path;
    }

    /**
     * The base sequence for the given path.
     *
     * @param path the list of vertexes that make up the path
     * @param expandSource if true and if we encounter a source node, then expand (and reverse) the character sequence for that node
     * @return  non-null sequence of bases corresponding to the given path
     */
    @Ensures({"result != null"})
    public byte[] getBasesForPath(final List<MultiDeBruijnVertex> path, final boolean expandSource) {
        if ( path == null ) throw new IllegalArgumentException("Path cannot be null");

        final StringBuilder sb = new StringBuilder();
        for ( final MultiDeBruijnVertex v : path ) {
            if ( expandSource && isSource(v) ) {
                final String seq = v.getSequenceString();
                sb.append(new StringBuilder(seq).reverse().toString());
            } else {
                sb.append((char)v.getSuffix());
            }
        }

        return sb.toString().getBytes();
    }

    /**
     * Finds the index of the best extent of the prefix match between the provided paths, for dangling head merging.
     * Assumes that path1.length >= maxIndex and path2.length >= maxIndex.
     *
     * @param path1  the first path
     * @param path2  the second path
     * @param maxIndex the maximum index to traverse (not inclusive)
     * @return the index of the ideal prefix match or -1 if it cannot find one, must be less than maxIndex
     */
    protected int bestPrefixMatch(final byte[] path1, final byte[] path2, final int maxIndex) {
        final int maxMismatches = getMaxMismatches(maxIndex);
        int mismatches = 0;
        int index = 0;
        int lastGoodIndex = -1;
        while ( index < maxIndex ) {
            if ( path1[index] != path2[index] ) {
                if ( ++mismatches > maxMismatches )
                    return -1;
                lastGoodIndex = index;
            }
            index++;
        }
        // if we got here then we hit the max index
        return lastGoodIndex;
    }

    /**
     * Determine the maximum number of mismatches permitted on the branch.
     * Unless it's preset (e.g. by unit tests) it should be the length of the branch divided by the kmer size.
     *
     * @param lengthOfDanglingBranch  the length of the branch itself
     * @return positive integer
     */
    private int getMaxMismatches(final int lengthOfDanglingBranch) {
        return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : Math.max(1, (lengthOfDanglingBranch / kmerSize));
    }

    protected boolean extendDanglingPathAgainstReference(final DanglingChainMergeHelper danglingHeadMergeResult, final int numNodesToExtend) {

        final int indexOfLastDanglingNode = danglingHeadMergeResult.danglingPath.size() - 1;
        final int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
        if ( indexOfRefNodeToUse >= danglingHeadMergeResult.referencePath.size() )
            return false;

        final MultiDeBruijnVertex danglingSource = danglingHeadMergeResult.danglingPath.remove(indexOfLastDanglingNode);
        final StringBuilder sb = new StringBuilder();
        final byte[] refSourceSequence = danglingHeadMergeResult.referencePath.get(indexOfRefNodeToUse).getSequence();
        for ( int i = 0; i < numNodesToExtend; i++ )
            sb.append((char)refSourceSequence[i]);
        sb.append(danglingSource.getSequenceString());
        final byte[] sequenceToExtend = sb.toString().getBytes();

        // clean up the source and edge
        final MultiSampleEdge sourceEdge = outgoingEdgeOf(danglingSource);
        MultiDeBruijnVertex prevV = getEdgeTarget(sourceEdge);
        removeEdge(danglingSource, prevV);

        // extend the path
        for ( int i = numNodesToExtend; i > 0; i-- ) {
            final MultiDeBruijnVertex newV = new MultiDeBruijnVertex(Arrays.copyOfRange(sequenceToExtend, i, i+kmerSize));
            addVertex(newV);
            final MultiSampleEdge newE = addEdge(newV, prevV);
            newE.setMultiplicity(sourceEdge.getMultiplicity());
            danglingHeadMergeResult.danglingPath.add(newV);
            prevV = newV;
        }

        return true;
    }
}