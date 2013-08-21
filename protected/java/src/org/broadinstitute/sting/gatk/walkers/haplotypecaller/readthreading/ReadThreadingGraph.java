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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller.readthreading;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.KMerCounter;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.Kmer;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.smithwaterman.SWParameterSet;
import org.broadinstitute.sting.utils.smithwaterman.SmithWaterman;
import org.jgrapht.EdgeFactory;
import org.jgrapht.alg.CycleDetector;

import java.io.File;
import java.util.*;

public class ReadThreadingGraph extends BaseGraph<MultiDeBruijnVertex, MultiSampleEdge> {
    /**
     * Edge factory that encapsulates the numPruningSamples assembly parameter
     */
    private static class MyEdgeFactory implements EdgeFactory<MultiDeBruijnVertex, MultiSampleEdge> {
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

    private final static Logger logger = Logger.getLogger(ReadThreadingGraph.class);

    private final static String ANONYMOUS_SAMPLE = "XXX_UNNAMED_XXX";
    private final static boolean WRITE_GRAPH = false;
    private final static boolean DEBUG_NON_UNIQUE_CALC = false;

    private final static int MAX_CIGAR_COMPLEXITY = 3;
    private final static int MIN_DANGLING_TAIL_LENGTH = 5;  // SNP + 3 stabilizing nodes + the LCA

    /** for debugging info printing */
    private static int counter = 0;

    /**
     * Sequences added for read threading before we've actually built the graph
     */
    private final Map<String, List<SequenceForKmers>> pending = new LinkedHashMap<>();

    /**
     * A set of non-unique kmers that cannot be used as merge points in the graph
     */
    private Set<Kmer> nonUniqueKmers;

    /**
     * A map from kmers -> their corresponding vertex in the graph
     */
    private Map<Kmer, MultiDeBruijnVertex> uniqueKmers = new LinkedHashMap<>();

    /**
     *
     */
    final int kmerSize;
    final boolean debugGraphTransformations;
    final byte minBaseQualityToUseInAssembly;

    protected boolean increaseCountsBackwards = true;
    protected boolean increaseCountsThroughBranches = false; // this may increase the branches without bounds

    // --------------------------------------------------------------------------------
    // state variables, initialized in resetToInitialState()
    // --------------------------------------------------------------------------------
    private Kmer refSource;
    private boolean alreadyBuilt;

    public ReadThreadingGraph() {
        this(25, false, (byte)6, 1);
    }

    public ReadThreadingGraph(final int kmerSize) {
        this(kmerSize, false, (byte)6, 1);
    }

    /**
     * Create a new ReadThreadingAssembler using kmerSize for matching
     * @param kmerSize must be >= 1
     */
    protected ReadThreadingGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly, final int numPruningSamples) {
        super(kmerSize, new MyEdgeFactory(numPruningSamples));

        if ( kmerSize < 1 ) throw new IllegalArgumentException("bad minkKmerSize " + kmerSize);
        this.kmerSize = kmerSize;
        this.debugGraphTransformations = debugGraphTransformations;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;

        resetToInitialState();
    }

    /**
     * Reset this assembler to its initial state, so we can create another assembly with a different set of reads
     */
    private void resetToInitialState() {
        pending.clear();
        nonUniqueKmers = null;
        uniqueKmers.clear();
        refSource = null;
        alreadyBuilt = false;
    }

    /**
     * Add the all bases in sequence to the graph
     * @param sequence a non-null sequence
     * @param isRef is this the reference sequence?
     */
    protected void addSequence(final byte[] sequence, final boolean isRef) {
        addSequence("anonymous", sequence, null, isRef);
    }

    /**
     * Add all bases in sequence to this graph
     *
     * @see #addSequence(String, String, byte[], int, int, int[], boolean) for full information
     */
    public void addSequence(final String seqName, final byte[] sequence, final int[] counts, final boolean isRef) {
        addSequence(seqName, ANONYMOUS_SAMPLE, sequence, 0, sequence.length, counts, isRef);
    }

    /**
     * Add bases in sequence to this graph
     *
     * @param seqName a useful seqName for this read, for debugging purposes
     * @param sequence non-null sequence of bases
     * @param counts a vector of counts for each bases, indicating how many times that base was observed in the sequence.
     *               This allows us to support reduced reads in the ReadThreadingAssembler.  Can be null, meaning that
     *               each base is only observed once.  If not null, must have length == sequence.length.
     * @param start the first base offset in sequence that we should use for constructing the graph using this sequence, inclusive
     * @param stop the last base offset in sequence that we should use for constructing the graph using this sequence, exclusive
     * @param isRef is this the reference sequence.
     */
    public void addSequence(final String seqName, final String sampleName, final byte[] sequence, final int start, final int stop, final int[] counts, final boolean isRef) {
        // note that argument testing is taken care of in SequenceForKmers
        if ( alreadyBuilt ) throw new IllegalStateException("Graph already built");

        // get the list of sequences for this sample
        List<SequenceForKmers> sampleSequences = pending.get(sampleName);
        if ( sampleSequences == null ) { // need to create
            sampleSequences = new LinkedList<>();
            pending.put(sampleName, sampleSequences);
        }

        // add the new sequence to the list of sequences for sample
        sampleSequences.add(new SequenceForKmers(seqName, sequence, start, stop, counts, isRef));
    }

    /**
     * Return a count appropriate for a kmer starting at kmerStart in sequence for kmers
     *
     * @param seqForKmers a non-null sequence for kmers object
     * @param kmerStart the position where the kmer starts in sequence
     * @return a count for a kmer from start -> start + kmerSize in seqForKmers
     */
    private int getCountGivenKmerStart(final SequenceForKmers seqForKmers, final int kmerStart) {
        return seqForKmers.getCount(kmerStart + kmerSize - 1);
    }

    /**
     * Thread sequence seqForKmers through the current graph, updating the graph as appropriate
     * @param seqForKmers a non-null sequence
     */
    private void threadSequence(final SequenceForKmers seqForKmers) {
        final Pair<MultiDeBruijnVertex,Integer> startingInfo = findStart(seqForKmers);
        if ( startingInfo == null )
            return;

        final MultiDeBruijnVertex startingVertex = startingInfo.getFirst();
        final int uniqueStartPos = startingInfo.getSecond();

        // increase the counts of all edges incoming into the starting vertex supported by going back in sequence
        if ( increaseCountsBackwards )
            increaseCountsInMatchedKmers(seqForKmers, startingVertex, startingVertex.getSequence(), kmerSize - 2);

        if ( debugGraphTransformations ) startingVertex.addRead(seqForKmers.name);

        // keep track of information about the reference source
        if ( seqForKmers.isRef ) {
            if ( refSource != null ) throw new IllegalStateException("Found two refSources! prev: " + refSource + ", new: " + startingVertex);
            refSource = new Kmer(seqForKmers.sequence, seqForKmers.start, kmerSize);
        }

        // loop over all of the bases in sequence, extending the graph by one base at each point, as appropriate
        MultiDeBruijnVertex vertex = startingVertex;
        for ( int i = uniqueStartPos + 1; i <= seqForKmers.stop - kmerSize; i++ ) {
            final int count = getCountGivenKmerStart(seqForKmers, i);

            vertex = extendChainByOne(vertex, seqForKmers.sequence, i, count, seqForKmers.isRef);
            if ( debugGraphTransformations ) vertex.addRead(seqForKmers.name);
        }
    }

    /**
     * Class to keep track of the important dangling tail merging data
     */
    protected final class DanglingTailMergeResult {
        final List<MultiDeBruijnVertex> danglingPath, referencePath;
        final byte[] danglingPathString, referencePathString;
        final Cigar cigar;

        public DanglingTailMergeResult(final List<MultiDeBruijnVertex> danglingPath,
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
     * Attempt to attach vertex with out-degree == 0 to the graph
     *
     * @param vertex the vertex to recover
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return 1 if we successfully recovered the vertex and 0 otherwise
     */
    protected int recoverDanglingChain(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        if ( outDegreeOf(vertex) != 0 ) throw new IllegalStateException("Attempting to recover a dangling tail for " + vertex + " but it has out-degree > 0");

        // generate the CIGAR string from Smith-Waterman between the dangling tail and reference paths
        final DanglingTailMergeResult danglingTailMergeResult = generateCigarAgainstReferencePath(vertex, pruneFactor);

        // if the CIGAR is too complex (or couldn't be computed) then we do not allow the merge into the reference path
        if ( danglingTailMergeResult == null || ! cigarIsOkayToMerge(danglingTailMergeResult.cigar) )
            return 0;

        // merge
        return mergeDanglingTail(danglingTailMergeResult);
    }

    /**
     * Determine whether the provided cigar is okay to merge into the reference path
     *
     * @param cigar    the cigar to analyze
     * @return true if it's okay to merge, false otherwise
     */
    protected boolean cigarIsOkayToMerge(final Cigar cigar) {

        final List<CigarElement> elements = cigar.getCigarElements();
        final int numElements = elements.size();

        // don't allow more than a couple of different ops
        if ( numElements > MAX_CIGAR_COMPLEXITY )
            return false;

        // the last element must be an M
        if ( elements.get(numElements - 1).getOperator() != CigarOperator.M )
            return false;

        // TODO -- do we want to check whether the Ms mismatch too much also?

        return true;
    }

    /**
     * Actually merge the dangling tail if possible
     *
     * @param danglingTailMergeResult   the result from generating a Cigar for the dangling tail against the reference
     * @return 1 if merge was successful, 0 otherwise
     */
    protected int mergeDanglingTail(final DanglingTailMergeResult danglingTailMergeResult) {

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

        addEdge(danglingTailMergeResult.danglingPath.get(altIndexToMerge), danglingTailMergeResult.referencePath.get(refIndexToMerge), ((MyEdgeFactory)getEdgeFactory()).createEdge(false, 1));
        return 1;
    }

    /**
     * Generates the CIGAR string from the Smith-Waterman alignment of the dangling path (where the
     * provided vertex is the sink) and the reference path.
     *
     * @param vertex   the sink of the dangling tail
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return a SmithWaterman object which can be null if no proper alignment could be generated
     */
    protected DanglingTailMergeResult generateCigarAgainstReferencePath(final MultiDeBruijnVertex vertex, final int pruneFactor) {

        // find the lowest common ancestor path between vertex and the reference sink if available
        final List<MultiDeBruijnVertex> altPath = findPathToLowestCommonAncestorOfReference(vertex, pruneFactor);
        if ( altPath == null || isRefSource(altPath.get(0)) || altPath.size() < MIN_DANGLING_TAIL_LENGTH )
            return null;

        // now get the reference path from the LCA
        final List<MultiDeBruijnVertex> refPath = getReferencePath(altPath.get(0));

        // create the Smith-Waterman strings to use
        final byte[] refBases = getBasesForPath(refPath);
        final byte[] altBases = getBasesForPath(altPath);

        // run Smith-Waterman to determine the best alignment (and remove trailing deletions since they aren't interesting)
        final SmithWaterman alignment = new SWPairwiseAlignment(refBases, altBases, SWParameterSet.STANDARD_NGS, SWPairwiseAlignment.OVERHANG_STRATEGY.LEADING_INDEL);
        return new DanglingTailMergeResult(altPath, refPath, altBases, refBases, AlignmentUtils.removeTrailingDeletions(alignment.getCigar()));
    }

    /**
     * Finds the path upwards in the graph from this vertex to the reference sequence, including the lowest common ancestor vertex.
     * Note that nodes are excluded if their pruning weight is less than the pruning factor.
     *
     * @param vertex   the original vertex
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     * @return the path if it can be determined or null if this vertex either doesn't merge onto the reference path or
     *  has an ancestor with multiple incoming edges before hitting the reference path
     */
    protected List<MultiDeBruijnVertex> findPathToLowestCommonAncestorOfReference(final MultiDeBruijnVertex vertex, final int pruneFactor) {
        final LinkedList<MultiDeBruijnVertex> path = new LinkedList<>();

        MultiDeBruijnVertex v = vertex;
        while ( ! isReferenceNode(v) && inDegreeOf(v) == 1 ) {
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

        return isReferenceNode(v) ? path : null;
    }

    /**
     * Finds the path downwards in the graph from this vertex to the reference sink, including this vertex
     *
     * @param start   the reference vertex to start from
     * @return the path (non-null, non-empty)
     */
    protected List<MultiDeBruijnVertex> getReferencePath(final MultiDeBruijnVertex start) {
        if ( ! isReferenceNode(start) ) throw new IllegalArgumentException("Cannot construct the reference path from a vertex that is not on that path");

        final List<MultiDeBruijnVertex> path = new ArrayList<>();

        MultiDeBruijnVertex v = start;
        while ( v != null ) {
            path.add(v);
            v = getNextReferenceVertex(v);
        }

        return path;
    }

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    public void buildGraphIfNecessary() {
        if ( alreadyBuilt ) return;

        // determine the kmer size we'll use, and capture the set of nonUniques for that kmer size
        final NonUniqueResult result = determineKmerSizeAndNonUniques(kmerSize, kmerSize);
        nonUniqueKmers = result.nonUniques;

        if ( DEBUG_NON_UNIQUE_CALC ) {
            logger.info("using " + kmerSize + " kmer size for this assembly with the following non-uniques");
        }

        // go through the pending sequences, and add them to the graph
        for ( final List<SequenceForKmers> sequencesForSample : pending.values() ) {
            for ( final SequenceForKmers sequenceForKmers : sequencesForSample ) {
                threadSequence(sequenceForKmers);
                if ( WRITE_GRAPH ) printGraph(new File("threading." + counter++ + "." + sequenceForKmers.name.replace(" ", "_") + ".dot"), 0);
            }

            // flush the single sample edge values from the graph
            for ( final MultiSampleEdge e : edgeSet() ) e.flushSingleSampleMultiplicity();
        }

        // clear
        pending.clear();
        alreadyBuilt = true;
    }

    /**
     * @return true if the graph has cycles, false otherwise
     */
    public boolean hasCycles() {
        return new CycleDetector<>(this).detectCycles();
    }

    /**
     * Does the graph not have enough complexity?  We define low complexity as a situation where the number
     * of non-unique kmers is more than 20% of the total number of kmers.
     *
     * @return true if the graph has low complexity, false otherwise
     */
    public boolean isLowComplexity() {
        return nonUniqueKmers.size() * 4 > uniqueKmers.size();
    }

    /**
     * Try to recover dangling tails
     *
     * @param pruneFactor  the prune factor to use in ignoring chain pieces
     */
    public void recoverDanglingTails(final int pruneFactor) {
        if ( ! alreadyBuilt )  throw new IllegalStateException("recoverDanglingTails requires the graph be already built");

        int attempted = 0;
        int nRecovered = 0;
        for ( final MultiDeBruijnVertex v : vertexSet() ) {
            if ( outDegreeOf(v) == 0 && ! isRefNodeAndRefSink(v) ) {
                attempted++;
                nRecovered += recoverDanglingChain(v, pruneFactor);
            }
        }

        if ( debugGraphTransformations ) logger.info("Recovered " + nRecovered + " of " + attempted + " dangling tails");
    }

    /** structure that keeps track of the non-unique kmers for a given kmer size */
    private static class NonUniqueResult {
        final Set<Kmer> nonUniques;
        final int kmerSize;

        private NonUniqueResult(Set<Kmer> nonUniques, int kmerSize) {
            this.nonUniques = nonUniques;
            this.kmerSize = kmerSize;
        }
    }

    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param minKmerSize the minimum kmer size to consider when constructing the graph
     * @param maxKmerSize the maximum kmer size to consider
     * @return a non-null NonUniqueResult
     */
    protected NonUniqueResult determineKmerSizeAndNonUniques(final int minKmerSize, final int maxKmerSize) {
        final Collection<SequenceForKmers> withNonUniques = getAllPendingSequences();
        final Set<Kmer> nonUniqueKmers = new HashSet<Kmer>();

        // go through the sequences and determine which kmers aren't unique within each read
        int kmerSize = minKmerSize;
        for ( ; kmerSize <= maxKmerSize; kmerSize++) {
            // clear out set of non-unique kmers
            nonUniqueKmers.clear();

            // loop over all sequences that have non-unique kmers in them from the previous iterator
            final Iterator<SequenceForKmers> it = withNonUniques.iterator();
            while ( it.hasNext() ) {
                final SequenceForKmers sequenceForKmers = it.next();

                // determine the non-unique kmers for this sequence
                final Collection<Kmer> nonUniquesFromSeq = determineNonUniqueKmers(sequenceForKmers, kmerSize);
                if ( nonUniquesFromSeq.isEmpty() ) {
                    // remove this sequence from future consideration
                    it.remove();
                } else {
                    // keep track of the non-uniques for this kmerSize, and keep it in the list of sequences that have non-uniques
                    nonUniqueKmers.addAll(nonUniquesFromSeq);
                }
            }

            if ( nonUniqueKmers.isEmpty() )
                // this kmerSize produces no non-unique sequences, so go ahead and use it for our assembly
                break;
        }

        // necessary because the loop breaks with kmerSize = max + 1
        return new NonUniqueResult(nonUniqueKmers, Math.min(kmerSize, maxKmerSize));
    }

    /**
     * Get the collection of all sequences for kmers across all samples in no particular order
     * @return non-null Collection
     */
    private Collection<SequenceForKmers> getAllPendingSequences() {
        final LinkedList<SequenceForKmers> result = new LinkedList<SequenceForKmers>();
        for ( final List<SequenceForKmers> oneSampleWorth : pending.values() ) result.addAll(oneSampleWorth);
        return result;
    }

    /**
     * Get the collection of non-unique kmers from sequence for kmer size kmerSize
     * @param seqForKmers a sequence to get kmers from
     * @param kmerSize the size of the kmers
     * @return a non-null collection of non-unique kmers in sequence
     */
    private Collection<Kmer> determineNonUniqueKmers(final SequenceForKmers seqForKmers, final int kmerSize) {
        // count up occurrences of kmers within each read
        final KMerCounter counter = new KMerCounter(kmerSize);
        final int stopPosition = seqForKmers.stop - kmerSize;
        for ( int i = 0; i <= stopPosition; i++ ) {
            final Kmer kmer = new Kmer(seqForKmers.sequence, i, kmerSize);
            counter.addKmer(kmer, 1);
        }

        return counter.getKmersWithCountsAtLeast(2);
    }

    /**
     * Convert this kmer graph to a simple sequence graph.
     *
     * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
     * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
     *
     * @return a newly allocated SequenceGraph
     */
    // TODO -- should override base class method
    public SeqGraph convertToSequenceGraph() {
        buildGraphIfNecessary();

        final SeqGraph seqGraph = new SeqGraph(kmerSize);
        final Map<MultiDeBruijnVertex, SeqVertex> vertexMap = new HashMap<MultiDeBruijnVertex, SeqVertex>();

        // create all of the equivalent seq graph vertices
        for ( final MultiDeBruijnVertex dv : vertexSet() ) {
            final SeqVertex sv = new SeqVertex(dv.getAdditionalSequence(isSource(dv)));
            sv.setAdditionalInfo(dv.additionalInfo());
            vertexMap.put(dv, sv);
            seqGraph.addVertex(sv);
        }

        // walk through the nodes and connect them to their equivalent seq vertices
        for( final MultiSampleEdge e : edgeSet() ) {
            final SeqVertex seqInV = vertexMap.get(getEdgeSource(e));
            final SeqVertex seqOutV = vertexMap.get(getEdgeTarget(e));
            //logger.info("Adding edge " + seqInV + " -> " + seqOutV);
            seqGraph.addEdge(seqInV, seqOutV, new BaseEdge(e.isRef(), e.getMultiplicity()));
        }

        return seqGraph;
    }

    private void increaseCountsInMatchedKmers(final SequenceForKmers seqForKmers,
                                              final MultiDeBruijnVertex vertex,
                                              final byte[] originalKmer,
                                              final int offset) {
        if ( offset == -1 ) return;

        for ( final MultiSampleEdge edge : incomingEdgesOf(vertex) ) {
            final MultiDeBruijnVertex prev = getEdgeSource(edge);
            final byte suffix = prev.getSuffix();
            final byte seqBase = originalKmer[offset];
//            logger.warn(String.format("Increasing counts for %s -> %s via %s at %d with suffix %s vs. %s",
//                    prev, vertex, edge, offset, (char)suffix, (char)seqBase));
            if ( suffix == seqBase && (increaseCountsThroughBranches || inDegreeOf(vertex) == 1) ) {
                edge.incMultiplicity(seqForKmers.getCount(offset));
                increaseCountsInMatchedKmers(seqForKmers, prev, originalKmer, offset-1);
            }
        }
    }

    /**
     * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
     *
     * @param seqForKmers the sequence we want to thread into the graph
     * @return a pair of the starting vertex and its position in seqForKmer
     */
    private Pair<MultiDeBruijnVertex, Integer> findStart(final SequenceForKmers seqForKmers) {
        final int uniqueStartPos = seqForKmers.isRef ? 0 : findUniqueStartPosition(seqForKmers.sequence, seqForKmers.start, seqForKmers.stop);

        if ( uniqueStartPos == -1 )
            return null;

        return getOrCreateKmerVertex(seqForKmers.sequence, uniqueStartPos, true);
    }

    /**
     * Find a starting point in sequence that begins a unique kmer among all kmers in the graph
     * @param sequence the sequence of bases
     * @param start the first base to use in sequence
     * @param stop the last base to use in sequence
     * @return the index into sequence that begins a unique kmer of size kmerSize, or -1 if none could be found
     */
    private int findUniqueStartPosition(final byte[] sequence, final int start, final int stop) {
        for ( int i = start; i < stop - kmerSize; i++ ) {
            final Kmer kmer1 = new Kmer(sequence, i, kmerSize);
            if ( uniqueKmers.containsKey(kmer1) )
                return i;
        }
        return -1;
    }

    /**
     * Get the vertex for the kmer in sequence starting at start
     * @param sequence the sequence
     * @param start the position of the kmer start
     * @param allowRefSource if true, we will allow matches to the kmer that represents the reference starting kmer
     * @return a non-null vertex
     */
    private Pair<MultiDeBruijnVertex, Integer> getOrCreateKmerVertex(final byte[] sequence, final int start, final boolean allowRefSource) {
        final Kmer kmer = new Kmer(sequence, start, kmerSize);
        final MultiDeBruijnVertex vertex = getUniqueKmerVertex(kmer, allowRefSource);
        if ( vertex != null ) {
            return new Pair<>(vertex, start);
        } else {
            return new Pair<>(createVertex(kmer), start);
        }
    }

    /**
     * Get the unique vertex for kmer, or null if not possible.
     *
     * @param allowRefSource if true, we will allow kmer to match the reference source vertex
     * @return a vertex for kmer, or null if it's not unique
     */
    private MultiDeBruijnVertex getUniqueKmerVertex(final Kmer kmer, final boolean allowRefSource) {
        if ( ! allowRefSource && kmer.equals(refSource) ) return null;
        return uniqueKmers.get(kmer);
    }

    /**
     * Create a new vertex for kmer.  Add it to the uniqueKmers map if appropriate.
     *
     * kmer must not have a entry in unique kmers, or an error will be thrown
     *
     * @param kmer the kmer we want to create a vertex for
     * @return the non-null created vertex
     */
    private MultiDeBruijnVertex createVertex(final Kmer kmer) {
        final MultiDeBruijnVertex newVertex = new MultiDeBruijnVertex(kmer.bases());
        final int prevSize = vertexSet().size();
        addVertex(newVertex);

        // make sure we aren't adding duplicates (would be a bug)
        if ( vertexSet().size() != prevSize + 1) throw new IllegalStateException("Adding vertex " + newVertex + " to graph didn't increase the graph size");

        // add the vertex to the unique kmer map, if it is in fact unique
        if ( ! nonUniqueKmers.contains(kmer) && ! uniqueKmers.containsKey(kmer) ) // TODO -- not sure this last test is necessary
            uniqueKmers.put(kmer, newVertex);

        return newVertex;
    }

    /**
     * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
     * the graph one bp according to the bases in sequence.
     *
     * @param prevVertex a non-null vertex where sequence was last anchored in the graph
     * @param sequence the sequence we're threading through the graph
     * @param kmerStart the start of the current kmer in graph we'd like to add
     * @param count the number of observations of this kmer in graph (can be > 1 for reduced reads)
     * @param isRef is this the reference sequence?
     * @return a non-null vertex connecting prevVertex to in the graph based on sequence
     */
    private MultiDeBruijnVertex extendChainByOne(final MultiDeBruijnVertex prevVertex, final byte[] sequence, final int kmerStart, final int count, final boolean isRef) {
        final Set<MultiSampleEdge> outgoingEdges = outgoingEdgesOf(prevVertex);

        final int nextPos = kmerStart + kmerSize - 1;
        for ( final MultiSampleEdge outgoingEdge : outgoingEdges ) {
            final MultiDeBruijnVertex target = getEdgeTarget(outgoingEdge);
            if ( target.getSuffix() == sequence[nextPos] ) {
                // we've got a match in the chain, so simply increase the count of the edge by 1 and continue
                outgoingEdge.incMultiplicity(count);
                return target;
            }
        }

        // none of our outgoing edges had our unique suffix base, so we check for an opportunity to merge back in
        final Kmer kmer = new Kmer(sequence, kmerStart, kmerSize);
        final MultiDeBruijnVertex uniqueMergeVertex = getUniqueKmerVertex(kmer, false);

        if ( isRef && uniqueMergeVertex != null )
            throw new IllegalStateException("Found a unique vertex to merge into the reference graph " + prevVertex + " -> " + uniqueMergeVertex);

        // either use our unique merge vertex, or create a new one in the chain
        final MultiDeBruijnVertex nextVertex = uniqueMergeVertex == null ? createVertex(kmer) : uniqueMergeVertex;
        addEdge(prevVertex, nextVertex, ((MyEdgeFactory)getEdgeFactory()).createEdge(isRef, count));
        return nextVertex;
    }

    /**
     * Add the given read to the sequence graph.  Ultimately the read will get sent through addSequence(), but first
     * this method ensures we only use high quality bases and accounts for reduced reads, etc.
     *
     * @param read a non-null read
     */
    protected void addRead(final GATKSAMRecord read) {
        final byte[] sequence = read.getReadBases();
        final byte[] qualities = read.getBaseQualities();
        final int[] reducedReadCounts = read.getReducedReadCounts();  // will be null if read is not reduced

        int lastGood = -1; // the index of the last good base we've seen
        for( int end = 0; end <= sequence.length; end++ ) {
            if ( end == sequence.length || ! baseIsUsableForAssembly(sequence[end], qualities[end]) ) {
                // the first good base is at lastGood, can be -1 if last base was bad
                final int start = lastGood;
                // the stop base is end - 1 (if we're not at the end of the sequence)
                final int len = end - start;

                if ( start != -1 && len >= kmerSize ) {
                    // if the sequence is long enough to get some value out of, add it to the graph
                    final String name = read.getReadName() + "_" + start + "_" + end;
                    addSequence(name, read.getReadGroup().getSample(), read.getReadBases(), start, end, reducedReadCounts, false);
                }

                lastGood = -1; // reset the last good base
            } else if ( lastGood == -1 ) {
                lastGood = end; // we're at a good base, the last good one is us
            }
        }
    }

    /**
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    protected boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }

    /**
     * Get the set of non-unique kmers in this graph.  For debugging purposes
     * @return a non-null set of kmers
     */
    protected Set<Kmer> getNonUniqueKmers() {
        return nonUniqueKmers;
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{" +
                "kmerSize=" + kmerSize +
                '}';
    }
}