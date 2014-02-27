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

import java.util.ArrayList;
import java.util.Map;
import java.util.PriorityQueue;

/**
* General recursive sub-haplotype finder.
* <p>
*   Provides the k-best sub-haplotypes looking into the outgoing set of vertices (that contain at least one solution).
* </p>
* <p>
*  This is done efficiently by keeping an priority-queue on best subhaplotype solutions and pulling them on demand
*  as needed.
*  </p>
* <p>
*  Solutions are cached for repeated retrieval so that we save compute at vertices that share sub-haplotypes
*     (share descendant vertices). This aspect is controlled by {@link KBestSubHaplotypeFinder} that instantiate
*     a unique {@link KBestSubHaplotypeFinder} for each vertex in the graph that belongs to a valid path
*     between the source and sink node.
* </p>
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
class RecursiveSubHaplotypeFinder implements KBestSubHaplotypeFinder {

    private final SeqGraph graph;

    private final SeqVertex vertex;

    private final Map<BaseEdge,KBestSubHaplotypeFinder> children;

    private boolean childrenWereProcessed = false;

    /**
     * Holds the number of possible paths from this source node vertex to the sink vertex.
     *
     * <p>Updated by {@link #processChildrenIfNeeded()}</p>
     */
    private int possibleHaplotypeCount;

    /**
     * Holds the best {@code i} paths to the sink so far calculated where {@code i+1} is the length of rankedResults.
     *
     * <p>As more results are requested the array will grow. All positions and solutions are calculated up to {@code i}</p>.
     */
    private ArrayList<KBestHaplotype> rankedResults;

    /**
     * Priority queue with best sub-haplotype solutions that haven't been calculated and cached on {@link #rankedResults} yet.
     */
    private PriorityQueue<ChildKBestSubHaplotype> nextChildrenKBestHaplotypePath;

    /**
     * Creates a recursive sub-haplotype finder give the target graph, first vertex and all possible outgoing edges
     *  with the corresponding sub-sub-haplotype finders.
     *
     * <p>For efficiency shake, it will not verify the content of {@code children} map; i.e. that indeed all keys
     * are outgoing edges from {@code vertex} on {@code graph} and that the value sub-haplotype resolver have as
     * the first vertex the adjacent vertex through that key edge.</p>
     *
     * @param graph the search graph.
     * @param vertex first vertex for all sub-haplotype solutions provided by this finder
     * @param children map from outgoing edge to the corresponding sub-sub-haplotype finder.
     */
    public RecursiveSubHaplotypeFinder(final SeqGraph graph, final SeqVertex vertex,
                                       final Map<BaseEdge, KBestSubHaplotypeFinder> children) {
        if (vertex == null) throw new IllegalArgumentException("the vertex provided cannot be null");
        if (graph == null) throw new IllegalArgumentException("the graph provided cannot be null");
        this.vertex = vertex;
        this.children = children;
        this.graph = graph;
    }

    @Override
    public int getCount() {
        processChildrenIfNeeded();
        return possibleHaplotypeCount;
    }

    /**
     * Process children and initialize structures if not done before.
     */
    private void processChildrenIfNeeded() {
        if (childrenWereProcessed) return;
        long possibleHaplotypeCount = 0;

        nextChildrenKBestHaplotypePath = new PriorityQueue<>(children.size());

        for (final Map.Entry<BaseEdge,KBestSubHaplotypeFinder> entry : children.entrySet()) {
            final KBestSubHaplotypeFinder child = entry.getValue();
            final BaseEdge edge = entry.getKey();
            final int childPossibleHaplotypePathCount = child.getCount();
            if (childPossibleHaplotypePathCount != 0) // paranoia check, should not happen at this point.
                nextChildrenKBestHaplotypePath.add(new ChildKBestSubHaplotype(-1,edge,child,0));
            possibleHaplotypeCount += childPossibleHaplotypePathCount;
        }

        // Just make sure we won't incur in overflow here for very large graphs; who is ever going to ask for more than 2G paths!!!)
        this.possibleHaplotypeCount = (int) Math.min(Integer.MAX_VALUE,possibleHaplotypeCount);

        // 10 is a bit arbitrary as it is difficult to anticipate what would be the number of requested
        // best sub-haplotypes for any node. It shouldn't be too large so that it does not waste space
        // but not too small so that there is no need to resize when just a few best solutions are requested.
        rankedResults = new ArrayList<>(Math.min(this.possibleHaplotypeCount,10));

        childrenWereProcessed = true;
    }

    @Override
    public KBestHaplotype getKBest(int k) {
        if (k < 0)
            throw new IllegalArgumentException("the rank requested cannot be negative");
        processChildrenIfNeeded();
        if (k >= possibleHaplotypeCount)
            throw new IllegalArgumentException("the rank requested cannot be equal or greater to the number of possible haplotypes");
        if (rankedResults.size() > k)
            return rankedResults.get(k);

        rankedResults.ensureCapacity(k+1);
        for (int i = rankedResults.size(); i <= k; i++) {
            // since k < possibleHaplotypeCount is guarantee no to be empty.
            if (nextChildrenKBestHaplotypePath.isEmpty())
                throw new IllegalStateException("what the heck " + k + " " + possibleHaplotypeCount);
            final ChildKBestSubHaplotype nextResult = nextChildrenKBestHaplotypePath.remove();
            nextResult.rank = i;
            rankedResults.add(nextResult);
            final int childRank = nextResult.subpath.rank();
            final KBestSubHaplotypeFinder child = nextResult.child;

            // if there is no further solution from the same child we cannot add another solution from that child.
            if (childRank + 1 >= nextResult.child.getCount())
                continue;
            nextChildrenKBestHaplotypePath.add(new ChildKBestSubHaplotype(-1,nextResult.edge, child, childRank + 1));
        }
        return rankedResults.get(k);
    }

    /**
     * Custom extension of the {@link KBestHaplotype} used for solutions generated by this class.
     */
    private class ChildKBestSubHaplotype extends KBestHaplotype implements Comparable<ChildKBestSubHaplotype>{
        private final int score;
        private int rank;
        private final KBestSubHaplotypeFinder child;
        private final BaseEdge edge;
        private final KBestHaplotype subpath;
        private final boolean isReference;

        public ChildKBestSubHaplotype(final int rank, final BaseEdge edge,
                                      final KBestSubHaplotypeFinder child, final int childRank) {
            this.child = child;
            this.edge = edge;
            this.rank = rank;
            this.subpath = child.getKBest(childRank);
            this.score = edge.getMultiplicity() + subpath.score();
            this.isReference = edge.isRef() && subpath.isReference();
        }

        @Override
        public SeqGraph graph() {
            return graph;
        }

        @Override
        public int compareTo(final ChildKBestSubHaplotype other) {
            if (other == null) throw new IllegalArgumentException("the other object cannot be null");
            return - Integer.compare(this.score,other.score);
        }

        @Override
        public int score() {
            return score;
        }

        @Override
        public int rank() {
            return rank;
        }


        @Override
        protected SeqVertex head() {
            return vertex;
        }

        @Override
        protected KBestHaplotype tail() {
            return subpath;
        }

        @Override
        public boolean isReference() {
            return isReference;
        }
    }
}
