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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.gatk.utils.collections.CountSet;
import org.broadinstitute.gatk.utils.collections.CountSet;
import org.broadinstitute.gatk.utils.collections.Pair;

import java.util.*;

/**
 * Encapsulates the graph traversals needed to find event-blocks.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class EventBlockFinder {

    private final HaplotypeGraph graph;

    private final Map<Pair<MultiDeBruijnVertex,MultiDeBruijnVertex>,EventBlock> eventBlockCache;

    /**
     * Constructs a new engine.
     *
     * @param graph the base haplotype graph to iterate over.
     */
    public EventBlockFinder(final HaplotypeGraph graph) {
        if (graph == null) throw new NullPointerException();
        this.graph = graph;
        eventBlockCache = new HashMap<>(20);
    }

    /**
     * Create a new traversal object based on a read anchoring.
     * @param anchoring
     * @return never {@code null}.
     */
    public Traversal traversal(final ReadAnchoring anchoring) {
        if (anchoring == null) throw new NullPointerException();
        return new Traversal(anchoring);
    }


    public class Traversal implements Iterable<EventBlock> {

        private final ReadAnchoring anchoring;

        private EventBlock lastEventBlock;


        private Traversal(final ReadAnchoring anchoring) {
            this.anchoring = anchoring;
            lastEventBlock = findLastEventBlock(anchoring);
        }

        @Override
        public java.util.Iterator<EventBlock> iterator() {
            return lastEventBlock == null ? Collections.EMPTY_SET.iterator() : new Iterator();
        }

        private class Iterator implements java.util.Iterator<EventBlock> {

            private MultiDeBruijnVertex currentVertex;

            private Iterator() {
                currentVertex = anchoring.leftAnchorVertex;
            }

            @Override
            public boolean hasNext() {
                return currentVertex != null;
            }

            @Override
            public EventBlock next() {
                final EventBlock result;
                if (currentVertex == null)
                    throw new NoSuchElementException("going beyond last event block");
                else if (currentVertex == lastEventBlock.getSource()) {
                    result = lastEventBlock;
                    currentVertex = null;
                } else {
                    final EventBlock candidate = findEventBlock(anchoring,false,currentVertex,lastEventBlock.getSource());
                    if (candidate == null) {
                        result = findEventBlock(anchoring,false,currentVertex,anchoring.rightAnchorVertex);
                        currentVertex = null;
                    } else {
                        result = candidate;
                        currentVertex = candidate.getSink();
                    }
                }
                return result;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        }
    }

    /**
     * Finds the last event block.
     * <p>
     * It can do it forward or backwards.
     * </p>
     *
     * @param anchoring target read anchoring information.
     * @return {@code null} if there is no event block, depending on {@code backwards} before or after current
     */
    private EventBlock findLastEventBlock(
            final ReadAnchoring anchoring) {
        return findEventBlock(anchoring,true,anchoring.leftAnchorVertex,anchoring.rightAnchorVertex);
    }

    /**
     * Finds an event block forward or backwards along the reference route.
     * @param anchoring the read anchoring information.
     * @param backwards true if the block should be constructed from right to left.
     * @param leftVertex the left vertex
     * @param rightVertex the right vertex
     * @return {@code null} if there is no such a event block between these coordinates.
     */
    private EventBlock findEventBlock(
            final ReadAnchoring anchoring, final boolean backwards,
            final MultiDeBruijnVertex leftVertex, final MultiDeBruijnVertex rightVertex) {

        MultiDeBruijnVertex currentVertex = backwards ? rightVertex : leftVertex;
        boolean foundEvent = false;
        final CountSet pathSizes = new CountSet(10); // typically more than enough.
        pathSizes.setTo(0);

        // Map between reference vertices where there is some expected open alternative path rejoining and the
        // predicted length of paths rejoining at that point counting from the beginning of the block.
        final Map<MultiDeBruijnVertex, CountSet> expectedAlternativePathRejoins = new HashMap<>(4);

        // Keeps record of possible left-clipping veritces; those that are located before any event path furcation
        // has been found. The value indicates the blockLength at the time we traverse that node.
        final Deque<Pair<MultiDeBruijnVertex, Integer>> possibleClippingPoints = new LinkedList<>();

        // We keep the distance from the beggining of the block (leftVertex).
        int blockLength = 0;
        while (currentVertex != null) {
            int openingDegree = backwards ? graph.outDegreeOf(currentVertex) : graph.inDegreeOf(currentVertex);
            if (openingDegree > 1) {
                final CountSet joiningPathLengths = expectedAlternativePathRejoins.remove(currentVertex);
                if (joiningPathLengths != null)
                    pathSizes.addAll(joiningPathLengths);
            }
            final boolean isValidBlockEnd = isValidBlockEnd(anchoring, currentVertex, expectedAlternativePathRejoins);
            if (foundEvent && isValidBlockEnd) // !gotcha we found a valid block end.
                break;
            else if (!foundEvent && isValidBlockEnd) // if no event has been found yet, still is a good clipping point.
                possibleClippingPoints.addLast(new Pair<>(currentVertex, blockLength));

            // We reached the end:
            if (currentVertex == (backwards ? leftVertex : rightVertex))
                break;

            // process next vertices, the next one on the reference and also possible start of alternative paths,
            // updates traversal structures accordingly.
            currentVertex = advanceOnReferencePath(anchoring, backwards, currentVertex, pathSizes, expectedAlternativePathRejoins);
            foundEvent |= expectedAlternativePathRejoins.size() > 0;
            pathSizes.incAll(1);
            blockLength++;
        }

        // we have not found an event, thus there is no block to report:
        if (!foundEvent)
            return null;

        // We try to clip off as much as we can from the beginning of the block before any event, but at least
        // leaving enough block length to meet the shortest path unless all paths have the same size (SNPs only)
        final int maxClipping = pathSizes.size() <= 1 ? blockLength : pathSizes.min();
        MultiDeBruijnVertex clippingEnd = backwards ? anchoring.rightAnchorVertex : anchoring.leftAnchorVertex;
        while (!possibleClippingPoints.isEmpty()) {
            final Pair<MultiDeBruijnVertex, Integer> candidate = possibleClippingPoints.removeLast();
            if (candidate.getSecond() <= maxClipping) {
                clippingEnd = candidate.getFirst();
                break;
            }
        }

        return resolveEventBlock(backwards ? new Pair<>(currentVertex, clippingEnd) : new Pair<>(clippingEnd, currentVertex));
    }

    /**
     * Gets or constructs a event-block through the cache.
     * @param borders the source and sink vertex pair for the requested event block.
     * @return never {@code null}
     */
    @Requires("borders != null && border.getFirst() != null && border.getSecond() != null")
    private EventBlock resolveEventBlock(final Pair<MultiDeBruijnVertex,MultiDeBruijnVertex> borders) {
        EventBlock result = eventBlockCache.get(borders);
        if (result == null)
            eventBlockCache.put(borders,result = new EventBlock(graph, borders.getFirst(),borders.getSecond()));
        return result;
    }

    /**
     * Move on vertex along the reference path checking for the presence of new opening alternative paths.
     *
     * @param anchoring anchoring information on the targeted read.
     * @param backwards whether we are extending the block backwards or forwards.
     * @param currentVertex the current vertex.
     * @param pathSizes current block path sizes.
     * @param expectedAlternativePathRejoins information about location of vertices along the reference path where open alternative paths will rejoin.
     * @return the next current-vertex, never {@code null} unless there is a bug.
     */
    private MultiDeBruijnVertex advanceOnReferencePath(final ReadAnchoring anchoring, final boolean backwards, final MultiDeBruijnVertex currentVertex, final CountSet pathSizes, final Map<MultiDeBruijnVertex, CountSet> expectedAlternativePathRejoins) {
        final Set<MultiSampleEdge> nextEdges = backwards ? graph.incomingEdgesOf(currentVertex) : graph.outgoingEdgesOf(currentVertex);
        MultiDeBruijnVertex nextReferenceVertex = null;
        for (final MultiSampleEdge e : nextEdges) {
            final MultiDeBruijnVertex nextVertex = backwards ? graph.getEdgeSource(e) : graph.getEdgeTarget(e);
            if (e.isRef())
                nextReferenceVertex = nextVertex;
            else {
                final CountSet pathSizesPlusOne = pathSizes.clone();
                pathSizesPlusOne.incAll(1);
                graph.calculateRejoins(nextVertex, expectedAlternativePathRejoins, anchoring.referenceWithinAnchorsMap.keySet(), pathSizesPlusOne, true, backwards);
            }
        }
        return nextReferenceVertex;
    }

    /**
     * Check whether the current vertex is a valid block end.
     *
     * @param anchoring reads anchoring information necessary to make the evaluation.
     * @param currentVertex target potential block end
     * @param expectedAlternativePathRejoins traversal states regarding open alternative paths.
     *
     * @return {@code true} iff so.
     */
    private boolean isValidBlockEnd(final ReadAnchoring anchoring, final MultiDeBruijnVertex currentVertex, final Map<MultiDeBruijnVertex, CountSet> expectedAlternativePathRejoins) {
        final boolean isUniqueKmer = anchoring.uniqueKmerOffsets.containsKey(currentVertex);
        final boolean isAnchorable = graph.getAnchorableVertices().contains(currentVertex) && isUniqueKmer && expectedAlternativePathRejoins.size() == 0;
        return isUniqueKmer && isAnchorable;
    }
}
