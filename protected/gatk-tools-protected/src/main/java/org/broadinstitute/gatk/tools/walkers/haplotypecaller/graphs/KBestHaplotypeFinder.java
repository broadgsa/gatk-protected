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

import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.jgrapht.alg.CycleDetector;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Efficient algorithm to obtain the list of best haplotypes given the {@link SeqGraph instace}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class KBestHaplotypeFinder extends AbstractList<KBestHaplotype> implements Iterable<KBestHaplotype>  {

    /**
     * The search graph.
     */
    private final SeqGraph graph;

    /**
     * Map of sub-haplotype finder by their source vertex.
     */
    protected Map<SeqVertex,KBestSubHaplotypeFinder> finderByVertex;

    /**
     * Possible haplotype sink vertices.
     */
    protected Set<SeqVertex> sinks;

    /**
     * Possible haplotype source vertices.
     */
    protected Set<SeqVertex> sources;

    /**
     * The top finder.
     *
     * <p>If there is only a single source vertex, its finder is the top finder. However when there
     * is more than one possible source, we create a composite finder that alternates between individual source vertices
     * for their best haplotypes.</p>
     */
    private final KBestSubHaplotypeFinder topFinder;

    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the seq-graph to search.
     * @param source the source vertex for all haplotypes.
     * @param sink sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code source} or {@code sink} is {@code null} or</li>
     *     <li>either {@code source} or {@code sink} is not a vertex in {@code graph}.</li>
     * </ul>
     */
    public KBestHaplotypeFinder(final SeqGraph graph, final SeqVertex source, final SeqVertex sink) {
        this(graph,Collections.singleton(source),Collections.singleton(sink));
    }

    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the seq-graph to search.
     * @param sources source vertices for all haplotypes.
     * @param sinks sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code sources} or {@code sinks} is {@code null} or</li>
     *     <li>any of {@code sources}' or any {@code sinks}' member is not a vertex in {@code graph}.</li>
     * </ul>
     */
    public KBestHaplotypeFinder(final SeqGraph graph, final Set<SeqVertex> sources, final Set<SeqVertex> sinks) {
        if (graph == null) throw new IllegalArgumentException("graph cannot be null");
        if (sources == null) throw new IllegalArgumentException("source cannot be null");
        if (sinks == null) throw new IllegalArgumentException("sink cannot be null");
        if (!graph.containsAllVertices(sources)) throw new IllegalArgumentException("source does not belong to the graph");
        if (!graph.containsAllVertices(sinks)) throw new IllegalArgumentException("sink does not belong to the graph");

        //TODO dealing with cycles here due to a bug in some of the graph transformations that produces cycles.
        //TODO Once that is solve, the if-else below should be substituted by a throw if there is any cycles,
        //TODO just the line commented out below if you want to trade early-bug-fail for speed.
        //this.graph = graph;
        this.graph = new CycleDetector<>(graph).detectCycles() ? removeCycles(graph,sources,sinks) : graph;

        finderByVertex = new HashMap<>(this.graph.vertexSet().size());
        this.sinks = sinks;
        this.sources = sources;
        if (sinks.size() == 0 || sources.size() == 0)
            topFinder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
        else if (sources.size() == 1)
            topFinder = createVertexFinder(sources.iterator().next());
        else
            topFinder = createAggregatedFinder();
    }

    /**
     * Constructs a new best haplotype finder.
     * <p>
     *     It will consider all source and sink vertex when looking for haplotypes.
     * </p>
     *
     * @param graph the seq-graph to search for the best haplotypes.
     */
    public KBestHaplotypeFinder(SeqGraph graph) {
        this(graph,graph.getSources(),graph.getSinks());
    }

    /**
     * Creates an aggregated recursive finder to try all possible source vertices.
     *
     * @return never {@code null}.
     */
    private KBestSubHaplotypeFinder createAggregatedFinder() {
        final List<KBestSubHaplotypeFinder> sourceFinders = new ArrayList<>(sources.size());
        for (final SeqVertex source : sources)
            sourceFinders.add(createVertexFinder(source));
        return new AggregatedSubHaplotypeFinder(sourceFinders);
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     *
     * @param original graph to modify.
     * @param sources considered source vertices.
     * @param sinks considered sink vertices.
     * @return never {@code null}.
     */
    private static SeqGraph removeCycles(final SeqGraph original, final Set<SeqVertex> sources, final Set<SeqVertex> sinks) {
        final Set<BaseEdge> edgesToRemove = new HashSet<>(original.edgeSet().size());
        final Set<SeqVertex> vertexToRemove = new HashSet<>(original.vertexSet().size());

        boolean foundSomePath = false;
        for (final SeqVertex source : sources)
            foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove,
                    vertexToRemove, new HashSet<SeqVertex>(original.vertexSet().size())) | foundSomePath;

        if (!foundSomePath)
            throw new IllegalStateException("could not find any path from the source vertex to the sink vertex after removing cycles: "
                    + Arrays.toString(sources.toArray()) + " => " + Arrays.toString(sinks.toArray()));

        if (edgesToRemove.isEmpty() && vertexToRemove.isEmpty())
            throw new IllegalStateException("cannot find a way to remove the cycles");

        final SeqGraph result = (SeqGraph) original.clone();
        result.removeAllEdges(edgesToRemove);
        result.removeAllVertices(vertexToRemove);
        return result;
    }

    /**
     * Recursive call that looks for edges and vertices that need to be removed to get rid of cycles.
     *
     * @param graph the original graph.
     * @param currentVertex current search vertex.
     * @param sinks considered sink vertices.
     * @param edgesToRemove collection  of edges that need to be removed in order to get rid of cycles.
     * @param verticesToRemove collection of vertices that can be removed.
     * @param parentVertices collection of vertices that preceded the {@code currentVertex}; i.e. the it can be
     *                       reached from those vertices using edges existing in {@code graph}.
     *
     * @return {@code true} to indicate that the some sink vertex is reachable by {@code currentVertex},
     *  {@code false} otherwise.
     */
    private static boolean findGuiltyVerticesAndEdgesToRemoveCycles(final SeqGraph graph,
                                                                    final SeqVertex currentVertex,
                                                                    final Set<SeqVertex> sinks,
                                                                    final Set<BaseEdge> edgesToRemove,
                                                                    final Set<SeqVertex> verticesToRemove,
                                                                    final Set<SeqVertex> parentVertices) {
        if (sinks.contains(currentVertex)) return true;

        final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
        boolean reachesSink = false;
        parentVertices.add(currentVertex);

        for (final BaseEdge edge : outgoingEdges) {
            final SeqVertex child = graph.getEdgeTarget(edge);
            if (parentVertices.contains(child))
                edgesToRemove.add(edge);
            else {
              final boolean childReachSink = findGuiltyVerticesAndEdgesToRemoveCycles(graph, child, sinks,
                      edgesToRemove, verticesToRemove, parentVertices);
              reachesSink = reachesSink || childReachSink;
            }
        }
        parentVertices.remove(currentVertex);
        if (!reachesSink) verticesToRemove.add(currentVertex);
        return reachesSink;
    }

    @Override
    public KBestHaplotype get(final int index) {
        if (index < 0 || index >= size())
            throw new IndexOutOfBoundsException();
        return topFinder.getKBest(index);
    }

    @Override
    public Iterator<KBestHaplotype> iterator() {
        return new Iterator<KBestHaplotype>() {
            private int nextK = 0;
            private final int maxK = topFinder.getCount();


            @Override
            public boolean hasNext() {
                return nextK < maxK;
            }

            @Override
            public KBestHaplotype next() {
                if (nextK >= maxK) throw new NoSuchElementException();
                return topFinder.getKBest(nextK++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public int size() {
        return topFinder.getCount();
    }

    /**
     * Returns an iterator on the first k best haplotypes.
     * <p>
     *     It might return less than k haplotypes if the total number of possible haplotypes is smaller.
     * </p>
     *
     * @param k the maximum number of haplotypes to return.
     * @return never {@code null}, but perhaps a iterator that return no haplotype.
     */
    public Iterator<KBestHaplotype> iterator(final int k) {

        return new Iterator<KBestHaplotype>() {
            private int nextK = 0;
            private final int maxK = Math.min(size(), k);

            @Override
            public boolean hasNext() {
                return nextK < maxK;
            }

            @Override
            public KBestHaplotype next() {
                if (nextK >= maxK) throw new NoSuchElementException();
                return topFinder.getKBest(nextK++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    /**
     * Creates a finder from a vertex.
     *
     * @param vertex the source vertex for the finder.
     *
     * @return never {@code null}, perhaps a finder that returns no haplotypes though.
     */
    protected KBestSubHaplotypeFinder createVertexFinder(final SeqVertex vertex) {
        KBestSubHaplotypeFinder finder = finderByVertex.get(vertex);
        if (finder == null) {
            if (sinks.contains(vertex))
                finder = new EmptyPathHaplotypeFinderNode(graph,vertex);
            else {
                final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(vertex);
                if (outgoingEdges.isEmpty())
                    finder = DeadEndKBestSubHaplotypeFinder.INSTANCE;
                else {
                    final Map<BaseEdge,KBestSubHaplotypeFinder> undeadChildren = createChildrenFinders(outgoingEdges);
                    finder = undeadChildren.isEmpty() ? DeadEndKBestSubHaplotypeFinder.INSTANCE :
                            new RecursiveSubHaplotypeFinder(graph,vertex,undeadChildren);
                }
            }
            finderByVertex.put(vertex, finder);
        }
        return finder;
    }

    /**
     * Creates finder for target vertices of a collection of edges.
     * <p>
     *     This peculiar signature is convenient for when we want to create finders for the children of a vertex.
     * </p>
     *
     * @param baseEdges target collection of edges.
     *
     * @return never {@code null}, perhaps an empty map if there is no children with valid paths to any sink for this
     *  finder.
     */
    private Map<BaseEdge, KBestSubHaplotypeFinder> createChildrenFinders(final Set<BaseEdge> baseEdges) {
        final Map<BaseEdge,KBestSubHaplotypeFinder> result = new LinkedHashMap<>(baseEdges.size());
        for (final BaseEdge edge : baseEdges) {
            final KBestSubHaplotypeFinder targetFinder = createVertexFinder(graph.getEdgeTarget(edge));
            if (targetFinder.getCount() == 0) continue;
            result.put(edge, targetFinder);
        }
        return result;
    }

    /**
     * Print a DOT representation of search graph.
     *
     * @param out character stream printer where to print the DOT representation to.
     *
     * @throws IllegalArgumentException if {@code out} is {@code null}.
     */
    public void printDOT(final PrintWriter out) {
        if (out == null)
            throw new IllegalArgumentException("the out writer cannot be null");
        out.println("digraph {");
        out.println("    rankdir = LR");
        out.println("    node [shape=box, margin=0.01]");
        out.println("    subgraph cluster_dummy { style = invis; x [label=\"\",shape=none,margin=0] }");
        final StringBuilder referenceCluster = new StringBuilder(1000);

        referenceCluster.append("    subgraph cluster_ref {\n");
        referenceCluster.append("        node [penwidth=2]\n");
        for (final KBestSubHaplotypeFinder finder : finderByVertex.values() ) {
            final String id = finder.id();
            final String line = String.format("    %s [label=<%s>]",id,finder.label());
            if (finder.isReference())
                referenceCluster.append("    ").append(line).append('\n');
            else
                out.println(line);
        }
        referenceCluster.append("    }");
        out.println(referenceCluster.toString());

        for (final KBestSubHaplotypeFinder finder : finderByVertex.values())
            for (final Pair<? extends KBestSubHaplotypeFinder,String> subFinderLabel : finder.subFinderLabels()) {
                final KBestSubHaplotypeFinder subFinder = subFinderLabel.getFirst();

                final String edgeLabel = subFinderLabel.getSecond();
                out.println(String.format("    %s -> %s [label=%s]",finder.id(),subFinder.id(),edgeLabel));
            }
        out.println("}");
    }

    /**
     * Print a DOT representation of search graph.
     *
     * @param file file where to print the DOT representation to.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws FileNotFoundException if {@code file} cannot be created or written.
     * @throws IllegalStateException if there was some trouble when writing the DOT representation.
     */
    public void printDOT(final File file) throws FileNotFoundException {
        if (file == null)
            throw new IllegalArgumentException("the output file cannot be null");
        final PrintWriter out = new PrintWriter(file);
        printDOT(out);
        if (out.checkError())
            throw new IllegalStateException("error occurred while writing k-best haplotype search graph into file '"
                    + file.getAbsolutePath() + "'");
        out.close();
    }

    /**
     * Print a DOT representation of search graph.
     *
     * @param fileName name of the file where to print the DOT representation to.
     *
     * @throws IllegalArgumentException if {@code fileName} is {@code null}.
     * @throws FileNotFoundException if no file named {@code fileName} cannot be created or written.
     * @throws IllegalStateException if there was some trouble when writing the DOT representation.
     */
    @SuppressWarnings("unused") // Available for debugging purposes.
    public void printDOTFile(final String fileName) throws FileNotFoundException {
        printDOT(new File(fileName));
    }

    /**
     * Get the score of a give sequence of bases
     *
     * @param bases the base sequence.
     *
     * @return {@link Double#NaN} if there is no score for the sequence, i.e. there is no such a haplotype accessible
     *   throw this finder.
     */
    public double score(final byte[] bases) {
        return topFinder.score(bases,0,bases.length);
    }

    /**
     * Get the score of a give sequence of bases
     *
     * @param haplotype the haplotype.
     *
     * @return {@link Double#NaN} if there is no score for the sequence, i.e. there is no such a haplotype accessible
     *   throw this finder.
     */
    public double score(final Haplotype haplotype) {
        return score(haplotype.getBases());
    }


    /**
     * Returns a unique list of haplotypes solutions.
     * <p>
     *     The result will not contain more than one haplotype with the same base sequence. The solution of the best
     *     score is returned.
     * </p>
     * <p>
     *     This makes sense when there are more than one possible path through the graph to create the same haplotype.
     * </p>
     * <p>
     *     The resulting list is sorted by the score with more likely haplotype search results first.
     * </p>
     *
     *
     * @return never {@code null}, perhaps an empty list.
     */
    public List<KBestHaplotype> unique() {
        final int requiredCapacity = size();
        final Set<Haplotype> haplotypes = new HashSet<>(requiredCapacity);
        final List<KBestHaplotype> result = new ArrayList<>(requiredCapacity);
        for (final KBestHaplotype kbh : this) {
            if (haplotypes.add(kbh.haplotype())) {
                result.add(kbh);
            }
        }
        return result;
    }
}
