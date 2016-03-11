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


import java.util.List;
import java.util.ListIterator;

/**
 * Represents a route or path through a graph.
 * <p>
 *     In contrast with a {@link Path}, a route keeps track of the
 * path taken at furcations in order to speed up some path comparisions like the
 * one implemented by {@link #isSuffix}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class Route<V extends BaseVertex, E extends BaseEdge> extends Path<V,E> {

    protected final Route<V,E> previousRouteWithLastVertexThatIsForkOrJoin;
    protected final boolean lastVertexIsForkOrJoin;

    /**
     * Create a zero length route with a start in a particular vertex:
     *
     * @param initialVertex the first vertex of the route.
     * @param graph the new route's graph.
     *
     * @throws IllegalArgumentException if {@code initialVertex} or {@code graph} are {@code null}.
     *      or if {@code initialVertex} does not belong to {@code graph}.
     */
    public Route(final V initialVertex, final BaseGraph<V, E> graph) {
        super(initialVertex, graph);
        previousRouteWithLastVertexThatIsForkOrJoin = null;
        lastVertexIsForkOrJoin = graph.inDegreeOf(initialVertex) > 1;
    }

    @Override
    public boolean equals(final Object other) {
        if (other == null) return false;
        if (other == this) return true;
        if (! (other instanceof Route)) return false;
        @SuppressWarnings("unchecked")
        final Route<V,E> otherRoute = (Route<V,E>) other;
        return otherRoute.length() == this.length() && isSuffix(otherRoute);
    }

    /**
     * Extends a route into a new instance.
     *
     * @param prefix the route to extend.
     * @param nextVertex the vertex to extend the route to.
     *
     * @throws IllegalArgumentException if {@code prefix} is {@code null} or {@code nextVertex} is {@code null}
     *   or {@code nextVertex} does not belong to {@code prefix}'s graph or there is no edge that in the graph
     *   that would connect {@code prefix}'s last vertex with {@code nextVertex} directly.
     */
    public Route(final Route<V,E> prefix, final V nextVertex) {
        this(prefix,resolveSuffixEdge(prefix,nextVertex));
    }


    /**
     * Extends a route into a new instance.
     *
     * @param prevVertex the vertex to extend the route to.
     * @param suffix the route to extend.
     *
     * @throws IllegalArgumentException if {@code suffix} is {@code null} or {@code prevVertex} is {@code null}
     *   or {@code prevVertex} does not belong to {@code suffix}'s graph or there is no edge that in the graph
     *   that would connect {@code suffix}'s first vertex with {@code prevVertex} directly.
     */
    public Route(final V prevVertex, final Route<V,E> suffix) {
        this(resolvePrefixEdge(prevVertex, suffix),suffix);
    }

    /**
     * Resolves the prefix edge as required by {@link Route(V,Route)}.
     */
    private static <V extends BaseVertex,E extends BaseEdge>  E resolvePrefixEdge(final V prevVertex, final Route<V, E> suffix) {
        if (prevVertex == null) throw new NullPointerException();
        if (!suffix.getGraph().containsVertex(prevVertex)) throw new IllegalArgumentException();
        final E result = suffix.getGraph().getEdge(prevVertex,suffix.getFirstVertex());
        if (result == null)
            throw new IllegalArgumentException("there is no such edge in the graph");
        return result;
    }

    /**
     * Resolves the suffix edge as required by {@link Route(Route,V)}
     */
    private static <V extends BaseVertex,E extends BaseEdge>  E resolveSuffixEdge(final Route<V,E> prefix, final V nextVertex) {
        if (nextVertex == null) throw new IllegalArgumentException();
        if (!prefix.getGraph().containsVertex(nextVertex)) throw new IllegalArgumentException();
        final E result = prefix.getGraph().getEdge(prefix.getLastVertex(),nextVertex);
        if (result == null)
            throw new IllegalArgumentException("there is no such edge in the graph");
        return result;
    }

    /**
     * Extends a route by prefixing an edge.
     *
     * @param initialEdge the extending edge.
     * @param suffix the original path.
     *
     * @throws IllegalArgumentException if {@code suffix} or {@code initialEdge} are {@code null}, or {@code initialEdge} is
     * not part of {@code suffix}'s graph, or {@code initialEdge} does not have as a target the first vertex in {@code suffix}.
     */
    public Route(final E initialEdge, final Route<V,E> suffix) {
        super(initialEdge,suffix);
        final V firstVertex = getFirstVertex();
        if(suffix.length() == 0) {
            lastVertexIsForkOrJoin = suffix.lastVertexIsForkOrJoin || graph.outDegreeOf(firstVertex) > 1;
            previousRouteWithLastVertexThatIsForkOrJoin = graph.inDegreeOf(firstVertex) > 1 ? new Route<>(firstVertex,graph) : null;
        } else {
            lastVertexIsForkOrJoin = suffix.lastVertexIsForkOrJoin;
            if (suffix.previousRouteWithLastVertexThatIsForkOrJoin != null)
                previousRouteWithLastVertexThatIsForkOrJoin = new Route<>(initialEdge,suffix.previousRouteWithLastVertexThatIsForkOrJoin);
            else
                previousRouteWithLastVertexThatIsForkOrJoin = graph.outDegreeOf(firstVertex) > 1 ?
                        new Route<>(new Route<>(firstVertex,graph),edgesInOrder.get(0)) :
                            graph.inDegreeOf(firstVertex) > 1 ? new Route<>(firstVertex,graph) : null;
        }
    }

    /**
     * Create copy of an existing route.
     * @param route the route to copy
     *
     * @throws NullPointerException if {@code route} is {@code null}.
     */
    protected Route(final Route<V, E> route) {
        super(route);
        lastVertexIsForkOrJoin = route.lastVertexIsForkOrJoin;
        previousRouteWithLastVertexThatIsForkOrJoin = route.previousRouteWithLastVertexThatIsForkOrJoin;
    }

    /**
     * Create a new Route extending another one with an edge
     *
     * @param route the route to extend.
     * @param edge the edge to extend the route with.
     *
     * @throws IllegalArgumentException if {@code route} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code route}'s graph, or {@code edge} does not have as a source the last vertex in {@code route}.
     */
    public Route(final Route<V, E> route, final E edge) {
        super(route, edge);
        lastVertexIsForkOrJoin = graph.outDegreeOf(route.lastVertex) > 1 || graph.inDegreeOf(lastVertex) > 1;
        previousRouteWithLastVertexThatIsForkOrJoin = route.lastVertexIsForkOrJoin ? route : route.previousRouteWithLastVertexThatIsForkOrJoin;
    }

    @Override
    public boolean isSuffix(final Path<V,E> other) {
        if (other == this)
            return true;
        else if (other == null)
            throw new IllegalArgumentException("other path must not be null");
        else if (getGraph() != other.getGraph())
            throw new IllegalArgumentException("other path must be part of the same graph");
        else if (other instanceof Route)
            return isRouteSuffix((Route<V,E>)other);
        else
            return super.isSuffix(other);
    }

    @Override
    public String toString() {
        return super.toString().replace("Path{", "Route{");
    }

    /**
     * Faster version when comparing with a route.
     */
    protected boolean isRouteSuffix(final Route<V,E> other) {
        if (other.getGraph() != this.getGraph())
            throw new IllegalArgumentException("you cannot compare routes on different graphs");
        else if (lastVertex != other.lastVertex)  // obvious case.
            return false;
        else if (this.previousRouteWithLastVertexThatIsForkOrJoin == null
                && other.previousRouteWithLastVertexThatIsForkOrJoin != null) // I am shorter or different path for sure.
            return false;
        else if (this.edgesInOrder.size() < other.edgesInOrder.size())  // I am shorter regardless of path, no way Jose!
            return false;
        else if (this.previousRouteWithLastVertexThatIsForkOrJoin == null || other.previousRouteWithLastVertexThatIsForkOrJoin == null) {
            final ListIterator<E> myEdges = edgesInOrder.listIterator(edgesInOrder.size());
            final ListIterator<E> otherEdges = other.edgesInOrder.listIterator(other.edgesInOrder.size());
            while (otherEdges.hasPrevious())
                if (myEdges.previous() != otherEdges.previous())
                    return false;
            return true;
        } else
            return (other.previousRouteWithLastVertexThatIsForkOrJoin == this.previousRouteWithLastVertexThatIsForkOrJoin)
                || (previousRouteWithLastVertexThatIsForkOrJoin.lastVertex == other.previousRouteWithLastVertexThatIsForkOrJoin.lastVertex
              && previousRouteWithLastVertexThatIsForkOrJoin.isRouteSuffix(other.previousRouteWithLastVertexThatIsForkOrJoin));
    }

    /**
     * Checks whether the last vertex in the route is a fork or a joining vertex.
     * @return {@code true} iff so.
     */
    public boolean lastVertexIsForkOrJoin() {
        return lastVertexIsForkOrJoin;
    }

    /**
     * Returns the longest prefix route that has as a last vertex a join or furcation vertex.
     *
     * @return never {@code null}.
     */
    public Route<V,E> getPrefixRouteWithLastVertexThatIsForkOrJoin() {
        return previousRouteWithLastVertexThatIsForkOrJoin;
    }



    /**
     * Splice out the first few vertices of the route.
     *
     * @param length how many vertices to splice out
     * @return a new route without those spliced vertices.
     *
     * @throws IllegalArgumentException if {@code length} is equal to the route's length or greater or if it is negative.
     * Notice that non-vertex route are no legal routes.
     */
    public Route<V,E> splicePrefix(final int length) {
        if (length == 0)
            return this;
        if (length >= length())
            throw new IllegalArgumentException("prefix slicing to long");
        if (length < 0)
            throw new IllegalArgumentException("prefix cannot be negative");

        final List<E> resultEdges = getEdges().subList(length,length());
        Route<V,E> result = new Route<>(graph.getEdgeSource(resultEdges.get(0)),graph);
        for (final E edge : resultEdges)
            result = new Route<>(result,edge);
        return result;
    }
}
