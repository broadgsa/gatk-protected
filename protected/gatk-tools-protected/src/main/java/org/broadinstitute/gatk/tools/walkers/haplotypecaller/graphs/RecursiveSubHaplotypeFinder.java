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

import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.Pair;

import java.util.*;

/**
* General recursive sub-haplotype finder.
* <p>
*   Provides the k-best sub-haplotypes from a vertex provided map between outgoing edges and its target finders
* </p>
* <p>
*  This is done efficiently by keeping an priority-queue on best sub-haplotype solutions and pulling them on demand
*  as needed.
* </p>
* <p>
*  Solutions are cached for repeated retrieval so that we save compute at vertices that share sub-haplotypes
*     (share descendant vertices). This aspect is controlled by {@link KBestSubHaplotypeFinder} that instantiate
*     a unique {@link KBestSubHaplotypeFinder} for each vertex in the graph that belongs to a valid path
*     between the source and sink node.
* </p>
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
class RecursiveSubHaplotypeFinder extends AggregatedSubHaplotypeFinder<RecursiveSubHaplotypeFinder.EdgeSubHaplotypeFinder> {


    private final SeqVertex vertex;
    private final boolean isReference;

    /**
     * Creates a recursive sub-haplotype finder give the target graph, first vertex and all possible outgoing edges
     *  with the corresponding sub-sub-haplotype finders.
     *
     * <p>For efficiency shake, it will not verify the content of {@code children} map; i.e. that indeed all keys
     * are outgoing edges from {@code vertex} on {@code graph} and that the value sub-haplotype resolver have as
     * the first vertex the adjacent vertex through that key edge.</p>
     *
     * @param vertex first vertex for all sub-haplotype solutions provided by this finder
     * @param children map from outgoing edge to the corresponding sub-sub-haplotype finder.
     */
    public RecursiveSubHaplotypeFinder(final SeqGraph graph, final SeqVertex vertex,
                                       final Map<BaseEdge, KBestSubHaplotypeFinder> children) {
        super(createChildFinderCollection(vertex, children));
        this.vertex = vertex;
        this.isReference = graph.isReferenceNode(vertex);
    }

    /**
     * Wraps the descendant vertices finders in order to take advantage of the {@link AggregatedSubHaplotypeFinder}
     * common code.
     * <p>
     * Automatically calibrates the edge score (cost) so that it takes into account the total across all edges.
     * </p>
     *
     * @param vertex the parent vertex.
     * @param finders the child vertices indexed by the connecting edge.
     * @return never {@code null} but potentially an empty collection if there is child returning some sub-haplotype
     *    solution.
     */
    private static Collection<EdgeSubHaplotypeFinder> createChildFinderCollection(final SeqVertex vertex,
                                                             final Map<BaseEdge,KBestSubHaplotypeFinder> finders) {
        if (finders == null) throw new IllegalArgumentException("the edge to child map cannot be null");
        final ArrayList<EdgeSubHaplotypeFinder> result = new ArrayList<>(finders.size());
        for (final Map.Entry<BaseEdge,KBestSubHaplotypeFinder> e : finders.entrySet()) {
            final EdgeSubHaplotypeFinder subFinder = new EdgeSubHaplotypeFinder(vertex,e.getKey(), e.getValue());
            if (subFinder.getCount() == 0) continue;
            result.add(subFinder);
        }
        if (result.size() == 0)
            return Collections.emptySet();
        else if (result.size() == 1) // no calibration needed, by default edgeScore is 0.
            return Collections.singleton(result.get(0));
        else {
            double totalEdgeMultiplicityAcrossEdges = 0;
            for (final EdgeSubHaplotypeFinder finder : result)
                totalEdgeMultiplicityAcrossEdges += Math.max(0.5,finder.edge.getMultiplicity());
            final double log10TotalEdgeMultiplicityAcrossEdges = Math.log10(totalEdgeMultiplicityAcrossEdges);
            for (final EdgeSubHaplotypeFinder finder : result)
                finder.calibrateEdgeScore(log10TotalEdgeMultiplicityAcrossEdges);
            return result;
        }
    }

    @Override
    public boolean isReference() {
        return isReference;
    }

    @Override
    public String label() {
        return vertex.getSequenceString();
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        final Set<Pair<? extends KBestSubHaplotypeFinder,String>> result = new LinkedHashSet<>(subFinders.size());
        for (final EdgeSubHaplotypeFinder subFinder : subFinders)
            result.add(new Pair<>(subFinder,simplifyZeros(String.format("%.4f", subFinder.edgeScore))));
        return result;
    }

    /**
     * Removes zeros decimal positions from edge-labels.
     *
     * @param edgeLabel the original label to reformat.
     * @return never {@code null}, the reformatted label.
     */
    private String simplifyZeros(final String edgeLabel) {
        if (edgeLabel.equals("0.000") || edgeLabel.equals("-0.000") )
            return "0.";
        int i = edgeLabel.length() - 1;
        while (edgeLabel.charAt(i) == '0')
            i--;
        return (i == edgeLabel.length() - 1) ? edgeLabel : edgeLabel.substring(0,i);
    }

    protected static class EdgeSubHaplotypeFinder implements KBestSubHaplotypeFinder {

        private final KBestSubHaplotypeFinder childFinder;

        private final SeqVertex vertex;

        private final BaseEdge edge;

        private double edgeScore = 0;

        private EdgeSubHaplotypeFinder(final SeqVertex vertex, final BaseEdge edge, final KBestSubHaplotypeFinder childFinder) {
            this.childFinder = childFinder;
            this.edge = edge;
            this.vertex = vertex;
            this.edgeScore = 0;
        }

        private void calibrateEdgeScore(final double log10TotalMultiplicityAcrossOutgoingEdges) {
            edgeScore = Math.log10(Math.max(edge.getMultiplicity(),0.5)) - log10TotalMultiplicityAcrossOutgoingEdges;
        }

        @Override
        public String id() {
            return childFinder.id();
        }

        @Override
        public String label() {
            return childFinder.label();
        }

        @Override
        public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
            return childFinder.subFinderLabels();
        }

        @Override
        public int getCount() {
            return childFinder.getCount();
        }

        @Override
        public KBestHaplotype getKBest(int k) {
            return new ChildKBestSubHaplotype(vertex,edge,childFinder.getKBest(k),edgeScore);
        }

        @Override
        public boolean isReference() {
            return childFinder.isReference();
        }

        @Override
        public double score(final byte[] bases, final int offset, final int length) {
            if (length == 0)
                return 0;
            final byte[] vertexSequence = vertex.getSequence();
            if (length < vertexSequence.length) // query is not long enough to have any score.
                return Double.NaN;
            else if (!Utils.equalRange(vertexSequence,0,bases,offset,vertexSequence.length))
                return Double.NaN;
            else
                return edgeScore + childFinder.score(bases,offset + vertexSequence.length,length - vertexSequence.length);
        }
    }

    @Override
    public String id() {
        return "v" + vertex.getId();
    }

    /**
     * Custom extension of the {@link KBestHaplotype} used for solutions generated by this class.
     *
     * <p>
     *     These by delegating on the encapsulated solution from outgoing edge's finder by adding
     *     the edge score and prefixing this outer finder
     *     source vertex.
     * </p>
     */
    private static class ChildKBestSubHaplotype extends KBestHaplotype {

        private final double score;
        private final KBestHaplotype child;
        private final SeqVertex vertex;
        private final boolean isReference;


        public ChildKBestSubHaplotype(final SeqVertex vertex, final BaseEdge edge, final KBestHaplotype child, final double edgeScore) {
            this.score = edgeScore + child.score();
            this.vertex = vertex;
            this.child = child;
            this.isReference = edge.isRef() && child.isReference();
        }

        @Override
        public SeqGraph graph() {
            return child.graph();
        }

        @Override
        public double score() {
            return score;
        }

        @Override
        public int rank() {
            return child.rank();
        }

        @Override
        protected SeqVertex head() {
            return vertex;
        }

        @Override
        protected KBestHaplotype tail() {
            return child;
        }

        @Override
        public boolean isReference() {
            return isReference;
        }
    }
}
