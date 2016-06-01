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

import com.google.java.contract.Requires;

import java.util.*;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V.  This replaces each
 * Vi with three nodes prefix, seq_i, and suffix connected in a simple chain.
 *
 * This operation can be performed in a very general case, without too much worry about the incoming
 * and outgoing edge structure of each Vi.  The partner algorithm SharedSequenceMerger can
 * put these pieces back together in a smart way that maximizes the sharing of nodes
 * while respecting complex connectivity.
 *
 * User: depristo
 * Date: 3/22/13
 * Time: 8:31 AM
 */
public class CommonSuffixSplitter {
    /**
     * Create a new graph that contains the vertices in toMerge with their shared suffix and prefix
     * sequences extracted out.
     *
     */
    public CommonSuffixSplitter() {}

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @param graph the graph containing the vertices in toMerge
     * @param v The bottom node whose incoming vertices we'd like to split
     * @return true if some useful splitting was done, false otherwise
     */
    public boolean split(final SeqGraph graph, final SeqVertex v) {
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( v == null ) throw new IllegalArgumentException("v cannot be null");
        if ( ! graph.vertexSet().contains(v) ) throw new IllegalArgumentException("graph doesn't contain vertex v " + v);

        final Collection<SeqVertex> toSplit = graph.incomingVerticesOf(v);
        if ( toSplit.size() < 2 )
            // Can only split at least 2 vertices
            return false;
        else if ( ! safeToSplit(graph, v, toSplit) ) {
            return false;
        } else {
            final SeqVertex suffixVTemplate = commonSuffix(toSplit);
            if ( suffixVTemplate.isEmpty() ) {
                return false;
            } else if ( wouldEliminateRefSource(graph, suffixVTemplate, toSplit) ) {
                return false;
            } else if ( allVerticesAreTheCommonSuffix(suffixVTemplate, toSplit) ) {
                return false;
            } else {
                final List<BaseEdge> edgesToRemove = new LinkedList<BaseEdge>();

//                graph.printGraph(new File("split.pre_" + v.getSequenceString() + "." + counter + ".dot"), 0);
                for ( final SeqVertex mid : toSplit ) {
                    // create my own copy of the suffix
                    final SeqVertex suffixV = new SeqVertex(suffixVTemplate.getSequence());
                    graph.addVertex(suffixV);
                    final SeqVertex prefixV = mid.withoutSuffix(suffixV.getSequence());
                    final BaseEdge out = graph.outgoingEdgeOf(mid);

                    final SeqVertex incomingTarget;
                    if ( prefixV == null ) {
                        // this node is entirely explained by suffix
                        incomingTarget = suffixV;
                    } else {
                        incomingTarget = prefixV;
                        graph.addVertex(prefixV);
                        graph.addEdge(prefixV, suffixV, new BaseEdge(out.isRef(), 1));
                        edgesToRemove.add(out);
                    }

                    graph.addEdge(suffixV, graph.getEdgeTarget(out), out.copy());

                    for ( final BaseEdge in : graph.incomingEdgesOf(mid) ) {
                        graph.addEdge(graph.getEdgeSource(in), incomingTarget, in.copy());
                        edgesToRemove.add(in);
                    }
                }

                graph.removeAllVertices(toSplit);
                graph.removeAllEdges(edgesToRemove);
//                graph.printGraph(new File("split.post_" + v.getSequenceString() + "." + counter++ + ".dot"), 0);

                return true;
            }
        }
    }

    /**
     * Would factoring out this suffix result in elimating the reference source vertex?
     * @param graph the graph
     * @param commonSuffix the common suffix of all toSplits
     * @param toSplits the list of vertices we're are trying to split
     * @return true if toSplit contains the reference source and this ref source has all and only the bases of commonSuffix
     */
    private boolean wouldEliminateRefSource(final SeqGraph graph, final SeqVertex commonSuffix, final Collection<SeqVertex> toSplits) {
        for ( final SeqVertex toSplit : toSplits ) {
            if ( graph.isRefSource(toSplit) )
                return toSplit.length() == commonSuffix.length();
        }
        return false;
    }

//    private static int counter = 0;

    /**
     * Would all vertices that we'd split just result in the common suffix?
     *
     * That is, suppose we have prefix nodes ABC and ABC.  After splitting all of the vertices would
     * just be ABC again, and we'd enter into an infinite loop.
     *
     * @param commonSuffix the common suffix of all vertices in toSplits
     * @param toSplits the collection of vertices we want to split
     * @return true if all of the vertices are equal to the common suffix
     */
    private boolean allVerticesAreTheCommonSuffix(final SeqVertex commonSuffix, final Collection<SeqVertex> toSplits) {
        for ( final SeqVertex toSplit : toSplits ) {
            if ( toSplit.length() != commonSuffix.length() )
                return false;
        }

        return true;
    }

    /**
     * Can we safely split up the vertices in toMerge?
     *
     * @param graph a graph
     * @param bot a vertex whose incoming vertices we want to split
     * @param toMerge the set of vertices we'd be splitting up
     * @return true if we can safely split up toMerge
     */
    private boolean safeToSplit(final SeqGraph graph, final SeqVertex bot, final Collection<SeqVertex> toMerge) {
        final Set<SeqVertex> outgoingOfBot = new HashSet<SeqVertex>(graph.outgoingVerticesOf(bot));
        for ( final SeqVertex m : toMerge ) {
            final Set<BaseEdge> outs = graph.outgoingEdgesOf(m);
            if ( m == bot || outs.size() != 1 || ! graph.outgoingVerticesOf(m).contains(bot) )
                // m == bot => don't allow self cycles in the graph
                return false;
            if ( outgoingOfBot.contains(m) )
                // forbid cycles from bottom -> mid
                return false;
        }

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
     * @return a single vertex that contains the common suffix of all middle vertices
     */
    @Requires("!middleVertices.isEmpty()")
    protected static SeqVertex commonSuffix(final Collection<SeqVertex> middleVertices) {
        final List<byte[]> kmers = GraphUtils.getKmers(middleVertices);
        final int min = GraphUtils.minKmerLength(kmers);
        final int suffixLen = GraphUtils.compSuffixLen(kmers, min);
        final byte[] kmer = kmers.get(0);
        final byte[] suffix = Arrays.copyOfRange(kmer, kmer.length - suffixLen, kmer.length);
        return new SeqVertex(suffix);
    }
}