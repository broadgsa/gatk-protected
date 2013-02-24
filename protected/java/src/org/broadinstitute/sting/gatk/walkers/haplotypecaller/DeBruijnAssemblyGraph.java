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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 2/6/13
 */

public class DeBruijnAssemblyGraph extends DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> {

    public DeBruijnAssemblyGraph() {
        super(DeBruijnEdge.class);
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    public boolean isReferenceNode( final DeBruijnVertex v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        for( final DeBruijnEdge e : edgesOf(v) ) {
            if( e.isRef() ) { return true; }
        }
        return false;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node
     */
    public boolean isSource( final DeBruijnVertex v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        return inDegreeOf(v) == 0;
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    @Ensures({"result != null"})
    public byte[] getAdditionalSequence( final DeBruijnVertex v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to pull sequence from a null vertex."); }
        return ( isSource(v) ? v.getSequence() : v.getSuffix() );
    }

    /**
     * @param e the edge to test
     * @return  true if this edge is a reference source edge
     */
    public boolean isRefSource( final DeBruijnEdge e ) {
        if( e == null ) { throw new IllegalArgumentException("Attempting to test a null edge."); }
        for( final DeBruijnEdge edgeToTest : incomingEdgesOf(getEdgeSource(e)) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    public boolean isRefSource( final DeBruijnVertex v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        for( final DeBruijnEdge edgeToTest : incomingEdgesOf(v) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @param e the edge to test
     * @return  true if this edge is a reference sink edge
     */
    public boolean isRefSink( final DeBruijnEdge e ) {
        if( e == null ) { throw new IllegalArgumentException("Attempting to test a null edge."); }
        for( final DeBruijnEdge edgeToTest : outgoingEdgesOf(getEdgeTarget(e)) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference sink
     */
    public boolean isRefSink( final DeBruijnVertex v ) {
        if( v == null ) { throw new IllegalArgumentException("Attempting to test a null vertex."); }
        for( final DeBruijnEdge edgeToTest : outgoingEdgesOf(v) ) {
            if( edgeToTest.isRef() ) { return false; }
        }
        return true;
    }

    /**
     * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public DeBruijnVertex getReferenceSourceVertex( ) {
        for( final DeBruijnVertex v : vertexSet() ) {
            if( isReferenceNode(v) && isRefSource(v) ) {
                return v;
            }
        }
        return null;
    }

    /**
     * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public DeBruijnVertex getReferenceSinkVertex( ) {
        for( final DeBruijnVertex v : vertexSet() ) {
            if( isReferenceNode(v) && isRefSink(v) ) {
                return v;
            }
        }
        return null;
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the next reference vertex if it exists
     */
    public DeBruijnVertex getNextReferenceVertex( final DeBruijnVertex v ) {
        if( v == null ) { return null; }
        for( final DeBruijnEdge edgeToTest : outgoingEdgesOf(v) ) {
            if( edgeToTest.isRef() ) {
                return getEdgeTarget(edgeToTest);
            }
        }
        return null;
    }

    /**
     * Traverse the graph and get the previous reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the previous reference vertex if it exists
     */
    public DeBruijnVertex getPrevReferenceVertex( final DeBruijnVertex v ) {
        if( v == null ) { return null; }
        for( final DeBruijnEdge edgeToTest : incomingEdgesOf(v) ) {
            if( isReferenceNode(getEdgeSource(edgeToTest)) ) {
                return getEdgeSource(edgeToTest);
            }
        }
        return null;
    }

    /**
     * Does a reference path exist between the two vertices?
     * @param fromVertex    from this vertex, can be null
     * @param toVertex      to this vertex, can be null
     * @return              true if a reference path exists in the graph between the two vertices
     */
    public boolean referencePathExists(final DeBruijnVertex fromVertex, final DeBruijnVertex toVertex) {
        DeBruijnVertex v = fromVertex;
        if( v == null ) {
            return false;
        }
        v = getNextReferenceVertex(v);
        if( v == null ) {
            return false;
        }
        while( !v.equals(toVertex) ) {
            v = getNextReferenceVertex(v);
            if( v == null ) {
                return false;
            }
        }
        return true;
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    public byte[] getReferenceBytes( final DeBruijnVertex fromVertex, final DeBruijnVertex toVertex, final boolean includeStart, final boolean includeStop ) {
        if( fromVertex == null ) { throw new IllegalArgumentException("Starting vertex in requested path cannot be null."); }
        if( toVertex == null ) { throw  new IllegalArgumentException("From vertex in requested path cannot be null."); }

        byte[] bytes = null;
        DeBruijnVertex v = fromVertex;
        if( includeStart ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        v = getNextReferenceVertex(v); // advance along the reference path
        while( v != null && !v.equals(toVertex) ) {
            bytes = ArrayUtils.addAll( bytes, getAdditionalSequence(v) );
            v = getNextReferenceVertex(v); // advance along the reference path
        }
        if( includeStop && v != null && v.equals(toVertex)) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        return bytes;
    }

    /**
     * Pull kmers out of the given long sequence and throw them on in the graph
     * @param sequence      byte array holding the sequence with which to build the assembly graph
     * @param KMER_LENGTH   the desired kmer length to use
     * @param isRef         if true the kmers added to the graph will have reference edges linking them
     */
    public void addSequenceToGraph( final byte[] sequence, final int KMER_LENGTH, final boolean isRef ) {
        if( sequence.length < KMER_LENGTH + 1 ) { throw new IllegalArgumentException("Provided sequence is too small for the given kmer length"); }
        final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
        for( int iii = 0; iii < kmersInSequence - 1; iii++ ) {
            addKmersToGraph(Arrays.copyOfRange(sequence, iii, iii + KMER_LENGTH), Arrays.copyOfRange(sequence, iii + 1, iii + 1 + KMER_LENGTH), isRef);
        }
    }

    /**
     * Add edge to assembly graph connecting the two kmers
     * @param kmer1 the source kmer for the edge
     * @param kmer2 the target kmer for the edge
     * @param isRef true if the added edge is a reference edge
     * @return      will return false if trying to add a reference edge which creates a cycle in the assembly graph
     */
    public boolean addKmersToGraph( final byte[] kmer1, final byte[] kmer2, final boolean isRef ) {
        if( kmer1 == null ) { throw new IllegalArgumentException("Attempting to add a null kmer to the graph."); }
        if( kmer2 == null ) { throw new IllegalArgumentException("Attempting to add a null kmer to the graph."); }
        if( kmer1.length != kmer2.length ) { throw new IllegalArgumentException("Attempting to add a kmers to the graph with different lengths."); }

        final int numVertexBefore = vertexSet().size();
        final DeBruijnVertex v1 = new DeBruijnVertex( kmer1, kmer1.length );
        addVertex(v1);
        final DeBruijnVertex v2 = new DeBruijnVertex( kmer2, kmer2.length );
        addVertex(v2);
        if( isRef && vertexSet().size() == numVertexBefore ) { return false; }

        final DeBruijnEdge targetEdge = getEdge(v1, v2);
        if ( targetEdge == null ) {
            addEdge(v1, v2, new DeBruijnEdge( isRef ));
        } else {
            if( isRef ) {
                targetEdge.setIsRef( true );
            }
            targetEdge.setMultiplicity(targetEdge.getMultiplicity() + 1);
        }
        return true;
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param GRAPH_WRITER  PrintStream to write to
     */
    public void printGraph( final PrintStream GRAPH_WRITER ) {
        if( GRAPH_WRITER == null ) { throw new IllegalArgumentException("PrintStream cannot be null."); }

        GRAPH_WRITER.println("digraph assembly {");
        for( final DeBruijnEdge edge : edgeSet() ) {
            GRAPH_WRITER.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [" + "label=\""+ edge.getMultiplicity() +"\"" + "];");
            if( edge.isRef() ) {
                GRAPH_WRITER.println("\t" + getEdgeSource(edge).toString() + " -> " + getEdgeTarget(edge).toString() + " [color=red];");
            }
        }
        for( final DeBruijnVertex v : vertexSet() ) {
            final String label = ( inDegreeOf(v) == 0 ? v.toString() : v.getSuffixString() );
            GRAPH_WRITER.println("\t" + v.toString() + " [label=\"" + label + "\"]");
        }
        GRAPH_WRITER.println("}");
    }
}
