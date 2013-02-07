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
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks, rpoplin
 * Date: Mar 14, 2011
 */

public class SimpleDeBruijnAssembler extends LocalAssemblyEngine {

    private static final int KMER_OVERLAP = 5; // the additional size of a valid chunk of sequence, used to string together k-mers
    private static final int NUM_BEST_PATHS_PER_KMER_GRAPH = 11;
    private static final byte MIN_QUALITY = (byte) 16;

    // Smith-Waterman parameters originally copied from IndelRealigner
    private static final double SW_MATCH = 5.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -22.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -1.2; //-1.0/.0;

    private final boolean DEBUG;
    private final PrintStream GRAPH_WRITER;
    private final List<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>> graphs = new ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>>();
    private final int MIN_KMER;

    private int PRUNE_FACTOR = 2;
    
    public SimpleDeBruijnAssembler( final boolean debug, final PrintStream graphWriter, final int minKmer ) {
        super();
        DEBUG = debug;
        GRAPH_WRITER = graphWriter;
        MIN_KMER = minKmer;
    }

    /**
     * Main entry point into the assembly engine. Build a set of deBruijn graphs out of the provided reference sequence and list of reads
     * @param activeRegion              ActiveRegion object holding the reads which are to be used during assembly
     * @param refHaplotype              reference haplotype object
     * @param fullReferenceWithPadding  byte array holding the reference sequence with padding
     * @param refLoc                    GenomeLoc object corresponding to the reference sequence with padding
     * @param PRUNE_FACTOR              prune kmers from the graph if their weight is <= this value
     * @param activeAllelesToGenotype   the alleles to inject into the haplotypes during GGA mode
     * @return                          a non-empty list of all the haplotypes that are produced during assembly
     */
    @Ensures({"result.contains(refHaplotype)"})
    public List<Haplotype> runLocalAssembly( final ActiveRegion activeRegion, final Haplotype refHaplotype, final byte[] fullReferenceWithPadding, final GenomeLoc refLoc, final int PRUNE_FACTOR, final List<VariantContext> activeAllelesToGenotype ) {
        if( activeRegion == null ) { throw new IllegalArgumentException("Assembly engine cannot be used with a null ActiveRegion."); }
        if( refHaplotype == null ) { throw new IllegalArgumentException("Reference haplotype cannot be null."); }
        if( fullReferenceWithPadding.length != refLoc.size() ) { throw new IllegalArgumentException("Reference bases and reference loc must be the same size."); }
        if( PRUNE_FACTOR < 0 ) { throw new IllegalArgumentException("Pruning factor cannot be negative"); }

        // set the pruning factor for this run of the assembly engine
        this.PRUNE_FACTOR = PRUNE_FACTOR;

        // create the graphs
        createDeBruijnGraphs( activeRegion.getReads(), refHaplotype );

        // clean up the graphs by pruning and merging
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            pruneGraph( graph, PRUNE_FACTOR );
            //eliminateNonRefPaths( graph );
            mergeNodes( graph );
        }

        // print the graphs if the appropriate debug option has been turned on
        if( GRAPH_WRITER != null ) {
            printGraphs();
        }

        // find the best paths in the graphs and return them as haplotypes
        return findBestPaths( refHaplotype, fullReferenceWithPadding, refLoc, activeAllelesToGenotype, activeRegion.getExtendedLoc() );
    }

    @Requires({"reads != null", "refHaplotype != null"})
    protected void createDeBruijnGraphs( final List<GATKSAMRecord> reads, final Haplotype refHaplotype ) {
        graphs.clear();

        final int maxKmer = refHaplotype.getBases().length;
        // create the graph for each possible kmer
        for( int kmer = MIN_KMER; kmer <= maxKmer; kmer += 6 ) {
            final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = createGraphFromSequences( reads, kmer, refHaplotype, DEBUG );
            if( graph != null ) { // graphs that fail during creation ( for example, because there are cycles in the reference graph ) will show up here as a null graph object
                graphs.add(graph);
            }
        }
    }

    @Requires({"graph != null"})
    protected static void mergeNodes( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
        boolean foundNodesToMerge = true;
        while( foundNodesToMerge ) {
            foundNodesToMerge = false;

            for( final DeBruijnEdge e : graph.edgeSet() ) {
                final DeBruijnVertex outgoingVertex = graph.getEdgeTarget(e);
                final DeBruijnVertex incomingVertex = graph.getEdgeSource(e);
                if( !outgoingVertex.equals(incomingVertex) && graph.inDegreeOf(outgoingVertex) == 1 && graph.outDegreeOf(incomingVertex) == 1) {
                    final Set<DeBruijnEdge> outEdges = graph.outgoingEdgesOf(outgoingVertex);
                    final Set<DeBruijnEdge> inEdges = graph.incomingEdgesOf(incomingVertex);
                    if( inEdges.size() == 1 && outEdges.size() == 1 ) {
                        inEdges.iterator().next().setMultiplicity( inEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() / 2 ) );
                        outEdges.iterator().next().setMultiplicity( outEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() / 2 ) );
                    } else if( inEdges.size() == 1 ) {
                        inEdges.iterator().next().setMultiplicity( inEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() - 1 ) );
                    } else if( outEdges.size() == 1 ) {
                        outEdges.iterator().next().setMultiplicity( outEdges.iterator().next().getMultiplicity() + ( e.getMultiplicity() - 1 ) );
                    }

                    final DeBruijnVertex addedVertex = new DeBruijnVertex( ArrayUtils.addAll(incomingVertex.getSequence(), outgoingVertex.getSuffix()), outgoingVertex.kmer );
                    graph.addVertex(addedVertex);
                    for( final DeBruijnEdge edge : outEdges ) {
                        graph.addEdge(addedVertex, graph.getEdgeTarget(edge), new DeBruijnEdge(edge.isRef(), edge.getMultiplicity()));
                    }
                    for( final DeBruijnEdge edge : inEdges ) {
                        graph.addEdge(graph.getEdgeSource(edge), addedVertex, new DeBruijnEdge(edge.isRef(), edge.getMultiplicity()));
                    }

                    graph.removeVertex( incomingVertex );
                    graph.removeVertex( outgoingVertex );
                    foundNodesToMerge = true;
                    break;
                }
            }
        }
    }

    protected static void pruneGraph( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final int pruneFactor ) {
        final List<DeBruijnEdge> edgesToRemove = new ArrayList<DeBruijnEdge>();
        for( final DeBruijnEdge e : graph.edgeSet() ) {
            if( e.getMultiplicity() <= pruneFactor && !e.isRef() ) { // remove non-reference edges with weight less than or equal to the pruning factor
                edgesToRemove.add(e);
            }
        }
        graph.removeAllEdges(edgesToRemove);

        // Run through the graph and clean up singular orphaned nodes
        final List<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
        for( final DeBruijnVertex v : graph.vertexSet() ) {
            if( graph.inDegreeOf(v) == 0 && graph.outDegreeOf(v) == 0 ) {
                verticesToRemove.add(v);
            }
        }
        graph.removeAllVertices(verticesToRemove);
    }

    protected static void eliminateNonRefPaths( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
        final List<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
        boolean done = false;
        while( !done ) {
            done = true;
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                if( graph.inDegreeOf(v) == 0 || graph.outDegreeOf(v) == 0 ) {
                    boolean isRefNode = false;
                    for( final DeBruijnEdge e : graph.edgesOf(v) ) {
                        if( e.isRef() ) {
                            isRefNode = true;
                            break;
                        }
                    }
                    if( !isRefNode ) {
                        done = false;
                        verticesToRemove.add(v);
                    }
                }
            }
            graph.removeAllVertices(verticesToRemove);
            verticesToRemove.clear();
        }
    }

    @Requires({"reads != null", "KMER_LENGTH > 0", "refHaplotype != null"})
    protected static DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> createGraphFromSequences( final List<GATKSAMRecord> reads, final int KMER_LENGTH, final Haplotype refHaplotype, final boolean DEBUG ) {

        final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        // First pull kmers from the reference haplotype and add them to the graph
        final byte[] refSequence = refHaplotype.getBases();
        if( refSequence.length >= KMER_LENGTH + KMER_OVERLAP ) {
            final int kmersInSequence = refSequence.length - KMER_LENGTH + 1;
            for( int iii = 0; iii < kmersInSequence - 1; iii++ ) {
                if( !addKmersToGraph(graph, Arrays.copyOfRange(refSequence, iii, iii + KMER_LENGTH), Arrays.copyOfRange(refSequence, iii + 1, iii + 1 + KMER_LENGTH), true) ) {
                    if( DEBUG ) {
                        System.out.println("Cycle detected in reference graph for kmer = " + KMER_LENGTH + " ...skipping");
                    }
                    return null;
                }
            }
        }

        // Next pull kmers out of every read and throw them on the graph
        for( final GATKSAMRecord read : reads ) {
            final byte[] sequence = read.getReadBases();
            final byte[] qualities = read.getBaseQualities();
            final byte[] reducedReadCounts = read.getReducedReadCounts();  // will be null if read is not reduced
            if( sequence.length > KMER_LENGTH + KMER_OVERLAP ) {
                final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
                for( int iii = 0; iii < kmersInSequence - 1; iii++ ) {                    
                    // if the qualities of all the bases in the kmers are high enough
                    boolean badKmer = false;
                    for( int jjj = iii; jjj < iii + KMER_LENGTH + 1; jjj++) {
                        if( qualities[jjj] < MIN_QUALITY ) {
                            badKmer = true;
                            break;
                        }
                    }
                    if( !badKmer ) {
                        int countNumber = 1;
                        if( read.isReducedRead() ) {
                            // compute mean number of reduced read counts in current kmer span
                            // precise rounding can make a difference with low consensus counts
                            countNumber = MathUtils.arrayMax(Arrays.copyOfRange(reducedReadCounts, iii, iii + KMER_LENGTH));
                        }

                        final byte[] kmer1 = Arrays.copyOfRange(sequence, iii, iii + KMER_LENGTH);
                        final byte[] kmer2 = Arrays.copyOfRange(sequence, iii + 1, iii + 1 + KMER_LENGTH);

                        for( int kkk=0; kkk < countNumber; kkk++ ) {
                            addKmersToGraph(graph, kmer1, kmer2, false);
                        }
                    }
                }
            }
        }
        return graph;
    }

    @Requires({"graph != null", "kmer1.length > 0", "kmer2.length > 0"})
    protected static boolean addKmersToGraph( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final byte[] kmer1, final byte[] kmer2, final boolean isRef ) {

        final int numVertexBefore = graph.vertexSet().size();
        final DeBruijnVertex v1 = new DeBruijnVertex( kmer1, kmer1.length );
        graph.addVertex(v1);
        final DeBruijnVertex v2 = new DeBruijnVertex( kmer2, kmer2.length );
        graph.addVertex(v2);
        if( isRef && graph.vertexSet().size() == numVertexBefore ) { return false; }

        final DeBruijnEdge targetEdge = graph.getEdge(v1, v2);
        if ( targetEdge == null ) {
            graph.addEdge(v1, v2, new DeBruijnEdge( isRef ));
        } else {
            if( isRef ) {
                targetEdge.setIsRef( true );
            }
            targetEdge.setMultiplicity(targetEdge.getMultiplicity() + 1);
        }
        return true;
    }

    protected void printGraphs() {
        int count = 0;
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            GRAPH_WRITER.println("digraph kmer" + count++ +" {");
            for( final DeBruijnEdge edge : graph.edgeSet() ) {
                if( edge.getMultiplicity() > PRUNE_FACTOR ) {
                    GRAPH_WRITER.println("\t" + graph.getEdgeSource(edge).toString() + " -> " + graph.getEdgeTarget(edge).toString() + " [" + (edge.getMultiplicity() <= PRUNE_FACTOR ? "style=dotted,color=grey" : "label=\""+ edge.getMultiplicity() +"\"") + "];");
                }
                if( edge.isRef() ) {
                    GRAPH_WRITER.println("\t" + graph.getEdgeSource(edge).toString() + " -> " + graph.getEdgeTarget(edge).toString() + " [color=red];");
                }
                if( !edge.isRef() && edge.getMultiplicity() <= PRUNE_FACTOR ) { System.out.println("Graph pruning warning!"); }
            }
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                final String label = ( graph.inDegreeOf(v) == 0 ? v.toString() : v.getSuffixString() );
                GRAPH_WRITER.println("\t" + v.toString() + " [label=\"" + label + "\"]");
            }
            GRAPH_WRITER.println("}");
        }
    }

    @Ensures({"result.contains(refHaplotype)"})
    private List<Haplotype> findBestPaths( final Haplotype refHaplotype, final byte[] fullReferenceWithPadding, final GenomeLoc refLoc, final List<VariantContext> activeAllelesToGenotype, final GenomeLoc activeRegionWindow ) {
        final List<Haplotype> returnHaplotypes = new ArrayList<Haplotype>();

        // add the reference haplotype separately from all the others
        final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( fullReferenceWithPadding, refHaplotype.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        refHaplotype.setAlignmentStartHapwrtRef( swConsensus.getAlignmentStart2wrt1() );
        refHaplotype.setCigar( swConsensus.getCigar() );
        if( !returnHaplotypes.add( refHaplotype ) ) {
            throw new ReviewedStingException("Unable to add reference haplotype during assembly: " + refHaplotype);
        }

        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();
        final int activeRegionStop = refHaplotype.getAlignmentStartHapwrtRef() + refHaplotype.getCigar().getReferenceLength();

        // for GGA mode, add the desired allele into the haplotype
        for( final VariantContext compVC : activeAllelesToGenotype ) {
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart());
                addHaplotype( insertedRefHaplotype, fullReferenceWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop, true );
            }
        }

        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            for ( final KBestPaths.Path path : KBestPaths.getKBestPaths(graph, NUM_BEST_PATHS_PER_KMER_GRAPH) ) {

                final Haplotype h = new Haplotype( path.getBases() );
                if( addHaplotype( h, fullReferenceWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop, false ) ) {

                    // for GGA mode, add the desired allele into the haplotype if it isn't already present
                    if( !activeAllelesToGenotype.isEmpty() ) {
                        final Map<Integer,VariantContext> eventMap = GenotypingEngine.generateVCsFromAlignment( h, h.getAlignmentStartHapwrtRef(), h.getCigar(), fullReferenceWithPadding, h.getBases(), refLoc, "HCassembly" ); // BUGBUG: need to put this function in a shared place
                        for( final VariantContext compVC : activeAllelesToGenotype ) { // for GGA mode, add the desired allele into the haplotype if it isn't already present
                            final VariantContext vcOnHaplotype = eventMap.get(compVC.getStart());

                            // This if statement used to additionally have:
                            //      "|| !vcOnHaplotype.hasSameAllelesAs(compVC)"
                            //  but that can lead to problems downstream when e.g. you are injecting a 1bp deletion onto
                            //  a haplotype that already contains a 1bp insertion (so practically it is reference but
                            //  falls into the bin for the 1bp deletion because we keep track of the artificial alleles).
                            if( vcOnHaplotype == null ) {
                                for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                                    addHaplotype( h.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart()), fullReferenceWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop, false );
                                }
                            }
                        }
                    }
                }
            }
        }

        if( DEBUG ) { 
            if( returnHaplotypes.size() > 1 ) {
                System.out.println("Found " + returnHaplotypes.size() + " candidate haplotypes to evaluate every read against.");
            } else {
                System.out.println("Found only the reference haplotype in the assembly graph.");
            }
            for( final Haplotype h : returnHaplotypes ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + h.getCigar() );
            }
        }

        return returnHaplotypes;
    }

    // this function is slated for removal when SWing is removed
    private boolean addHaplotype( final Haplotype haplotype, final byte[] ref, final List<Haplotype> haplotypeList, final int activeRegionStart, final int activeRegionStop, final boolean FORCE_INCLUSION_FOR_GGA_MODE ) {
        if( haplotype == null ) { return false; }

        final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, haplotype.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        haplotype.setAlignmentStartHapwrtRef( swConsensus.getAlignmentStart2wrt1() );
        haplotype.setCigar( AlignmentUtils.leftAlignIndel(swConsensus.getCigar(), ref, haplotype.getBases(), swConsensus.getAlignmentStart2wrt1(), 0, true) );

        if( swConsensus.getCigar().toString().contains("S") || swConsensus.getCigar().getReferenceLength() < 60 ) { // protect against SW failures
            return false;
        }

        final int hapStart = ReadUtils.getReadCoordinateForReferenceCoordinate( haplotype.getAlignmentStartHapwrtRef(), haplotype.getCigar(), activeRegionStart, ReadUtils.ClippingTail.LEFT_TAIL, true );
        int hapStop = ReadUtils.getReadCoordinateForReferenceCoordinate( haplotype.getAlignmentStartHapwrtRef(), haplotype.getCigar(), activeRegionStop, ReadUtils.ClippingTail.RIGHT_TAIL, true );
        if( hapStop == ReadUtils.CLIPPING_GOAL_NOT_REACHED && activeRegionStop == haplotype.getAlignmentStartHapwrtRef() + haplotype.getCigar().getReferenceLength() ) {
            hapStop = activeRegionStop; // contract for getReadCoordinateForReferenceCoordinate function says that if read ends at boundary then it is outside of the clipping goal
        }
        byte[] newHaplotypeBases;
        // extend partial haplotypes to contain the full active region sequence
        int leftBreakPoint = 0;
        int rightBreakPoint = 0;
        if( hapStart == ReadUtils.CLIPPING_GOAL_NOT_REACHED && hapStop == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            newHaplotypeBases = ArrayUtils.addAll( ArrayUtils.addAll( ArrayUtils.subarray(ref, activeRegionStart, swConsensus.getAlignmentStart2wrt1()),
                                                   haplotype.getBases()),
                                                   ArrayUtils.subarray(ref, swConsensus.getAlignmentStart2wrt1() + swConsensus.getCigar().getReferenceLength(), activeRegionStop) );
            leftBreakPoint = swConsensus.getAlignmentStart2wrt1() - activeRegionStart;
            rightBreakPoint = leftBreakPoint + haplotype.getBases().length;
            //newHaplotypeBases = haplotype.getBases();
            //return false; // piece of haplotype isn't anchored within the active region so don't build a haplotype out of it
        } else if( hapStart == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            //return false;
            newHaplotypeBases = ArrayUtils.addAll( ArrayUtils.subarray(ref, activeRegionStart, swConsensus.getAlignmentStart2wrt1()), ArrayUtils.subarray(haplotype.getBases(), 0, hapStop) );
            //newHaplotypeBases = ArrayUtils.subarray(haplotype.getBases(), 0, hapStop);
            leftBreakPoint = swConsensus.getAlignmentStart2wrt1() - activeRegionStart;
        } else if( hapStop == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            //return false;
            newHaplotypeBases = ArrayUtils.addAll( ArrayUtils.subarray(haplotype.getBases(), hapStart, haplotype.getBases().length), ArrayUtils.subarray(ref, swConsensus.getAlignmentStart2wrt1() + swConsensus.getCigar().getReferenceLength(), activeRegionStop) );
            //newHaplotypeBases = ArrayUtils.subarray(haplotype.getBases(), hapStart, haplotype.getBases().length);
            rightBreakPoint = haplotype.getBases().length - hapStart;
        } else {
            newHaplotypeBases = ArrayUtils.subarray(haplotype.getBases(), hapStart, hapStop);
        }

        final Haplotype h = new Haplotype( newHaplotypeBases );
        final SWPairwiseAlignment swConsensus2 = new SWPairwiseAlignment( ref, h.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );

        h.setAlignmentStartHapwrtRef( swConsensus2.getAlignmentStart2wrt1() );
        h.setCigar( AlignmentUtils.leftAlignIndel(swConsensus2.getCigar(), ref, h.getBases(), swConsensus2.getAlignmentStart2wrt1(), 0, true) );
        if ( haplotype.isArtificialHaplotype() ) {
            h.setArtificialEvent(haplotype.getArtificialEvent());
        }
        h.leftBreakPoint = leftBreakPoint;
        h.rightBreakPoint = rightBreakPoint;
        if( swConsensus2.getCigar().toString().contains("S") || swConsensus2.getCigar().getReferenceLength() != activeRegionStop - activeRegionStart ) { // protect against SW failures
            return false;
        }

        if( FORCE_INCLUSION_FOR_GGA_MODE || !haplotypeList.contains(h) ) {
            haplotypeList.add(h);
            return true;
        } else {
            return false;
        }
    }
}