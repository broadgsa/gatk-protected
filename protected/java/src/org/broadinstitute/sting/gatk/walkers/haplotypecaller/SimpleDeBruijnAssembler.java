package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
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
    private static final byte MIN_QUALITY = (byte) 17;

    // Smith-Waterman parameters originally copied from IndelRealigner
    private static final double SW_MATCH = 5.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -22.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -1.2; //-1.0/.0;

    private final boolean DEBUG;
    private final PrintStream GRAPH_WRITER;
    private final ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>> graphs = new ArrayList<DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>>();

    private int PRUNE_FACTOR = 1;
    
    public SimpleDeBruijnAssembler( final boolean debug, final PrintStream graphWriter ) {
        super();
        DEBUG = debug;
        GRAPH_WRITER = graphWriter;
    }

    public ArrayList<Haplotype> runLocalAssembly( final ActiveRegion activeRegion, final Haplotype refHaplotype, final byte[] fullReferenceWithPadding, final GenomeLoc refLoc, final int PRUNE_FACTOR, final ArrayList<VariantContext> activeAllelesToGenotype ) {
        this.PRUNE_FACTOR = PRUNE_FACTOR;

        // create the graphs
        createDeBruijnGraphs( activeRegion.getReads(), refHaplotype );

        // clean up the graphs by pruning and merging
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            pruneGraph( graph, PRUNE_FACTOR );
            //eliminateNonRefPaths( graph );
            mergeNodes( graph );
        }

        if( GRAPH_WRITER != null ) {
            printGraphs();
        }

        // find the best paths in the graphs
        return findBestPaths( refHaplotype, fullReferenceWithPadding, refLoc, activeAllelesToGenotype, activeRegion.getExtendedLoc() );
    }

    private void createDeBruijnGraphs( final ArrayList<GATKSAMRecord> reads, final Haplotype refHaplotype ) {
        graphs.clear();

        // create the graph
        for( int kmer = 31; kmer <= 75; kmer += 6 ) {
            final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);
            if( createGraphFromSequences( graph, reads, kmer, refHaplotype, DEBUG ) ) {
                graphs.add(graph);
            }
        }
    }

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
                        graph.addEdge(addedVertex, graph.getEdgeTarget(edge), new DeBruijnEdge(edge.getIsRef(), edge.getMultiplicity()));
                    }
                    for( final DeBruijnEdge edge : inEdges ) {
                        graph.addEdge(graph.getEdgeSource(edge), addedVertex, new DeBruijnEdge(edge.getIsRef(), edge.getMultiplicity()));
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
        final ArrayList<DeBruijnEdge> edgesToRemove = new ArrayList<DeBruijnEdge>();
        for( final DeBruijnEdge e : graph.edgeSet() ) {
            if( e.getMultiplicity() <= pruneFactor && !e.getIsRef() ) { // remove non-reference edges with weight less than or equal to the pruning factor
                edgesToRemove.add(e);
            }
        }
        graph.removeAllEdges(edgesToRemove);

        // Run through the graph and clean up singular orphaned nodes
        final ArrayList<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
        for( final DeBruijnVertex v : graph.vertexSet() ) {
            if( graph.inDegreeOf(v) == 0 && graph.outDegreeOf(v) == 0 ) {
                verticesToRemove.add(v);
            }
        }
        graph.removeAllVertices(verticesToRemove);
    }

    protected static void eliminateNonRefPaths( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph ) {
        final ArrayList<DeBruijnVertex> verticesToRemove = new ArrayList<DeBruijnVertex>();
        boolean done = false;
        while( !done ) {
            done = true;
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                if( graph.inDegreeOf(v) == 0 || graph.outDegreeOf(v) == 0 ) {
                    boolean isRefNode = false;
                    for( final DeBruijnEdge e : graph.edgesOf(v) ) {
                        if( e.getIsRef() ) {
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

    private static boolean createGraphFromSequences( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, final ArrayList<GATKSAMRecord> reads, final int KMER_LENGTH, final Haplotype refHaplotype, final boolean DEBUG ) {
        final byte[] refSequence = refHaplotype.getBases();
        if( refSequence.length >= KMER_LENGTH + KMER_OVERLAP ) {
            final int kmersInSequence = refSequence.length - KMER_LENGTH + 1;
            for (int i = 0; i < kmersInSequence - 1; i++) {
                // get the kmers
                final byte[] kmer1 = new byte[KMER_LENGTH];
                System.arraycopy(refSequence, i, kmer1, 0, KMER_LENGTH);
                final byte[] kmer2 = new byte[KMER_LENGTH];
                System.arraycopy(refSequence, i+1, kmer2, 0, KMER_LENGTH);
                if( !addKmersToGraph(graph, kmer1, kmer2, true) ) {
                    if( DEBUG ) {
                        System.out.println("Cycle detected in reference graph for kmer = " + KMER_LENGTH + " ...skipping");
                    }
                    return false;
                }
            }
        }

        for( final GATKSAMRecord read : reads ) {
            final byte[] sequence = read.getReadBases();
            final byte[] qualities = read.getBaseQualities();
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
                        // get the kmers
                        final byte[] kmer1 = new byte[KMER_LENGTH];
                        System.arraycopy(sequence, iii, kmer1, 0, KMER_LENGTH);
                        final byte[] kmer2 = new byte[KMER_LENGTH];
                        System.arraycopy(sequence, iii+1, kmer2, 0, KMER_LENGTH);

                        addKmersToGraph(graph, kmer1, kmer2, false);
                    }
                }
            }
        }
        return true;
    }

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

    private void printGraphs() {
        int count = 0;
        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            GRAPH_WRITER.println("digraph kmer" + count++ +" {");
            for( final DeBruijnEdge edge : graph.edgeSet() ) {
                if( edge.getMultiplicity() > PRUNE_FACTOR ) {
                    GRAPH_WRITER.println("\t" + graph.getEdgeSource(edge).toString() + " -> " + graph.getEdgeTarget(edge).toString() + " [" + (edge.getMultiplicity() <= PRUNE_FACTOR ? "style=dotted,color=grey" : "label=\""+ edge.getMultiplicity() +"\"") + "];");
                }
                if( edge.getIsRef() ) {
                    GRAPH_WRITER.println("\t" + graph.getEdgeSource(edge).toString() + " -> " + graph.getEdgeTarget(edge).toString() + " [color=red];");
                }
                if( !edge.getIsRef() && edge.getMultiplicity() <= PRUNE_FACTOR ) { System.out.println("Graph pruning warning!"); }
            }
            for( final DeBruijnVertex v : graph.vertexSet() ) {
                final String label = ( graph.inDegreeOf(v) == 0 ? v.toString() : v.getSuffixString() );
                GRAPH_WRITER.println("\t" + v.toString() + " [label=\"" + label + "\"]");
            }
            GRAPH_WRITER.println("}");
        }
    }

    @Ensures({"result.contains(refHaplotype)"})
    private ArrayList<Haplotype> findBestPaths( final Haplotype refHaplotype, final byte[] fullReferenceWithPadding, final GenomeLoc refLoc, final ArrayList<VariantContext> activeAllelesToGenotype, final GenomeLoc activeRegionWindow ) {
        final ArrayList<Haplotype> returnHaplotypes = new ArrayList<Haplotype>();

        // add the reference haplotype separately from all the others
        final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( fullReferenceWithPadding, refHaplotype.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        refHaplotype.setAlignmentStartHapwrtRef( swConsensus.getAlignmentStart2wrt1() );
        refHaplotype.setCigar( swConsensus.getCigar() );
        if( !returnHaplotypes.add( refHaplotype ) ) {
            throw new ReviewedStingException("Unable to add reference haplotype during assembly: " + refHaplotype);
        }

        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();
        final int activeRegionStop = refHaplotype.getAlignmentStartHapwrtRef() + refHaplotype.getCigar().getReferenceLength();

        for( final VariantContext compVC : activeAllelesToGenotype ) { // for GGA mode, add the desired allele into the haplotype
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart());
                if( !addHaplotype( insertedRefHaplotype, fullReferenceWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop ) ) {
                    return returnHaplotypes;
                    //throw new ReviewedStingException("Unable to add reference+allele haplotype during GGA-enabled assembly: " + insertedRefHaplotype);
                }
            }
        }

        for( final DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph : graphs ) {
            for ( final KBestPaths.Path path : KBestPaths.getKBestPaths(graph, NUM_BEST_PATHS_PER_KMER_GRAPH) ) {
                final Haplotype h = new Haplotype( path.getBases( graph ), path.getScore() );
                if( addHaplotype( h, fullReferenceWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop ) ) {
                    if( !activeAllelesToGenotype.isEmpty() ) { // for GGA mode, add the desired allele into the haplotype if it isn't already present
                        final HashMap<Integer,VariantContext> eventMap = GenotypingEngine.generateVCsFromAlignment( h.getAlignmentStartHapwrtRef(), h.getCigar(), fullReferenceWithPadding, h.getBases(), refLoc, "HCassembly", 0 ); // BUGBUG: need to put this function in a shared place
                        for( final VariantContext compVC : activeAllelesToGenotype ) { // for GGA mode, add the desired allele into the haplotype if it isn't already present
                            final VariantContext vcOnHaplotype = eventMap.get(compVC.getStart());
                            if( vcOnHaplotype == null || !vcOnHaplotype.hasSameAllelesAs(compVC) ) {
                                for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                                    addHaplotype( h.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart()), fullReferenceWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop );
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

    private boolean addHaplotype( final Haplotype haplotype, final byte[] ref, final ArrayList<Haplotype> haplotypeList, final int activeRegionStart, final int activeRegionStop ) {
        //final int sizeOfActiveRegion = activeRegionStop - activeRegionStart;
        final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, haplotype.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        haplotype.setAlignmentStartHapwrtRef( swConsensus.getAlignmentStart2wrt1() );
        haplotype.setCigar( AlignmentUtils.leftAlignIndel(swConsensus.getCigar(), ref, haplotype.getBases(), swConsensus.getAlignmentStart2wrt1(), 0) );

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
        h.setCigar( AlignmentUtils.leftAlignIndel(swConsensus2.getCigar(), ref, h.getBases(), swConsensus2.getAlignmentStart2wrt1(), 0) );
        h.leftBreakPoint = leftBreakPoint;
        h.rightBreakPoint = rightBreakPoint;
        if( swConsensus2.getCigar().toString().contains("S") || swConsensus2.getCigar().getReferenceLength() != activeRegionStop - activeRegionStart ) { // protect against SW failures
            return false;
        }

        if( !haplotypeList.contains(h) ) {
            haplotypeList.add(h);
            return true;
        } else {
            return false;
        }
    }
}