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
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks, rpoplin
 * Date: Mar 14, 2011
 */

public class DeBruijnAssembler extends LocalAssemblyEngine {
    private final static Logger logger = Logger.getLogger(DeBruijnAssembler.class);

    private static final int KMER_OVERLAP = 5; // the additional size of a valid chunk of sequence, used to string together k-mers

    // TODO -- this number is very low, and limits our ability to explore low-frequnecy variants.  It should
    // TODO -- be increased to a large number of eliminated altogether when moving to the bubble caller where
    // TODO -- we are no longer considering a combinatorial number of haplotypes as the number of bubbles increases
    private static final int NUM_BEST_PATHS_PER_KMER_GRAPH = 25;
    private static final int GRAPH_KMER_STEP = 6;

    // Smith-Waterman parameters originally copied from IndelRealigner, only used during GGA mode
    private static final double SW_MATCH = 5.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -22.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -1.2; //-1.0/.0;

    private final boolean debug;
    private final boolean debugGraphTransformations;
    private final int minKmer;

    private final int onlyBuildKmersOfThisSizeWhenDebuggingGraphAlgorithms;


    protected DeBruijnAssembler() {
        this(false, -1, 11);
    }

    public DeBruijnAssembler(final boolean debug,
                             final int debugGraphTransformations,
                             final int minKmer) {
        super();
        this.debug = debug;
        this.debugGraphTransformations = debugGraphTransformations > 0;
        this.onlyBuildKmersOfThisSizeWhenDebuggingGraphAlgorithms = debugGraphTransformations;
        this.minKmer = minKmer;
    }

    /**
     * Main entry point into the assembly engine. Build a set of deBruijn graphs out of the provided reference sequence and list of reads
     * @param activeRegion              ActiveRegion object holding the reads which are to be used during assembly
     * @param refHaplotype              reference haplotype object
     * @param fullReferenceWithPadding  byte array holding the reference sequence with padding
     * @param refLoc                    GenomeLoc object corresponding to the reference sequence with padding
     * @param activeAllelesToGenotype   the alleles to inject into the haplotypes during GGA mode
     * @return                          a non-empty list of all the haplotypes that are produced during assembly
     */
    @Ensures({"result.contains(refHaplotype)"})
    public List<Haplotype> runLocalAssembly( final ActiveRegion activeRegion, final Haplotype refHaplotype, final byte[] fullReferenceWithPadding, final GenomeLoc refLoc, final List<VariantContext> activeAllelesToGenotype ) {
        if( activeRegion == null ) { throw new IllegalArgumentException("Assembly engine cannot be used with a null ActiveRegion."); }
        if( refHaplotype == null ) { throw new IllegalArgumentException("Reference haplotype cannot be null."); }
        if( fullReferenceWithPadding.length != refLoc.size() ) { throw new IllegalArgumentException("Reference bases and reference loc must be the same size."); }
        if( pruneFactor < 0 ) { throw new IllegalArgumentException("Pruning factor cannot be negative"); }

        // create the graphs
        final List<SeqGraph> graphs = createDeBruijnGraphs( activeRegion.getReads(), refHaplotype );

        // print the graphs if the appropriate debug option has been turned on
        if( graphWriter != null ) {
            printGraphs(graphs);
        }

        // find the best paths in the graphs and return them as haplotypes
        return findBestPaths( graphs, refHaplotype, fullReferenceWithPadding, refLoc, activeAllelesToGenotype, activeRegion.getExtendedLoc() );
    }

    @Requires({"reads != null", "refHaplotype != null"})
    protected List<SeqGraph> createDeBruijnGraphs( final List<GATKSAMRecord> reads, final Haplotype refHaplotype ) {
        final List<SeqGraph> graphs = new LinkedList<SeqGraph>();

        final int maxKmer = ReadUtils.getMaxReadLength(reads) - KMER_OVERLAP - 1;
        if( maxKmer < minKmer) {
            // Reads are too small for assembly so don't try to create any assembly graphs
            return Collections.emptyList();
        }
        // create the graph for each possible kmer
        for( int kmer = maxKmer; kmer >= minKmer; kmer -= GRAPH_KMER_STEP ) {
            if ( debugGraphTransformations && kmer > onlyBuildKmersOfThisSizeWhenDebuggingGraphAlgorithms)
                continue;

            if ( debug ) logger.info("Creating de Bruijn graph for " + kmer + " kmer using " + reads.size() + " reads");
            DeBruijnGraph graph = createGraphFromSequences( reads, kmer, refHaplotype);
            if( graph != null ) { // graphs that fail during creation ( for example, because there are cycles in the reference graph ) will show up here as a null graph object
                // do a series of steps to clean up the raw assembly graph to make it analysis-ready
                if ( debugGraphTransformations ) graph.printGraph(new File("unpruned.dot"), pruneFactor);

                if ( shouldErrorCorrectKmers() ) {
                    throw new UserException("Error correction no longer supported because of the " +
                            "incredibly naive way this was implemented.  The command line argument remains because some" +
                            " future subsystem will actually go and error correct the reads");
                }

                final SeqGraph seqGraph = toSeqGraph(graph);

                if ( seqGraph != null ) { // if the graph contains interesting variation from the reference
                    sanityCheckReferenceGraph(seqGraph, refHaplotype);
                    graphs.add(seqGraph);

                    if ( debugGraphTransformations ) // we only want to use one graph size
                        break;
                }
            }

        }

        return graphs;
    }

    private SeqGraph toSeqGraph(final DeBruijnGraph deBruijnGraph) {
        final SeqGraph seqGraph = deBruijnGraph.convertToSequenceGraph();
        if ( debugGraphTransformations ) seqGraph.printGraph(new File("sequenceGraph.1.dot"), pruneFactor);

        // TODO -- we need to come up with a consistent pruning algorithm.  The current pruning algorithm
        // TODO -- works well but it doesn't differentiate between an isolated chain that doesn't connect
        // TODO -- to anything from one that's actuall has good support along the chain but just happens
        // TODO -- to have a connection in the middle that has weight of < pruneFactor.  Ultimately
        // TODO -- the pruning algorithm really should be an error correction algorithm that knows more
        // TODO -- about the structure of the data and can differeniate between an infrequent path but
        // TODO -- without evidence against it (such as occurs when a region is hard to get any reads through)
        // TODO -- from a error with lots of weight going along another similar path
        // the very first thing we need to do is zip up the graph, or pruneGraph will be too aggressive
        seqGraph.zipLinearChains();
        if ( debugGraphTransformations ) seqGraph.printGraph(new File("sequenceGraph.2.zipped.dot"), pruneFactor);

        // now go through and prune the graph, removing vertices no longer connected to the reference chain
        // IMPORTANT: pruning must occur before we call simplifyGraph, as simplifyGraph adds 0 weight
        // edges to maintain graph connectivity.
        seqGraph.pruneGraph(pruneFactor);
        seqGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();

        if ( debugGraphTransformations ) seqGraph.printGraph(new File("sequenceGraph.3.pruned.dot"), pruneFactor);
        seqGraph.simplifyGraph();
        if ( debugGraphTransformations ) seqGraph.printGraph(new File("sequenceGraph.4.merged.dot"), pruneFactor);

        // The graph has degenerated in some way, so the reference source and/or sink cannot be id'd.  Can
        // happen in cases where for example the reference somehow manages to acquire a cycle, or
        // where the entire assembly collapses back into the reference sequence.
        if ( seqGraph.getReferenceSourceVertex() == null || seqGraph.getReferenceSinkVertex() == null )
            return null;

        seqGraph.removePathsNotConnectedToRef();
        seqGraph.simplifyGraph();
        if ( seqGraph.vertexSet().size() == 1 ) {
            // we've prefectly assembled into a single reference haplotype, add a empty seq vertex to stop
            // the code from blowing up.
            // TODO -- ref properties should really be on the vertices, not the graph itself
            final SeqVertex complete = seqGraph.vertexSet().iterator().next();
            final SeqVertex dummy = new SeqVertex("");
            seqGraph.addVertex(dummy);
            seqGraph.addEdge(complete, dummy, new BaseEdge(true, 0));
        }
        if ( debugGraphTransformations ) seqGraph.printGraph(new File("sequenceGraph.5.final.dot"), pruneFactor);

        return seqGraph;
    }

    protected <T extends BaseVertex> void sanityCheckReferenceGraph(final BaseGraph<T> graph, final Haplotype refHaplotype) {
        if( graph.getReferenceSourceVertex() == null ) {
            throw new IllegalStateException("All reference graphs must have a reference source vertex.");
        }
        if( graph.getReferenceSinkVertex() == null ) {
            throw new IllegalStateException("All reference graphs must have a reference sink vertex.");
        }
        if( !Arrays.equals(graph.getReferenceBytes(graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex(), true, true), refHaplotype.getBases()) ) {
            throw new IllegalStateException("Mismatch between the reference haplotype and the reference assembly graph path." +
                    " graph = " + new String(graph.getReferenceBytes(graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex(), true, true)) +
                    " haplotype = " + new String(refHaplotype.getBases())
            );
        }
    }

    @Requires({"reads != null", "kmerLength > 0", "refHaplotype != null"})
    protected DeBruijnGraph createGraphFromSequences( final List<GATKSAMRecord> reads, final int kmerLength, final Haplotype refHaplotype ) {
        final DeBruijnGraph graph = new DeBruijnGraph(kmerLength);

        // First pull kmers from the reference haplotype and add them to the graph
        if ( ! addReferenceKmersToGraph(graph, refHaplotype.getBases()) )
            // something went wrong, so abort right now with a null graph
            return null;

        // now go through the graph already seeded with the reference sequence and add the read kmers to it
        if ( ! addReadKmersToGraph(graph, reads) )
            // some problem was detected adding the reads to the graph, return null to indicate we failed
            return null;

        graph.cleanNonRefPaths();
        return graph;
    }

    /**
     * Add the high-quality kmers from the reads to the graph
     *
     * @param graph a graph to add the read kmers to
     * @param reads a non-null list of reads whose kmers we want to add to the graph
     * @return true if we successfully added the read kmers to the graph without corrupting it in some way
     */
    protected boolean addReadKmersToGraph(final DeBruijnGraph graph, final List<GATKSAMRecord> reads) {
        final int kmerLength = graph.getKmerSize();

        // Next pull kmers out of every read and throw them on the graph
        for( final GATKSAMRecord read : reads ) {
            final byte[] sequence = read.getReadBases();
            final byte[] qualities = read.getBaseQualities();
            final byte[] reducedReadCounts = read.getReducedReadCounts();  // will be null if read is not reduced
            if( sequence.length > kmerLength + KMER_OVERLAP ) {
                final int kmersInSequence = sequence.length - kmerLength + 1;
                for( int iii = 0; iii < kmersInSequence - 1; iii++ ) {
                    // if the qualities of all the bases in the kmers are high enough
                    boolean badKmer = false;
                    for( int jjj = iii; jjj < iii + kmerLength + 1; jjj++) {
                        if( qualities[jjj] < minBaseQualityToUseInAssembly ) {
                            badKmer = true;
                            break;
                        }
                    }
                    if( !badKmer ) {
                        // how many observations of this kmer have we seen?  A normal read counts for 1, but
                        // a reduced read might imply a higher multiplicity for our the edge
                        int countNumber = 1;
                        if( read.isReducedRead() ) {
                            // compute mean number of reduced read counts in current kmer span
                            // precise rounding can make a difference with low consensus counts
                            // TODO -- optimization: should extend arrayMax function to take start stop values
                            countNumber = MathUtils.arrayMax(Arrays.copyOfRange(reducedReadCounts, iii, iii + kmerLength));
                        }

                        graph.addKmerPairFromSeqToGraph(sequence, iii, false, countNumber);
                    }
                }
            }
        }

        // always returns true now, but it's possible that we'd add reads and decide we don't like the graph in some way
        return true;
    }

    /**
     * Add the kmers from the reference sequence to the DeBruijnGraph
     *
     * @param graph the graph to add the reference kmers to. Must be empty
     * @param refSequence the reference sequence from which we'll get our kmers
     * @return true if we succeeded in creating a good graph from the reference sequence, false otherwise
     */
    protected boolean addReferenceKmersToGraph(final DeBruijnGraph graph, final byte[] refSequence) {
        if ( graph == null ) throw new IllegalArgumentException("graph cannot be null");
        if ( graph.vertexSet().size() != 0 ) throw new IllegalArgumentException("Reference sequences must be added before any other vertices, but got a graph with " + graph.vertexSet().size() + " vertices in it already: " + graph);
        if ( refSequence == null ) throw new IllegalArgumentException("refSequence cannot be null");


        final int kmerLength = graph.getKmerSize();
        if( refSequence.length < kmerLength + KMER_OVERLAP ) {
            // not enough reference sequence to build a kmer graph of this length, return null
            return false;
        }

        final int kmersInSequence = refSequence.length - kmerLength + 1;
        for( int iii = 0; iii < kmersInSequence - 1; iii++ ) {
            graph.addKmerPairFromSeqToGraph(refSequence, iii, true, 1);
        }

        // we expect that every kmer in the sequence is unique, so that the graph has exactly kmersInSequence vertices
        if ( graph.vertexSet().size() != kmersInSequence ) {
            if( debug ) logger.info("Cycle detected in reference graph for kmer = " + kmerLength + " ...skipping");
            return false;
        }

        return true;
    }

    protected void printGraphs(final List<SeqGraph> graphs) {
        final int writeFirstGraphWithSizeSmallerThan = 50;

        graphWriter.println("digraph assemblyGraphs {");
        for( final SeqGraph graph : graphs ) {
            if ( debugGraphTransformations && graph.getKmerSize() >= writeFirstGraphWithSizeSmallerThan ) {
                logger.info("Skipping writing of graph with kmersize " + graph.getKmerSize());
                continue;
            }

            graph.printGraph(graphWriter, false, pruneFactor);

            if ( debugGraphTransformations )
                break;
        }

        graphWriter.println("}");
    }

    @Requires({"refWithPadding.length > refHaplotype.getBases().length", "refLoc.containsP(activeRegionWindow)"})
    @Ensures({"result.contains(refHaplotype)"})
    private List<Haplotype> findBestPaths( final List<SeqGraph> graphs, final Haplotype refHaplotype, final byte[] refWithPadding, final GenomeLoc refLoc, final List<VariantContext> activeAllelesToGenotype, final GenomeLoc activeRegionWindow ) {

        // add the reference haplotype separately from all the others to ensure that it is present in the list of haplotypes
        // TODO -- this use of an array with contains lower may be a performance problem returning in an O(N^2) algorithm
        final List<Haplotype> returnHaplotypes = new ArrayList<Haplotype>();
        refHaplotype.setAlignmentStartHapwrtRef(activeRegionWindow.getStart() - refLoc.getStart());
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);
        returnHaplotypes.add( refHaplotype );

        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();
        final int activeRegionStop = refHaplotype.getAlignmentStartHapwrtRef() + refHaplotype.getCigar().getReferenceLength();

        // for GGA mode, add the desired allele into the haplotype
        for( final VariantContext compVC : activeAllelesToGenotype ) {
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart());
                addHaplotypeForGGA( insertedRefHaplotype, refWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop, true );
            }
        }

        for( final SeqGraph graph : graphs ) {
            for ( final Path<SeqVertex> path : new KBestPaths<SeqVertex>().getKBestPaths(graph, NUM_BEST_PATHS_PER_KMER_GRAPH) ) {
//                logger.info("Found path " + path);
                Haplotype h = new Haplotype( path.getBases() );
                if( !returnHaplotypes.contains(h) ) {
                    final Cigar cigar = path.calculateCigar();
                    if( cigar.isEmpty() ) {
                        throw new IllegalStateException("Smith-Waterman alignment failure. Cigar = " + cigar + " with reference length " + cigar.getReferenceLength() + " but expecting reference length of " + refHaplotype.getCigar().getReferenceLength());
                    } else if ( pathIsTooDivergentFromReference(cigar) || cigar.getReferenceLength() < 60 ) { // N cigar elements means that a bubble was too divergent from the reference so skip over this path
                        continue;
                    } else if( cigar.getReferenceLength() != refHaplotype.getCigar().getReferenceLength() ) { // SW failure
                        throw new IllegalStateException("Smith-Waterman alignment failure. Cigar = " + cigar + " with reference length " + cigar.getReferenceLength() + " but expecting reference length of " + refHaplotype.getCigar().getReferenceLength());
                    }
                    h.setCigar(cigar);

                    // extend partial haplotypes which are anchored in the reference to include the full active region
                    h = extendPartialHaplotype(h, activeRegionStart, refWithPadding);
                    final Cigar leftAlignedCigar = leftAlignCigarSequentially(AlignmentUtils.consolidateCigar(h.getCigar()), refWithPadding, h.getBases(), activeRegionStart, 0);
                    if( leftAlignedCigar.getReferenceLength() != refHaplotype.getCigar().getReferenceLength() ) { // left alignment failure
                        continue;
                    }
                    if( !returnHaplotypes.contains(h) ) {
                        h.setAlignmentStartHapwrtRef(activeRegionStart);
                        h.setCigar(leftAlignedCigar);
                        h.setScore(path.getScore());
                        returnHaplotypes.add(h);

                        if ( debug )
                            logger.info("Adding haplotype " + h.getCigar() + " from debruijn graph with kmer " + graph.getKmerSize());

                        // for GGA mode, add the desired allele into the haplotype if it isn't already present
                        if( !activeAllelesToGenotype.isEmpty() ) {
                            final Map<Integer,VariantContext> eventMap = GenotypingEngine.generateVCsFromAlignment( h, refWithPadding, refLoc, "HCassembly" ); // BUGBUG: need to put this function in a shared place
                            for( final VariantContext compVC : activeAllelesToGenotype ) { // for GGA mode, add the desired allele into the haplotype if it isn't already present
                                final VariantContext vcOnHaplotype = eventMap.get(compVC.getStart());

                                // This if statement used to additionally have:
                                //      "|| !vcOnHaplotype.hasSameAllelesAs(compVC)"
                                //  but that can lead to problems downstream when e.g. you are injecting a 1bp deletion onto
                                //  a haplotype that already contains a 1bp insertion (so practically it is reference but
                                //  falls into the bin for the 1bp deletion because we keep track of the artificial alleles).
                                if( vcOnHaplotype == null ) {
                                    for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                                        addHaplotypeForGGA( h.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart()), refWithPadding, returnHaplotypes, activeRegionStart, activeRegionStop, false );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // add genome locs to the haplotypes
        for ( final Haplotype h : returnHaplotypes ) h.setGenomeLocation(activeRegionWindow);

        if ( returnHaplotypes.size() < returnHaplotypes.size() )
            logger.info("Found " + returnHaplotypes.size() + " candidate haplotypes of " + returnHaplotypes.size() + " possible combinations to evaluate every read against at " + refLoc);

        if( debug ) {
            if( returnHaplotypes.size() > 1 ) {
                logger.info("Found " + returnHaplotypes.size() + " candidate haplotypes of " + returnHaplotypes.size() + " possible combinations to evaluate every read against.");
            } else {
                logger.info("Found only the reference haplotype in the assembly graph.");
            }
            for( final Haplotype h : returnHaplotypes ) {
                logger.info( h.toString() );
                logger.info( "> Cigar = " + h.getCigar() + " : " + h.getCigar().getReferenceLength() + " score " + h.getScore() );
            }
        }

        return returnHaplotypes;
    }

    /**
     * Extend partial haplotypes which are anchored in the reference to include the full active region
     * @param haplotype             the haplotype to extend
     * @param activeRegionStart     the place where the active region starts in the ref byte array
     * @param refWithPadding        the full reference byte array with padding which encompasses the active region
     * @return                      a haplotype fully extended to encompass the active region
     */
    @Requires({"haplotype != null", "activeRegionStart >= 0", "refWithPadding != null", "refWithPadding.length > 0"})
    @Ensures({"result != null", "result.getCigar() != null"})
    private Haplotype extendPartialHaplotype( final Haplotype haplotype, final int activeRegionStart, final byte[] refWithPadding ) {
        final Cigar cigar = haplotype.getCigar();
        final Cigar newCigar = new Cigar();
        byte[] newHaplotypeBases = haplotype.getBases();
        int refPos = activeRegionStart;
        int hapPos = 0;
        for( int iii = 0; iii < cigar.getCigarElements().size(); iii++ ) {
            final CigarElement ce = cigar.getCigarElement(iii);
            switch (ce.getOperator()) {
                case M:
                    refPos += ce.getLength();
                    hapPos += ce.getLength();
                    newCigar.add(ce);
                    break;
                case I:
                    hapPos += ce.getLength();
                    newCigar.add(ce);
                    break;
                case D:
                    if( iii == 0 || iii == cigar.getCigarElements().size() - 1 ) {
                        newHaplotypeBases = ArrayUtils.addAll( Arrays.copyOfRange(newHaplotypeBases, 0, hapPos),
                                ArrayUtils.addAll(Arrays.copyOfRange(refWithPadding, refPos, refPos + ce.getLength()),
                                Arrays.copyOfRange(newHaplotypeBases, hapPos, newHaplotypeBases.length)));
                        hapPos += ce.getLength();
                        refPos += ce.getLength();
                        newCigar.add(new CigarElement(ce.getLength(), CigarOperator.M));
                    } else {
                        refPos += ce.getLength();
                        newCigar.add(ce);
                    }
                    break;
                default:
                    throw new IllegalStateException("Unsupported cigar operator detected: " + ce.getOperator());
            }
        }
        final Haplotype returnHaplotype = new Haplotype(newHaplotypeBases, haplotype.isReference());
        returnHaplotype.setCigar( newCigar );
        return returnHaplotype;
    }

    /**
     * We use CigarOperator.N as the signal that an incomplete or too divergent bubble was found during bubble traversal
     * @param c the cigar to test
     * @return  true if we should skip over this path
     */
    @Requires("c != null")
    private boolean pathIsTooDivergentFromReference( final Cigar c ) {
        for( final CigarElement ce : c.getCigarElements() ) {
            if( ce.getOperator().equals(CigarOperator.N) ) {
                return true;
            }
        }
        return false;
    }

    /**
     * Left align the given cigar sequentially. This is needed because AlignmentUtils doesn't accept cigars with more than one indel in them.
     * This is a target of future work to incorporate and generalize into AlignmentUtils for use by others.
     * @param cigar     the cigar to left align
     * @param refSeq    the reference byte array
     * @param readSeq   the read byte array
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return          the left-aligned cigar
     */
    @Ensures({"cigar != null", "refSeq != null", "readSeq != null", "refIndex >= 0", "readIndex >= 0"})
    protected Cigar leftAlignCigarSequentially(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        final Cigar cigarToReturn = new Cigar();
        Cigar cigarToAlign = new Cigar();
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            final CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                cigarToAlign.add(ce);
                for( final CigarElement toAdd : AlignmentUtils.leftAlignIndel(cigarToAlign, refSeq, readSeq, refIndex, readIndex, false).getCigarElements() ) {
                    cigarToReturn.add(toAdd);
                }
                refIndex += cigarToAlign.getReferenceLength();
                readIndex += cigarToAlign.getReadLength();
                cigarToAlign = new Cigar();
            } else {
                cigarToAlign.add(ce);
            }
        }
        if( !cigarToAlign.isEmpty() ) {
            for( final CigarElement toAdd : cigarToAlign.getCigarElements() ) {
                cigarToReturn.add(toAdd);
            }
        }
        return cigarToReturn;
    }

    /**
     * Take a haplotype which was generated by injecting an allele into a string of bases and run SW against the reference to determine the variants on the haplotype.
     * Unfortunately since this haplotype didn't come from the assembly graph you can't straightforwardly use the bubble traversal algorithm to get this information.
     * This is a target for future work as we rewrite the HaplotypeCaller to be more bubble-caller based.
     * @param haplotype                     the candidate haplotype
     * @param ref                           the reference bases to align against
     * @param haplotypeList                 the current list of haplotypes
     * @param activeRegionStart             the start of the active region in the reference byte array
     * @param activeRegionStop              the stop of the active region in the reference byte array
     * @param FORCE_INCLUSION_FOR_GGA_MODE  if true will include in the list even if it already exists
     * @return                              true if the candidate haplotype was successfully incorporated into the haplotype list
     */
    @Requires({"ref != null", "ref.length >= activeRegionStop - activeRegionStart"})
    private boolean addHaplotypeForGGA( final Haplotype haplotype, final byte[] ref, final List<Haplotype> haplotypeList, final int activeRegionStart, final int activeRegionStop, final boolean FORCE_INCLUSION_FOR_GGA_MODE ) {
        if( haplotype == null ) { return false; }

        final SWPairwiseAlignment swConsensus = new SWPairwiseAlignment( ref, haplotype.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        haplotype.setAlignmentStartHapwrtRef( swConsensus.getAlignmentStart2wrt1() );

        if( swConsensus.getCigar().toString().contains("S") || swConsensus.getCigar().getReferenceLength() < 60 || swConsensus.getAlignmentStart2wrt1() < 0 ) { // protect against unhelpful haplotype alignments
            return false;
        }

        haplotype.setCigar( AlignmentUtils.leftAlignIndel(swConsensus.getCigar(), ref, haplotype.getBases(), swConsensus.getAlignmentStart2wrt1(), 0, true) );

        final int hapStart = ReadUtils.getReadCoordinateForReferenceCoordinate(haplotype.getAlignmentStartHapwrtRef(), haplotype.getCigar(), activeRegionStart, ReadUtils.ClippingTail.LEFT_TAIL, true);
        int hapStop = ReadUtils.getReadCoordinateForReferenceCoordinate( haplotype.getAlignmentStartHapwrtRef(), haplotype.getCigar(), activeRegionStop, ReadUtils.ClippingTail.RIGHT_TAIL, true );
        if( hapStop == ReadUtils.CLIPPING_GOAL_NOT_REACHED && activeRegionStop == haplotype.getAlignmentStartHapwrtRef() + haplotype.getCigar().getReferenceLength() ) {
            hapStop = activeRegionStop; // contract for getReadCoordinateForReferenceCoordinate function says that if read ends at boundary then it is outside of the clipping goal
        }
        byte[] newHaplotypeBases;
        // extend partial haplotypes to contain the full active region sequence
        if( hapStart == ReadUtils.CLIPPING_GOAL_NOT_REACHED && hapStop == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            newHaplotypeBases = ArrayUtils.addAll( ArrayUtils.addAll( ArrayUtils.subarray(ref, activeRegionStart, swConsensus.getAlignmentStart2wrt1()),
                    haplotype.getBases()),
                    ArrayUtils.subarray(ref, swConsensus.getAlignmentStart2wrt1() + swConsensus.getCigar().getReferenceLength(), activeRegionStop) );
        } else if( hapStart == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            newHaplotypeBases = ArrayUtils.addAll( ArrayUtils.subarray(ref, activeRegionStart, swConsensus.getAlignmentStart2wrt1()), ArrayUtils.subarray(haplotype.getBases(), 0, hapStop) );
        } else if( hapStop == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            newHaplotypeBases = ArrayUtils.addAll( ArrayUtils.subarray(haplotype.getBases(), hapStart, haplotype.getBases().length), ArrayUtils.subarray(ref, swConsensus.getAlignmentStart2wrt1() + swConsensus.getCigar().getReferenceLength(), activeRegionStop) );
        } else {
            newHaplotypeBases = ArrayUtils.subarray(haplotype.getBases(), hapStart, hapStop);
        }

        final Haplotype h = new Haplotype( newHaplotypeBases );
        final SWPairwiseAlignment swConsensus2 = new SWPairwiseAlignment( ref, h.getBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );

        h.setAlignmentStartHapwrtRef( swConsensus2.getAlignmentStart2wrt1() );
        if ( haplotype.isArtificialHaplotype() ) {
            h.setArtificialEvent(haplotype.getArtificialEvent());
        }
        if( swConsensus2.getCigar().toString().contains("S") || swConsensus2.getCigar().getReferenceLength() != activeRegionStop - activeRegionStart || swConsensus2.getAlignmentStart2wrt1() < 0 ) { // protect against unhelpful haplotype alignments
            return false;
        }

        h.setCigar( AlignmentUtils.leftAlignIndel(swConsensus2.getCigar(), ref, h.getBases(), swConsensus2.getAlignmentStart2wrt1(), 0, true) );

        if( FORCE_INCLUSION_FOR_GGA_MODE || !haplotypeList.contains(h) ) {
            haplotypeList.add(h);
            return true;
        } else {
            return false;
        }
    }
}