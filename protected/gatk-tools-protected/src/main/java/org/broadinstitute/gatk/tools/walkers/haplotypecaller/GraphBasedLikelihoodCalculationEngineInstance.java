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

import htsjdk.variant.variantcontext.Allele;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.genotyper.AlleleList;
import org.broadinstitute.gatk.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.Route;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.CountSet;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pairhmm.FlexibleHMM;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Fast pseudo-likelihood calculation engine based on the assembly haplotype graph.
 *
 * <p>
 *     An instance is good for active region. {@link GraphBasedLikelihoodCalculationEngine} instance them on demand
 *     as requested by the {@code HaplotypeCaller} code.
 * </p>
 */
public class GraphBasedLikelihoodCalculationEngineInstance {

    private final static Logger logger = Logger.getLogger(GraphBasedLikelihoodCalculationEngineInstance.class);


    /**
     * Unified kmer size used for the Haplotype graph.
     */
    protected final int kmerSize;

    /**
     * Reference to the haplotype graph.
     */
    protected final HaplotypeGraph haplotypeGraph;

    /**
     * Haplotypes included in the haplotype graph.
     */
    private final List<Haplotype> haplotypes;

    /**
     * Whether there is some variation present in the haplotype assembly.
     */
    private final boolean hasVariation;


    /**
     * Counts of reads that anchoread somewhere.
     *
     * <p>Used for debugging purposes</p>
     */
    private int anchoredReads = 0;

    /**
     * Count of reads that didn't anchor anywere.
     *
     * <p>Used for debugging purposes</p>
     */
    private int nonAnchoredReads = 0;

    /**
     * Pair-hmm implementation to use to calculate read likelihoods.
     */
    private final FlexibleHMM hmm;

    /**
     * Holds the log10 probability of passing from a indel to a match.
     */
    private final double indelToMatchTransitionLog10Probability;

    /**
     * Maximum likelihood difference between the reference haplotype and the best alternative haplotype.
     *
     * <p>If the difference is greater for a read, the reference haplotype likelihood is increase in order to not go
     *  beyond this limit</p>
     */
    protected final double log10globalReadMismappingRate;

    protected final EventBlockFinder eventBlockSearchEngine;


    /**
     * Constructs a new engine based on the results of the assembly.
     *
     * @param assemblyResultSet                           assembly-result set
     * @param hmm                           fast-hmm implementation to use.
     * @param log10globalReadMismappingRate maximum cost for the reference haplotype vs the best alternative available.
     * @param heterogeneousKmerSizeResolution multi-kmersize dataset resolution.
     * @throws NullPointerException     if any argument is null.
     * @throws IllegalArgumentException if log10globalReadMismappingRate >= 0.
     */
    public GraphBasedLikelihoodCalculationEngineInstance(final AssemblyResultSet assemblyResultSet, final FlexibleHMM hmm, final double log10globalReadMismappingRate, final HeterogeneousKmerSizeResolution heterogeneousKmerSizeResolution) {
        if (heterogeneousKmerSizeResolution  == null) throw new NullPointerException("the kmerSize resolution cannot be null");
        if (assemblyResultSet == null) throw new NullPointerException("the assembly result set cannot be null");
        if (hmm == null) throw new NullPointerException("the fast-hmm component cannot be null");
        if (log10globalReadMismappingRate >= 0)
            throw new IllegalArgumentException("the global reading mismapping rate cannot be positive or zero");

        this.hmm = hmm;
        this.indelToMatchTransitionLog10Probability = QualityUtils.qualToProbLog10(hmm.getGapExtensionPenalty());
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;

        haplotypes = new ArrayList<>(assemblyResultSet.getHaplotypeList());
        Collections.sort(haplotypes, Haplotype.ALPHANUMERICAL_COMPARATOR);

        // make sure that kmerSize is not bigger than the smallest haplotype. It can well happen when there are cycles and kmerSize inflates.
        final Haplotype referenceHaplotype = assemblyResultSet.getReferenceHaplotype();
        int minHaplotypeLength = referenceHaplotype.length();
        for (final Haplotype h : haplotypes)
            if (minHaplotypeLength > h.length())
                minHaplotypeLength = h.length();

        // Determine the kmerSize to use for the unified haplotype assembly graph

        kmerSize = Math.min(minHaplotypeLength,
                heterogeneousKmerSizeResolution.useMaximum() ? assemblyResultSet.getMaximumKmerSize() : assemblyResultSet.getMinimumKmerSize());

        haplotypeGraph = new HaplotypeGraph(kmerSize,haplotypes);


        if (haplotypeGraph.hasCycles())
            Utils.warnUser(logger, "cycle caused at merging haplotypes with different kmerSizes: active region " + assemblyResultSet.getRegionForGenotyping() + " will be skipped");

        //TODO haplpotypeGraph.getReferenceSourceVertex() == null
        //TODO Is a quick patch to ignore cases where the trimming has rendered kmerSize so big that is bigger than the haplotype
        //TODO and reduction to the minimum haplotype size result in no unique kmers.
        //TODO the actual solution: we need to impose a maximum trimming at least for Graph-based HC runs as we are loosing
        //TODO a bit of sensitivity as trimming results in lack of unique kmers.
        if (haplotypeGraph.hasCycles() || haplotypeGraph.getReferenceHaplotype() == null) {
            hasVariation = false;
            eventBlockSearchEngine = null;
            return;
        }

        haplotypeGraph.mergeCommonChains();
        //TODO recover dangling ends. Did not work the last time I tried but may be worth to retry.
        //haplotypeGraph.recoverDanglingTails(-1);
        logger.debug("using haplotype graph with kmerSize " + haplotypeGraph.getKmerSize());

        hasVariation = !haplotypeGraph.hasCycles() && haplotypeGraph.getHaplotypes().size() > 1;

        eventBlockSearchEngine = new EventBlockFinder(haplotypeGraph);
    }

    /**
     * Determines whether based on result from assembly and the relevant user options we can reuse th existing
     *
     * @param assemblyResultSet assembly result set.
     * @param kmerSize intended kmerSize for the haplotype graph.
     * @param heterogeneousKmerSizeResolution user instruction as to how to resolve situation where we have haplotypes comming from different kmer sizes.
     * @return {@code true} iff we can reuse an existing read-threading graph with that kmerSize in the assembly result set.
     */
    @SuppressWarnings("unused")
    private static boolean canReuseReadThreadingGraphAsHaplotypeGraph(final AssemblyResultSet assemblyResultSet, final int kmerSize, final HeterogeneousKmerSizeResolution heterogeneousKmerSizeResolution) {
        return !assemblyResultSet.wasTrimmed() && (!assemblyResultSet.hasMultipleKmerSizes() || heterogeneousKmerSizeResolution.combinesKmerSizes()) &&
                assemblyResultSet.getUniqueReadThreadingGraph(kmerSize) != null;
    }

    /**
     * Checks whether the underlying haplotype graph assembly contains any variation worth analyzing.
     *
     * @return {@code true} iff so.
     */
    @SuppressWarnings("unused")
    public boolean hasVariation() {
        return hasVariation;
    }

    /**
     * Calculates the likelihood of reads across many samples evaluated against haplotypes resulting from the
     * active region assembly process.
     *
     * @param haplotypes to evaluate.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @throws NullPointerException if either parameter is {@code null}.
     *
     * @return never {@code null}, and with at least one entry for input sample (keys in {@code perSampleReadList}.
     *    The value maps can be potentially empty though.
     */
    public ReadLikelihoods<Haplotype> computeReadLikelihoods(final List<Haplotype> haplotypes, final SampleList samples,
            final Map<String, List<GATKSAMRecord>> perSampleReadList) {
        // General preparation on the input haplotypes:
        final List<Haplotype> sortedHaplotypes = new ArrayList<>(haplotypes);
        Collections.sort(sortedHaplotypes, Haplotype.ALPHANUMERICAL_COMPARATOR);
        final AlleleList<Haplotype> alleles = new IndexedAlleleList<>(sortedHaplotypes);
        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, alleles, perSampleReadList);

        // The actual work:
        final int sampleCount = result.sampleCount();
        for (int s = 0; s < sampleCount; s++) {
            final List<GATKSAMRecord> sampleReads = result.sampleReads(s);

            // Get the cost/likelihood of each read at relevant subpaths on the tree:
            final Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> costsByEndingVertex = calculatePathCostsByRead(sampleReads);
            // Create the resulting per-read maps:
            calculatePerReadAlleleLikelihoodMap(costsByEndingVertex,result.sampleMatrix(s) );
        }
        result.normalizeLikelihoods(true,log10globalReadMismappingRate);
        logger.debug("Likelihood analysis summary: reads anchored " + anchoredReads + "/" + (anchoredReads + nonAnchoredReads) + "");
        return result;
    }


    /**
     * Prints a graph into a dot file.
     *
     * @param fileName name of the output file.
     */
    public void printGraph(final String fileName) {
        haplotypeGraph.printGraph(fileName);
    }

    /**
     * Returns the kmerSize the engine is using to match read vs graph kmers thus reducing computation.
     *
     * @return greater than 0.
     */
    public int getKmerSize() {
        return kmerSize;
    }

    /**
     * Tells whether the underlying haplotype graph contained cycles.
     *
     * @return {@code true} iff so.
     */
    @SuppressWarnings("unused")
    public boolean hasCycles() {
        return haplotypeGraph.hasCycles();
    }


    /**
     * Builds the result per-read allele likelihood map.
     *
     * @param costsEndingByVertex Read vs haplotype graph sub-paths cost indexed by ending vertex.
     * @param likelihoods matrix where to set the likelihoods where the first index in the haplotype's and the second
     *                    the read.
     */
    protected void calculatePerReadAlleleLikelihoodMap(final Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> costsEndingByVertex,
                                                       final ReadLikelihoods.Matrix<Haplotype> likelihoods) {
        final int alleleCount = likelihoods.alleleCount();
        for (int h = 0; h < alleleCount; h++)
            calculatePerReadAlleleLikelihoodMapHaplotypeProcessing(h, likelihoods, costsEndingByVertex);
    }

    /**
     * Work done per haplotype to build the result per-read allele likelihood map.
     * <p/>
     * <p>
     * Basically for each haplotype we go through its path in the graph collecting all the read cost that we find
     * on the way. For each read present we add up all its cost resulting in a single value per read, i.e. its
     * "likelihood".
     * </p>
     *
     * @param haplotypeIndex                 the target haplotype index in the {@code likelihoods} matrix.
     * @param likelihoods               matrix of likelihoods.
     * @param costsEndingByVertex       read costs assorted by their end vertex.
     */
    private void calculatePerReadAlleleLikelihoodMapHaplotypeProcessing(final int haplotypeIndex,
                                                                        final ReadLikelihoods.Matrix<Haplotype> likelihoods,
                                                                        final Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> costsEndingByVertex) {
        final Haplotype haplotype = likelihoods.alleleAt(haplotypeIndex);
        final HaplotypeRoute haplotypeRoute = haplotypeGraph.getHaplotypeRoute(haplotype);
        final Set<MultiDeBruijnVertex> haplotypeVertices = haplotypeRoute.vertexSet();
        final Map<GATKSAMRecord, ReadCost> readCostByRead = new HashMap<>();
        final Set<MultiDeBruijnVertex> visitedVertices = new HashSet<>(haplotypeVertices.size());
        final List<MultiSampleEdge> edgeList = haplotypeRoute.getEdges();

        MultiDeBruijnVertex currentVertex = haplotypeRoute.getFirstVertex();
        Route<MultiDeBruijnVertex, MultiSampleEdge> pathSoFar = new Route<>(currentVertex, haplotypeGraph);
        final Iterator<MultiSampleEdge> edgeIterator = edgeList.iterator();

        while (true) {
            visitedVertices.add(currentVertex);
            final Set<ReadSegmentCost> finishingAtElementCostSet = costsEndingByVertex.get(currentVertex);
            updateReadCosts(readCostByRead, visitedVertices, pathSoFar, finishingAtElementCostSet);
            if (!edgeIterator.hasNext()) break;
            final MultiSampleEdge nextEdge = edgeIterator.next();
            pathSoFar = new Route<>(pathSoFar, nextEdge);
            currentVertex = pathSoFar.getLastVertex();
        }

        int readIndex = 0;
        for (final GATKSAMRecord read : likelihoods.reads()) {
            final ReadCost rc = readCostByRead.get(read);
            //if (rc != null)
            likelihoods.set(haplotypeIndex,readIndex,rc == null ? Double.NEGATIVE_INFINITY : rc.getCost());
            readIndex++;
        }
    }

    /**
     * Update the read cost based on the path cost found at a vertex.
     *
     * @param readCosts                 collection of read costs so far
     * @param visitedVertices           visited vertices collection.
     * @param pathSoFar                 the haplotype path visited so far.
     * @param finishingAtElementCostSet collection of path cost to process
     */
    private void updateReadCosts(final Map<GATKSAMRecord, ReadCost> readCosts,
                                 final Set<MultiDeBruijnVertex> visitedVertices,
                                 final Route<MultiDeBruijnVertex, MultiSampleEdge> pathSoFar,
                                 final Set<ReadSegmentCost> finishingAtElementCostSet) {
        if (finishingAtElementCostSet != null) {
            for (final ReadSegmentCost pc : finishingAtElementCostSet) {
                if (!visitedVertices.contains(pc.path.getFirstVertex()))
                    continue;
                if (!pathSoFar.isSuffix(pc.path))
                    continue;
                ReadCost rc = readCosts.get(pc.read);
                if (rc == null)
                    readCosts.put(pc.read, rc = new ReadCost(pc.read,indelToMatchTransitionLog10Probability));
                rc.addCost(pc.getCost());
            }
        }
    }

    /**
     * Likelihood penalty for unreported haplotype vs read likelihood with respect to the worst reported one.
     */
    private static final int UNREPORTED_HAPLOTYPE_LIKELIHOOD_PENALTY = -3;

    /**
     * Re-scales all haplotype vs read likelihoods so that for read, the best haplotype, hash likelihood 0.
     *
     * @param alleleVersions map between input haplotypes and output alleles.
     * @param result where to change the likelihoods.
     * @param mayNeedAdjustment set of read that might need adjustment. Others might be ignored.
     * @param maxAlternative map from each read and the maximum alternative haplotype likelihood.
     */
    @SuppressWarnings("unused")
    private void makeLikelihoodAdjustment(final Map<Haplotype, Allele> alleleVersions,
                                          final PerReadAlleleLikelihoodMap result,
                                          final Set<GATKSAMRecord> mayNeedAdjustment,
                                          final Map<GATKSAMRecord, Double> maxAlternative) {
        final Map<GATKSAMRecord, Map<Allele, Double>> map = result.getLikelihoodReadMap();

        for (final GATKSAMRecord read : mayNeedAdjustment) {
            final Map<Allele, Double> existingLikelihoods = map.get(read);
            if (existingLikelihoods != null) {
                double worstRelativeLikelihood = 0;
                double bestRelativeLikelihood = Double.NEGATIVE_INFINITY;
                for (final Map.Entry<Allele, Double> entry : map.get(read).entrySet()) {
                    final double candidateRelativeLikelihood = entry.getValue();
                    if (candidateRelativeLikelihood > bestRelativeLikelihood) {
                        bestRelativeLikelihood = candidateRelativeLikelihood;
                    }
                    if (!Double.isInfinite(candidateRelativeLikelihood) && worstRelativeLikelihood > candidateRelativeLikelihood)
                        worstRelativeLikelihood = candidateRelativeLikelihood;
                }
                if (Double.isInfinite(bestRelativeLikelihood)) {
                    bestRelativeLikelihood = 0;
                    worstRelativeLikelihood = 0;
                } else {
                    worstRelativeLikelihood += UNREPORTED_HAPLOTYPE_LIKELIHOOD_PENALTY;
                }
                final double bestLikelihood = 0.0; // the best becomes zero.
                maxAlternative.put(read, bestLikelihood);
                for (final Map.Entry<Haplotype, Allele> entry : alleleVersions.entrySet()) {
                    final Allele a = entry.getValue();
                    final Double relativeLikelihoodO = existingLikelihoods.get(a);
                    final double relativeLikelihood = relativeLikelihoodO == null || Double.isInfinite(relativeLikelihoodO) ? worstRelativeLikelihood : relativeLikelihoodO;
                    final double likelihood = relativeLikelihood - bestRelativeLikelihood + bestLikelihood;
                    if (likelihood > 0)
                        throw new IllegalStateException("Likelihood larger than 1 with read " + read.getReadName());
                    existingLikelihoods.put(a, likelihood);
                }
            }
        }
    }

    /**
     * Calculates path costs for a set of reads.
     * <p/>
     * <p>
     * The resulting map has one entry per read, where the read is the key and the value list of path-cost sets.
     * Each element in that list corresponds to an event block. Each path cost in one of those sets indicate the
     * likelihood (cost) of traversing a possible path across the event block using that read.
     * </p>
     *
     * @param reads reads to analyze.
     * @return never {@code null}.
     */
    protected Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> calculatePathCostsByRead(
            final List<GATKSAMRecord> reads) {
        final Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> result = new HashMap<>(reads.size());
        if (!hasVariation)
            return Collections.emptyMap();
        for (final GATKSAMRecord r : reads) {
            calculatePathCostsByRead(r, result);
        }
        return result;
    }

    /**
     * Calculates path cost for a single read.
     *
     * @param read              target read.
     * @param result       map where to add the result.
     */
    private void calculatePathCostsByRead(final GATKSAMRecord read,
                                          final Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> result) {

        final ReadAnchoring anchoring = new ReadAnchoring(read,haplotypeGraph);
        // cannot anchor so go the tradition pair-hmm way.
        hmm.loadRead(read);
        if (!anchoring.isAnchoredSomewhere()) {
            defaultToRegularPairHMM(anchoring, result);
            nonAnchoredReads++;
        } else {
            calculateReadSegmentCosts(anchoring, hmm, result);
            if (!anchoring.isPerfectAnchoring()) danglingEndPathCosts(anchoring, hmm, result);
            anchoredReads++;
        }
    }

    /**
     * Calculates read vs haplotype likelihoods using the classic PairHMM approach.
     * <p/>
     * <p>
     * It basically compares the read with each haplotype full path without short cuts.
     * </p>
     *
     * @param anchoring      anchoring information of the read.
     * @param destination where to leave the results indexed by ending veretex.
     */
    private void defaultToRegularPairHMM(final ReadAnchoring anchoring, final Map<MultiDeBruijnVertex,Set<ReadSegmentCost>> destination) {

        for (final Map.Entry<Haplotype, HaplotypeRoute> entry : haplotypeGraph.getHaplotypeRouteMap().entrySet()) {
            if (entry.getValue() == null) continue;
            final byte[] haplotypeBases = entry.getKey().getBases();
            hmm.loadHaplotypeBases(haplotypeBases);
            final double cost = hmm.calculateLocalLikelihood(0, anchoring.read.getReadLength(), 0, haplotypeBases.length, false);
            final ReadSegmentCost readSegmentCost = new ReadSegmentCost(anchoring.read, entry.getValue(), cost);
            addReadSegmentCost(destination, readSegmentCost);
        }
    }

    /**
     * Add a new read-segment-cost to an ending vertex indexed map.
     * @param destination where to add the read-segment-cost.
     * @param cost the read-segment-cost to add.
     */
    private void addReadSegmentCost(final Map<MultiDeBruijnVertex,Set<ReadSegmentCost>> destination, final ReadSegmentCost cost) {
        final MultiDeBruijnVertex endVertex = cost.path.getLastVertex();
        Set<ReadSegmentCost> vpcSet = destination.get(endVertex);
        if (vpcSet == null)
            destination.put(endVertex, vpcSet = new HashSet<>(10));
        vpcSet.add(cost);
    }

    /**
     * Calculate the likelihood cost of path section of a read across the graph.
     * <p/>
     * <p>
     * Given a read, its anchors and other unique kmer mapable to the reference path we can divide the graph
     * into event blocks: a set of one or more variations and the possible path across that block.
     * </p>
     * <p/>
     * <p>
     * The result value will have  one element fo reach block. Each element is the set of all path costs (likelihoods)
     * to traverse the block using all possible paths (different haplotypes).
     * </p>
     * <p/>
     * <p>
     * The current implementation has some added complexity in order to avoid a situation in where the last part
     * of the anchored section of the read is thrown out. We first determine the last event block boundaries and we
     * make sure that we won't run over its left limit when covering for earlier event blocks.
     * </p>
     *
     * @param anchoring        target read graph anchoring information.
     * @param hmm                       the pair-hmm calculation engine. It must have been loaded with the same {@code read} already.
     * @param destination where to add the costs.
     */
    private void calculateReadSegmentCosts(final ReadAnchoring anchoring, final FlexibleHMM hmm, final Map<MultiDeBruijnVertex, Set<ReadSegmentCost>> destination) {

        final EventBlockFinder.Traversal traversal = eventBlockSearchEngine.traversal(anchoring);

        for (final EventBlock eventBlock : traversal) {

          //  final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> acrossBlockPaths =
          //          calculateAllPathsBetweenVertices(anchoring,
          //                  eventBlock.getSource(), eventBlock.getSink());//eventBlock.getRoutesAcross();

            final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> acrossBlockPaths = eventBlock.getRoutesAcross();

            int leftBlockBoundaryIndex = anchoring.uniqueKmerOffsets.get(eventBlock.getSource());
            int rightBlockBoundaryIndex = anchoring.uniqueKmerOffsets.get(eventBlock.getSink());
            calculateCostForPathSet(anchoring.read, acrossBlockPaths, hmm, leftBlockBoundaryIndex, rightBlockBoundaryIndex, true, false, null, null, destination);
        }
    }

    /**
     * Calculate path cost for a set of paths across a event block.
     *
     * @param read the target read.
     * @param acrossBlockPaths event block paths to evaluate.
     * @param hmm pair-hmm engine to use to calculate likelihoods.
     * @param beforeBlockReadOffset kmer offset on the read for the vertex kmer before the block.
     * @param afterBlockReadOffset  kmer offset on the read for the vertex kmer after the block.
     * @param doClipping whether to perform any clipping in order to save cpu time.
     * @param prependVertex if not null, the end cost path with be prepended with this vertex.
     * @param appendVertex if not null, the end cost path will be appended with this vertex.
     * @param includePathEnds whether to include or exclude the vertices at the very end or beginning of the paths.
     */
    private void calculateCostForPathSet(
            final GATKSAMRecord read, final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> acrossBlockPaths,
            final FlexibleHMM hmm, final int beforeBlockReadOffset, final int afterBlockReadOffset,
            final boolean doClipping, final boolean includePathEnds,
            final MultiDeBruijnVertex prependVertex,
            final MultiDeBruijnVertex appendVertex,
            final Map<MultiDeBruijnVertex,Set<ReadSegmentCost>> destination) {


        final Set<ReadSegmentCost> readSegmentCosts = new TreeSet<>(ReadSegmentComparator.INSTANCE);

        final int readStart = beforeBlockReadOffset + kmerSize;
        final int readEnd = Math.max(readStart, afterBlockReadOffset + kmerSize - 1);
        final byte[][] pathBases = new byte[acrossBlockPaths.size()][];
        final CountSet pathSizes = new CountSet(acrossBlockPaths.size());
        int nextPath = 0;

        // Complete the read segment cost with the corresponding path bases
        for (final Route<MultiDeBruijnVertex, MultiSampleEdge> p : acrossBlockPaths) {
            final ReadSegmentCost readSegmentCost = new ReadSegmentCost(read, p, 0);
            pathBases[nextPath++] = readSegmentCost.bases = eventBlockPathBases(p, includePathEnds);
            pathSizes.add(readSegmentCost.bases.length);
            readSegmentCosts.add(readSegmentCost);
        }

        // Add the read 'path size'.
        pathSizes.add(readEnd - readStart);

        final byte[] readBases = hmm.getReadBases();

        // Perform right clipping of bases that are common to all paths and read.
        int rightClipping = !doClipping ? 0 : calculateRightClipping(readEnd, pathBases, readBases,pathSizes);

        // Calculate the costs.
        for (final ReadSegmentCost readSegmentCost : readSegmentCosts) {
            hmm.loadHaplotypeBases(readSegmentCost.bases);
            readSegmentCost.setCost(hmm.calculateLocalLikelihood(Math.max(0, readStart), readEnd - rightClipping, 0, readSegmentCost.bases.length - rightClipping, false));
            if (prependVertex != null)
                readSegmentCost.path = new Route<>(prependVertex,readSegmentCost.path);
            if (appendVertex != null)
                readSegmentCost.path = new Route<>(readSegmentCost.path,appendVertex);
            addReadSegmentCost(destination,readSegmentCost);
        }
    }

    /**
     * Determines how much we can clip away from the right side of a set of path without loosing accuracy when comparing
     * likelihood vs the read.
     *
     * @param readEnd     exclusive position right after the last one of the region considered.
     * @param pathBases   bases of possible path in the same event block.
     * @param readBases   full length read bases.
     * @param pathSizes   path size set.
     *
     * @return 0 or greater.
     */
    private int calculateRightClipping(final int readEnd, final byte[][] pathBases,
                                       final byte[] readBases, final CountSet pathSizes) {
        final int maxClipping = pathSizes.size() > 1 ? 0 : Math.min(pathSizes.min(), kmerSize - 1);
        int rightClipping = 0;
        while (rightClipping < maxClipping) {
            final byte readBase = readBases[readEnd - rightClipping - 1];
            boolean dontGoFurther = false;
            for (int i = 0; !dontGoFurther && i < pathBases.length; i++)
                if (pathBases[i][pathBases[i].length - rightClipping - 1] != readBase)
                    dontGoFurther = true;
            if (dontGoFurther)
                break;
            rightClipping++;
        }
        return rightClipping;
    }

    /**
     * Calculates a graph path bases.
     * <p/>
     * <p>
     * When the path starts on a source vertex, all its sequence is considered as part of the path bases. For regular
     * vertices start only the suffix (last) base is considered.
     * </p>
     *
     * @param path            the targeted path.
     * @param includePathEnds whether the bases included in the first and last vertex of the path should be included or excluded.
     * @return never {@code null} but perhaps a zero-length base array if the final requested path length is zero.
     */
    //TODO this method could be moved to the Path class, but require consider how to make the API more concise.
    private byte[] eventBlockPathBases(final Path<MultiDeBruijnVertex, MultiSampleEdge> path,
                                       final boolean includePathEnds) {
        // We first calculate the size of the return.
        final List<MultiDeBruijnVertex> vertices = path.getVertices();
        final boolean pathStartsAtSource = haplotypeGraph.isSource(path.getFirstVertex());
        final int resultLength = includePathEnds
                ? vertices.size() + (pathStartsAtSource ? path.getFirstVertex().getSequence().length - 1 : 0)
                : vertices.size() - 2;
        // Trivial empty return cases:
        if (resultLength <= 0)
            return new byte[0];
        final byte[] result = new byte[resultLength];
        if (result.length == 0) {
            return result;
        }
        // General return cases:
        final ListIterator<MultiDeBruijnVertex> it = vertices.listIterator(includePathEnds ? 0 : 1); // skip the vertex (exclusive)
        for (int i = 0; i < resultLength; i++) { // i < resultLength implicitly skips the last vertex (exclusive).
            final MultiDeBruijnVertex vertex = it.next();
            if (i == 0 && includePathEnds && pathStartsAtSource) {
                System.arraycopy(vertex.getSequence(), 0, result, 0, kmerSize);
                i = kmerSize - 1;
            } else
                result[i] = vertex.getSuffix();
        }
        return result;
    }

    /**
     * Calculate the path cost of dangling ends.
     * <p/>
     * <p>
     * A dangling end is the section of the read that falls before the left anchor or after the right anchor.
     * </p>
     *
     * @param anchoring anchoring information of the read vs the haplotype assembly graph.
     * @param hmm       the PairHMM engine to use to calculate likelihoods.
     * @param destination cost destination.
     */
    private void danglingEndPathCosts(final ReadAnchoring anchoring, final FlexibleHMM hmm, final Map<MultiDeBruijnVertex,Set<ReadSegmentCost>> destination) {
        if (anchoring.leftAnchorIndex > 0 || anchoring.leftAnchorIndex == 0
                && anchoring.leftAnchorVertex.hasAmbiguousSequence())
            leftDanglingEndPathCosts(anchoring, hmm,destination);

        if (anchoring.rightAnchorIndex < anchoring.read.getReadLength() - kmerSize)
            rightDanglingEndPathCosts(anchoring, hmm, destination);
    }

    /**
     * Generates all relevant right dangling end path costs.
     *
     * @param anchoring the anchoring information for the read under analysis.
     * @param hmm       pair-hmm implementation to use to calculate likelihoods. It is assumed to be loaded with
     *                  the same read as {@code anchoring} refers to.
     * @param destination where the place the resulting read-segment-costs.
     */
    private void rightDanglingEndPathCosts(final ReadAnchoring anchoring, final FlexibleHMM hmm,
                                           final Map<MultiDeBruijnVertex,Set<ReadSegmentCost>> destination) {
        final int readStart = anchoring.rightAnchorIndex;
        final int readEnd = anchoring.read.getReadLength() - kmerSize + 1;
        final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> haplotypeRoutes =
                extendsHaplotypeRoutesForwards(anchoring.rightAnchorVertex);
        if (haplotypeRoutes.size() >= 2)
            calculateCostForPathSet(anchoring.read,
                    haplotypeRoutes, hmm, readStart, readEnd, false, true,anchoring.rightAnchorVertex,null,destination);

    }

    /**
     * Generates all relevant left dangling end path costs.
     *
     * @param anchoring the anchoring information for the read under analysis.
     * @param hmm       pair-hmm implementation to use to calculate likelihoods. It is assumed to be loaded with
     *                  the same read as {@code anchoring} refers to.
     * @param destination where the place the resulting read-segment-costs.
     */
    private void leftDanglingEndPathCosts(final ReadAnchoring anchoring, final FlexibleHMM hmm,
                                          final Map<MultiDeBruijnVertex,Set<ReadSegmentCost>> destination) {
        final int readStart = -kmerSize;
        final int readEnd = anchoring.leftAnchorIndex;
        final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> haplotypeRoutes =
                extendsHaplotypeRoutesBackwards(anchoring.leftAnchorVertex);
        if (haplotypeRoutes.size() >= 2) // if there is just one haplotype route there is no relevant variation in the dangling end.
            calculateCostForPathSet(anchoring.read, haplotypeRoutes, hmm,
                    readStart, readEnd, false, true, null, anchoring.leftAnchorVertex, destination);
    }

    /**
     * Construct haplotype routes prefixes to an anchor vertex.
     * <p/>
     * <p>
     * The output should contain a route for each haplotype that includes the input anchor vertex.
     * This route would be the prefix of the haplotype that finishes at that vertex.
     * </p>
     *
     * @param anchorVertex the target anchor vertex.
     * @return never {@code null}.
     */
    private Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> extendsHaplotypeRoutesBackwards(
            final MultiDeBruijnVertex anchorVertex) {
        final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> result = new HashSet<>(haplotypes.size());
        for (final MultiDeBruijnVertex parent : haplotypeGraph.incomingVerticesOf(anchorVertex))
            extendsHaplotypeRoutesFrom(parent, result, false);
        return result;
    }

    /**
     * Construct haplotype routes suffix from an anchor vertex.
     * <p/>
     * <p>
     * The output should contain a route for each haplotype that includes the input anchor vertex.
     * This route would be the suffix of the haplotype that starts at that vertex.
     * </p>
     *
     * @param anchorVertex the target anchor vertex.
     * @return never {@code null}.
     */
    private Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> extendsHaplotypeRoutesForwards(
            final MultiDeBruijnVertex anchorVertex) {
        final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> result = new HashSet<>(haplotypes.size());
        for (final MultiDeBruijnVertex parent : haplotypeGraph.outgoingVerticesOf(anchorVertex))
            extendsHaplotypeRoutesFrom(parent, result, true);
        return result;
    }

    /**
     * Extends from a vertex considering path furcations that are part of some valid haplotype
     * <p/>
     * <p>
     * In other words, it will ignore subpaths that are not valid part of an assembled haplotype.
     * </p>
     *
     * @param start   start seed vertex.
     * @param result  destination for found extensions.
     * @param forward whether to traverse edges forward or backwards.
     */
    private void extendsHaplotypeRoutesFrom(final MultiDeBruijnVertex start, final Set<Route<MultiDeBruijnVertex, MultiSampleEdge>> result, final boolean forward) {
        final Set<HaplotypeRoute> validHaplotypeRoutes = haplotypeGraph.getEnclosingHaplotypeRoutes(start);
        if (validHaplotypeRoutes.size() == 0) return;
        final Deque<Pair<Route<MultiDeBruijnVertex, MultiSampleEdge>, Set<HaplotypeRoute>>> queue = new LinkedList<>();
        queue.add(new Pair<>(new Route<>(start, haplotypeGraph), validHaplotypeRoutes));
        while (!queue.isEmpty()) {
            final Pair<Route<MultiDeBruijnVertex, MultiSampleEdge>, Set<HaplotypeRoute>> current = queue.remove();
            final Route<MultiDeBruijnVertex, MultiSampleEdge> path = current.getFirst();
            final MultiDeBruijnVertex vertex = forward ? path.getLastVertex() : path.getFirstVertex();
            final Set<HaplotypeRoute> validRoutes = current.getSecond();
            for (final HaplotypeRoute hr : validRoutes) {
                final MultiDeBruijnVertex routeEndVertex = forward ? hr.getLastVertex() : hr.getFirstVertex();
                if (vertex.equals(routeEndVertex)) {
                    result.add(path);
                    break;
                }
            }
            final Set<MultiDeBruijnVertex> nextVertices = forward ? haplotypeGraph.outgoingVerticesOf(vertex) :
                    haplotypeGraph.incomingVerticesOf(vertex);
            for (final MultiDeBruijnVertex candidate : nextVertices) {
                extendsHaplotypeRoutesFrom$ProcessCandidateExtendingVertex(forward, queue, path, validRoutes, candidate);
            }
        }
    }

    /**
     * Check on an candidate vertice to exted a path.
     *
     * <p>
     *     This method updates the traversal queue accordingly.
     * </p>
     *
     * @param forward whether the extension is forward, or backwards.
     * @param queue queue with open paths yet to be explored.
     * @param path path extension to evaluate.
     * @param validRoutes collection of valid haplotype routes used to discard non-informative extensions.
     * @param candidate the candidate extending vertex.
     */
    private void extendsHaplotypeRoutesFrom$ProcessCandidateExtendingVertex(
            final boolean forward,
            final Deque<Pair<Route<MultiDeBruijnVertex, MultiSampleEdge>, Set<HaplotypeRoute>>> queue,
            final Route<MultiDeBruijnVertex, MultiSampleEdge> path,
            final Set<HaplotypeRoute> validRoutes, final MultiDeBruijnVertex candidate) {
        final Set<HaplotypeRoute> parentValidHaplotypes = haplotypeGraph.getEnclosingHaplotypeRoutes(candidate);
        switch (parentValidHaplotypes.size()) {
            case 0:
                return;
            case 1:
                if (validRoutes.containsAll(parentValidHaplotypes))
                    queue.add(new Pair<>(forward ? new Route<>(path, candidate) : new Route<>(candidate, path), parentValidHaplotypes));
                else
                    return;
                break;
            default:
                if (parentValidHaplotypes.size() == validRoutes.size() && parentValidHaplotypes.containsAll(validRoutes)) {
                    queue.add(new Pair<>(forward ? new Route<>(path, candidate) : new Route<>(candidate, path), parentValidHaplotypes));
                } else {
                    final Set<HaplotypeRoute> newValidHaplotypeRoutes = new HashSet<>(validRoutes.size());
                    for (final HaplotypeRoute hr : validRoutes)
                        if (parentValidHaplotypes.contains(hr))
                            newValidHaplotypeRoutes.add(hr);
                    if (newValidHaplotypeRoutes.size() == 0)
                        return;
                    queue.add(new Pair<>(forward ? new Route<>(path, candidate) : new Route<>(candidate, path), newValidHaplotypeRoutes));
                }
        }
    }

    public List<Haplotype> getHaplotypeList() {
        return new ArrayList<>(haplotypeGraph.getHaplotypes());
    }

    /**
     * Returns the haplotype graph associated with this instance.
     * @return never {@code null}
     */
    public HaplotypeGraph getHaplotypeGraph() {
        return haplotypeGraph;
    }
}