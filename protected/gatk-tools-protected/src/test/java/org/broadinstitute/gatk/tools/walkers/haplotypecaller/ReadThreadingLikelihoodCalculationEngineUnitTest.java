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
import org.broadinstitute.gatk.utils.genotyper.SampleListUtils;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pairhmm.*;
import org.broadinstitute.gatk.utils.sam.ClippedGATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: valentin
 * Date: 8/4/13
 * Time: 10:20 PM
 * To change this template use File | Settings | File Templates.
 */
@Test(enabled=false)
public class ReadThreadingLikelihoodCalculationEngineUnitTest extends ActiveRegionTestDataSetUnitTest {


//    private static FastHMM hmm = new MLLog10PairHMM((byte)10); // new FastLoglessPairHMM((byte)10);

    private static FlexibleHMM hmm = new FastLoglessPairHMM((byte)10);

    @Test(dataProvider="activeRegionTestDataSets",enabled=false)
    public void testActiveRegionsDataSet(final ActiveRegionTestDataSet as, final int kmerSize, final int readLength, final String variation, final int readCount, final int regionSize, final byte bq, final byte iq, final byte dq) {
        super.testActiveRegionsDataSet(as,kmerSize,readLength,variation,readCount,regionSize,bq,iq,dq);
    }
        /** How many missing read record are tolerated in the graph based approach. For example a read is missed
         * if it does not map to the reference path with at least two kmers in non overlapping positions. This constant
         * indictes the proportion of reads reacords that we can miss with respect to all possible
         */
    private static final double READ_SKIP_TOLERANCE = 0.01;

    //final PairHMMLikelihoodCalculationEngine fullPairHMM = new PairHMMLikelihoodCalculationEngine((byte)10, false,
    //        PairHMM.HMM_IMPLEMENTATION.LOGLESS_CACHING, -3);
    final PairHMMLikelihoodCalculationEngine fullPairHMM = new PairHMMLikelihoodCalculationEngine((byte)10,
            PairHMM.HMM_IMPLEMENTATION.LOGLESS_CACHING, PairHMM.HMM_SUB_IMPLEMENTATION.UNVECTORIZED, true, Double.NEGATIVE_INFINITY,
            true, PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.NONE);

    // When using likelihoods it should be around 0.05 since
    // When using maximum-likelihoods it can be as low as 0.00001
    private static final double SIGNIFICANT_LnLK_RATIO_DIFF_FRACTION = hmm instanceof FastLoglessPairHMM ? 0.1 : 0.00001;


    // Some case is kind of expected to have differences between PairHMM and GraphBased Flexible PairHMM.
    // Is therefore difficult to test for the to give similar results in a unit test....
    // This is left for for example Integration tests like GraphBasedVsLoglessAccuracyIntegrationTest.
    // There code herein is maintain around for historical purposes, but disabled.
    @Test(dataProvider="readLikekihoodRatioTestData",enabled=false)
    public void testReadLikelihoodRatios(final ActiveRegionTestDataSet ds, final GATKSAMRecord read, final Allele a1,
                                         final Allele a2, final PerReadAlleleLikelihoodMap loglessLks,
                                         final PerReadAlleleLikelihoodMap graphLks, final List<Civar.ElementOffset> readEventOffsets, final List<Civar.ElementOffset> firstAlleleCivar, final List<Civar.ElementOffset> secondAlleleCivar ) {

        checkForLongEventsThatMightCauseFailures(readEventOffsets, firstAlleleCivar, secondAlleleCivar);
        final Map<Allele,Double> logless = loglessLks.getLikelihoodReadMap().get(read);
        final Map<Allele,Double> graph = graphLks.getLikelihoodReadMap().get(read);
        final double loglessA1Lk = logless.get(a1);
        final double loglessA2Lk = logless.get(a2);
        if (graph == null)
            throw new SkipException("no likelihoods produced for this read using the graph method: Lla1= " + loglessA1Lk + " Lla2= " + loglessA2Lk + "LlDiff=" + (loglessA2Lk - loglessA1Lk)  );

        final Double graphA1Lk = graph.get(a1);
        final Double graphA2Lk = graph.get(a2);
        if (graphA1Lk == null)
            throw new SkipException("no likelihoods produced for this read in the first haplotype: Lla1= " + loglessA1Lk + " Lla2= " + loglessA2Lk + "LlDiff=" + (loglessA2Lk - loglessA1Lk)  );
        if (graphA2Lk == null)
            throw new SkipException("no likelihoods produced for this read in the second haplotype:  Lla1= " + loglessA1Lk + " Lla2= " + loglessA2Lk + "LlDiff=" + (loglessA2Lk - loglessA1Lk)  );

        final double loglessDiff = loglessA1Lk - loglessA2Lk;
        final double graphDiff = graphA1Lk - graphA2Lk;
        final double epsilon = calculateEpsilon(graphDiff,loglessDiff);
        Assert.assertEquals(graphDiff,loglessDiff,epsilon,String.format("Delta(%f,%f) = %f > %f",graphDiff,loglessDiff,Math.abs(graphDiff - loglessDiff),epsilon));

    }

    private double calculateEpsilon(final double graphDiff, final double loglessDiff) {
        if (hmm instanceof FastLoglessPairHMM)
            return Math.max(0.01,Math.max(Math.abs(loglessDiff),Math.abs(graphDiff)) * SIGNIFICANT_LnLK_RATIO_DIFF_FRACTION);
        else
            return SIGNIFICANT_LnLK_RATIO_DIFF_FRACTION;
    }

    private static final double MIN_READ_ACROSS_SIZE_FOR_INDEL_EVENTS = 0.8; // 50% percent.
    private static final double MIN_LARGE_INDEL = 4;

    private void checkForLongEventsThatMightCauseFailures(final List<Civar.ElementOffset> read, final List<Civar.ElementOffset> a1, final List<Civar.ElementOffset> a2) {

        int sequenceLength = Math.max(a1.get(a1.size() - 1).templateTo, a2.get(a2.size() - 1).templateTo) + 1;

        boolean tai1 = thereAreIndels(a1);
        boolean tai2 = thereAreIndels(a2);
        boolean tair = thereAreIndels(read);
        boolean thereAreIndels = tai1 || tai2 || tair;
        if (!thereAreIndels) return;

        final boolean[] inserts = new boolean[sequenceLength];
        final boolean[] deletions = new boolean[sequenceLength];
        final int[] range = new int[2];

        int refStart = Integer.MAX_VALUE;
        int refEnd = -1;

        for (final Civar.ElementOffset ce : read) {
            if (refStart > ce.templateFrom)
                refStart = ce.templateFrom;
            if (refEnd < ce.templateTo)
                refEnd = ce.templateTo;
            switch (ce.element.operator()) {
                case DELETION:
                    deletions[ce.templateFrom] = deletions[ce.templateTo] = true;
                    break;
                case INSERTION:
                    inserts[ce.templateFrom] = inserts[ce.templateTo] = true;
                    break;
                case MATCH:
                    break;
            }
        }

        range[0] = refStart;
        range[1] = refEnd;

        checkForLongEventsThatMightCauseFailures_allele(refStart,refEnd,inserts,deletions,a1);
        checkForLongEventsThatMightCauseFailures_allele(refStart,refEnd,inserts,deletions,a2);
    }

    private void checkForLongEventsThatMightCauseFailures_allele(final int refStart, final int refEnd, final boolean[] inserts, final boolean[] deletions, final List<Civar.ElementOffset> a1) {
        for (final Civar.ElementOffset ce : a1) {
            if (ce.templateFrom <= refStart) continue;
            if (ce.templateTo >= refEnd) continue;
            int size;
            switch (ce.element.operator()) {
                case DELETION:
                    size = ce.templateTo - ce.templateFrom;
                    if (deletions[ce.templateFrom] || deletions[ce.templateTo]) continue;
                    break;
                case INSERTION:
                    size = ce.sequenceTo - ce.sequenceFrom;
                    if (inserts[ce.templateFrom] || inserts[ce.templateTo]) continue;
                    break;
                default:
                    continue;
            }
            int minMargin = (int) Math.ceil(size * MIN_READ_ACROSS_SIZE_FOR_INDEL_EVENTS);
            if (ce.templateFrom - refStart < minMargin)
                throw new SkipException("Large Indel");
            if (refEnd - ce.templateTo < minMargin)
                throw new SkipException("Large Indel");
        }

    }

    private boolean thereAreIndels(final List<Civar.ElementOffset> a1) {
        for (final Civar.ElementOffset ce : a1)  {
             switch (ce.element.operator()) {
                 case DELETION:
                     if (ce.templateTo - ce.templateFrom >= MIN_LARGE_INDEL) return true;
                     break;
                 case INSERTION:
                     if (ce.sequenceTo - ce.sequenceFrom >= MIN_LARGE_INDEL) return true;
                     break;
             }
        }
        return false;
    }


    @DataProvider(name="readLikekihoodRatioTestData")
    public Iterator<Object[]> readLikelihoodRatioTestDataSets() {
           final Iterator<Object[]> activeRegionTestDataSetIterator = super.activeRegionTestDataSets();
           return new java.util.Iterator<Object[]>() {

               public static final boolean INTRODUCE_READ_ERRORS = true;

               private List<Pair<Allele,Allele>> allelePairs;
               private Iterator<Pair<Allele,Allele>> allelePairsIt;
               private Iterator<GATKSAMRecord> readIt;
               private GATKSAMRecord read;
               private Iterator<List<Civar.ElementOffset>> civarEventOffsetsIt;
               private List<Civar.ElementOffset> civarEventOffsets;
               private ActiveRegionTestDataSet dataSet;
               private GraphBasedLikelihoodCalculationEngineInstance graphEngine;
               private PerReadAlleleLikelihoodMap graphLks;
               private PerReadAlleleLikelihoodMap loglessLks;
               private Map<Allele,Civar> civarByAllele;
               private String reference;

               @Override
               public boolean hasNext() {
                   return activeRegionTestDataSetIterator.hasNext()  || (readIt != null && readIt.hasNext()) || (allelePairsIt != null && allelePairsIt.hasNext());
               }

               @Override
               public Object[] next() {
                   if (allelePairsIt != null && allelePairsIt.hasNext()) {
                       final Pair<Allele,Allele> allelePair = allelePairsIt.next();
                       return new Object[] { dataSet, read, allelePair.getFirst(), allelePair.getSecond(), loglessLks, graphLks, civarEventOffsets, civarByAllele.get(allelePair.getFirst()).eventOffsets(reference,0,Integer.MAX_VALUE), civarByAllele.get(allelePair.getSecond()).eventOffsets(reference,0,Integer.MAX_VALUE)};
                   }
                   if (readIt != null && readIt.hasNext()) {
                        allelePairsIt = allelePairs.iterator();
                        final Pair<Allele,Allele> allelePair = allelePairsIt.next();
                        return new Object[] {dataSet, read = readIt.next(), allelePair.getFirst(), allelePair.getSecond(), loglessLks, graphLks, civarEventOffsets = civarEventOffsetsIt.next(), civarByAllele.get(allelePair.getFirst()).eventOffsets(reference,0,Integer.MAX_VALUE), civarByAllele.get(allelePair.getSecond()).eventOffsets(reference,0,Integer.MAX_VALUE) };
                   }
                   final Object[] params = activeRegionTestDataSetIterator.next();
                   dataSet = (ActiveRegionTestDataSet) params[0];
                   if (INTRODUCE_READ_ERRORS) dataSet.introduceErrors(new Random(13));
                   graphEngine = new GraphBasedLikelihoodCalculationEngineInstance(dataSet.assemblyResultSet(),hmm,Double.NEGATIVE_INFINITY, HeterogeneousKmerSizeResolution.COMBO_MAX);
                   graphLks = graphEngine.computeReadLikelihoods(dataSet.haplotypeList(), SampleListUtils.singletonList("anonymous"),Collections.singletonMap("anonymous",dataSet.readList())).toPerReadAlleleLikelihoodMap(0);

                   // clip reads at the anchors.
                   final Map<GATKSAMRecord,GATKSAMRecord> clippedReads = anchorClippedReads(graphEngine.getHaplotypeGraph(),dataSet.readList());
                   final List<GATKSAMRecord> clippedReadList = new ArrayList<>(dataSet.readList().size());

                   for (final GATKSAMRecord r : dataSet.readList()) {
                       clippedReadList.add(clippedReads.containsKey(r) ? clippedReads.get(r) : r);
                   }

                   loglessLks = fullPairHMM.computeReadLikelihoods(dataSet.assemblyResultSet(),SampleListUtils.singletonList("anonymous"),Collections.singletonMap("anonymous",clippedReadList)).toPerReadAlleleLikelihoodMap(0);

                   // Change clipped by unclipped in the resulting likelihood map.
                   for (final GATKSAMRecord r : clippedReads.keySet()) {
                       loglessLks.getLikelihoodReadMap().put(r,loglessLks.getLikelihoodReadMap().remove(clippedReads.get(r)));
                   }
                   final List<Haplotype> haplotypes = dataSet.haplotypeList();
                   final Map<Haplotype,Allele> alleleByHaplotype = new HashMap<>(haplotypes.size());
                   final Map<String,Civar> civarBySequence = new HashMap<>(haplotypes.size());
                   final Map<String,Haplotype> haplotypeBySequence = new HashMap<>(haplotypes.size());
                   civarByAllele = new HashMap<>(haplotypes.size());
                   final List<Civar> unrolledCivars = dataSet.unrolledCivars();
                   for (int i = 0; i < haplotypes.size(); i++)  {
                       final Haplotype h = haplotypes.get(i);
                       haplotypeBySequence.put(h.getBaseString(),h);
                       civarBySequence.put(h.getBaseString(),unrolledCivars.get(i));
                   }
                   for (final Allele a : loglessLks.getAllelesSet()) {
                       alleleByHaplotype.put(haplotypeBySequence.get(a.getBaseString()),a);
                       civarByAllele.put(a,civarBySequence.get(a.getBaseString()));
                   }
                   allelePairs = new ArrayList<>(haplotypes.size() * 2);
                   final Haplotype[] haplotypeArray = haplotypes.toArray(new Haplotype[haplotypes.size()]);
                   for (int i = 0; i < haplotypeArray.length; i++)
                       for (int j = i + 1; j < haplotypeArray.length; j++)
                           allelePairs.add(new Pair<>(alleleByHaplotype.get(haplotypeArray[i]),alleleByHaplotype.get(haplotypeArray[j])));
                   allelePairsIt = allelePairs.iterator();
                   readIt = dataSet.readList().iterator();
                   final Pair<Allele,Allele> allelePair = allelePairsIt.next();
                   civarEventOffsetsIt = dataSet.readEventOffsetList().iterator();
                   reference = dataSet.getReference();
                   return new Object[] { dataSet , read = readIt.next(), allelePair.getFirst(), allelePair.getSecond(), loglessLks, graphLks, civarEventOffsets = civarEventOffsetsIt.next(), civarByAllele.get(allelePair.getFirst()).eventOffsets(reference,0,Integer.MAX_VALUE), civarByAllele.get(allelePair.getSecond()).eventOffsets(reference,0,Integer.MAX_VALUE)};
               }

               @Override
               public void remove() {
                   throw new UnsupportedOperationException();
               }
           };
    }


    /**
     * Returns the reads clipped at their anchors.
     *
     * @param reads target reads.
     * @return never {@code null}.
     */
    protected Map<GATKSAMRecord, GATKSAMRecord> anchorClippedReads(final HaplotypeGraph haplotypeGraph, final List<GATKSAMRecord> reads) {
        final Map<GATKSAMRecord, GATKSAMRecord> result = new HashMap<>(reads.size());
        for (final GATKSAMRecord r : reads) {
            final ReadAnchoring anchoring = new ReadAnchoring(r,haplotypeGraph);
            if (anchoring.isAnchoredSomewhere())
                continue;
            final int start = anchoring.leftAnchorIndex;
            final int end = anchoring.rightAnchorIndex + haplotypeGraph.getKmerSize();
            final GATKSAMRecord clipped = new ClippedGATKSAMRecord(r, start, end);
            result.put(r, clipped);
        }
        return result;
    }



}
