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

package org.broadinstitute.gatk.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.GeneralPloidyGenotypeLikelihoods;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

public class GeneralPloidyExactAFCalculator extends ExactAFCalculator {

    protected GeneralPloidyExactAFCalculator() {
    }

    @Override
    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        return subsetAlleles(vc, defaultPloidy, allelesToUse, false);
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        combineSinglePools(vc.getGenotypes(), defaultPloidy, vc.getNAlleles(), log10AlleleFrequencyPriors);
        return getResultFromFinalState(vc, log10AlleleFrequencyPriors, stateTracker);
    }

    /**
     * Simple wrapper class to hold values of combined pool likelihoods.
     * For fast hashing and fast retrieval, there's a hash map that shadows main list.
     *
     */
    static class CombinedPoolLikelihoods {
        private LinkedList<ExactACset> alleleCountSetList;
        private HashMap<ExactACcounts,ExactACset> conformationMap;
        private double maxLikelihood;


        public CombinedPoolLikelihoods() {
            // final int numElements = GenotypeLikelihoods.numLikelihoods();
            alleleCountSetList = new LinkedList<>();
            conformationMap = new HashMap<>();
            maxLikelihood = Double.NEGATIVE_INFINITY;
        }

        public void add(ExactACset set) {
            alleleCountSetList.add(set);
            conformationMap.put(set.getACcounts(), set);
            final double likelihood = set.getLog10Likelihoods()[0];

            if (likelihood > maxLikelihood )
                maxLikelihood = likelihood;

        }

        public boolean hasConformation(int[] ac) {
            return conformationMap.containsKey(new ExactACcounts(ac));

        }

        public double getLikelihoodOfConformation(int[] ac) {
            return conformationMap.get(new ExactACcounts(ac)).getLog10Likelihoods()[0];
        }

        public double getGLOfACZero() {
            return alleleCountSetList.get(0).getLog10Likelihoods()[0]; // AC 0 is always at beginning of list
        }

        public int getLength() {
            return alleleCountSetList.size();
        }
    }


    @Override
    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        final int numOriginalAltAlleles = likelihoodSums.length;
        final GenotypesContext genotypes = vc.getGenotypes();
        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            if (!genotype.hasPL())
                continue;
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (MathUtils.sum(gls) >= GATKVariantContextUtils.SUM_GL_THRESH_NOCALL)
                continue;

            final int PLindexOfBestGL = MathUtils.maxElementIndex(gls);

            final double bestToHomRefDiffGL = PLindexOfBestGL == PL_INDEX_OF_HOM_REF ? 0.0 : gls[PLindexOfBestGL] - gls[PL_INDEX_OF_HOM_REF];
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;

            final int[] acCount = GeneralPloidyGenotypeLikelihoods.getAlleleCountFromPLIndex(1 + numOriginalAltAlleles, ploidy, PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++)
                if (acCount[k] > 0 )
                    likelihoodSums[k-1].sum += acCount[k] * bestToHomRefDiffGL;
        }
    }

    /**
     * Simple non-optimized version that combines GLs from several pools and produces global AF distribution.
     * @param GLs                              Inputs genotypes context with per-pool GLs
     * @param numAlleles                       Number of alternate alleles
     * @param log10AlleleFrequencyPriors       Frequency priors
     */
    protected void combineSinglePools(final GenotypesContext GLs,
                                      final int defaultPloidy,
                                      final int numAlleles,
                                      final double[] log10AlleleFrequencyPriors) {

        // Combine each pool incrementally - likelihoods will be renormalized at each step

        // first element: zero ploidy, e.g. trivial degenerate distribution
        final int numAltAlleles = numAlleles - 1;
        final int[] zeroCounts = new int[numAlleles];
        final ExactACset set = new ExactACset(1, new ExactACcounts(zeroCounts));
        set.getLog10Likelihoods()[0] = 0.0;
        final StateTracker stateTracker = getStateTracker(false,numAltAlleles);
        int combinedPloidy = 0;
        CombinedPoolLikelihoods combinedPoolLikelihoods = new CombinedPoolLikelihoods();
        combinedPoolLikelihoods.add(set);

        for (final Genotype genotype : GLs.iterateInSampleNameOrder()) {
            // recover gls and check if they qualify.
            if (!genotype.hasPL())
                continue;
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (MathUtils.sum(gls) >= GATKVariantContextUtils.SUM_GL_THRESH_NOCALL)
                continue;
            stateTracker.reset();
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy < 1 ? defaultPloidy : declaredPloidy;
            // they do qualify so we proceed.
            combinedPoolLikelihoods = fastCombineMultiallelicPool(combinedPoolLikelihoods, gls,
                    combinedPloidy, ploidy, numAlleles, log10AlleleFrequencyPriors, stateTracker);
            combinedPloidy = ploidy + combinedPloidy; // total number of chromosomes in combinedLikelihoods
        }
        if (combinedPloidy == 0)
            stateTracker.setLog10LikelihoodOfAFzero(0.0);
    }

    private CombinedPoolLikelihoods fastCombineMultiallelicPool(final CombinedPoolLikelihoods originalPool,
                                                               double[] newGL,
                                                               int originalPloidy,
                                                               int newGLPloidy,
                                                               int numAlleles,
                                                               final double[] log10AlleleFrequencyPriors,
                                                               final StateTracker stateTracker) {
        final LinkedList<ExactACset> ACqueue = new LinkedList<>();
        // mapping of ExactACset indexes to the objects
        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<>();
        final CombinedPoolLikelihoods newPool = new CombinedPoolLikelihoods();

        // add AC=0 to the queue
        final int[] zeroCounts = new int[numAlleles];
        final int newPloidy = originalPloidy + newGLPloidy;
        zeroCounts[0] = newPloidy;

        ExactACset zeroSet = new ExactACset(1, new ExactACcounts(zeroCounts));

        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.getACcounts(), zeroSet);

        // keep processing while we have AC conformations that need to be calculated
        while ( !ACqueue.isEmpty() ) {
            stateTracker.incNEvaluations();
            // compute log10Likelihoods
            final ExactACset ACset = ACqueue.remove();

            calculateACConformationAndUpdateQueue(ACset, newPool, originalPool, newGL, log10AlleleFrequencyPriors, originalPloidy, newGLPloidy, ACqueue, indexesToACset, stateTracker);

            // clean up memory
            indexesToACset.remove(ACset.getACcounts());
        }
        return newPool;
    }

    // todo - refactor, function almost identical except for log10LofK computation in GeneralPloidyGenotypeLikelihoods
    /**
     *
     * @param set                       ExactACset holding conformation to be computed
     * @param newPool                   New pool likelihood holder
     * @param originalPool              Original likelihood holder
     * @param newGL                     New pool GL vector to combine
     * @param log10AlleleFrequencyPriors Prior object
     * @param originalPloidy             Total ploidy of original combined pool
     * @param newGLPloidy                Ploidy of GL vector
     * @param ACqueue                    Queue of conformations to compute
     * @param indexesToACset             AC indices of objects in queue
     * @return                           max log likelihood
     */
    private double calculateACConformationAndUpdateQueue(final ExactACset set,
                                                         final CombinedPoolLikelihoods newPool,
                                                         final CombinedPoolLikelihoods originalPool,
                                                         final double[] newGL,
                                                         final double[] log10AlleleFrequencyPriors,
                                                         final int originalPloidy,
                                                         final int newGLPloidy,
                                                         final LinkedList<ExactACset> ACqueue,
                                                         final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                                         final StateTracker stateTracker) {

        // compute likelihood in "set" of new set based on original likelihoods
        final int numAlleles = set.getACcounts().getCounts().length;
        final int newPloidy = set.getACsum();
        final double log10LofK = computeLofK(set, originalPool, newGL, log10AlleleFrequencyPriors, numAlleles, originalPloidy, newGLPloidy, stateTracker);


        // add to new pool
        if (!Double.isInfinite(log10LofK))
            newPool.add(set);

        if ( stateTracker.abort(log10LofK, set.getACcounts(), true, true) )
            return log10LofK;

        // iterate over higher frequencies if possible
        // by convention, ACcounts contained in set have full vector of possible pool ac counts including ref count.
        // so, if first element is zero, it automatically means we have no wiggle since we're in a corner of the conformation space
        final int ACwiggle = set.getACcounts().getCounts()[0];
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;


        // add conformations for other cases
        for ( int allele = 1; allele < numAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // is this a valid conformation?
            int altSum = (int)MathUtils.sum(ACcountsClone) - ACcountsClone[0];
            ACcountsClone[0] = newPloidy - altSum;
            if (ACcountsClone[0] < 0)
                continue;


            GeneralPloidyGenotypeLikelihoods.updateACset(ACcountsClone, ACqueue, indexesToACset);
        }


        return log10LofK;
    }


    /**
     * Compute likelihood of a particular AC conformation and update AFresult object
     * @param set                     Set of AC counts to compute
     * @param firstGLs                  Original pool likelihoods before combining
     * @param secondGL                  New GL vector with additional pool
     * @param log10AlleleFrequencyPriors     Allele frequency priors
     * @param numAlleles                Number of alleles (including ref)
     * @param ploidy1                   Ploidy of original pool (combined)
     * @param ploidy2                   Ploidy of new pool
     * @return                          log-likelihood of requested conformation
     */
    private double computeLofK(final ExactACset set,
                               final CombinedPoolLikelihoods firstGLs,
                               final double[] secondGL,
                               final double[] log10AlleleFrequencyPriors,
                               final int numAlleles, final int ploidy1, final int ploidy2, final StateTracker stateTracker) {

        final int newPloidy = ploidy1 + ploidy2;

        // sanity check
        int totalAltK = set.getACsum();
        if (newPloidy != totalAltK)
            throw new ReviewedGATKException("BUG: inconsistent sizes of set.getACsum and passed ploidy values");

        totalAltK -= set.getACcounts().getCounts()[0];
        // totalAltK has sum of alt alleles of conformation now


        // special case for k = 0 over all k
        if ( totalAltK == 0 ) {   // all-ref case
            final double log10Lof0 = firstGLs.getGLOfACZero() + secondGL[HOM_REF_INDEX];
            set.getLog10Likelihoods()[0] = log10Lof0;
            stateTracker.setLog10LikelihoodOfAFzero(log10Lof0);
            stateTracker.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return log10Lof0;

        }   else {

            // initialize result with denominator
            // ExactACset holds by convention the conformation of all alleles, and the sum of all allele count is just the ploidy.
            // To compute n!/k1!k2!k3!... we need to compute first n!/(k2!k3!...) and then further divide by k1! where k1=ploidy-sum_k_i

            int[] currentCount = set.getACcounts().getCounts();
            double denom =  -MathUtils.log10MultinomialCoefficient(newPloidy, currentCount);

            // for current conformation, get all possible ways to break vector K into two components G1 and G2
            final GeneralPloidyGenotypeLikelihoods.SumIterator innerIterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(numAlleles,ploidy2);
            set.getLog10Likelihoods()[0] = Double.NEGATIVE_INFINITY;
            while (innerIterator.hasNext()) {
                // check if breaking current conformation into g1 and g2 is feasible.
                final int[] acCount2 = innerIterator.getCurrentVector();
                final int[] acCount1 = MathUtils.vectorDiff(currentCount, acCount2);
                final int idx2 = innerIterator.getLinearIndex();
                // see if conformation is valid and if original pool had this conformation
                // for conformation to be valid, all elements of g2 have to be <= elements of current AC set
                if (isValidConformation(acCount1,ploidy1) && firstGLs.hasConformation(acCount1)) {
                    final double gl2 = secondGL[idx2];
                    if (!Double.isInfinite(gl2)) {
                        final double firstGL = firstGLs.getLikelihoodOfConformation(acCount1);
                        final double num1 = MathUtils.log10MultinomialCoefficient(ploidy1, acCount1);
                        final double num2 = MathUtils.log10MultinomialCoefficient(ploidy2, acCount2);
                        final double sum = firstGL + gl2 + num1 + num2;

                        set.getLog10Likelihoods()[0] = MathUtils.approximateLog10SumLog10(set.getLog10Likelihoods()[0], sum);
                    }
                }
                innerIterator.next();
            }

            set.getLog10Likelihoods()[0] += denom;
        }

        double log10LofK = set.getLog10Likelihoods()[0];

        // update the MLE if necessary
        final int altCounts[] = Arrays.copyOfRange(set.getACcounts().getCounts(),1, set.getACcounts().getCounts().length);
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        stateTracker.updateMLEifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        // apply the priors over each alternate allele
        for (final int ACcount : altCounts ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        stateTracker.updateMAPifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        return log10LofK;
    }

    /**
     * Small helper routine - is a particular AC conformation vector valid? ie are all elements non-negative and sum to ploidy?
     * @param set                            AC conformation vector
     * @param ploidy                         Ploidy of set
     * @return                               Valid conformation
     */
    private static boolean isValidConformation(final int[] set, final int ploidy) {
        int sum=0;
        for (final int ac: set) {
            if (ac < 0)
                return false;
            sum += ac;

        }

        return (sum == ploidy);
    }

    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PLs, ADs and SACs, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information
     *                                          for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @return                                  GenotypesContext with new PLs, SACs and AD.
     */
    @Override
    public GenotypesContext subsetAlleles(final VariantContext vc, final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {

        final GenotypesContext result = GenotypesContext.create();

        // Subset genotypes for each sample
        for (final Genotype g : vc.getGenotypes()) // If it really needs to process order by sample name do so.
            result.add(subsetGenotypeAlleles(g, allelesToUse, vc, defaultPloidy, assignGenotypes));
        return GATKVariantContextUtils.fixADFromSubsettedAlleles(result, vc, allelesToUse);
    }

    /**
     * From a given genotype, extract a given subset of alleles and update genotype PLs and SACs.
     * @param g                                 genotype to subset
     * @param allelesToUse                      alleles to subset
     * @param vc                                variant context with alleles and genotypes
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information for a sample.
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @return                                  Genotypes with new PLs and SACs
     */
    private Genotype subsetGenotypeAlleles(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc, final int defaultPloidy,
                                           boolean assignGenotypes) {
        final int ploidy = g.getPloidy() <= 0 ? defaultPloidy : g.getPloidy();
        if (!g.hasLikelihoods())
            return GenotypeBuilder.create(g.getSampleName(),GATKVariantContextUtils.noCallAlleles(ploidy));
        else {
            // subset likelihood alleles
            final double[] newLikelihoods = subsetLikelihoodAlleles(g, allelesToUse, vc, ploidy);
            if (MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL)
                return GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy));
            else  // just now we would care about newSACs
                return subsetGenotypeAllelesWithLikelihoods(g, allelesToUse, vc, ploidy, assignGenotypes, newLikelihoods);
        }
    }

    /**
     * From a given genotype, extract a given subset of alleles and return the new PLs
     * @param g                                 genotype to subset
     * @param allelesToUse                      alleles to subset
     * @param vc                                variant context with alleles and genotypes
     * @param ploidy                            number of chromosomes
     * @return                                  the subsetted PLs
     */
    private double[] subsetLikelihoodAlleles(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc, final int ploidy){

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // create the new likelihoods array from the alleles we are allowed to use
        final double[] originalLikelihoods = g.getLikelihoods().getAsVector();

        if ( numOriginalAltAlleles != numNewAltAlleles ) {
            // might need to re-normalize the new likelihoods
            return MathUtils.normalizeFromLog10(GeneralPloidyGenotypeLikelihoods.subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse),
                    false, true);
        }
        else
            return originalLikelihoods;
    }

    /**
     * From a given genotype, subset the PLs and SACs
     * @param g                                 genotype to subset
     * @param allelesToUse                      alleles to subset
     * @param vc                                variant context with alleles and genotypes
     * @param ploidy                            number of chromosomes
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @param newLikelihoods                    the PL values
     * @return genotype with the subsetted PLsL and SACs
     */
    private Genotype subsetGenotypeAllelesWithLikelihoods(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc, int ploidy,
                                                          final boolean assignGenotypes, final double[] newLikelihoods) {

        final GenotypeBuilder gb = new GenotypeBuilder(g);
        final String sampleName = g.getSampleName();

        // add likelihoods
        gb.PL(newLikelihoods);

        // get and add subsetted SACs
        final int[] newSACs = subsetSACAlleles(g, allelesToUse, vc);
        if (newSACs != null)
            gb.attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, newSACs);
        if (assignGenotypes)
            assignGenotype(gb, vc, sampleName, newLikelihoods, allelesToUse, ploidy);
        else
            gb.alleles(GATKVariantContextUtils.noCallAlleles(ploidy));

        return gb.make();
    }

    /**
     * From a given genotype, extract a given subset of alleles and return the new SACs
     * @param g                             genotype to subset
     * @param allelesToUse                  alleles to subset
     * @param vc                            variant context with alleles and genotypes
     * @return                              the subsetted SACs
     */
    private int[] subsetSACAlleles(final Genotype g, final List<Allele> allelesToUse, final VariantContext vc){

        if ( !g.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY) )
            return null;

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;
        final List<Integer> sacIndexesToUse = numOriginalAltAlleles == numNewAltAlleles ? null : GATKVariantContextUtils.determineSACIndexesToUse(vc, allelesToUse);

        return GATKVariantContextUtils.makeNewSACs(g, sacIndexesToUse);
    }

    /**
     * Assign genotypes (GTs) to the samples in the VariantContext greedily based on the PLs
     *
     * @param gb                   the GenotypeBuilder to modify
     * @param vc                   the VariantContext
     * @param sampleName           the sample name
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes       Number of chromosomes per pool
     */
    private void assignGenotype(final GenotypeBuilder gb,
                                final VariantContext vc,
                                final String sampleName,
                                final double[] newLikelihoods,
                                final List<Allele> allelesToUse,
                                final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(numChromosomes,allelesToUse.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);

        gb.alleles(alleleCounts.asAlleleList(allelesToUse));

        removePLsIfMaxNumPLValuesExceeded(gb, vc, sampleName, newLikelihoods);

        // TODO - deprecated so what is the appropriate method to call?
        if ( numNewAltAlleles > 0 )
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
    }
}
