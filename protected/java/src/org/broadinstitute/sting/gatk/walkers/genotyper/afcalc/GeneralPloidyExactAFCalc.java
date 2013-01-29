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

package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.gatk.walkers.genotyper.GeneralPloidyGenotypeLikelihoods;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

public class GeneralPloidyExactAFCalc extends ExactAFCalc {
    static final int MAX_LENGTH_FOR_POOL_PL_LOGGING = 10; // if PL vectors longer than this # of elements, don't log them

    private final int ploidy;
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6
    private final static boolean VERBOSE = false;

    protected GeneralPloidyExactAFCalc(final int nSamples, final int maxAltAlleles, final int ploidy) {
        super(nSamples, maxAltAlleles, ploidy);
        this.ploidy = ploidy;
    }

    @Override
    protected VariantContext reduceScope(VariantContext vc) {
        // don't try to genotype too many alternate alleles
        if ( vc.getAlternateAlleles().size() > getMaxAltAlleles()) {
            logger.warn("this tool is currently set to genotype at most " + getMaxAltAlleles() + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart() + " has " + (vc.getAlternateAlleles().size()) + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

            final List<Allele> alleles = new ArrayList<Allele>(getMaxAltAlleles() + 1);
            alleles.add(vc.getReference());
            alleles.addAll(chooseMostLikelyAlternateAlleles(vc, getMaxAltAlleles(), ploidy));

            VariantContextBuilder builder = new VariantContextBuilder(vc);
            builder.alleles(alleles);
            builder.genotypes(subsetAlleles(vc, alleles, false, ploidy));
            return builder.make();
        } else {
            return vc;
        }
    }

    @Override
    public AFCalcResult computeLog10PNonRef(final VariantContext vc, final double[] log10AlleleFrequencyPriors) {
        combineSinglePools(vc.getGenotypes(), vc.getNAlleles(), ploidy, log10AlleleFrequencyPriors);
        return getResultFromFinalState(vc, log10AlleleFrequencyPriors);
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
            alleleCountSetList = new LinkedList<ExactACset>();
            conformationMap = new HashMap<ExactACcounts,ExactACset>();
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

    /**
     *
     * Chooses N most likely alleles in a set of pools (samples) based on GL sum over alt alleles
     * @param vc                          Input variant context
     * @param numAllelesToChoose          Number of alleles to choose
     * @param ploidy                      Ploidy per pool
     * @return                            list of numAllelesToChoose most likely alleles
     */

    private static final int PL_INDEX_OF_HOM_REF = 0;
    private static List<Allele> chooseMostLikelyAlternateAlleles(VariantContext vc, int numAllelesToChoose, int ploidy) {
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ )
            likelihoodSums[i] = new LikelihoodSum(vc.getAlternateAllele(i));

        // based on the GLs, find the alternate alleles with the most probability; sum the GLs for the most likely genotype
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes(), false);
        for ( final double[] likelihoods : GLs ) {

            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            final int[] acCount = GeneralPloidyGenotypeLikelihoods.getAlleleCountFromPLIndex(1 + numOriginalAltAlleles, ploidy, PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++) {
                if (acCount[k] > 0)
                    likelihoodSums[k-1].sum += acCount[k] * (likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF]);

            }
        }

        // sort them by probability mass and choose the best ones
        Collections.sort(Arrays.asList(likelihoodSums));
        final ArrayList<Allele> bestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( int i = 0; i < numAllelesToChoose; i++ )
            bestAlleles.add(likelihoodSums[i].allele);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<Allele>(numAllelesToChoose);
        for ( Allele allele : vc.getAlternateAlleles() ) {
            if ( bestAlleles.contains(allele) )
                orderedBestAlleles.add(allele);
        }

        return orderedBestAlleles;
    }


    /**
     * Simple non-optimized version that combines GLs from several pools and produces global AF distribution.
     * @param GLs                              Inputs genotypes context with per-pool GLs
     * @param numAlleles                       Number of alternate alleles
     * @param ploidyPerPool                    Number of samples per pool
     * @param log10AlleleFrequencyPriors       Frequency priors
     */
    protected void combineSinglePools(final GenotypesContext GLs,
                                      final int numAlleles,
                                      final int ploidyPerPool,
                                      final double[] log10AlleleFrequencyPriors) {

        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs, true);


        int combinedPloidy = 0;

        // Combine each pool incrementally - likelihoods will be renormalized at each step
         CombinedPoolLikelihoods combinedPoolLikelihoods = new CombinedPoolLikelihoods();

        // first element: zero ploidy, e.g. trivial degenerate distribution
        final int[] zeroCounts = new int[numAlleles];
        final ExactACset set = new ExactACset(1, new ExactACcounts(zeroCounts));
        set.getLog10Likelihoods()[0] = 0.0;

        combinedPoolLikelihoods.add(set);

        if ( genotypeLikelihoods.size() <= 1 ) {
            // no meaningful GLs at all, just set the tracker to non poly values
            getStateTracker().reset(); // just mimic-ing call below
            getStateTracker().setLog10LikelihoodOfAFzero(0.0);
        } else {
            for (int p=1; p<genotypeLikelihoods.size(); p++) {
                getStateTracker().reset(); // TODO -- why is this here?  It makes it hard to track the n evaluation
                combinedPoolLikelihoods = fastCombineMultiallelicPool(combinedPoolLikelihoods, genotypeLikelihoods.get(p),
                        combinedPloidy, ploidyPerPool, numAlleles, log10AlleleFrequencyPriors);
                combinedPloidy = ploidyPerPool + combinedPloidy; // total number of chromosomes in combinedLikelihoods
            }
        }
    }

    public CombinedPoolLikelihoods fastCombineMultiallelicPool(final CombinedPoolLikelihoods originalPool,
                                                               double[] newGL,
                                                               int originalPloidy,
                                                               int newGLPloidy,
                                                               int numAlleles,
                                                               final double[] log10AlleleFrequencyPriors) {
        final LinkedList<ExactACset> ACqueue = new LinkedList<ExactACset>();
        // mapping of ExactACset indexes to the objects
        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<ExactACcounts, ExactACset>();
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
            getStateTracker().incNEvaluations();
            // compute log10Likelihoods
            final ExactACset ACset = ACqueue.remove();
            final double log10LofKs = calculateACConformationAndUpdateQueue(ACset, newPool, originalPool, newGL, log10AlleleFrequencyPriors, originalPloidy, newGLPloidy, ACqueue, indexesToACset);

            // clean up memory
            indexesToACset.remove(ACset.getACcounts());
            if ( VERBOSE )
                System.out.printf(" *** removing used set=%s%n", ACset.getACcounts());

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
                                                         final HashMap<ExactACcounts, ExactACset> indexesToACset) {

        // compute likeihood in "set" of new set based on original likelihoods
        final int numAlleles = set.getACcounts().getCounts().length;
        final int newPloidy = set.getACsum();
        final double log10LofK = computeLofK(set, originalPool, newGL, log10AlleleFrequencyPriors, numAlleles, originalPloidy, newGLPloidy);


        // add to new pool
        if (!Double.isInfinite(log10LofK))
            newPool.add(set);

        // TODO -- change false to true this correct line when the implementation of this model is optimized (it's too slow now to handle this fix)
        if ( getStateTracker().abort(log10LofK, set.getACcounts(), false) ) {
            return log10LofK;
        }

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


//    /**
//     * Naive combiner of two multiallelic pools - number of alt alleles must be the same.
//     * Math is generalization of biallelic combiner.
//     *
//     * For vector K representing an allele count conformation,
//     * Pr(D | AC = K) = Sum_G Pr(D|AC1 = G) Pr (D|AC2=K-G) * F(G,K)
//     * where F(G,K) = choose(m1,[g0 g1 ...])*choose(m2,[...]) / choose(m1+m2,[k1 k2 ...])
//     * @param originalPool                    First log-likelihood pool GL vector
//     * @param yy                    Second pool GL vector
//     * @param ploidy1               Ploidy of first pool (# of chromosomes in it)
//     * @param ploidy2               Ploidy of second pool
//     * @param numAlleles            Number of alleles
//     * @param log10AlleleFrequencyPriors Array of biallelic priors
//     * @param resultTracker                Af calculation result object
//     */
//    public static void combineMultiallelicPoolNaively(CombinedPoolLikelihoods originalPool, double[] yy, int ploidy1, int ploidy2, int numAlleles,
//                                                      final double[] log10AlleleFrequencyPriors,
//                                                      final AFCalcResultTracker resultTracker) {
///*
//        final int dim1 = GenotypeLikelihoods.numLikelihoods(numAlleles, ploidy1);
//        final int dim2 = GenotypeLikelihoods.numLikelihoods(numAlleles, ploidy2);
//
//        if (dim1 != originalPool.getLength() || dim2 != yy.length)
//            throw new ReviewedStingException("BUG: Inconsistent vector length");
//
//        if (ploidy2 == 0)
//            return;
//
//        final int newPloidy = ploidy1 + ploidy2;
//
//        // Say L1(K) = Pr(D|AC1=K) * choose(m1,K)
//        // and L2(K) = Pr(D|AC2=K) * choose(m2,K)
//        GeneralPloidyGenotypeLikelihoods.SumIterator firstIterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(numAlleles,ploidy1);
//        final double[] x = originalPool.getLikelihoodsAsVector(true);
//        while(firstIterator.hasNext()) {
//            x[firstIterator.getLinearIndex()] += MathUtils.log10MultinomialCoefficient(ploidy1,firstIterator.getCurrentVector());
//            firstIterator.next();
//        }
//
//        GeneralPloidyGenotypeLikelihoods.SumIterator secondIterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(numAlleles,ploidy2);
//        final double[] y = yy.clone();
//        while(secondIterator.hasNext()) {
//            y[secondIterator.getLinearIndex()] += MathUtils.log10MultinomialCoefficient(ploidy2,secondIterator.getCurrentVector());
//            secondIterator.next();
//        }
//
//        // initialize output to -log10(choose(m1+m2,[k1 k2...])
//        final int outputDim = GenotypeLikelihoods.numLikelihoods(numAlleles, newPloidy);
//        final GeneralPloidyGenotypeLikelihoods.SumIterator outputIterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(numAlleles,newPloidy);
//
//
//        // Now, result(K) =  logSum_G (L1(G)+L2(K-G)) where G are all possible vectors that sum UP to K
//        while(outputIterator.hasNext()) {
//            final ExactACset set = new ExactACset(1, new ExactACcounts(outputIterator.getCurrentAltVector()));
//            double likelihood = computeLofK(set, x,y, log10AlleleFrequencyPriors, numAlleles, ploidy1, ploidy2, result);
//
//            originalPool.add(likelihood, set, outputIterator.getLinearIndex());
//            outputIterator.next();
//        }
//*/
//    }

    /**
     * Compute likelihood of a particular AC conformation and update AFresult object
     * @param set                     Set of AC counts to compute
     * @param firstGLs                  Original pool likelihoods before combining
     * @param secondGL                  New GL vector with additional pool
     * @param log10AlleleFrequencyPriors     Allele frequency priors
     * @param numAlleles                Number of alleles (including ref)
     * @param ploidy1                   Ploidy of original pool (combined)
     * @param ploidy2                   Ploidy of new pool
     * @return                          log-likehood of requested conformation
     */
    private double computeLofK(final ExactACset set,
                               final CombinedPoolLikelihoods firstGLs,
                               final double[] secondGL,
                               final double[] log10AlleleFrequencyPriors,
                               final int numAlleles, final int ploidy1, final int ploidy2) {

        final int newPloidy = ploidy1 + ploidy2;

        // sanity check
        int totalAltK = set.getACsum();
        if (newPloidy != totalAltK)
            throw new ReviewedStingException("BUG: inconsistent sizes of set.getACsum and passed ploidy values");

        totalAltK -= set.getACcounts().getCounts()[0];
        // totalAltK has sum of alt alleles of conformation now


        // special case for k = 0 over all k
        if ( totalAltK == 0 ) {   // all-ref case
            final double log10Lof0 = firstGLs.getGLOfACZero() + secondGL[HOM_REF_INDEX];
            set.getLog10Likelihoods()[0] = log10Lof0;

            getStateTracker().setLog10LikelihoodOfAFzero(log10Lof0);
            getStateTracker().setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
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
        getStateTracker().updateMLEifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        // apply the priors over each alternate allele
        for (final int ACcount : altCounts ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }
        // TODO -- GUILLERMO THIS CODE MAY PRODUCE POSITIVE LIKELIHOODS OR -INFINITY
        getStateTracker().updateMAPifNeeded(Math.max(log10LofK, -Double.MAX_VALUE), altCounts);

        return log10LofK;
    }

    /**
     * Small helper routine - is a particular AC conformationv vector valid? ie are all elements non-negative and sum to ploidy?
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
     * including updating the PL's, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @param ploidy                            number of chromosomes per sample (pool)
     * @return                                  GenotypesContext with new PLs
     */
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes,
                                          final int ploidy) {
        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();
        List<Allele> NO_CALL_ALLELES = new ArrayList<Allele>(ploidy);

        for (int k=0; k < ploidy; k++)
            NO_CALL_ALLELES.add(Allele.NO_CALL);

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;


        // create the new genotypes
        for ( int k = 0; k < oldGTs.size(); k++ ) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;

            // Optimization: if # of new alt alleles = 0 (pure ref call), keep original likelihoods so we skip normalization
            // and subsetting
            if ( numOriginalAltAlleles == numNewAltAlleles || numNewAltAlleles == 0) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = GeneralPloidyGenotypeLikelihoods.subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse);

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
            }
            else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                if ( numNewAltAlleles == 0 )
                    gb.noPL();
                else
                    gb.PL(newLikelihoods);

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL )
                    gb.alleles(NO_CALL_ALLELES);
                else
                    assignGenotype(gb, newLikelihoods, allelesToUse, ploidy);
                newGTs.add(gb.make());
            }
        }

        return newGTs;

    }

    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     *
     * @return genotype
     */
    private void assignGenotype(final GenotypeBuilder gb,
                                final double[] newLikelihoods,
                                final List<Allele> allelesToUse,
                                final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;



        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);

        final int[] mlAlleleCount = GeneralPloidyGenotypeLikelihoods.getAlleleCountFromPLIndex(allelesToUse.size(), numChromosomes, PLindex);
        final ArrayList<Double> alleleFreqs = new ArrayList<Double>();
        final ArrayList<Integer> alleleCounts = new ArrayList<Integer>();


        for (int k=1; k < mlAlleleCount.length; k++) {
            alleleCounts.add(mlAlleleCount[k]);
            final double freq = (double)mlAlleleCount[k] / (double)numChromosomes;
            alleleFreqs.add(freq);

        }

        // per-pool logging of AC and AF
        gb.attribute(VCFConstants.MLE_PER_SAMPLE_ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
        gb.attribute(VCFConstants.MLE_PER_SAMPLE_ALLELE_FRACTION_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);

        // remove PLs if necessary
        if (newLikelihoods.length > MAX_LENGTH_FOR_POOL_PL_LOGGING)
            gb.noPL();

        ArrayList<Allele> myAlleles = new ArrayList<Allele>();

        // add list of called ML genotypes to alleles list
        // TODO - too unwieldy?
        int idx = 0;
        for (int mlind = 0; mlind < mlAlleleCount.length; mlind++) {
            for (int k=0; k < mlAlleleCount[mlind]; k++)
                myAlleles.add(idx++,allelesToUse.get(mlind));
        }
        gb.alleles(myAlleles);

        if ( numNewAltAlleles > 0 )
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
    }

}
