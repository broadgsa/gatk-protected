/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class PoolAFCalculationModel extends AlleleFrequencyCalculationModel {
    static final int MAX_LENGTH_FOR_POOL_PL_LOGGING = 10; // if PL vectors longer than this # of elements, don't log them
    final protected UnifiedArgumentCollection UAC;

    private final int ploidy;
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6
    private final static boolean VERBOSE = false;

    protected PoolAFCalculationModel(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        ploidy = UAC.samplePloidy;
        this.UAC = UAC;

    }

    public List<Allele> getLog10PNonRef(final VariantContext vc,
                                        final double[] log10AlleleFrequencyPriors,
                                        final AlleleFrequencyCalculationResult result) {

        GenotypesContext GLs = vc.getGenotypes();
        List<Allele> alleles = vc.getAlleles();

        // don't try to genotype too many alternate alleles
        if ( vc.getAlternateAlleles().size() > MAX_ALTERNATE_ALLELES_TO_GENOTYPE ) {
            logger.warn("this tool is currently set to genotype at most " + MAX_ALTERNATE_ALLELES_TO_GENOTYPE + " alternate alleles in a given context, but the context at " + vc.getChr() + ":" + vc.getStart() + " has " + (vc.getAlternateAlleles().size()) + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

            alleles = new ArrayList<Allele>(MAX_ALTERNATE_ALLELES_TO_GENOTYPE + 1);
            alleles.add(vc.getReference());
            alleles.addAll(chooseMostLikelyAlternateAlleles(vc, MAX_ALTERNATE_ALLELES_TO_GENOTYPE, ploidy));


            GLs = subsetAlleles(vc, alleles, false, ploidy);
        }

        combineSinglePools(GLs, alleles.size(), ploidy, log10AlleleFrequencyPriors, result);

        return alleles;
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
            conformationMap.put(set.ACcounts, set);
            final double likelihood = set.log10Likelihoods[0];

            if (likelihood > maxLikelihood )
                maxLikelihood = likelihood;

        }

        public boolean hasConformation(int[] ac) {
            return conformationMap.containsKey(new ExactACcounts(ac));

        }

        public double getLikelihoodOfConformation(int[] ac) {
            return conformationMap.get(new ExactACcounts(ac)).log10Likelihoods[0];
        }

        public double getGLOfACZero() {
            return alleleCountSetList.get(0).log10Likelihoods[0]; // AC 0 is always at beginning of list
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

    private static List<Allele> chooseMostLikelyAlternateAlleles(VariantContext vc, int numAllelesToChoose, int ploidy) {
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ )
            likelihoodSums[i] = new LikelihoodSum(vc.getAlternateAllele(i));

        // based on the GLs, find the alternate alleles with the most probability; sum the GLs for the most likely genotype
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes());
        for ( final double[] likelihoods : GLs ) {

            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            final int[] acCount = PoolGenotypeLikelihoods.getAlleleCountFromPLIndex(1+numOriginalAltAlleles,ploidy,PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++) {
                if (acCount[k] > 0)
                    likelihoodSums[k-1].sum += likelihoods[PLindexOfBestGL];

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
     * @param result                           object to fill with output values
     */
    protected static void combineSinglePools(final GenotypesContext GLs,
                                             final int numAlleles,
                                             final int ploidyPerPool,
                                             final double[] log10AlleleFrequencyPriors,
                                             final AlleleFrequencyCalculationResult result) {

        final ArrayList<double[]> genotypeLikelihoods = getGLs(GLs);


        int combinedPloidy = 0;

        // Combine each pool incrementally - likelihoods will be renormalized at each step
         CombinedPoolLikelihoods combinedPoolLikelihoods = new CombinedPoolLikelihoods();

        // first element: zero ploidy, e.g. trivial degenerate distribution
        final int[] zeroCounts = new int[numAlleles];
        final ExactACset set = new ExactACset(1, new ExactACcounts(zeroCounts));
        set.log10Likelihoods[0] = 0.0;

        combinedPoolLikelihoods.add(set);
        for (int p=1; p<genotypeLikelihoods.size(); p++) {
            result.reset();
            combinedPoolLikelihoods = fastCombineMultiallelicPool(combinedPoolLikelihoods, genotypeLikelihoods.get(p), combinedPloidy, ploidyPerPool,
                    numAlleles, log10AlleleFrequencyPriors, result);
            combinedPloidy = ploidyPerPool + combinedPloidy; // total number of chromosomes in combinedLikelihoods
        }
    }

    public static CombinedPoolLikelihoods fastCombineMultiallelicPool(final CombinedPoolLikelihoods originalPool, double[] newGL, int originalPloidy, int newGLPloidy, int numAlleles,
                                                                      final double[] log10AlleleFrequencyPriors,
                                                                      final AlleleFrequencyCalculationResult result) {



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
        indexesToACset.put(zeroSet.ACcounts, zeroSet);

        // keep processing while we have AC conformations that need to be calculated
        double maxLog10L = Double.NEGATIVE_INFINITY;
        while ( !ACqueue.isEmpty() ) {
            // compute log10Likelihoods
            final ExactACset ACset = ACqueue.remove();
            final double log10LofKs = calculateACConformationAndUpdateQueue(ACset, newPool, originalPool, newGL, log10AlleleFrequencyPriors, originalPloidy, newGLPloidy, result, maxLog10L, ACqueue, indexesToACset);
            maxLog10L = Math.max(maxLog10L, log10LofKs);
            // clean up memory
            indexesToACset.remove(ACset.ACcounts);
            if ( VERBOSE )
                System.out.printf(" *** removing used set=%s%n", ACset.ACcounts);

        }
        return newPool;
    }

    // todo - refactor, function almost identical except for log10LofK computation in PoolGenotypeLikelihoods
    /**
     *
     * @param set                       ExactACset holding conformation to be computed
     * @param newPool                   New pool likelihood holder
     * @param originalPool              Original likelihood holder
     * @param newGL                     New pool GL vector to combine
     * @param log10AlleleFrequencyPriors Prior object
     * @param originalPloidy             Total ploidy of original combined pool
     * @param newGLPloidy                Ploidy of GL vector
     * @param result                     AFResult object
     * @param maxLog10L                  max likelihood observed so far
     * @param ACqueue                    Queue of conformations to compute
     * @param indexesToACset             AC indices of objects in queue
     * @return                           max log likelihood
     */
    private static double calculateACConformationAndUpdateQueue(final ExactACset set,
                                                                final CombinedPoolLikelihoods newPool,
                                                                final CombinedPoolLikelihoods originalPool,
                                                                final double[] newGL,
                                                                final double[] log10AlleleFrequencyPriors,
                                                                final int originalPloidy,
                                                                final int newGLPloidy,
                                                                final AlleleFrequencyCalculationResult result,
                                                                final double  maxLog10L,
                                                                final LinkedList<ExactACset> ACqueue,
                                                                final HashMap<ExactACcounts, ExactACset> indexesToACset) {

        // compute likeihood in "set" of new set based on original likelihoods
        final int numAlleles = set.ACcounts.counts.length;
        final int newPloidy = set.getACsum();
        final double log10LofK = computeLofK(set, originalPool, newGL, log10AlleleFrequencyPriors, numAlleles, originalPloidy, newGLPloidy, result);


        // add to new pool
        if (!Double.isInfinite(log10LofK))
            newPool.add(set);

        if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
            if ( VERBOSE )
                System.out.printf(" *** breaking early set=%s log10L=%.2f maxLog10L=%.2f%n", set.ACcounts, log10LofK, maxLog10L);
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        // by convention, ACcounts contained in set have full vector of possible pool ac counts including ref count.
        // so, if first element is zero, it automatically means we have no wiggle since we're in a corner of the conformation space
        final int ACwiggle = set.ACcounts.counts[0];
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;


        // add conformations for other cases
        for ( int allele = 1; allele < numAlleles; allele++ ) {
            final int[] ACcountsClone = set.ACcounts.getCounts().clone();
            ACcountsClone[allele]++;
            // is this a valid conformation?
            int altSum = (int)MathUtils.sum(ACcountsClone) - ACcountsClone[0];
            ACcountsClone[0] = newPloidy - altSum;
            if (ACcountsClone[0] < 0)
                continue;


            PoolGenotypeLikelihoods.updateACset(ACcountsClone, ACqueue, indexesToACset);
        }


        return log10LofK;
    }


    /**
     * Naive combiner of two multiallelic pools - number of alt alleles must be the same.
     * Math is generalization of biallelic combiner.
     *
     * For vector K representing an allele count conformation,
     * Pr(D | AC = K) = Sum_G Pr(D|AC1 = G) Pr (D|AC2=K-G) * F(G,K)
     * where F(G,K) = choose(m1,[g0 g1 ...])*choose(m2,[...]) / choose(m1+m2,[k1 k2 ...])
     * @param originalPool                    First log-likelihood pool GL vector
     * @param yy                    Second pool GL vector
     * @param ploidy1               Ploidy of first pool (# of chromosomes in it)
     * @param ploidy2               Ploidy of second pool
     * @param numAlleles            Number of alleles
     * @param log10AlleleFrequencyPriors Array of biallelic priors
     * @param result                Af calculation result object                  
     */
    public static void combineMultiallelicPoolNaively(CombinedPoolLikelihoods originalPool, double[] yy, int ploidy1, int ploidy2, int numAlleles,
                                                      final double[] log10AlleleFrequencyPriors,
                                                      final AlleleFrequencyCalculationResult result) {
/*
        final int dim1 = GenotypeLikelihoods.numLikelihoods(numAlleles, ploidy1);
        final int dim2 = GenotypeLikelihoods.numLikelihoods(numAlleles, ploidy2);

        if (dim1 != originalPool.getLength() || dim2 != yy.length)
            throw new ReviewedStingException("BUG: Inconsistent vector length");

        if (ploidy2 == 0)
            return;

        final int newPloidy = ploidy1 + ploidy2;

        // Say L1(K) = Pr(D|AC1=K) * choose(m1,K)
        // and L2(K) = Pr(D|AC2=K) * choose(m2,K)
        PoolGenotypeLikelihoods.SumIterator firstIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,ploidy1);
        final double[] x = originalPool.getLikelihoodsAsVector(true);
        while(firstIterator.hasNext()) {
            x[firstIterator.getLinearIndex()] += MathUtils.log10MultinomialCoefficient(ploidy1,firstIterator.getCurrentVector());
            firstIterator.next();
        }

        PoolGenotypeLikelihoods.SumIterator secondIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,ploidy2);
        final double[] y = yy.clone();
        while(secondIterator.hasNext()) {
            y[secondIterator.getLinearIndex()] += MathUtils.log10MultinomialCoefficient(ploidy2,secondIterator.getCurrentVector());
            secondIterator.next();
        }

        // initialize output to -log10(choose(m1+m2,[k1 k2...])
        final int outputDim = GenotypeLikelihoods.numLikelihoods(numAlleles, newPloidy);
        final PoolGenotypeLikelihoods.SumIterator outputIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,newPloidy);


        // Now, result(K) =  logSum_G (L1(G)+L2(K-G)) where G are all possible vectors that sum UP to K
        while(outputIterator.hasNext()) {
            final ExactACset set = new ExactACset(1, new ExactACcounts(outputIterator.getCurrentAltVector()));
            double likelihood = computeLofK(set, x,y, log10AlleleFrequencyPriors, numAlleles, ploidy1, ploidy2, result);

            originalPool.add(likelihood, set, outputIterator.getLinearIndex());
            outputIterator.next();
        }
*/
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
     * @param result                    AFResult object
     * @return                          log-likehood of requested conformation
     */
    private static double computeLofK(final ExactACset set,
                                      final CombinedPoolLikelihoods firstGLs,
                                      final double[] secondGL,
                                      final double[] log10AlleleFrequencyPriors,
                                      final int numAlleles, final int ploidy1, final int ploidy2,
                                      final AlleleFrequencyCalculationResult result) {

        final int newPloidy = ploidy1 + ploidy2;

        // sanity check
        int totalAltK = set.getACsum();
        if (newPloidy != totalAltK)
            throw new ReviewedStingException("BUG: inconsistent sizes of set.getACsum and passed ploidy values");

        totalAltK -= set.ACcounts.counts[0];
        // totalAltK has sum of alt alleles of conformation now


        // special case for k = 0 over all k
        if ( totalAltK == 0 ) {   // all-ref case
            final double log10Lof0 = firstGLs.getGLOfACZero() + secondGL[HOM_REF_INDEX];
            set.log10Likelihoods[0] = log10Lof0;

            result.setLog10LikelihoodOfAFzero(log10Lof0);
            result.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);

        }   else {

            // initialize result with denominator
            // ExactACset holds by convention the conformation of all alleles, and the sum of all allele count is just the ploidy.
            // To compute n!/k1!k2!k3!... we need to compute first n!/(k2!k3!...) and then further divide by k1! where k1=ploidy-sum_k_i

            int[] currentCount = set.ACcounts.getCounts();
            double denom =  -MathUtils.log10MultinomialCoefficient(newPloidy, currentCount);

            // for current conformation, get all possible ways to break vector K into two components G1 and G2
            final PoolGenotypeLikelihoods.SumIterator innerIterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles,ploidy2);
            set.log10Likelihoods[0] = Double.NEGATIVE_INFINITY;
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

                        set.log10Likelihoods[0] = MathUtils.approximateLog10SumLog10(set.log10Likelihoods[0], sum);
                    }
                }
                innerIterator.next();
            }

            set.log10Likelihoods[0] += denom;
        }

        double log10LofK = set.log10Likelihoods[0];

        // update the MLE if necessary
        final int altCounts[] = Arrays.copyOfRange(set.ACcounts.counts,1, set.ACcounts.counts.length);
        result.updateMLEifNeeded(log10LofK, altCounts);

        // apply the priors over each alternate allele
        for (final int ACcount : altCounts ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }
        result.updateMAPifNeeded(log10LofK, altCounts);

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
     * Combines naively two biallelic pools (of arbitrary size).
     * For two pools of size m1 and m2, we can compute the combined likelihood as:
     *   Pr(D|AC=k) = Sum_{j=0}^k Pr(D|AC1=j) Pr(D|AC2=k-j) * choose(m1,j)*choose(m2,k-j)/choose(m1+m2,k)
     * @param originalPool              Pool likelihood vector, x[k] = Pr(AC_i = k) for alt allele i
     * @param newPLVector               Second GL vector
     * @param ploidy1               Ploidy of first pool (# of chromosomes in it)
     * @param ploidy2               Ploidy of second pool
     * @param log10AlleleFrequencyPriors Array of biallelic priors
     * @param result                Af calculation result object
     * @return                Combined likelihood vector
     */
    public static ProbabilityVector combineBiallelicPoolsNaively(final ProbabilityVector originalPool, final double[] newPLVector,
                                                                 final int ploidy1, final int ploidy2, final double[] log10AlleleFrequencyPriors,
                                                                 final AlleleFrequencyCalculationResult result) {

        final int newPloidy = ploidy1 + ploidy2;

        final double[] combinedLikelihoods = new double[1+newPloidy];

        /** Pre-fill result array and incorporate weights into input vectors
         * Say L1(k) = Pr(D|AC1=k) * choose(m1,k)
         * and L2(k) = Pr(D|AC2=k) * choose(m2,k)
         * equation reduces to
         * Pr(D|AC=k) = 1/choose(m1+m2,k) * Sum_{j=0}^k L1(k) L2(k-j)
         * which is just plain convolution of L1 and L2 (with pre-existing vector)
         */

        // intialize result vector to -infinity
        Arrays.fill(combinedLikelihoods,Double.NEGATIVE_INFINITY);

        final double[] x = Arrays.copyOf(originalPool.getProbabilityVector(),1+ploidy1);
        for (int k=originalPool.getProbabilityVector().length; k< x.length; k++)
            x[k] = Double.NEGATIVE_INFINITY;

        final double[] y = newPLVector.clone();


        final double log10Lof0 = x[0]+y[0];
        result.setLog10LikelihoodOfAFzero(log10Lof0);
        result.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);

        double maxElement = log10Lof0;
        int maxElementIdx = 0;
        int[] alleleCounts = new int[1];
        for (int k= originalPool.getMinVal() ; k <= newPloidy; k++) {
            double[] acc = new double[k+1];
            Arrays.fill(acc,Double.NEGATIVE_INFINITY);
            double innerMax = Double.NEGATIVE_INFINITY;

            for (int j=0; j <=k; j++) {
                double x1,y1;


                if (k-j>=0 && k-j < y.length)
                    y1 = y[k-j] + MathUtils.log10BinomialCoefficient(ploidy2,k-j);
                else
                    continue;

                if (j < x.length)
                    x1 = x[j] + MathUtils.log10BinomialCoefficient(ploidy1,j);
                else
                    continue;

                if (Double.isInfinite(x1) || Double.isInfinite(y1))
                    continue;
                acc[j] = x1 + y1;
                if (acc[j] > innerMax)
                    innerMax = acc[j];
                else if (acc[j] < innerMax - MAX_LOG10_ERROR_TO_STOP_EARLY)
                    break;
            }
            combinedLikelihoods[k] = MathUtils.log10sumLog10(acc) - MathUtils.log10BinomialCoefficient(newPloidy,k);
            maxElementIdx = k;
            double maxDiff = combinedLikelihoods[k] - maxElement;
            if (maxDiff > 0)
                maxElement = combinedLikelihoods[k];
            else if (maxDiff < maxElement - MAX_LOG10_ERROR_TO_STOP_EARLY) {
                break;
            }

            alleleCounts[0] = k;
            result.updateMLEifNeeded(combinedLikelihoods[k],alleleCounts);
            result.updateMAPifNeeded(combinedLikelihoods[k] + log10AlleleFrequencyPriors[k],alleleCounts);


        }


        return new ProbabilityVector(MathUtils.normalizeFromLog10(Arrays.copyOf(combinedLikelihoods,maxElementIdx+1),false, true));
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
            if ( numOriginalAltAlleles == numNewAltAlleles) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = PoolGenotypeLikelihoods.subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(),allelesToUse);

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( MathUtils.sum(newLikelihoods) > VariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
            }
            else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                if ( numNewAltAlleles == 0 )
                    gb.noPL();
                else
                    gb.PL(newLikelihoods);

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || MathUtils.sum(newLikelihoods) > VariantContextUtils.SUM_GL_THRESH_NOCALL )
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
    private static void assignGenotype(final GenotypeBuilder gb,
                                       final double[] newLikelihoods,
                                       final List<Allele> allelesToUse,
                                       final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;



        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);

        final int[] mlAlleleCount = PoolGenotypeLikelihoods.getAlleleCountFromPLIndex(allelesToUse.size(), numChromosomes, PLindex);
        final ArrayList<Double> alleleFreqs = new ArrayList<Double>();
        final ArrayList<Integer> alleleCounts = new ArrayList<Integer>();


        for (int k=1; k < mlAlleleCount.length; k++) {
            alleleCounts.add(mlAlleleCount[k]);
            final double freq = (double)mlAlleleCount[k] / (double)numChromosomes;
            alleleFreqs.add(freq);

        }

        // per-pool logging of AC and AF
        gb.attribute(VCFConstants.MLE_ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
        gb.attribute(VCFConstants.MLE_ALLELE_FREQUENCY_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);

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
