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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Computes the conditional bi-allelic exact results
 *
 * Suppose vc contains 2 alt allele: A* with C and T.  This function first computes:
 *
 * (1) P(D | AF_c > 0 && AF_t == *) [i.e., T can be anything]
 *
 * it then computes the conditional probability on AF_c == 0:
 *
 * (2) P(D | AF_t > 0 && AF_c == 0)
 *
 * Thinking about this visually, we have the following likelihood matrix where each cell is
 * the P(D | AF_c == i && AF_t == j):
 *
 *     0 AF_c > 0
 *    -----------------
 * 0  |  |
 *    |--|-------------
 * a  |  |
 * f  |  |
 * _  |  |
 * t  |  |
 * >  |  |
 * 0  |  |
 *
 * What we really want to know how
 *
 * (3) P(D | AF_c == 0 & AF_t == 0)
 *
 * compares with
 *
 * (4) P(D | AF_c > 0 || AF_t > 0)
 *
 * This is effectively asking for the value in the upper left vs. the sum of all cells.
 *
 * This class implements the conditional likelihoods summation for any number of alt
 * alleles, where each alt allele has its EXACT probability of segregating calculated by
 * reducing each alt B into the case XB and computing P(D | AF_b > 0 ) as follows:
 *
 * Suppose we have for a A/B/C site the following GLs:
 *
 * AA AB BB AC BC CC
 *
 * and we want to get the bi-allelic GLs for X/B, where X is everything not B
 *
 * XX = AA + AC + CC (since X = A or C)
 * XB = AB + BC
 * BB = BB
 *
 * The posterior of the site being a variant site is calculated using
 * the likelihood of the AF whe all alternatives are collapsed to be zero.
 */
 public class IndependentAllelesDiploidExactAFCalculator extends DiploidExactAFCalculator {

    private final static int[] BIALLELIC_NON_INFORMATIVE_PLS = new int[]{0,0,0};
    private final static List<Allele> BIALLELIC_NOCALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Sorts AFCalcResults by their posteriors of AF > 0, so the
     */
    private final static class CompareAFCalculatorResultsByPNonRef implements Comparator<AFCalculationResult> {
        @Override
        public int compare(AFCalculationResult o1, AFCalculationResult o2) {
            return -1 * Double.compare(o1.getLog10PosteriorOfAFGT0(), o2.getLog10PosteriorOfAFGT0());
        }
    }

    private final static CompareAFCalculatorResultsByPNonRef compareAFCalcResultsByPNonRef = new CompareAFCalculatorResultsByPNonRef();

    /**
     * The AFCalc model we are using to do the bi-allelic computation
     */
    final AFCalculator biAlleleExactModel;

    protected IndependentAllelesDiploidExactAFCalculator() {
        super();
        biAlleleExactModel = new ReferenceDiploidExactAFCalculator();
    }

    @Override
    public AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                            final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        final List<AFCalculationResult> independentResultTrackers = computeAlleleIndependentExact(vc, defaultPloidy, log10AlleleFrequencyPriors);

        if ( independentResultTrackers.size() == 0 )
            throw new IllegalStateException("Independent alleles model returned an empty list of results at VC " + vc);

        if ( independentResultTrackers.size() == 1 ) {
            // fast path for the very common bi-allelic use case
            return independentResultTrackers.get(0);
        } else {
            final AFCalculationResult combinedAltAllelesResult = combineAltAlleleIndependentExact(vc,defaultPloidy,log10AlleleFrequencyPriors);
            // we are a multi-allelic, so we need to actually combine the results
            final List<AFCalculationResult> withMultiAllelicPriors = applyMultiAllelicPriors(independentResultTrackers);
            return combineIndependentPNonRefs(vc, withMultiAllelicPriors, combinedAltAllelesResult);
        }
    }

    private AFCalculationResult combineAltAlleleIndependentExact(final VariantContext vc, int defaultPloidy, double[] log10AlleleFrequencyPriors) {
        final VariantContext combinedAltAllelesVariantContext = makeCombinedAltAllelesVariantContext(vc);
        final AFCalculationResult resultTracker = biAlleleExactModel.getLog10PNonRef(combinedAltAllelesVariantContext, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
        return resultTracker;
    }

    private VariantContext makeCombinedAltAllelesVariantContext(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;

        if ( nAltAlleles == 1 )
            return vc;
        else {
            final VariantContextBuilder vcb = new VariantContextBuilder(vc);
            final Allele reference = vcb.getAlleles().get(0);
            vcb.alleles(Arrays.asList(reference, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE));
            final int genotypeCount = GenotypeLikelihoodCalculators.genotypeCount(2, vc.getNAlleles());
            final double[] hetLikelihoods = new double[vc.getNAlleles() - 1];
            final double[] homAltLikelihoods = new double[genotypeCount - hetLikelihoods.length - 1];
            final double[] newLikelihoods = new double[3];
            final List<Genotype> newGenotypes = new ArrayList<>(vc.getNSamples());
            for (final Genotype oldGenotype : vc.getGenotypes()) {
                final GenotypeBuilder gb = new GenotypeBuilder(oldGenotype);
                final List<Allele> oldAlleles = oldGenotype.getAlleles();
                if (oldAlleles != null) {
                    final List<Allele> newAlleles = new ArrayList<>(oldAlleles.size());
                    for (int i = 0; i < oldAlleles.size(); i++) {
                        final Allele oldAllele = oldAlleles.get(i);
                        if (oldAllele.isReference())
                            newAlleles.add(reference);
                        else if (oldAllele.isNoCall())
                            newAlleles.add(Allele.NO_CALL);
                        else
                            newAlleles.add(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                    }
                    gb.alleles(newAlleles);
                }
                if (combineAltAlleleLikelihoods(oldGenotype, genotypeCount, newLikelihoods, hetLikelihoods, homAltLikelihoods))
                    gb.PL(newLikelihoods);
                newGenotypes.add(gb.make());
            }
            return vcb.genotypes(newGenotypes).make();
        }
    }

    /**
     * Compute the conditional exact AFCalcResult for each allele in vc independently, returning
     * the result of each, in order of the alt alleles in VC
     *
     * @param vc the VariantContext we want to analyze, with at least 1 alt allele
     * @param log10AlleleFrequencyPriors the priors
     * @return a list of the AFCalcResults for each bi-allelic sub context of vc
     */
    @Requires({"vc != null", "log10AlleleFrequencyPriors != null"})
    @Ensures("goodIndependentResult(vc, result)")
    protected final List<AFCalculationResult> computeAlleleIndependentExact(final VariantContext vc, final int defaultPloidy,
                                                                     final double[] log10AlleleFrequencyPriors) {
        final List<AFCalculationResult> results = new LinkedList<AFCalculationResult>();

        for ( final VariantContext subvc : makeAlleleConditionalContexts(vc) ) {
            final AFCalculationResult resultTracker = biAlleleExactModel.getLog10PNonRef(subvc, defaultPloidy, vc.getNAlleles() - 1, log10AlleleFrequencyPriors);
            results.add(resultTracker);
        }

        return results;
    }

    /**
     * Returns the bi-allelic variant context for each alt allele in vc with bi-allelic likelihoods, in order
     *
     * @param vc the variant context to split.  Must have n.alt.alleles > 1
     * @return a bi-allelic variant context for each alt allele in vc
     */
    @Requires({"vc != null", "vc.getNAlleles() > 1"})
    @Ensures("result.size() == vc.getNAlleles() - 1")
    protected final List<VariantContext> makeAlleleConditionalContexts(final VariantContext vc) {
        final int nAltAlleles = vc.getNAlleles() - 1;

        if ( nAltAlleles == 1 ) {
            // fast path for bi-allelic case.
            return Collections.singletonList(vc);
        } else {
            // go through the work of ripping up the VC into its biallelic components
            final List<VariantContext> vcs = new LinkedList<VariantContext>();

            for ( int altI = 0; altI < nAltAlleles; altI++ ) {
                vcs.add(biallelicCombinedGLs(vc, altI + 1));
            }

            return vcs;
        }
    }

    /**
     * Create a single bi-allelic variant context from rootVC with alt allele with index altAlleleIndex
     *
     * @param rootVC the root (potentially multi-allelic) variant context
     * @param altAlleleIndex index of the alt allele, from 0 == first alt allele
     * @return a bi-allelic variant context based on rootVC
     */
    @Requires({"rootVC.getNAlleles() > 1", "altAlleleIndex < rootVC.getNAlleles()"})
    @Ensures({"result.isBiallelic()"})
    protected final VariantContext biallelicCombinedGLs(final VariantContext rootVC, final int altAlleleIndex) {
        if ( rootVC.isBiallelic() ) {
            return rootVC;
        } else {
            final int nAlts = rootVC.getNAlleles() - 1;
            final List<Genotype> biallelicGenotypes = new ArrayList<Genotype>(rootVC.getNSamples());
            for ( final Genotype g : rootVC.getGenotypes() )
                biallelicGenotypes.add(combineGLsPrecise(g, altAlleleIndex, nAlts));

            final VariantContextBuilder vcb = new VariantContextBuilder(rootVC);
            final Allele altAllele = rootVC.getAlternateAllele(altAlleleIndex - 1);
            vcb.alleles(Arrays.asList(rootVC.getReference(), altAllele));
            vcb.genotypes(biallelicGenotypes);
            return vcb.make();
        }
    }

    private static final double PHRED_2_LOG10_COEFF = -.1;

    /**
     * Returns a new Genotype with the PLs of the multi-allelic original reduced to a bi-allelic case.
     *
     * <p>Uses the log-sum-exp trick in order to work well with very low PLs</p>
     *
     * <p>This is handled in the following way:</p>
     *
     * <p>Suppose we have for a A/B/C site the following GLs:</p>
     *
     * <p>AA AB BB AC BC CC</p>
     *
     * <p>and we want to get the bi-allelic GLs for X/B, where X is everything not B</p>
     *
     * <p>XX = AA + AC + CC (since X = A or C)<br/>
     * XB = AB + BC                           <br/>
     * BB = BB                                     <br/>
     * </p>
     * <p>
     *     This implementation use the log sum trick in order to avoid numeric inestability.
     * </p>
     *
     * @param original the original multi-allelic genotype
     * @param altIndex the index of the alt allele we wish to keep in the bialleic case -- with ref == 0
     * @param nAlts the total number of alt alleles
     * @return a new biallelic genotype with appropriate PLs
     */
    @Requires({"original.hasLikelihoods()"})
    @Ensures({"result.hasLikelihoods()", "result.getPL().length == 3"})
    protected Genotype combineGLsPrecise(final Genotype original, final int altIndex, final int nAlts ) {

        if ( original.isNonInformative() )
            return new GenotypeBuilder(original).PL(BIALLELIC_NON_INFORMATIVE_PLS).alleles(BIALLELIC_NOCALL).make();

        if ( altIndex < 1 || altIndex > nAlts ) throw new IllegalStateException("altIndex must be between 1 and nAlts " + nAlts);

        final int[] pls = original.getPL();

        final int nAlleles = nAlts + 1;

        final int plCount = pls.length;

        double BB = 0;
        final double[] XBvalues = new double[nAlleles - 1];
        final double[] XXvalues = new double[plCount - nAlleles];

        int xbOffset = 0;
        int xxOffset = 0;
        for ( int index = 0; index < plCount; index++ ) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(index);
            int i = pair.alleleIndex1;
            int j = pair.alleleIndex2;
            if (i == j) {
              if (i == altIndex) BB = PHRED_2_LOG10_COEFF * pls[index]; else XXvalues[xxOffset++] = PHRED_2_LOG10_COEFF * pls[index];
            } else if (i == altIndex || j == altIndex)
              XBvalues[xbOffset++] = PHRED_2_LOG10_COEFF * pls[index];
            else
              XXvalues[xxOffset++] = PHRED_2_LOG10_COEFF * pls[index];
        }

        final double XB = MathUtils.log10sumLog10(XBvalues);
        final double XX = MathUtils.log10sumLog10(XXvalues);

        final double[] GLs = new double[] { XX, XB, BB};
        return new GenotypeBuilder(original).PL(GLs).alleles(BIALLELIC_NOCALL).make();
    }

    protected final List<AFCalculationResult> applyMultiAllelicPriors(final List<AFCalculationResult> conditionalPNonRefResults) {
        final ArrayList<AFCalculationResult> sorted = new ArrayList<AFCalculationResult>(conditionalPNonRefResults);

        // sort the results, so the most likely allele is first
        Collections.sort(sorted, compareAFCalcResultsByPNonRef);

        double lastPosteriorGt0 = sorted.get(0).getLog10PosteriorOfAFGT0();
        final double log10SingleAllelePriorOfAFGt0 = conditionalPNonRefResults.get(0).getLog10PriorOfAFGT0();

        for ( int i = 0; i < sorted.size(); i++ ) {
            if ( sorted.get(i).getLog10PosteriorOfAFGT0() > lastPosteriorGt0 )
                throw new IllegalStateException("pNonRefResults not sorted: lastPosteriorGt0 " + lastPosteriorGt0 + " but current is " + sorted.get(i).getLog10PosteriorOfAFGT0());

                final double log10PriorAFGt0 = (i + 1) * log10SingleAllelePriorOfAFGt0;
            final double log10PriorAFEq0 = Math.log10(1 - Math.pow(10, log10PriorAFGt0));
            final double[] thetaTONPriors = new double[] { log10PriorAFEq0, log10PriorAFGt0 };

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            sorted.set(i, sorted.get(i).withNewPriors(MathUtils.normalizeFromLog10(thetaTONPriors, true)));
        }

        return sorted;
    }

    /**
     * Take the independent estimates of pNonRef for each alt allele and combine them into a single result
     *
     * Given n independent calculations for each of n alternate alleles create a single
     * combined AFCalcResult with:
     *
     * priors for AF == 0 equal to theta^N for the nth least likely allele
     * posteriors that reflect the combined chance that any alleles are segregating and corresponding
     * likelihoods
     * combined MLEs in the order of the alt alleles in vc
     *
     * @param sortedResultsWithThetaNPriors the pNonRef result for each allele independently
     */
    protected AFCalculationResult combineIndependentPNonRefs(final VariantContext vc,
                                                      final List<AFCalculationResult> sortedResultsWithThetaNPriors,
                                                      final AFCalculationResult combinedAltAllelesResult) {


        int nEvaluations = 0;
        final int nAltAlleles = sortedResultsWithThetaNPriors.size();
        final int[] alleleCountsOfMLE = new int[nAltAlleles];
        final Map<Allele, Double> log10pRefByAllele = new HashMap<Allele, Double>(nAltAlleles);

        // the sum of the log10 posteriors for AF == 0 and AF > 0 to determine joint probs

        for ( final AFCalculationResult sortedResultWithThetaNPriors : sortedResultsWithThetaNPriors ) {
            final Allele altAllele = sortedResultWithThetaNPriors.getAllelesUsedInGenotyping().get(1);
            final int altI = vc.getAlleles().indexOf(altAllele) - 1;

            // MLE of altI allele is simply the MLE of this allele in altAlleles
            alleleCountsOfMLE[altI] = sortedResultWithThetaNPriors.getAlleleCountAtMLE(altAllele);

            // bind pNonRef for allele to the posterior value of the AF > 0 with the new adjusted prior
            log10pRefByAllele.put(altAllele, sortedResultWithThetaNPriors.getLog10PosteriorOfAFEq0());

            // trivial -- update the number of evaluations
            nEvaluations += sortedResultWithThetaNPriors.nEvaluations;
        }

        return new IndependentAlleleAFCalculationResult(alleleCountsOfMLE, nEvaluations, vc.getAlleles(),
                // necessary to ensure all values < 0
                MathUtils.normalizeFromLog10(new double[] { combinedAltAllelesResult.getLog10LikelihoodOfAFEq0(), combinedAltAllelesResult.getLog10LikelihoodOfAFGT0() }, true),
                // priors incorporate multiple alt alleles, must be normalized
                MathUtils.normalizeFromLog10(new double[] { combinedAltAllelesResult.getLog10PriorOfAFEq0(), combinedAltAllelesResult.getLog10PriorOfAFGT0() }, true),
                log10pRefByAllele, sortedResultsWithThetaNPriors);
    }

    private boolean combineAltAlleleLikelihoods(final Genotype g, final int plMaxIndex, final double[] dest,
                                                final double[] hetLikelihoods, final double[] homAltLikelihoods) {

        final int[] pls = g.getPL();
        if (pls == null)
            return false;
        int hetNextIndex = 0;
        int homAltNextIndex = 0;
        for (int plIndex = 1; plIndex < plMaxIndex; plIndex++) {
            final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(plIndex);
            if (alleles.alleleIndex1 == 0 || alleles.alleleIndex2 == 0)
                hetLikelihoods[hetNextIndex++] = pls[plIndex] * PHRED_2_LOG10_COEFF;
            else
                homAltLikelihoods[homAltNextIndex++] = pls[plIndex] * PHRED_2_LOG10_COEFF;
        }
        dest[0] = pls[0] * PHRED_2_LOG10_COEFF;
        dest[1] = MathUtils.approximateLog10SumLog10(hetLikelihoods);
        dest[2] = MathUtils.approximateLog10SumLog10(homAltLikelihoods);
        return true;
    }
}
