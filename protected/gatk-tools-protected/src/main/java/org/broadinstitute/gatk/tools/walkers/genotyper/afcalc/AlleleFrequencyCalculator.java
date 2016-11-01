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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.gatk.utils.Dirichlet;
import org.broadinstitute.gatk.utils.IndexRange;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFrequencyCalculator extends AFCalculator {
    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();
    private static final double THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE = 0.1;
    private static final int HOM_REF_GENOTYPE_INDEX = 0;

    private final double refPseudocount;
    private final double snpPseudocount;
    private final double indelPseudocount;
    private final int defaultPloidy;

    // this doesn't use the exact model, so number of evaluations is irrelevant
    private static final int DUMMY_N_EVALUATIONS = 1;


    public AlleleFrequencyCalculator(final double refPseudocount, final double snpPseudocount, final double indelPseudocount, final int defaultPloidy) {
        this.refPseudocount = refPseudocount;
        this.snpPseudocount = snpPseudocount;
        this.indelPseudocount = indelPseudocount;
        this.defaultPloidy = defaultPloidy;
    }

    public AFCalculationResult getLog10PNonRef(final VariantContext vc) {
        // maxAltAlleles is not used by getLog10PNonRef, so don't worry about the 0
        return getLog10PNonRef(vc, defaultPloidy, 0, null);
    }
    //TODO: this should be a class of static methods once the old AFCalculator is gone.
    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @param refSnpIndelPseudocounts a total hack.  A length-3 vector containing Dirichlet prior pseudocounts to
     *                                be given to ref, alt SNP, and alt indel alleles.  Hack won't be necessary when we destroy the old AF calculators
     * @return result (for programming convenience)
     */
    @Override
    public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles, final double[] refSnpIndelPseudocounts) {
        Utils.nonNull(vc, "vc is null");
        final int numAlleles = vc.getNAlleles();
        final List<Allele> alleles = vc.getAlleles();
        Utils.validateArg( numAlleles > 1, "VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);

        final double[] priorPseudocounts = alleles.stream()
                .mapToDouble(a -> a.isReference() ? refPseudocount : (a.length() > 1 ? snpPseudocount : indelPseudocount)).toArray();

        double[] alleleCounts = new double[numAlleles];
        final double flatLog10AlleleFrequency = -MathUtils.Log10Cache.get(numAlleles); // log10(1/numAlleles)
        double[] log10AlleleFrequencies = new IndexRange(0, numAlleles).mapToDouble(n -> flatLog10AlleleFrequency);
        double alleleCountsMaximumDifference = Double.POSITIVE_INFINITY;

        while (alleleCountsMaximumDifference > THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE) {
            final double[] newAlleleCounts = effectiveAlleleCounts(vc, log10AlleleFrequencies);
            alleleCountsMaximumDifference = Arrays.stream(MathArrays.ebeSubtract(alleleCounts, newAlleleCounts)).map(Math::abs).max().getAsDouble();
            alleleCounts = newAlleleCounts;
            final double[] posteriorPseudocounts = MathArrays.ebeAdd(priorPseudocounts, alleleCounts);

            // first iteration uses flat prior in order to avoid local minimum where the prior + no pseudocounts gives such a low
            // effective allele frequency that it overwhelms the genotype likelihood of a real variant
            // basically, we want a chance to get non-zero pseudocounts before using a prior that's biased against a variant
            log10AlleleFrequencies = new Dirichlet(posteriorPseudocounts).log10MeanWeights();
        }

        double[] log10POfZeroCountsByAllele = new double[numAlleles];
        double log10PNoVariant = 0;

        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods()) {
                continue;
            }
            final int ploidy = g.getPloidy() == 0 ? defaultPloidy : g.getPloidy();
            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, numAlleles);

            final double[] log10GenotypePosteriors = log10NormalizedGenotypePosteriors(g, glCalc, log10AlleleFrequencies);

            //the total probability
            log10PNoVariant += log10GenotypePosteriors[HOM_REF_GENOTYPE_INDEX];

            // per allele non-log space probabilities of zero counts for this sample
            // for each allele calculate the total probability of genotypes containing at least one copy of the allele
            final double[] log10ProbabilityOfNonZeroAltAlleles = new double[numAlleles];
            Arrays.fill(log10ProbabilityOfNonZeroAltAlleles, Double.NEGATIVE_INFINITY);

            for (int genotype = 0; genotype < glCalc.genotypeCount(); genotype++) {
                final double log10GenotypePosterior = log10GenotypePosteriors[genotype];
                glCalc.genotypeAlleleCountsAt(genotype).forEachAlleleIndexAndCount((alleleIndex, count) ->
                        log10ProbabilityOfNonZeroAltAlleles[alleleIndex] =
                                MathUtils.log10SumLog10(log10ProbabilityOfNonZeroAltAlleles[alleleIndex], log10GenotypePosterior));
            }

            for (int allele = 0; allele < numAlleles; allele++) {
                // if prob of non hom ref == 1 up to numerical precision, short-circuit to avoid NaN
                if (log10ProbabilityOfNonZeroAltAlleles[allele] >= 0) {
                    log10POfZeroCountsByAllele[allele] = Double.NEGATIVE_INFINITY;
                } else {
                    log10POfZeroCountsByAllele[allele] += MathUtils.log10OneMinusPow10(log10ProbabilityOfNonZeroAltAlleles[allele]);
                }
            }
        }

        // unfortunately AFCalculationResult expects integers for the MLE.  We really should emit the EM no-integer values
        // which are valuable (eg in CombineGVCFs) as the sufficient statistics of the Dirichlet posterior on allele frequencies
        final int[] integerAlleleCounts = Arrays.stream(alleleCounts).mapToInt(x -> (int) Math.round(x)).toArray();
        final int[] integerAltAlleleCounts = Arrays.copyOfRange(integerAlleleCounts, 1, numAlleles);

        //skip the ref allele (index 0)
        final Map<Allele, Double> log10PRefByAllele = IntStream.range(1, numAlleles).boxed()
                .collect(Collectors.toMap(alleles::get, a -> log10POfZeroCountsByAllele[a]));

        // we compute posteriors here and don't have the same prior that AFCalculationResult expects.  Therefore, we
        // give it our posterior as its "likelihood" along with a flat dummy prior
        final double[] dummyFlatPrior = {-1e-10, -1e-10};   //TODO: HACK must be negative for AFCalcResult
        final double[] log10PosteriorOfNoVariantYesVariant = {log10PNoVariant, MathUtils.log10OneMinusPow10(log10PNoVariant)};

        return new AFCalculationResult(integerAltAlleleCounts, DUMMY_N_EVALUATIONS, alleles, log10PosteriorOfNoVariantYesVariant, dummyFlatPrior, log10PRefByAllele);
    }

    // effectiveAlleleCounts[allele a] = SUM_{genotypes g} (posterior_probability(g) * num_copies of a in g), which we denote as SUM [n_g p_g]
    // for numerical stability we will do this in log space:
    // count = SUM 10^(log (n_g p_g)) = SUM 10^(log n_g + log p_g)
    // thanks to the log-sum-exp trick this lets us work with log posteriors alone
    private double[] effectiveAlleleCounts(final VariantContext vc, final double[] log10AlleleFrequencies) {
        final int numAlleles = vc.getNAlleles();
        Utils.validateArg(numAlleles == log10AlleleFrequencies.length, "number of alleles inconsistent");
        final double[] log10Result = new double[numAlleles];
        Arrays.fill(log10Result, Double.NEGATIVE_INFINITY);
        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods()) {
                continue;
            }
            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(g.getPloidy(), numAlleles);
            final double[] log10GenotypePosteriors = log10NormalizedGenotypePosteriors(g, glCalc, log10AlleleFrequencies);

            new IndexRange(0, glCalc.genotypeCount()).forEach(genotypeIndex ->
                    glCalc.genotypeAlleleCountsAt(genotypeIndex).forEachAlleleIndexAndCount((alleleIndex, count) ->
                            log10Result[alleleIndex] = MathUtils.log10SumLog10(log10Result[alleleIndex], log10GenotypePosteriors[genotypeIndex] + MathUtils.Log10Cache.get(count))));
        }
        return MathUtils.applyToArrayInPlace(log10Result, x -> Math.pow(10.0, x));
    }

    private static double[] log10NormalizedGenotypePosteriors(final Genotype g, final GenotypeLikelihoodCalculator glCalc, final double[] log10AlleleFrequencies) {
        final double[] log10Likelihoods = g.getLikelihoods().getAsVector();
        final double[] log10Posteriors = new IndexRange(0, glCalc.genotypeCount()).mapToDouble(genotypeIndex -> {
            final GenotypeAlleleCounts gac = glCalc.genotypeAlleleCountsAt(genotypeIndex);
            return gac.log10CombinationCount() + log10Likelihoods[genotypeIndex]
                    + gac.sumOverAlleleIndicesAndCounts((index, count) -> count * log10AlleleFrequencies[index]);
        });
        return MathUtils.normalizeFromLog10(log10Posteriors, true);
    }

    @Override   //Note: unused
    protected AFCalculationResult getResultFromFinalState(final VariantContext vc, final double[] priors, final StateTracker st) { return null; }

    @Override//Note: unused
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                      final double[] priors, final StateTracker st) { return null; }

    @Override   //Note: unused
    protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) { return null; }

    @Override   //trivial implementation -- new AFCalculator can handle multiallelics so we're not afraid
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        return vc;
    }

    @Override   //also trivial
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                                   final int defaultPloidy,
                                                   final List<Allele> allelesToUse,
                                                   final boolean assignGenotypes) {
        return vc.getGenotypes();
    }


}
