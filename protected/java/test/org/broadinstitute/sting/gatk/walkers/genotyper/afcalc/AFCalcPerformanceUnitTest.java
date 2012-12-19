package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AFCalcPerformanceUnitTest extends BaseTest {
    @DataProvider(name = "ScalingTests")
    public Object[][] makepolyTestProviderLotsOfAlleles() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // list of all high-quality models in the system
        final List<AFCalcFactory.Calculation> biAllelicModels = Arrays.asList(
                AFCalcFactory.Calculation.EXACT_INDEPENDENT,
                AFCalcFactory.Calculation.EXACT_REFERENCE);

        final List<AFCalcFactory.Calculation> multiAllelicModels = Arrays.asList(
                AFCalcFactory.Calculation.EXACT_INDEPENDENT);

//        for ( final int nonTypePLs : Arrays.asList(100) ) {
//            for ( final int nSamples : Arrays.asList(10000) ) {
//                final List<Integer> alleleCounts = Arrays.asList(50);
//                for ( final int nAltAlleles : Arrays.asList(1) ) {
        for ( final int nonTypePLs : Arrays.asList(100) ) {
            for ( final int nSamples : Arrays.asList(100, 1000) ) {
                final List<Integer> alleleCounts = Arrays.asList(0, 1, 2, 3, 4, 5, 10, 50, 500);
                    for ( final int nAltAlleles : Arrays.asList(1, 2, 3) ) {
                    final List<AFCalcFactory.Calculation> models = nAltAlleles > 1 ? multiAllelicModels : biAllelicModels;
                    for ( final AFCalcFactory.Calculation model : models ) {
                        for ( final List<Integer> ACs : Utils.makePermutations(alleleCounts, nAltAlleles, true) ) {
                            if ( MathUtils.sum(ACs) < nSamples * 2 ) {
                                final AFCalcTestBuilder testBuilder
                                        = new AFCalcTestBuilder(nSamples, nAltAlleles, model, AFCalcTestBuilder.PriorType.human);
                                tests.add(new Object[]{testBuilder, ACs, nonTypePLs});
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private Pair<Integer, Integer> estNumberOfEvaluations(final AFCalcTestBuilder testBuilder, final VariantContext vc, final int nonTypePL) {
        final int evalOverhead = 2; // 2
        final int maxEvalsPerSamplePerAC = 3;

        int minEvals = 0, maxEvals = 0;

        for ( final Allele alt : vc.getAlternateAlleles() ) {
            final int AC = vc.getCalledChrCount(alt);
            minEvals += AC + evalOverhead; // everyone is hom-var
            maxEvals += AC * maxEvalsPerSamplePerAC + 10;
        }

        return new Pair<Integer, Integer>(minEvals, maxEvals);
    }

    @Test(dataProvider = "ScalingTests")
    private void testScaling(final AFCalcTestBuilder testBuilder, final List<Integer> ACs, final int nonTypePL) {
        final AFCalc calc = testBuilder.makeModel();
        final double[] priors = testBuilder.makePriors();
        final VariantContext vc = testBuilder.makeACTest(ACs, 0, nonTypePL);
        final AFCalcResult result = calc.getLog10PNonRef(vc, priors);
        final Pair<Integer, Integer> expectedNEvaluation = estNumberOfEvaluations(testBuilder, vc, nonTypePL);
        final int minEvals = expectedNEvaluation.getFirst();
        final int maxEvals = expectedNEvaluation.getSecond();

        logger.warn(" min " + minEvals + " obs " + result.getnEvaluations() + " max " + maxEvals + " for test " + testBuilder + " sum(ACs)=" + (int)MathUtils.sum(ACs));

        Assert.assertTrue(result.getnEvaluations() >= minEvals,
                "Actual number of evaluations " + result.getnEvaluations() + " < min number of evals " + minEvals);
        Assert.assertTrue(result.getnEvaluations() <= maxEvals,
                "Actual number of evaluations " + result.getnEvaluations() + " > max number of evals " + minEvals);
    }
}