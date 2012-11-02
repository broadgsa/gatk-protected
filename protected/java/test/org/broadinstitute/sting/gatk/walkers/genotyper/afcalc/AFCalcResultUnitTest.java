package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AFCalcResultUnitTest extends BaseTest {
    private static class MyTest {
        final double[] Ls, expectedPosteriors;

        private MyTest(double[] ls, double[] expectedPosteriors) {
            Ls = ls;
            this.expectedPosteriors = expectedPosteriors;
        }

        @Override
        public String toString() {
            return "Ls [" + Utils.join(",", Ls) + "] expectedPosteriors [" + Utils.join(",", expectedPosteriors) + "]";
        }
    }

    @DataProvider(name = "TestComputePosteriors")
    public Object[][] makeTestCombineGLs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new MyTest(log10Even, log10Even)});

        for ( double L0 = -1e9; L0 < 0.0; L0 /= 10.0 ) {
            for ( double L1 = -1e2; L1 < 0.0; L1 /= 100.0 ) {
                final double[] input = new double[]{L0, L1};
                final double[] expected = MathUtils.normalizeFromLog10(input, true);
                tests.add(new Object[]{new MyTest(input, expected)});
            }
        }

        for ( double bigBadL = -1e50; bigBadL < -1e200; bigBadL *= 10 ) {
            // test that a huge bad likelihood remains, even with a massive better result
            for ( final double betterL : Arrays.asList(-1000.0, -100.0, -10.0, -1.0, -0.1, -0.01, -0.001, 0.0)) {
                tests.add(new Object[]{new MyTest(new double[]{bigBadL, betterL}, new double[]{bigBadL, 0.0})});
                tests.add(new Object[]{new MyTest(new double[]{betterL, bigBadL}, new double[]{0.0, bigBadL})});
            }
        }

        // test that a modest bad likelihood with an ~0.0 value doesn't get lost
        for ( final double badL : Arrays.asList(-10000.0, -1000.0, -100.0, -10.0)) {
            tests.add(new Object[]{new MyTest(new double[]{badL, -1e-9}, new double[]{badL, 0.0})});
            tests.add(new Object[]{new MyTest(new double[]{-1e-9, badL}, new double[]{0.0, badL})});
        }

        // test that a non-ref site gets reasonable posteriors with an ~0.0 value doesn't get lost
        for ( final double nonRefL : Arrays.asList(-100.0, -50.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0)) {
            tests.add(new Object[]{new MyTest(new double[]{0.0, nonRefL}, new double[]{0.0, nonRefL})});
        }

        return tests.toArray(new Object[][]{});
    }


    final static double[] log10Even = MathUtils.normalizeFromLog10(new double[]{0.5, 0.5}, true);
    final static Allele C = Allele.create("C");
    final static List<Allele> alleles = Arrays.asList(Allele.create("A", true), C);

    @Test(enabled = true, dataProvider = "TestComputePosteriors")
    private void testComputingPosteriors(final MyTest data) {
        final AFCalcResult result = new AFCalcResult(new int[]{0}, 1, alleles, data.Ls, log10Even, Collections.singletonMap(C, -1.0));

        Assert.assertEquals(result.getLog10PosteriorOfAFEq0(), data.expectedPosteriors[0], 1e-3, "AF = 0 not expected");
        Assert.assertEquals(result.getLog10PosteriorOfAFGT0(), data.expectedPosteriors[1], 1e-3, "AF > 0 not expected");

        final double[] actualPosteriors = new double[]{result.getLog10PosteriorOfAFEq0(), result.getLog10PosteriorOfAFGT0()};
        Assert.assertEquals(MathUtils.sumLog10(actualPosteriors), 1.0, 1e-3, "Posteriors don't sum to 1 with 1e-3 precision");
    }
}
