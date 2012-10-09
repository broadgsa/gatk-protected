package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


// SEE  private/R/pls.R if you want the truth output for these tests
public class IndependentAllelesDiploidExactAFCalcUnitTest extends BaseTest {
    @DataProvider(name = "TestCombineGLs")
    public Object[][] makeTestCombineGLs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{1, 1, makePL( 0, 10, 20), makePL( 0, 10, 20)});
        tests.add(new Object[]{1, 1, makePL(10,  0, 20), makePL(10,  0, 20)});
        tests.add(new Object[]{1, 1, makePL(20, 10,  0), makePL(20, 10,  0)});

        // AA AB BB AC BC CC => AA AB+BC CC
        tests.add(new Object[]{1, 2, makePL( 0, 10, 20, 30, 40, 50), makePL(0, 10, 20)});
        tests.add(new Object[]{2, 2, makePL( 0, 10, 20, 30, 40, 50), makePL(0, 30, 50)});

        tests.add(new Object[]{1, 2, makePL( 0, 10, 10, 10, 10, 10), makePL(0, 8, 11)});
        tests.add(new Object[]{2, 2, makePL( 0, 10, 10, 10, 10, 10), makePL(0, 8, 11)});

        tests.add(new Object[]{1, 2, makePL( 0, 1, 2, 3, 4, 5), makePL(0, 2, 5)});
        tests.add(new Object[]{2, 2, makePL( 0, 1, 2, 3, 4, 5), makePL(0, 4, 9)});

        tests.add(new Object[]{1, 2, makePL(  0, 50, 50, 50, 50, 50), makePL( 0, 47, 50)});
        tests.add(new Object[]{2, 2, makePL(  0, 50, 50, 50, 50, 50), makePL( 0, 47, 50)});

        tests.add(new Object[]{1, 2, makePL( 50,  0, 50, 50, 50, 50), makePL(45, 0, 50)});
        tests.add(new Object[]{2, 2, makePL( 50,  0, 50, 50, 50, 50), makePL( 0, 47, 50)});

        tests.add(new Object[]{1, 2, makePL( 50, 50, 0, 50, 50, 50), makePL(45, 47,  0)});
        tests.add(new Object[]{2, 2, makePL( 50, 50, 0, 50, 50, 50), makePL( 0, 47, 50)});

        tests.add(new Object[]{1, 2, makePL( 50, 50, 50,  0, 50, 50), makePL(0, 47, 50)});
        tests.add(new Object[]{2, 2, makePL( 50, 50, 50,  0, 50, 50), makePL(45, 0, 50)});

        tests.add(new Object[]{1, 2, makePL( 50, 50, 50, 50, 0, 50), makePL(45, 0, 50)});
        tests.add(new Object[]{2, 2, makePL( 50, 50, 50, 50, 0, 50), makePL(45, 0, 50)});

        tests.add(new Object[]{1, 2, makePL( 50, 50, 50, 50, 50,  0), makePL(0, 47, 50)});
        tests.add(new Object[]{2, 2, makePL( 50, 50, 50, 50, 50,  0), makePL(45, 47, 0)});

        return tests.toArray(new Object[][]{});
    }

    private Genotype makePL(final int ... PLs) {
        return ExactAFCalculationModelUnitTest.makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), PLs);
    }

    @Test(enabled = true, dataProvider = "TestCombineGLs")
    private void testCombineGLs(final int altIndex, final int nAlts, final Genotype testg, final Genotype expected) {
        final IndependentAllelesDiploidExactAFCalc calc = new IndependentAllelesDiploidExactAFCalc(1, 4);
        final Genotype combined = calc.combineGLs(testg, altIndex, nAlts);

        Assert.assertEquals(combined.getPL(), expected.getPL(),
                "Combined PLs " + Utils.join(",", combined.getPL()) + " != expected " + Utils.join(",", expected.getPL()));
    }
}