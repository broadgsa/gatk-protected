package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


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

    @DataProvider(name = "TestCombineGLsWithDrops")
    public Object[][] makeTestCombineGLsWithDrops() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Set<Integer> noDrops = Collections.emptySet();
        final Set<Integer> drop1 = Collections.singleton(1);
        final Set<Integer> drop2 = Collections.singleton(2);

        // AA AB BB AC BC CC
        //   drop1 (B): AA AC CC
        //   drop2 (C): AA AB BB
        tests.add(new Object[]{1, 2, makePL( 0, 1, 2, 3, 4, 5), makePL(0, 2, 5), noDrops});
        tests.add(new Object[]{2, 2, makePL( 0, 1, 2, 3, 4, 5), makePL(0, 4, 9), noDrops});
        tests.add(new Object[]{1, 2, makePL( 0, 1, 2, 3, 4, 5), makePL(0, 1, 2), drop2});
        tests.add(new Object[]{2, 2, makePL( 0, 1, 2, 3, 4, 5), makePL(0, 3, 5), drop1});

        tests.add(new Object[]{1, 2, makePL( 5, 4, 3, 2, 1, 0), makePL(0, 2, 6), noDrops});
        tests.add(new Object[]{2, 2, makePL( 5, 4, 3, 2, 1, 0), makePL(1, 0, 2), noDrops});
        tests.add(new Object[]{1, 2, makePL( 5, 4, 3, 2, 1, 0), makePL(2, 1, 0), drop2});
        tests.add(new Object[]{2, 2, makePL( 5, 4, 3, 2, 1, 0), makePL(5, 2, 0), drop1});

        tests.add(new Object[]{1, 2, makePL(10,10,10,10,10, 0), makePL( 0, 8,11), noDrops});
        tests.add(new Object[]{2, 2, makePL(10,10,10,10,10, 0), makePL( 5, 7, 0), noDrops});
        tests.add(new Object[]{1, 2, makePL(10,10,10,10,10, 0), makePL( 0, 0, 0), drop2});
        tests.add(new Object[]{2, 2, makePL(10,10,10,10,10, 0), makePL(10,10, 0), drop1});

        return tests.toArray(new Object[][]{});
    }

    private Genotype makePL(final int ... PLs) {
        return AFCalcUnitTest.makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), PLs);
    }

    @Test(enabled = true, dataProvider = "TestCombineGLs")
    private void testCombineGLs(final int altIndex, final int nAlts, final Genotype testg, final Genotype expected) {
        testCombineGLsWithDrops(altIndex, nAlts, testg, expected, Collections.<Integer>emptySet());
    }

    @Test(enabled = true, dataProvider = "TestCombineGLsWithDrops")
    private void testCombineGLsWithDrops(final int altIndex, final int nAlts, final Genotype testg, final Genotype expected, Set<Integer> allelesToDrop) {
        final IndependentAllelesDiploidExactAFCalc calc = (IndependentAllelesDiploidExactAFCalc)AFCalcFactory.createAFCalc(AFCalcFactory.Calculation.EXACT_INDEPENDENT, 1, 4);
        final Genotype combined = calc.combineGLs(testg, altIndex, allelesToDrop, nAlts);

        Assert.assertEquals(combined.getPL(), expected.getPL(),
                "Combined PLs " + Utils.join(",", combined.getPL()) + " != expected " + Utils.join(",", expected.getPL()));
    }


    static Allele A = Allele.create("A", true);
    static Allele C = Allele.create("C");
    static Allele G = Allele.create("G");

    @DataProvider(name = "TestMakeAlleleConditionalContexts")
    public Object[][] makeTestMakeAlleleConditionalContexts() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final VariantContextBuilder root = new VariantContextBuilder("x", "1", 1, 1, Arrays.asList(A));
        final VariantContextBuilder vcAC = new VariantContextBuilder(root).alleles(Arrays.asList(A, C));
        final VariantContextBuilder vcAG = new VariantContextBuilder(root).alleles(Arrays.asList(A, G));
        final VariantContextBuilder vcACG = new VariantContextBuilder(root).alleles(Arrays.asList(A, C, G));
        final VariantContextBuilder vcAGC = new VariantContextBuilder(root).alleles(Arrays.asList(A, G, C));

        final Genotype gACG = makePL( 0, 1, 2, 3, 4, 5);
        final Genotype gAGC = makePL( 0, 4, 5, 1, 3, 2);
        final Genotype gACcombined = makePL(0, 2, 5);
        final Genotype gAGcombined = makePL(0, 4, 9);
        final Genotype gACdropped  = makePL(0, 1, 2);
        final Genotype gAGdropped  = makePL(0, 3, 5);

        // biallelic
        tests.add(new Object[]{vcAC.genotypes(gACcombined).make(), Arrays.asList(vcAC.genotypes(gACcombined).make())});

        // tri-allelic
        tests.add(new Object[]{vcACG.genotypes(gACG).make(), Arrays.asList(vcAC.genotypes(gACcombined).make(), vcAG.genotypes(gAGdropped).make())});
        tests.add(new Object[]{vcAGC.genotypes(gAGC).make(), Arrays.asList(vcAG.genotypes(gAGcombined).make(), vcAC.genotypes(gACdropped).make())});

        return tests.toArray(new Object[][]{});
    }


    @Test(enabled = true, dataProvider = "TestMakeAlleleConditionalContexts")
    private void testMakeAlleleConditionalContexts(final VariantContext vc, final List<VariantContext> expectedVCs) {
        final IndependentAllelesDiploidExactAFCalc calc = (IndependentAllelesDiploidExactAFCalc)AFCalcFactory.createAFCalc(AFCalcFactory.Calculation.EXACT_INDEPENDENT, 1, 4);
        final List<VariantContext> biAllelicVCs = calc.makeAlleleConditionalContexts(vc);

        Assert.assertEquals(biAllelicVCs.size(), expectedVCs.size());

        for ( int i = 0; i < biAllelicVCs.size(); i++ ) {
            final VariantContext actual = biAllelicVCs.get(i);
            final VariantContext expected = expectedVCs.get(i);
            Assert.assertEquals(actual.getAlleles(), expected.getAlleles());

            for ( int j = 0; j < actual.getNSamples(); j++ )
                Assert.assertEquals(actual.getGenotype(j).getPL(), expected.getGenotype(j).getPL());
        }
    }

}