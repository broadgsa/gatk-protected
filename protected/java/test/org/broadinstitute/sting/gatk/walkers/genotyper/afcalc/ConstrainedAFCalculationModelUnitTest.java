package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class ConstrainedAFCalculationModelUnitTest extends BaseTest {
    static Allele A = Allele.create("A", true);
    static Allele C = Allele.create("C");
    static Allele G = Allele.create("G");

    protected static Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        return AFCalcUnitTest.makePL(expectedGT, pls);
    }

    @DataProvider(name = "MaxACsToVisit")
    public Object[][] makeMaxACsToVisit() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int nSamples = 10;

        for (int nNonInformative = 0; nNonInformative < nSamples - 1; nNonInformative++ ) {
            final int nChrom = (nSamples - nNonInformative) * 2;
            for ( int i = 0; i < nChrom; i++ ) {
                // bi-allelic
                tests.add(new Object[]{nSamples, Arrays.asList(i), nNonInformative, AFCalcFactory.Calculation.EXACT_CONSTRAINED});

                // tri-allelic
                for ( int j = 0; j < (nChrom - i); j++)
                    tests.add(new Object[]{nSamples, Arrays.asList(i, j), nNonInformative, AFCalcFactory.Calculation.EXACT_CONSTRAINED});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "MaxACsToVisit")
    public void testMaxACsToVisit(final int nSamples, final List<Integer> requestedACs, final int nNonInformative, final AFCalcFactory.Calculation modelType) {
        final int nAlts = requestedACs.size();
        final AFCalcTestBuilder testBuilder
                = new AFCalcTestBuilder(nSamples, nAlts, modelType,
                AFCalcTestBuilder.PriorType.human);

        final VariantContext vc = testBuilder.makeACTest(requestedACs, nNonInformative, 100);
        final int[] maxACsToVisit = ((ConstrainedDiploidExactAFCalc)testBuilder.makeModel()).computeMaxACs(vc);

        testExpectedACs(vc, maxACsToVisit);
    }

    private void testExpectedACs(final VariantContext vc, final int[] maxACsToVisit) {
        // this is necessary because cannot ensure that the tester gives us back the
        // requested ACs due to rounding errors
        final List<Integer> ACs = new ArrayList<Integer>();
        for ( final Allele a : vc.getAlternateAlleles() )
            ACs.add(vc.getCalledChrCount(a));

        for ( int i = 0; i < maxACsToVisit.length; i++ ) {
            Assert.assertEquals(maxACsToVisit[i], (int)ACs.get(i), "Maximum AC computed wasn't equal to the max possible in the construction for alt allele " + i);
        }
    }

    @DataProvider(name = "MaxACsGenotypes")
    public Object[][] makeMaxACsForGenotype() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);
        final List<Allele> CC = Arrays.asList(C, C);
        final List<Allele> AG = Arrays.asList(A, G);
        final List<Allele> GG = Arrays.asList(G, G);
        final List<Allele> CG = Arrays.asList(C, G);

        final VariantContext vc2 = new VariantContextBuilder("x","1", 1, 1, Arrays.asList(A, C)).make();
        final VariantContext vc3 = new VariantContextBuilder("x","1", 1, 1, Arrays.asList(A, C, G)).make();

        tests.add(new Object[]{vc2, makePL(AA, 0, 10, 10)});
        tests.add(new Object[]{vc2, makePL(AC, 10, 0, 10)});
        tests.add(new Object[]{vc2, makePL(CC, 10, 10, 0)});

        // make sure non-informative => 0
        tests.add(new Object[]{vc2, makePL(AA, 0, 0, 0)});
        tests.add(new Object[]{vc3, makePL(AA, 0, 0, 0, 0, 0, 0)});

        // multi-allelics
        tests.add(new Object[]{vc3, makePL(AG, 10, 10, 10, 0, 10, 10)});
        tests.add(new Object[]{vc3, makePL(CG, 10, 10, 10, 10, 0, 10)});
        tests.add(new Object[]{vc3, makePL(GG, 10, 10, 10, 10, 10, 0)});

        // deal with non-informatives third alleles
        tests.add(new Object[]{vc3, makePL(AC, 10, 0, 10,  0, 0, 10)});
        tests.add(new Object[]{vc3, makePL(AC, 10, 0, 10, 10, 0, 10)});
        tests.add(new Object[]{vc3, makePL(AC, 10, 0, 10, 10, 0,  0)});
        tests.add(new Object[]{vc3, makePL(AC, 10, 0, 10,  0, 0,  0)});
        tests.add(new Object[]{vc3, makePL(CC, 10, 10, 0,  0, 0, 10)});
        tests.add(new Object[]{vc3, makePL(CC, 10, 10, 0, 10, 0, 10)});
        tests.add(new Object[]{vc3, makePL(CC, 10, 10, 0, 10, 0,  0)});
        tests.add(new Object[]{vc3, makePL(CC, 10, 10, 0,  0, 0,  0)});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "MaxACsGenotypes")
    private void testMakeACByGenotype(final VariantContext vcRoot, final Genotype g) {
        final VariantContext vc = new VariantContextBuilder(vcRoot).genotypes(g).make();

        final AFCalcTestBuilder testBuilder
                = new AFCalcTestBuilder(1, vc.getNAlleles()-1, AFCalcFactory.Calculation.EXACT_CONSTRAINED,
                AFCalcTestBuilder.PriorType.human);

        final int[] maxACsToVisit = ((ConstrainedDiploidExactAFCalc)testBuilder.makeModel()).computeMaxACs(vc);

        testExpectedACs(vc, maxACsToVisit);
    }
}