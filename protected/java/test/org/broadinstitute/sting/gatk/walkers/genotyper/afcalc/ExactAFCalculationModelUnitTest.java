package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class ExactAFCalculationModelUnitTest extends BaseTest {
    static Allele A = Allele.create("A", true);
    static Allele C = Allele.create("C");
    static Allele G = Allele.create("G");

    static int sampleNameCounter = 0;
    static Genotype AA1, AB1, BB1, NON_INFORMATIVE1;
    static Genotype AA2, AB2, AC2, BB2, BC2, CC2, NON_INFORMATIVE2;
    final double[] FLAT_3SAMPLE_PRIORS = MathUtils.normalizeFromLog10(new double[2*3+1], true);  // flat priors
    final private static boolean INCLUDE_BIALLELIC = true;
    final private static boolean INCLUDE_TRIALLELIC = true;
    final private static boolean Guillermo_FIXME = false; // TODO -- can only be enabled when GdA fixes bug

    @BeforeSuite
    public void before() {
        AA1 = makePL(Arrays.asList(A, A), 0, 20, 20);
        AB1 = makePL(Arrays.asList(A, C), 20, 0, 20);
        BB1 = makePL(Arrays.asList(C, C), 20, 20, 0);
        NON_INFORMATIVE1 = makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0);

        AA2 = makePL(Arrays.asList(A, A), 0, 20, 20, 20, 20, 20);
        AB2 = makePL(Arrays.asList(A, C), 20, 0, 20, 20, 20, 20);
        BB2 = makePL(Arrays.asList(C, C), 20, 20, 0, 20, 20, 20);
        AC2 = makePL(Arrays.asList(A, G), 20, 20, 20, 0, 20, 20);
        BC2 = makePL(Arrays.asList(C, G), 20, 20, 20, 20, 0, 20);
        CC2 = makePL(Arrays.asList(G, G), 20, 20, 20, 20, 20, 0);
        NON_INFORMATIVE2 = makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), 0, 0, 0, 0, 0, 0);
    }

    protected static Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }

    private class GetGLsTest extends TestDataProvider {
        GenotypesContext GLs;
        int numAltAlleles;
        final ExactAFCalc calc;
        final int[] expectedACs;
        final double[] priors;
        final String priorName;

        private GetGLsTest(final ExactAFCalc calc, int numAltAlleles, List<Genotype> arg, final double[] priors, final String priorName) {
            super(GetGLsTest.class);
            GLs = GenotypesContext.create(new ArrayList<Genotype>(arg));
            this.numAltAlleles = numAltAlleles;
            this.calc = calc;
            this.priors = priors;
            this.priorName = priorName;

            expectedACs = new int[numAltAlleles+1];
            for ( int alleleI = 0; alleleI < expectedACs.length; alleleI++ ) {
                expectedACs[alleleI] = 0;
                final Allele allele = getAlleles().get(alleleI);
                for ( Genotype g : arg ) {
                    expectedACs[alleleI] += Collections.frequency(g.getAlleles(), allele);
                }
            }
        }

        public AFCalcResult execute() {
            return getCalc().getLog10PNonRef(getVC(), getPriors());
        }

        public AFCalcResult executeRef() {
            final ExactAFCalc ref = new ReferenceDiploidExactAFCalc(getCalc().nSamples, getCalc().getMaxAltAlleles());
            return ref.getLog10PNonRef(getVC(), getPriors());
        }

        public double[] getPriors() {
            return priors;
        }

        public ExactAFCalc getCalc() {
            return calc;
        }

        public VariantContext getVC() {
            VariantContextBuilder builder = new VariantContextBuilder("test", "1", 1, 1, getAlleles());
            builder.genotypes(GLs);
            return builder.make();
        }

        public List<Allele> getAlleles() {
            return Arrays.asList(Allele.create("A", true),
                    Allele.create("C"),
                    Allele.create("G"),
                    Allele.create("T")).subList(0, numAltAlleles+1);
        }

        public int getExpectedAltAC(final int alleleI) {
            return expectedACs[alleleI+1];
        }

        public String toString() {
            return String.format("%s model=%s prior=%s input=%s", super.toString(), calc.getClass().getSimpleName(),
                    priorName, GLs.size() > 5 ? String.format("%d samples", GLs.size()) : GLs);
        }
    }

    @DataProvider(name = "wellFormedGLs")
    public Object[][] createSimpleGLsData() {
        final List<Genotype> biAllelicSamples = Arrays.asList(AA1, AB1, BB1);
        final List<Genotype> triAllelicSamples = Arrays.asList(AA2, AB2, BB2, AC2, BC2, CC2);

        for ( final int nSamples : Arrays.asList(1, 2, 3, 4) ) {
            final ExactAFCalc diploidCalc = new ReferenceDiploidExactAFCalc(nSamples, 4);
//            final ExactAFCalc optDiploidCalc = new ConstrainedDiploidExactAFCalc(nSamples, 4);
            final ExactAFCalc generalCalc = new GeneralPloidyExactAFCalc(nSamples, 4, 2);
            final ExactAFCalc indCalc = new IndependentAllelesDiploidExactAFCalc(nSamples, 4);

            final int nPriorValues = 2*nSamples+1;
            final double[] flatPriors = MathUtils.normalizeFromLog10(new double[nPriorValues], true);  // flat priors
            final double[] humanPriors = new double[nPriorValues];
            UnifiedGenotyperEngine.computeAlleleFrequencyPriors(nPriorValues - 1, humanPriors, 0.001);

            for ( final double[] priors : Arrays.asList(flatPriors, humanPriors) ) { // , humanPriors) ) {
                for ( ExactAFCalc model : Arrays.asList(diploidCalc, generalCalc, indCalc) ) {
                    final String priorName = priors == humanPriors ? "human" : "flat";

                    // bi-allelic
                    if ( INCLUDE_BIALLELIC && nSamples <= biAllelicSamples.size() )
                        for ( List<Genotype> genotypes : Utils.makePermutations(biAllelicSamples, nSamples, true) )
                            new GetGLsTest(model, 1, genotypes, priors, priorName);

                    // tri-allelic
                    if ( INCLUDE_TRIALLELIC && ( ! priorName.equals("human") || model != generalCalc || Guillermo_FIXME ) )
                        for ( List<Genotype> genotypes : Utils.makePermutations(triAllelicSamples, nSamples, true) )
                            new GetGLsTest(model, 2, genotypes, priors, priorName);
                }
            }
        }

        return GetGLsTest.getTests(GetGLsTest.class);
    }

    private static class NonInformativeData {
        final Genotype nonInformative;
        final List<Genotype> called;
        final int nAltAlleles;

        private NonInformativeData(List<Genotype> called, Genotype nonInformative, int nAltAlleles) {
            this.called = called;
            this.nonInformative = nonInformative;
            this.nAltAlleles = nAltAlleles;
        }
    }

    @DataProvider(name = "GLsWithNonInformative")
    public Object[][] makeGLsWithNonInformative() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<NonInformativeData> nonInformativeTests = new LinkedList<NonInformativeData>();
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB1), NON_INFORMATIVE1, 1));
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB2), NON_INFORMATIVE2, 2));
        nonInformativeTests.add(new NonInformativeData(Arrays.asList(AB2, BC2), NON_INFORMATIVE2, 2));

        for ( final int nNonInformative : Arrays.asList(1, 10, 100) ) {
            for ( final NonInformativeData testData : nonInformativeTests ) {
                final List<Genotype> samples = new ArrayList<Genotype>();
                samples.addAll(testData.called);
                samples.addAll(Collections.nCopies(nNonInformative, testData.nonInformative));

                final int nSamples = samples.size();
                final ExactAFCalc diploidCalc = new ReferenceDiploidExactAFCalc(nSamples, 4);
//                final ExactAFCalc optDiploidCalc = new ConstrainedDiploidExactAFCalc(nSamples, 4);
                final ExactAFCalc generalCalc = new GeneralPloidyExactAFCalc(nSamples, 4, 2);
                final ExactAFCalc indCalc = new IndependentAllelesDiploidExactAFCalc(nSamples, 4);

                final double[] priors = MathUtils.normalizeFromLog10(new double[2*nSamples+1], true);  // flat priors

                for ( ExactAFCalc model : Arrays.asList(diploidCalc, generalCalc, indCalc) ) {
                    final GetGLsTest onlyInformative = new GetGLsTest(model, testData.nAltAlleles, testData.called, priors, "flat");

                    for ( int rotation = 0; rotation < nSamples; rotation++ ) {
                        Collections.rotate(samples, 1);
                        final GetGLsTest withNonInformative = new GetGLsTest(model, testData.nAltAlleles, samples, priors, "flat");
                        tests.add(new Object[]{onlyInformative, withNonInformative});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "wellFormedGLs")
    public void testGLs(GetGLsTest cfg) {
        testResultSimple(cfg);
    }

    @Test(enabled = true, dataProvider = "GLsWithNonInformative", dependsOnMethods = "testGLs")
    public void testGLsWithNonInformative(GetGLsTest onlyInformative, GetGLsTest withNonInformative) {
        final AFCalcResult expected = onlyInformative.execute();
        final AFCalcResult actual = withNonInformative.execute();

        testResultSimple(withNonInformative);
        compareAFCalcResults(actual, expected, onlyInformative.getCalc());
    }

    private void testResultSimple(final GetGLsTest cfg) {
        final AFCalcResult refResultTracker = cfg.executeRef();
        final AFCalcResult resultTracker = cfg.execute();

        compareAFCalcResults(resultTracker, refResultTracker, cfg.getCalc());

//        final int minNumberOfEvaluations = cfg.getVC().getCalledChrCount();
//        Assert.assertTrue(result.getnEvaluations() >= minNumberOfEvaluations,
//                "Number of evaluations " + result.getnEvaluations() + " must be at least " + minNumberOfEvaluations);
        Assert.assertNotNull(resultTracker.getAllelesUsedInGenotyping());
        Assert.assertTrue(cfg.getAlleles().containsAll(resultTracker.getAllelesUsedInGenotyping()), "Result object has alleles not in our initial allele list");

        for ( int altAlleleI = 0; altAlleleI < cfg.numAltAlleles; altAlleleI++ ) {
            int expectedAlleleCount = cfg.getExpectedAltAC(altAlleleI);
            int calcAC_MLE = resultTracker.getAlleleCountsOfMLE()[altAlleleI];

            final Allele allele = cfg.getAlleles().get(altAlleleI+1);
            Assert.assertEquals(calcAC_MLE, expectedAlleleCount, "MLE AC not equal to expected AC for allele " + allele);
        }

        // TODO
        // TODO -- enable when we understand the contract between AC_MAP and pNonRef
        // TODO
//        final int AC_MAP = (int)MathUtils.sum(result.getAlleleCountsOfMAP());
//        if ( AC_MAP > 0  ) {
//            Assert.assertTrue(result.getNormalizedPosteriorOfAFzero() < 0.50, "MAP AC " + AC_MAP + " > 0 but we had posterior AF = 0 > 0.5 of " + result.getNormalizedPosteriorOfAFzero());
//        } else {
//            Assert.assertTrue(result.getNormalizedPosteriorOfAFzero() > 0.50, "MAP AC " + AC_MAP + " == 0 but we had posterior AF = 0 < 0.5 of " + result.getNormalizedPosteriorOfAFzero());
//        }
    }

    private void compareAFCalcResults(final AFCalcResult actual, final AFCalcResult expected, final ExactAFCalc calc) {
        final double TOLERANCE = 1;  // TODO -- tighten up tolerances

        Assert.assertEquals(actual.getLog10PriorOfAFEq0(), expected.getLog10PriorOfAFEq0(), TOLERANCE);
        Assert.assertEquals(actual.getLog10PriorOfAFGT0(), expected.getLog10PriorOfAFGT0(), TOLERANCE);
        Assert.assertEquals(actual.getLog10LikelihoodOfAFEq0(), expected.getLog10LikelihoodOfAFEq0(), TOLERANCE);
        Assert.assertEquals(actual.getLog10LikelihoodOfAFGT0(), expected.getLog10LikelihoodOfAFGT0(), TOLERANCE);
        Assert.assertEquals(actual.getLog10PosteriorOfAFEq0(), expected.getLog10PosteriorOfAFEq0(), TOLERANCE);
        Assert.assertEquals(actual.getLog10PosteriorOfAFGT0(), expected.getLog10PosteriorOfAFGT0(), TOLERANCE);
        Assert.assertEquals(actual.getAlleleCountsOfMLE(), expected.getAlleleCountsOfMLE());
        Assert.assertEquals(actual.getAllelesUsedInGenotyping(), expected.getAllelesUsedInGenotyping());

        for ( final Allele a : expected.getAllelesUsedInGenotyping() ) {
            if ( ! a.isReference() ) {
                Assert.assertEquals(actual.getAlleleCountAtMLE(a), expected.getAlleleCountAtMLE(a));
                if ( ! ( calc instanceof GeneralPloidyExactAFCalc ) )
                    // TODO -- delete when general ploidy works properly with multi-allelics
                    Assert.assertEquals(actual.isPolymorphic(a, 0.0), expected.isPolymorphic(a, 0.0));
            }
        }
    }

    @Test(enabled = true, dataProvider = "Models")
    public void testLargeGLs(final ExactAFCalc calc) {
        final Genotype BB = makePL(Arrays.asList(C, C), 20000000, 20000000, 0);
        GetGLsTest cfg = new GetGLsTest(calc, 1, Arrays.asList(BB, BB, BB), FLAT_3SAMPLE_PRIORS, "flat");

        final AFCalcResult resultTracker = cfg.execute();

        int calculatedAlleleCount = resultTracker.getAlleleCountsOfMLE()[0];
        Assert.assertEquals(calculatedAlleleCount, 6);
    }

    @Test(enabled = true, dataProvider = "Models")
    public void testMismatchedGLs(final ExactAFCalc calc) {
        final Genotype AB = makePL(Arrays.asList(A, C), 2000, 0, 2000, 2000, 2000, 2000);
        final Genotype AC = makePL(Arrays.asList(A, G), 100, 100, 100, 0, 100, 100);
        GetGLsTest cfg = new GetGLsTest(calc, 2, Arrays.asList(AB, AC), FLAT_3SAMPLE_PRIORS, "flat");

        final AFCalcResult resultTracker = cfg.execute();

        Assert.assertEquals(resultTracker.getAlleleCountsOfMLE()[0], 1);
        Assert.assertEquals(resultTracker.getAlleleCountsOfMLE()[1], 1);
    }

    // --------------------------------------------------------------------------------
    //
    // Code to test that the pNonRef value is meaningful
    //
    // --------------------------------------------------------------------------------

    private static class PNonRefData {
        final Genotype g;
        final double pNonRef, tolerance;
        final boolean canScale;
        final List<ExactAFCalculationTestBuilder.ModelType> badModels;
        final VariantContext vc;

        private PNonRefData(final VariantContext vc, Genotype g, double pNonRef, double tolerance, final boolean canScale) {
            this(vc, g, pNonRef, tolerance, canScale, Collections.<ExactAFCalculationTestBuilder.ModelType>emptyList());
        }

        private PNonRefData(final VariantContext vc, Genotype g, double pNonRef, double tolerance, final boolean canScale, final List<ExactAFCalculationTestBuilder.ModelType> badModels) {
            this.g = g;
            this.pNonRef = pNonRef;
            this.tolerance = tolerance;
            this.canScale = canScale;
            this.badModels = badModels;
            this.vc = vc;
        }

        public PNonRefData scale(final int scaleFactor) {
            if ( canScale ) {
                final int[] PLs = new int[g.getPL().length];
                for ( int i = 0; i < PLs.length; i++ ) PLs[i] = g.getPL()[i] * ((int)Math.log10(scaleFactor)+1);
                final Genotype scaledG = new GenotypeBuilder(g).PL(PLs).make();
                final double scaledPNonRef = pNonRef < 0.5 ? pNonRef / scaleFactor : 1 - ((1-pNonRef) / scaleFactor);
                return new PNonRefData(vc, scaledG, scaledPNonRef, tolerance, true);
            } else {
                return this;
            }
        }
    }

    @DataProvider(name = "PNonRef")
    public Object[][] makePNonRefTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);
        final List<Allele> CC = Arrays.asList(C, C);
        final List<Allele> AG = Arrays.asList(A, G);
        final List<Allele> GG = Arrays.asList(G, G);
        final List<Allele> CG = Arrays.asList(C, G);

        final VariantContext vc2 = new VariantContextBuilder("x","1", 1, 1, Arrays.asList(A, C)).make();
        final VariantContext vc3 = new VariantContextBuilder("x","1", 1, 1, Arrays.asList(A, C, G)).make();
        final ExactAFCalculationTestBuilder.PriorType priorType = ExactAFCalculationTestBuilder.PriorType.flat;

        final List<ExactAFCalculationTestBuilder.ModelType> constrainedModel = Arrays.asList(ExactAFCalculationTestBuilder.ModelType.ConstrainedDiploidExact);

        final double TOLERANCE = 0.5;

        final List<PNonRefData> initialPNonRefData = Arrays.asList(
                // bi-allelic sites
                new PNonRefData(vc2, makePL(AA, 0, 10, 10), 0.1666667, TOLERANCE, true),
                new PNonRefData(vc2, makePL(AA, 0,  1, 10), 0.4721084, TOLERANCE, false, constrainedModel),
                new PNonRefData(vc2, makePL(AA, 0,  1,  1), 0.6136992, TOLERANCE, false, constrainedModel),
                new PNonRefData(vc2, makePL(AA, 0,  5,  5), 0.3874259, TOLERANCE, false, constrainedModel),
                new PNonRefData(vc2, makePL(AC, 10, 0, 10), 0.9166667, TOLERANCE, true),
                new PNonRefData(vc2, makePL(CC, 10, 10, 0), 0.9166667, TOLERANCE, true),

                // tri-allelic sites -- cannot scale because of the naivety of our scaling estimator
                new PNonRefData(vc3, makePL(AA, 0, 10, 10, 10, 10, 10), 0.3023255813953489, TOLERANCE * 2, false), // more tolerance because constrained model is a bit inaccurate
                new PNonRefData(vc3, makePL(AC, 10, 0, 10, 10, 10, 10), 0.9166667, TOLERANCE, false),
                new PNonRefData(vc3, makePL(CC, 10, 10, 0, 10, 10, 10), 0.9166667, TOLERANCE, false),
                new PNonRefData(vc3, makePL(AG, 10, 10, 10, 0, 10, 10), 0.9166667, TOLERANCE, false),
                new PNonRefData(vc3, makePL(CG, 10, 10, 10, 10, 0, 10), 0.80, TOLERANCE, false),
                new PNonRefData(vc3, makePL(GG, 10, 10, 10, 10, 10, 0), 0.9166667, TOLERANCE, false)
        );

        for ( ExactAFCalculationTestBuilder.ModelType modelType : ExactAFCalculationTestBuilder.ModelType.values() ) {
            for ( int nNonInformative = 0; nNonInformative < 3; nNonInformative++ ) {
                for ( final PNonRefData rootData : initialPNonRefData ) {
                    for ( int plScale = 1; plScale <= 100000; plScale *= 10 ) {
                        if ( ! rootData.badModels.contains(modelType) && (plScale == 1 || rootData.canScale) ) {
                            final PNonRefData data = rootData.scale(plScale);
                            tests.add(new Object[]{data.vc, modelType, priorType, Arrays.asList(data.g), data.pNonRef, data.tolerance, nNonInformative});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "PNonRef")
    private void testPNonRef(final VariantContext vcRoot,
                             ExactAFCalculationTestBuilder.ModelType modelType,
                             ExactAFCalculationTestBuilder.PriorType priorType,
                             final List<Genotype> genotypes,
                             final double expectedPNonRef,
                             final double tolerance,
                             final int nNonInformative) {
        final ExactAFCalculationTestBuilder testBuilder
                = new ExactAFCalculationTestBuilder(1, vcRoot.getNAlleles()-1, modelType, priorType);

        final VariantContextBuilder vcb = new VariantContextBuilder(vcRoot);
        vcb.genotypes(genotypes);

        final AFCalcResult resultTracker = testBuilder.makeModel().getLog10PNonRef(vcb.make(), testBuilder.makePriors());

        Assert.assertEquals(resultTracker.getLog10PosteriorOfAFGT0(), Math.log10(expectedPNonRef), tolerance,
                "Actual pNonRef not within tolerance " + tolerance + " of expected");
    }

    // --------------------------------------------------------------------------------
    //
    // Test priors
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "Models")
    public Object[][] makeModels() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new ReferenceDiploidExactAFCalc(2, 4)});
//        tests.add(new Object[]{new ConstrainedDiploidExactAFCalc(2, 4)});
//        tests.add(new Object[]{new GeneralPloidyExactAFCalc(2, 4, 2)});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "Models")
    public void testBiallelicPriors(final ExactAFCalc model) {
        final int REF_PL = 10;
        final Genotype AB = makePL(Arrays.asList(A,C), REF_PL, 0, 10000);

        for ( int log10NonRefPrior = 1; log10NonRefPrior < 10*REF_PL; log10NonRefPrior += 1 ) {
            final double refPrior = 1 - QualityUtils.qualToErrorProb(log10NonRefPrior);
            final double[] priors = MathUtils.toLog10(new double[]{refPrior, (1-refPrior) / 2, (1-refPrior) / 2});
            GetGLsTest cfg = new GetGLsTest(model, 1, Arrays.asList(AB), priors, "pNonRef" + log10NonRefPrior);
            final AFCalcResult resultTracker = cfg.execute();
            final int actualAC = resultTracker.getAlleleCountsOfMLE()[0];

            final double pRefWithPrior = AB.getLikelihoods().getAsVector()[0] + priors[0];
            final double pHetWithPrior = AB.getLikelihoods().getAsVector()[1] + priors[1];
            final double nonRefPost = Math.pow(10, pHetWithPrior) / (Math.pow(10, pRefWithPrior) + Math.pow(10, pHetWithPrior));

            if ( nonRefPost < 0.1 )
                Assert.assertTrue(resultTracker.isPolymorphic(C, -1));

            final int expectedMLEAC = 1; // the MLE is independent of the prior
            Assert.assertEquals(actualAC, expectedMLEAC,
                    "actual AC with priors " + log10NonRefPrior + " not expected "
                            + expectedMLEAC + " priors " + Utils.join(",", priors));
        }
    }

    @Test(enabled = false, dataProvider = "Models")
    public void testTriallelicPriors(final ExactAFCalc model) {
        // TODO
        // TODO
        // TODO THIS SEEMS TO ID A BUG IN THE EXACT MODEL FOR MULTI-ALLELICS, AS THE
        // TODO SECOND ALLELE ISN'T HAVING A SQUARED PRIOR.  TALK TO ERIC AND CONFIRM
        // TODO
        // TODO
        final int REF_PL_AB = 10, REF_PL_AC = 20; // first AC goes, then AB
        final Genotype AB = makePL(Arrays.asList(A,C), REF_PL_AB, 0, 10000, 10000, 10000);
        final Genotype AC = makePL(Arrays.asList(A, G), REF_PL_AC, 10000, 10000, 0, 10000, 10000);

        for ( int log10NonRefPrior = 1; log10NonRefPrior < 100*REF_PL_AC; log10NonRefPrior += 1 ) {
            final double refPrior = 1 - QualityUtils.qualToErrorProb(log10NonRefPrior);
            final double nonRefPrior = (1-refPrior) / 2;
            final double[] priors = MathUtils.toLog10(new double[]{refPrior, nonRefPrior, nonRefPrior, nonRefPrior, nonRefPrior, nonRefPrior});
            GetGLsTest cfg = new GetGLsTest(model, 2, Arrays.asList(AB, AC), priors, "pNonRef" + log10NonRefPrior);
            final AFCalcResult resultTracker = cfg.execute();
            final int actualAC_AB = resultTracker.getAlleleCountsOfMLE()[0];

            final double pRefABWithPrior = AB.getLikelihoods().getAsVector()[0] + priors[0];
            final double pHetABWithPrior = AB.getLikelihoods().getAsVector()[1] + priors[1];
            final int expectedAC_AB = pRefABWithPrior <= pHetABWithPrior ? 1 : 0;
            Assert.assertEquals(actualAC_AB, expectedAC_AB,
                    "actual AC with priors " + log10NonRefPrior + " not expected "
                            + expectedAC_AB + " priors " + Utils.join(",", priors));

            final double nonRefPriorSecondAllele = Math.pow(nonRefPrior, 2);
            final double refPriorSecondAllele = 1 - nonRefPriorSecondAllele;
            final int actualAC_AC = resultTracker.getAlleleCountsOfMLE()[1];
            final double pRefACWithPrior = AB.getLikelihoods().getAsVector()[0] + Math.log10(refPriorSecondAllele);
            final double pHetACWithPrior = AC.getLikelihoods().getAsVector()[3] + Math.log10(nonRefPriorSecondAllele);
            final int expectedAC_AC = pRefACWithPrior <= pHetACWithPrior ? 1 : 0;
            Assert.assertEquals(actualAC_AC, expectedAC_AC,
                    "actual AC with priors " + log10NonRefPrior + " not expected "
                            + expectedAC_AC + " priors " + Utils.join(",", priors));
        }
    }

    @DataProvider(name = "MaxACsToVisit")
    public Object[][] makeMaxACsToVisit() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int nSamples = 10;
        final ExactAFCalculationTestBuilder.ModelType modelType = ExactAFCalculationTestBuilder.ModelType.ConstrainedDiploidExact;

        for (int nNonInformative = 0; nNonInformative < nSamples - 1; nNonInformative++ ) {
            final int nChrom = (nSamples - nNonInformative) * 2;
            for ( int i = 0; i < nChrom; i++ ) {
                // bi-allelic
                tests.add(new Object[]{nSamples, Arrays.asList(i), nNonInformative, modelType});

                // tri-allelic
                for ( int j = 0; j < (nChrom - i); j++)
                    tests.add(new Object[]{nSamples, Arrays.asList(i, j), nNonInformative, modelType});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "MaxACsToVisit")
    public void testMaxACsToVisit(final int nSamples, final List<Integer> requestedACs, final int nNonInformative, final ExactAFCalculationTestBuilder.ModelType modelType) {
        final int nAlts = requestedACs.size();
        final ExactAFCalculationTestBuilder testBuilder
                = new ExactAFCalculationTestBuilder(nSamples, nAlts, modelType,
                ExactAFCalculationTestBuilder.PriorType.human);

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

        final ExactAFCalculationTestBuilder.ModelType modelType = ExactAFCalculationTestBuilder.ModelType.ConstrainedDiploidExact;
        final ExactAFCalculationTestBuilder testBuilder
                = new ExactAFCalculationTestBuilder(1, vc.getNAlleles()-1, modelType,
                ExactAFCalculationTestBuilder.PriorType.human);
        final int[] maxACsToVisit = ((ConstrainedDiploidExactAFCalc)testBuilder.makeModel()).computeMaxACs(vc);
        testExpectedACs(vc, maxACsToVisit);
    }
}