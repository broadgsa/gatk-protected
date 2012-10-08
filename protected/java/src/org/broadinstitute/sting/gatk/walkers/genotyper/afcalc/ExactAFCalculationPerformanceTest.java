package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.TTCCLayout;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 10/2/12
 * Time: 10:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class ExactAFCalculationPerformanceTest {
    final static Logger logger = Logger.getLogger(ExactAFCalculationPerformanceTest.class);

    private static abstract class Analysis {
        final GATKReport report;

        public Analysis(final String name, final List<String> columns) {
            report = GATKReport.newSimpleReport(name, columns);
        }

        public abstract void run(final ExactAFCalculationTestBuilder testBuilder,
                                 final List<Object> coreColumns);

        public String getName() {
            return getTable().getTableName();
        }

        public GATKReportTable getTable() {
            return report.getTables().iterator().next();
        }
    }

    private static class AnalyzeByACAndPL extends Analysis {
        public AnalyzeByACAndPL(final List<String> columns) {
            super("AnalyzeByACAndPL", Utils.append(columns, "non.type.pls", "ac", "n.alt.seg", "other.ac"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final ExactAFCalc calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                for ( int[] ACs : makeACs(testBuilder.numAltAlleles, testBuilder.nSamples*2) ) {
                    final VariantContext vc = testBuilder.makeACTest(ACs, 0, nonTypePL);

                    timer.start();
                    final AFCalcResultTracker resultTracker = calc.getLog10PNonRef(vc, priors);
                    final long runtime = timer.getElapsedTimeNano();

                    int otherAC = 0;
                    int nAltSeg = 0;
                    for ( int i = 0; i < ACs.length; i++ ) {
                        nAltSeg += ACs[i] > 0 ? 1 : 0;
                        if ( i > 0 ) otherAC += ACs[i];
                    }

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, resultTracker.getnEvaluations(), nonTypePL, ACs[0], nAltSeg, otherAC));
                    report.addRowList(columns);
                }
            }
        }

        private List<int[]> makeACs(final int nAltAlleles, final int nChrom) {
            if ( nAltAlleles > 2 ) throw new IllegalArgumentException("nAltAlleles must be < 3");

            final List<int[]> ACs = new LinkedList<int[]>();

            final List<Integer> ACsToTry = MathUtils.log10LinearRange(0, nChrom, 0.1); //Arrays.asList(0, 1, 2, 3, 6, 10, 20, 40, 60, 100, 200, 400, 600, 1000, 2000, 4000, 6000, 10000, 100000);

            for ( int i : ACsToTry ) {
                if ( i < nChrom ) {
                    if ( nAltAlleles == 1 ) {
                        ACs.add(new int[]{i});
                    } else if ( nAltAlleles == 2 ) {
                        for ( int j : ACsToTry ) {
                            if ( j < nChrom - i )
                                ACs.add(new int[]{i, j});
                        }
                    } else {
                        throw new IllegalStateException("cannot get here");
                    }
                }
            }

            return ACs;
        }
    }

    private static class AnalyzeBySingletonPosition extends Analysis {
        public AnalyzeBySingletonPosition(final List<String> columns) {
            super("AnalyzeBySingletonPosition", Utils.append(columns, "non.type.pls", "position.of.singleton"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final ExactAFCalc calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                final int[] ac = new int[testBuilder.numAltAlleles];
                ac[0] = 1;
                final VariantContext vc = testBuilder.makeACTest(ac, 0, nonTypePL);

                for ( final int position : MathUtils.log10LinearRange(0, vc.getNSamples(), 0.1) ) {
                    final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                    final List<Genotype> genotypes = new ArrayList<Genotype>(vc.getGenotypes());
                    Collections.rotate(genotypes, position);
                    vcb.genotypes(genotypes);

                    timer.start();
                    final AFCalcResultTracker resultTracker = calc.getLog10PNonRef(vcb.make(), priors);
                    final long runtime = timer.getElapsedTimeNano();

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, resultTracker.getnEvaluations(), nonTypePL, position));
                    report.addRowList(columns);
                }
            }
        }
    }

    private static class AnalyzeByNonInformative extends Analysis {
        public AnalyzeByNonInformative(final List<String> columns) {
            super("AnalyzeByNonInformative", Utils.append(columns, "non.type.pls", "n.non.informative"));
        }

        public void run(final ExactAFCalculationTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final ExactAFCalc calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                final int[] ac = new int[testBuilder.numAltAlleles];
                ac[0] = 1;

                for ( int nNonInformative = 0; nNonInformative < testBuilder.nSamples; nNonInformative++ ) {
                    final VariantContext vc = testBuilder.makeACTest(ac, nNonInformative, nonTypePL);

                    timer.start();
                    final AFCalcResultTracker resultTracker = calc.getLog10PNonRef(vc, priors);
                    final long runtime = timer.getElapsedTimeNano();

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, resultTracker.getnEvaluations(), nonTypePL, nNonInformative));
                    report.addRowList(columns);
                }
            }
        }
    }

    private static class ModelParams {
        final ExactAFCalculationTestBuilder.ModelType modelType;
        final int maxBiNSamples, maxTriNSamples;

        private ModelParams(ExactAFCalculationTestBuilder.ModelType modelType, int maxBiNSamples, int maxTriNSamples) {
            this.modelType = modelType;
            this.maxBiNSamples = maxBiNSamples;
            this.maxTriNSamples = maxTriNSamples;
        }

        public boolean meetsConstraints(final int nAltAlleles, final int nSamples) {
            if ( nAltAlleles == 1 )
                return nSamples <= maxBiNSamples;
            else if ( nAltAlleles == 2 )
                return nSamples <= maxTriNSamples;
            else
                throw new IllegalStateException("Unexpected number of alt alleles " + nAltAlleles);
        }
    }

    public enum Operation {
        ANALYZE,
        SINGLE
    }
    public static void main(final String[] args) throws Exception {
        final TTCCLayout layout = new TTCCLayout();
        layout.setThreadPrinting(false);
        layout.setCategoryPrefixing(false);
        layout.setContextPrinting(false);
        logger.addAppender(new ConsoleAppender(layout));

        final Operation op = Operation.valueOf(args[0]);

        switch ( op ) {
            case ANALYZE: analyze(args); break;
            case SINGLE: profileBig(args); break;
            default: throw new IllegalAccessException("unknown operation " + op);
        }
    }

    private static void profileBig(final String[] args) throws Exception {
        final int nSamples = Integer.valueOf(args[1]);
        final int ac = Integer.valueOf(args[2]);

        final ExactAFCalculationTestBuilder testBuilder = new ExactAFCalculationTestBuilder(nSamples, 1,
                ExactAFCalculationTestBuilder.ModelType.IndependentDiploidExact,
                ExactAFCalculationTestBuilder.PriorType.human);

        final VariantContext vc = testBuilder.makeACTest(new int[]{ac}, 0, 100);

        final SimpleTimer timer = new SimpleTimer().start();
        final AFCalcResultTracker resultTracker = testBuilder.makeModel().getLog10PNonRef(vc, testBuilder.makePriors());
        final long runtime = timer.getElapsedTimeNano();
        logger.info("result " + resultTracker.getNormalizedPosteriorOfAFGTZero());
        logger.info("runtime " + runtime);
    }

    private static void analyze(final String[] args) throws Exception {
        final List<String> coreColumns = Arrays.asList("iteration", "n.alt.alleles", "n.samples",
                "exact.model", "prior.type", "runtime", "n.evaluations");

        final PrintStream out = new PrintStream(new FileOutputStream(args[1]));

        final List<ModelParams> modelParams = Arrays.asList(
                new ModelParams(ExactAFCalculationTestBuilder.ModelType.ReferenceDiploidExact, 10000, 10),
//                new ModelParams(ExactAFCalculationTestBuilder.ModelType.GeneralExact, 100, 10),
                new ModelParams(ExactAFCalculationTestBuilder.ModelType.ConstrainedDiploidExact, 10000, 100),
                new ModelParams(ExactAFCalculationTestBuilder.ModelType.IndependentDiploidExact, 10000, 1000));

        final boolean ONLY_HUMAN_PRIORS = false;
        final List<ExactAFCalculationTestBuilder.PriorType> priorTypes = ONLY_HUMAN_PRIORS
                ? Arrays.asList(ExactAFCalculationTestBuilder.PriorType.values())
                : Arrays.asList(ExactAFCalculationTestBuilder.PriorType.human);

        final List<Analysis> analyzes = new ArrayList<Analysis>();
        analyzes.add(new AnalyzeByACAndPL(coreColumns));
        analyzes.add(new AnalyzeBySingletonPosition(coreColumns));
        //analyzes.add(new AnalyzeByNonInformative(coreColumns));

        for ( int iteration = 0; iteration < 1; iteration++ ) {
            for ( final int nAltAlleles : Arrays.asList(1, 2) ) {
                for ( final int nSamples : Arrays.asList(1, 10, 100, 1000, 10000) ) {
                    for ( final ModelParams modelToRun : modelParams) {
                        if ( modelToRun.meetsConstraints(nAltAlleles, nSamples) ) {
                            for ( final ExactAFCalculationTestBuilder.PriorType priorType : priorTypes ) {
                                final ExactAFCalculationTestBuilder testBuilder
                                        = new ExactAFCalculationTestBuilder(nSamples, nAltAlleles, modelToRun.modelType, priorType);

                                for ( final Analysis analysis : analyzes ) {
                                    logger.info(Utils.join("\t", Arrays.asList(iteration, nAltAlleles, nSamples, modelToRun.modelType, priorType, analysis.getName())));
                                    final List<?> values = Arrays.asList(iteration, nAltAlleles, nSamples, modelToRun.modelType, priorType);
                                    analysis.run(testBuilder, (List<Object>)values);
                                }
                            }
                        }
                    }
                }
            }
        }

        final GATKReport report = new GATKReport();
        for ( final Analysis analysis : analyzes )
            report.addTable(analysis.getTable());
        report.print(out);
        out.close();
    }
}