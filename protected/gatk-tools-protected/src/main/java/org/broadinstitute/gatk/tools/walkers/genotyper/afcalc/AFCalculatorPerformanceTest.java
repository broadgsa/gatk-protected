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

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.TTCCLayout;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.SimpleTimer;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;

import java.io.*;
import java.util.*;

/**
 * A simple GATK utility (i.e, runs from command-line) for assessing the performance of
 * the exact model
 */
public class AFCalculatorPerformanceTest {
    final static Logger logger = Logger.getLogger(AFCalculatorPerformanceTest.class);

    private static abstract class Analysis {
        final GATKReport report;

        public Analysis(final String name, final List<String> columns) {
            report = GATKReport.newSimpleReport(name, columns);
        }

        public abstract void run(final AFCalculatorTestBuilder testBuilder,
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

        public void run(final AFCalculatorTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final AFCalculator calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                for ( int[] ACs : makeACs(testBuilder.numAltAlleles, testBuilder.nSamples*2) ) {
                    final VariantContext vc = testBuilder.makeACTest(ACs, 0, nonTypePL);

                    timer.start();
                    final AFCalculationResult resultTracker = calc.getLog10PNonRef(vc, HomoSapiensConstants.DEFAULT_PLOIDY, testBuilder.numAltAlleles, priors);
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

        public void run(final AFCalculatorTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final AFCalculator calc = testBuilder.makeModel();
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
                    final AFCalculationResult resultTracker = calc.getLog10PNonRef(vcb.make(), HomoSapiensConstants.DEFAULT_PLOIDY, testBuilder.numAltAlleles, priors);
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

        public void run(final AFCalculatorTestBuilder testBuilder, final List<Object> coreValues) {
            final SimpleTimer timer = new SimpleTimer();

            for ( final int nonTypePL : Arrays.asList(100) ) {
                final AFCalculator calc = testBuilder.makeModel();
                final double[] priors = testBuilder.makePriors();

                final int[] ac = new int[testBuilder.numAltAlleles];
                ac[0] = 1;

                for ( int nNonInformative = 0; nNonInformative < testBuilder.nSamples; nNonInformative++ ) {
                    final VariantContext vc = testBuilder.makeACTest(ac, nNonInformative, nonTypePL);

                    timer.start();
                    final AFCalculationResult resultTracker = calc.getLog10PNonRef(vc, HomoSapiensConstants.DEFAULT_PLOIDY, testBuilder.numAltAlleles, priors);
                    final long runtime = timer.getElapsedTimeNano();

                    final List<Object> columns = new LinkedList<Object>(coreValues);
                    columns.addAll(Arrays.asList(runtime, resultTracker.getnEvaluations(), nonTypePL, nNonInformative));
                    report.addRowList(columns);
                }
            }
        }
    }

    private static class ModelParams {
        final AFCalculatorImplementation modelType;
        final int maxBiNSamples, maxTriNSamples;

        private ModelParams(AFCalculatorImplementation modelType, int maxBiNSamples, int maxTriNSamples) {
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
        SINGLE,
        EXACT_LOG
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
            case EXACT_LOG: exactLog(args); break;
            default: throw new IllegalAccessException("unknown operation " + op);
        }
    }

    private static void exactLog(final String[] args) throws Exception {
        final File ref = new File(args[1]);
        final File exactLogFile = new File(args[2]);
        final List<Integer> startsToUse = new LinkedList<Integer>();

        for ( int i = 3; i < args.length; i++ )
            startsToUse.add(Integer.valueOf(args[i]));

        final CachingIndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(ref);
        final GenomeLocParser parser = new GenomeLocParser(seq);
        final BufferedReader reader = new BufferedReader(new FileReader(exactLogFile));
        final List<ExactCallLogger.ExactCall> loggedCalls = ExactCallLogger.readExactLog(reader, startsToUse, parser);

        for ( final ExactCallLogger.ExactCall call : loggedCalls ) {
            final AFCalculatorTestBuilder testBuilder = new AFCalculatorTestBuilder(call.vc.getNSamples(), 1,
                    AFCalculatorImplementation.EXACT_INDEPENDENT,
                    AFCalculatorTestBuilder.PriorType.human);
            logger.info(call);
            final SimpleTimer timer = new SimpleTimer().start();
            final AFCalculationResult result = testBuilder.makeModel().getLog10PNonRef(call.vc, HomoSapiensConstants.DEFAULT_PLOIDY, testBuilder.numAltAlleles,testBuilder.makePriors());
            final long newNanoTime = timer.getElapsedTimeNano();
            if ( call.originalCall.anyPolymorphic(-1) || result.anyPolymorphic(-1) ) {
                logger.info("**** ONE IS POLY");
            }
            logger.info("\t\t getLog10PosteriorOfAFGT0: " + call.originalCall.getLog10PosteriorOfAFGT0() + " vs " + result.getLog10PosteriorOfAFGT0());
            final double speedup = call.runtime / (1.0 * newNanoTime);
            logger.info("\t\t runtime:                  " + call.runtime + " vs " + newNanoTime + " speedup " + String.format("%.2f", speedup) + "x");
            for ( final Allele a : call.originalCall.getAllelesUsedInGenotyping() ) {
                if ( a.isNonReference() ) {
                    final String warningmeMLE = call.originalCall.getAlleleCountAtMLE(a) != result.getAlleleCountAtMLE(a) ? " DANGER-MLE-DIFFERENT" : "";
                    logger.info("\t\t   MLE       " + a + ":            " + call.originalCall.getAlleleCountAtMLE(a) + " vs " + result.getAlleleCountAtMLE(a) + warningmeMLE);
                    final String warningmePost = call.originalCall.getLog10PosteriorOfAFEq0ForAllele(a) == 0 && result.getLog10PosteriorOfAFEq0ForAllele(a) < -10 ? " DANGER-POSTERIORS-DIFFERENT" : "";
                    logger.info("\t\t   Posterior " + a + ":            " + call.originalCall.getLog10PosteriorOfAFEq0ForAllele(a) + " vs " + result.getLog10PosteriorOfAFEq0ForAllele(a) + warningmePost);
                }
            }
        }
    }

    private static void profileBig(final String[] args) throws Exception {
        final int nSamples = Integer.valueOf(args[1]);
        final int ac = Integer.valueOf(args[2]);

        final AFCalculatorTestBuilder testBuilder = new AFCalculatorTestBuilder(nSamples, 1,
                AFCalculatorImplementation.EXACT_INDEPENDENT,
                AFCalculatorTestBuilder.PriorType.human);

        final VariantContext vc = testBuilder.makeACTest(new int[]{ac}, 0, 100);

        final SimpleTimer timer = new SimpleTimer().start();
        final AFCalculationResult resultTracker = testBuilder.makeModel().getLog10PNonRef(vc, HomoSapiensConstants.DEFAULT_PLOIDY, testBuilder.numAltAlleles, testBuilder.makePriors());
        final long runtime = timer.getElapsedTimeNano();
        logger.info("result " + resultTracker.getLog10PosteriorOfAFGT0());
        logger.info("runtime " + runtime);
    }

    private static void analyze(final String[] args) throws Exception {
        final List<String> coreColumns = Arrays.asList("iteration", "n.alt.alleles", "n.samples",
                "exact.model", "prior.type", "runtime", "n.evaluations");

        final PrintStream out = new PrintStream(new FileOutputStream(args[1]));

        final List<ModelParams> modelParams = Arrays.asList(
                new ModelParams(AFCalculatorImplementation.EXACT_REFERENCE, 10000, 10),
//                new ModelParams(AFCalcTestBuilder.ModelType.GeneralExact, 100, 10),
                new ModelParams(AFCalculatorImplementation.EXACT_INDEPENDENT, 10000, 1000));

        final boolean ONLY_HUMAN_PRIORS = false;
        final List<AFCalculatorTestBuilder.PriorType> priorTypes = ONLY_HUMAN_PRIORS
                ? Arrays.asList(AFCalculatorTestBuilder.PriorType.values())
                : Arrays.asList(AFCalculatorTestBuilder.PriorType.human);

        final List<Analysis> analyzes = new ArrayList<Analysis>();
        analyzes.add(new AnalyzeByACAndPL(coreColumns));
        analyzes.add(new AnalyzeBySingletonPosition(coreColumns));
        //analyzes.add(new AnalyzeByNonInformative(coreColumns));

        for ( int iteration = 0; iteration < 1; iteration++ ) {
            for ( final int nAltAlleles : Arrays.asList(1, 2) ) {
                for ( final int nSamples : Arrays.asList(1, 10, 100, 1000, 10000) ) {
                    for ( final ModelParams modelToRun : modelParams) {
                        if ( modelToRun.meetsConstraints(nAltAlleles, nSamples) ) {
                            for ( final AFCalculatorTestBuilder.PriorType priorType : priorTypes ) {
                                final AFCalculatorTestBuilder testBuilder
                                        = new AFCalculatorTestBuilder(nSamples, nAltAlleles, modelToRun.modelType, priorType);

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