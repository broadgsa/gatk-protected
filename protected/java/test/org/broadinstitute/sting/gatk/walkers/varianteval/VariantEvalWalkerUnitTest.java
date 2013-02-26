/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.varianteval;


// the imports for unit testing.

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class VariantEvalWalkerUnitTest extends BaseTest {
    VariantEval VEwalker;
    VariantContext eval;


    @BeforeMethod
    public void init() {
        VEwalker = new VariantEval();
        eval = new VariantContextBuilder("x", "chr1", 1, 1, Collections.singleton(Allele.create("A", true))).make();
    }

    // --------------------------------------------------------------------------------
    //
    // Test stratifications / evaluations
    //
    // --------------------------------------------------------------------------------

    private class StratifiedEvalTestProvider extends TestDataProvider {
        final List<VariantStratifier> stratificationObjects = new ArrayList<VariantStratifier>();
        final Set<Class<? extends VariantEvaluator>> evaluationObjects = new HashSet<Class<? extends VariantEvaluator>>();
        final List<Integer> expectedCounts;
        final int maxI;

        /**
         *
         * @param maxI test integers from 1 ... maxI
         * @param expectedCounts the expected number of integers from 1 ... maxI divisible by each combination, in order, of allStates
         * @param allStates all stratification tests, in order
         */
        public StratifiedEvalTestProvider(int maxI,
                                          final List<Integer> expectedCounts,
                                          final List<Integer> ... allStates) {
            super(StratifiedEvalTestProvider.class);

            this.maxI = maxI;
            this.expectedCounts = expectedCounts;
            this.evaluationObjects.add(CounterEval.class);

            String stateName = "";
            for ( List<Integer> states : allStates ) {
                stratificationObjects.add(new IntegerStratifier(states));
                stateName = stateName + Utils.join(",", states) + " ";
            }

            setName(String.format("maxI=%d expectedCounts=%s states=%s", maxI, Utils.join(",", expectedCounts), stateName));
        }
    }

    /**
     * Test stratifier -> holds a list of integers, and the states are if the integer value of evalName is divisable
     * by that number
     */
    public static class IntegerStratifier extends VariantStratifier {
        final List<Integer> integers;

        private IntegerStratifier(final List<Integer> integers) {
            this.integers = integers;
            initialize();
        }

        @Override
        public void initialize() {
            states.addAll(integers);
        }

        @Override
        public List<Object> getRelevantStates(final ReferenceContext ref, final RefMetaDataTracker tracker, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName) {
            int i = Integer.valueOf(evalName); // a terrible hack, but we can now provide accessible states
            List<Object> states = new ArrayList<Object>();
            for ( int state : integers )
                if ( i % state == 0 )
                    states.add(state);
            return states;
        }
    }

    /**
     * Test evaluator -> just counts the number of calls to update1
     */
    public static class CounterEval extends VariantEvaluator {
        public int count = 0;

        @Override public int getComparisonOrder() { return 1; }

        @Override
        public void update1(final VariantContext eval, final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
            count++;
        }

        @Override
        public boolean supportsCombine() {
            return true;
        }

        @Override
        public void combine(final VariantEvaluator other) {
            this.count += ((CounterEval)other).count;
        }
    }

    private void initialize(StratifiedEvalTestProvider cfg) {
        VEwalker.createStratificationStates(cfg.stratificationObjects, cfg.evaluationObjects);

        final RefMetaDataTracker tracker = new RefMetaDataTracker();
        final ReferenceContext ref = null;
        final VariantContext comp = null;
        final String compName = null, sampleName = null;

        // increment eval counts for each stratification of divisors of i from from 1...maxI
        for ( int i = 1; i <= cfg.maxI; i++ ) {
            final String evalName = String.valueOf(i); // terrible hack to stratify by divisor
            for ( EvaluationContext nec : VEwalker.getEvaluationContexts(tracker, ref, eval, evalName, comp, compName, sampleName) ) {
                synchronized (nec) {
                    nec.apply(tracker, ref, null, comp, eval);
                }
            }
        }
    }

    @DataProvider(name = "StratifiedEvalTestProvider")
    public Object[][] makeStratifiedEvalTestProvider() {

        new StratifiedEvalTestProvider(4, // test 1, 2, 3, 4
                Arrays.asList(4, 2), //  4 divisible by 1, 2 by 2
                Arrays.asList(1, 2));

        new StratifiedEvalTestProvider(6, // test 1, 2, 3, 4, 5, 6
                Arrays.asList(6, 3, 2), //  6 divisible by 1, 3 by 2, 2 by 3
                Arrays.asList(1, 2, 3));

        // test that some states can be empty -- does this work in VE?
        new StratifiedEvalTestProvider(6,
                Arrays.asList(3, 2),
                Arrays.asList(2, 3));

        // test a single stratification
        new StratifiedEvalTestProvider(6,
                Arrays.asList(3),
                Arrays.asList(2));

        // test a meaningless state
        new StratifiedEvalTestProvider(4, // test 1, 2, 3, 4
                Arrays.asList(4, 2), //  4 divisible by 1, 2 by 2
                Arrays.asList(1, 2), Arrays.asList(1));

        // test a adding a state that divides space in half
        new StratifiedEvalTestProvider(4,
                Arrays.asList(2, 2),
                Arrays.asList(1, 2), Arrays.asList(2));

        // test pairs of strats
        new StratifiedEvalTestProvider(12,
                Arrays.asList(4, 3, 2, 3),
                Arrays.asList(1, 2), Arrays.asList(3, 4));

        return StratifiedEvalTestProvider.getTests(StratifiedEvalTestProvider.class);
    }

    /**
     * Ensures that counting and stratifications all are working properly by iterating
     * over integers 1...cfg.N and stratify according to cfg, and that the counts in
     * each bin are as expected.
     *
     * @param cfg
     */
    @Test(dataProvider = "StratifiedEvalTestProvider")
    public void testBasicOperation(StratifiedEvalTestProvider cfg) {
        initialize(cfg);
        checkStratificationCountsAreExpected(VEwalker.stratManager, cfg.expectedCounts);
    }

    private final void checkStratificationCountsAreExpected(final StratificationManager<VariantStratifier, EvaluationContext> manager,
                                                            final List<Integer> expectedCounts) {
        for ( int key = 0; key < manager.size(); key++ ) {
            final String stratStateString = manager.getStratsAndStatesStringForKey(key);
            final EvaluationContext nec = manager.get(key);

            for ( final VariantEvaluator ve : nec.getVariantEvaluators() ) {
                // test for count here
                final CounterEval counterEval = (CounterEval)ve;
                final int expected = expectedCounts.get(key);
                Assert.assertEquals(counterEval.count, expected, "Count seen of " + counterEval.count + " not expected " + expected + " at " + stratStateString);
            }
        }
    }

    /**
     * A derived test on testBasicOperation that checks that combining stratifications
     * works as expected by ensuring the results are the same when the remapped
     * strats are the identity map (A -> A, B -> B, etc)
     */
    @Test(dataProvider = "StratifiedEvalTestProvider", dependsOnMethods = {"testBasicOperation"})
    public void testIdentityCombine(StratifiedEvalTestProvider cfg) {
        for ( int i = 0; i < cfg.stratificationObjects.size(); i++ ) {
            initialize(cfg);
            final VariantStratifier toReplace = cfg.stratificationObjects.get(i);
            final VariantStratifier newStrat = cfg.stratificationObjects.get(i);
            final Map<Object, Object> remappedStates = Utils.makeIdentityFunctionMap(newStrat.getAllStates());
            StratificationManager<VariantStratifier, EvaluationContext> combined =
                    VEwalker.stratManager.combineStrats(toReplace, newStrat, EvaluationContext.COMBINER, remappedStates);
            checkStratificationCountsAreExpected(combined, cfg.expectedCounts);
        }
    }

//    /**
//     * A derived test on testBasicOperation that checks that combining stratifications
//     * works as expected. We look into cfg, and if there are multiple states we create
//     * dynamically create a combinations of the stratifications, and ensure that the
//     * combined results are as we expected.
//     */
//    @Test(dataProvider = "StratifiedEvalTestProvider", dependsOnMethods = {"testBasicOperation"})
//    public void testCombinedEachStrat(StratifiedEvalTestProvider cfg) {
//        for ( int i = 0; i < cfg.stratificationObjects.size(); i++ ) {
//            initialize(cfg);
//            final VariantStratifier toReplace = cfg.stratificationObjects.get(i);
//
//            // TODO -- replace this code with something that combines values in strat
//            final VariantStratifier newStrat = cfg.stratificationObjects.get(i);
//            final Map<Object, Object> remappedStates = Utils.makeIdentityFunctionMap(newStrat.getAllStates());
//            final List<Integer> expected = cfg.expectedCounts;
//
//            StratificationManager<VariantStratifier, EvaluationContext> combined =
//                    VEwalker.stratManager.combineStrats(toReplace, newStrat, EvaluationContext.COMBINER, remappedStates);
//            checkStratificationCountsAreExpected(combined, expected);
//        }
//    }
}