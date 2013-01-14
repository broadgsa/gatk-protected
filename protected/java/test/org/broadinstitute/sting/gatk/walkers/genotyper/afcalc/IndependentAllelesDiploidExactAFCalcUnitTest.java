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

package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
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

    private Genotype makePL(final int ... PLs) {
        return AFCalcUnitTest.makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), PLs);
    }

    @Test(enabled = true, dataProvider = "TestCombineGLs")
    private void testCombineGLs(final int altIndex, final int nAlts, final Genotype testg, final Genotype expected) {
        final IndependentAllelesDiploidExactAFCalc calc = (IndependentAllelesDiploidExactAFCalc)AFCalcFactory.createAFCalc(AFCalcFactory.Calculation.EXACT_INDEPENDENT, 1, 4);
        final Genotype combined = calc.combineGLs(testg, altIndex, nAlts);

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
        final Genotype gACcombined2 = makePL(0, 1, 4);
        final Genotype gAGcombined = makePL(0, 4, 9);

        // biallelic
        tests.add(new Object[]{vcAC.genotypes(gACcombined).make(), Arrays.asList(vcAC.genotypes(gACcombined).make())});

        // tri-allelic
        tests.add(new Object[]{vcACG.genotypes(gACG).make(), Arrays.asList(vcAC.genotypes(gACcombined).make(), vcAG.genotypes(gAGcombined).make())});
        tests.add(new Object[]{vcAGC.genotypes(gAGC).make(), Arrays.asList(vcAG.genotypes(gAGcombined).make(), vcAC.genotypes(gACcombined2).make())});

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
                Assert.assertEquals(actual.getGenotype(j).getPL(), expected.getGenotype(j).getPL(),
                        "expected PLs " + Utils.join(",", expected.getGenotype(j).getPL()) + " not equal to actual " + Utils.join(",", actual.getGenotype(j).getPL()));
        }
    }


    @DataProvider(name = "ThetaNTests")
    public Object[][] makeThetaNTests() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Double> log10LAlleles = Arrays.asList(0.0, -1.0, -2.0, -3.0, -4.0);

        for ( final double log10pRef : Arrays.asList(-1, -2, -3) ) {
            for ( final int ploidy : Arrays.asList(1, 2, 3, 4) ) {
                for ( List<Double> permutations : Utils.makePermutations(log10LAlleles, ploidy, true)) {
                    tests.add(new Object[]{permutations, Math.pow(10, log10pRef)});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ThetaNTests")
    public void testThetaNTests(final List<Double> log10LAlleles, final double pRef) {
        // biallelic
        final double[] rawPriors = MathUtils.toLog10(new double[]{pRef, 1-pRef});

        final double log10pNonRef = Math.log10(1-pRef);

        final List<AFCalcResult> originalPriors = new LinkedList<AFCalcResult>();
        final List<Double> pNonRefN = new LinkedList<Double>();
        for ( int i = 0; i < log10LAlleles.size(); i++ ) {
            final double log10LAllele1 = log10LAlleles.get(i);
            final double[] L1 = MathUtils.normalizeFromLog10(new double[]{log10LAllele1, 0.0}, true);
            final AFCalcResult result1 = new AFCalcResult(new int[]{1}, 1, Arrays.asList(A, C), L1, rawPriors, Collections.singletonMap(C, -10000.0));
            originalPriors.add(result1);
            pNonRefN.add(log10pNonRef*(i+1));
        }

        final IndependentAllelesDiploidExactAFCalc calc = (IndependentAllelesDiploidExactAFCalc)AFCalcFactory.createAFCalc(AFCalcFactory.Calculation.EXACT_INDEPENDENT, 1, 2);
        final List<AFCalcResult> thetaNPriors = calc.applyMultiAllelicPriors(originalPriors);

        double prevPosterior = 0.0;
        for ( int i = 0; i < log10LAlleles.size(); i++ ) {
            final AFCalcResult thetaN = thetaNPriors.get(i);
            AFCalcResult orig = null;
            for ( final AFCalcResult x : originalPriors )
                if ( x.getAllelesUsedInGenotyping().equals(thetaN.getAllelesUsedInGenotyping()))
                    orig = x;

            Assert.assertNotNull(orig, "couldn't find original AFCalc");

            Assert.assertEquals(orig.getLog10PriorOfAFGT0(), log10pNonRef, 1e-6);
            Assert.assertEquals(thetaN.getLog10PriorOfAFGT0(), pNonRefN.get(i), 1e-6);

            Assert.assertTrue(orig.getLog10PosteriorOfAFGT0() <= prevPosterior, "AFCalc results should be sorted but " + prevPosterior + " is > original posterior " + orig.getLog10PosteriorOfAFGT0());
            prevPosterior = orig.getLog10PosteriorOfAFGT0();
        }
    }
}
