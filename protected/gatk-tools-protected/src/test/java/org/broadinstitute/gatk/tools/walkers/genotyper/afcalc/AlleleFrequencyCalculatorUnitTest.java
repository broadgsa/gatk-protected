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

import htsjdk.variant.variantcontext.*;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 7/28/16.
 */
public class AlleleFrequencyCalculatorUnitTest extends BaseTest {
    private static final double EPS = 1.0e-8;
    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();

    private static final Allele A = Allele.create("A", true);
    private static final Allele B = Allele.create("C");
    private static final Allele C = Allele.create("G");
    private static final Allele indel1 = Allele.create("AA");

    private static final int HAPLOID = 1;
    private static final int DIPLOID = 2;
    private static final int TRIPLOID = 3;

    private static final int BIALLELIC = 2;
    private static final int TRIALLELIC = 3;

    private static final int EXTREMELY_CONFIDENT_PL = 1000;
    private static final int FAIRLY_CONFIDENT_PL = 20;
    private static final int LOW_CONFIDENCE_PL = 10;

    private static final int DEFAULT_PLOIDY = 2;

    private static int sampleNameCounter = 0;

    @Test
    public void testSymmetries() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B,C);
        final Genotype AA = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,2}, FAIRLY_CONFIDENT_PL);
        final Genotype BB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {1,2}, FAIRLY_CONFIDENT_PL);
        final Genotype CC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {2,2}, FAIRLY_CONFIDENT_PL);
        final Genotype AB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,1,1}, FAIRLY_CONFIDENT_PL);
        final Genotype AC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,2,1}, FAIRLY_CONFIDENT_PL);

        final Genotype BBB = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {1,3}, FAIRLY_CONFIDENT_PL);
        final Genotype CCC = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {2,3}, FAIRLY_CONFIDENT_PL);

        // make pairs of VCs tht differ only by B <--> C
        final List<Pair<VariantContext, VariantContext>> switchBWithCPairs = Arrays.asList(
                new Pair<>(makeVC(alleles, AA, BB), makeVC(alleles, AA, CC)),
                new Pair<>(makeVC(alleles, AA, AB), makeVC(alleles, AA, AC)),
                new Pair<>(makeVC(alleles, AB, AB), makeVC(alleles, AC, AC)),
                new Pair<>(makeVC(alleles, AA, AA, BB), makeVC(alleles, AA, AA, CC)),
                new Pair<>(makeVC(alleles, AA, AB, AB), makeVC(alleles, AA, AC, AC)),
                new Pair<>(makeVC(alleles, AA, BBB), makeVC(alleles, AA, CCC))
        );
        for (final Pair<VariantContext, VariantContext> pair : switchBWithCPairs) {
            final VariantContext vc1 = pair.getFirst();
            final VariantContext vc2 = pair.getSecond();
            final AFCalculationResult result1 = afCalc.getLog10PNonRef(vc1);
            final AFCalculationResult result2 = afCalc.getLog10PNonRef(vc2);
            Assert.assertEquals(result1.getLog10PosteriorOfAFEq0(), result2.getLog10PosteriorOfAFEq0(), EPS);
            Assert.assertEquals(result1.getLog10PosteriorOfAFEq0ForAllele(B), result2.getLog10PosteriorOfAFEq0ForAllele(C), EPS);
            Assert.assertEquals(result1.getLog10PosteriorOfAFEq0ForAllele(C), result2.getLog10PosteriorOfAFEq0ForAllele(B), EPS);
        }
    }

    @Test
    public void testMLECounts() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 1, 1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B,C);
        final Genotype AA = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,2}, FAIRLY_CONFIDENT_PL);
        final Genotype BB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {1,2}, FAIRLY_CONFIDENT_PL);
        final Genotype AB = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,1,1}, FAIRLY_CONFIDENT_PL);
        final Genotype AC = genotypeWithObviousCall(DIPLOID, TRIALLELIC, new int[] {0,1,2,1}, FAIRLY_CONFIDENT_PL);

        final Genotype BBB = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {1,3}, FAIRLY_CONFIDENT_PL);
        final Genotype CCC = genotypeWithObviousCall(TRIPLOID, TRIALLELIC, new int[] {2,3}, FAIRLY_CONFIDENT_PL);

        final List<Pair<VariantContext, int[]>> vcWithExpectedCounts = Arrays.asList(
                new Pair<>(makeVC(alleles, AA, BB), new int[] {2,0}),
                new Pair<>(makeVC(alleles, AA, AB), new int[] {1,0}),
                new Pair<>(makeVC(alleles, AB, AB), new int[] {2,0}),
                new Pair<>(makeVC(alleles, AA, AA, BB), new int[] {2,0}),
                new Pair<>(makeVC(alleles, AA, AB, AB), new int[] {2,0}),
                new Pair<>(makeVC(alleles, AA, BBB), new int[] {3,0}),
                new Pair<>(makeVC(alleles, AA, BBB, CCC), new int[] {3,3}),
                new Pair<>(makeVC(alleles, AA, AB, AC), new int[] {1,1}),
                new Pair<>(makeVC(alleles, AA, AB, AC, BBB, CCC), new int[] {4,4})

        );
        for (final Pair<VariantContext, int[]> pair : vcWithExpectedCounts) {
            final VariantContext vc = pair.getFirst();
            final int[] expected = pair.getSecond();
            final int[] actual = afCalc.getLog10PNonRef(vc).getAlleleCountsOfMLE();
            Assert.assertEquals(Arrays.asList(expected), Arrays.asList(actual));
        }
    }

    // many samples with low confidence should yield a non-zero MLE, in contrast to the old exact model
    @Test
    public void testManySamplesWithLowConfidence() {
        // prior corresponding to 1000 observations of ref, 1 of a SNP
        // for this test, we want many pseudocounts in the prior because the new AF calculator learns the allele frequency
        // and we don't want the complication of the posterior being differetn from the prior
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1000, 1, 1, DEFAULT_PLOIDY);    //prior corresponding to 1000 observations of ref, 1 of a SNP
        final List<Allele> alleles = Arrays.asList(A,B);

        // for FAIRLY_CONFIDENT_PL = 20, this genotype has about 100 times greater likelihood to be het than hom ref
        // with our prior giving 1000 times as much weight to ref, this implies a 1 in 5 chance of each sample having a copy of the alt allele
        // (that is, 100/1000 times the combinatorial factor of 2).  Thus the MLE for up to 2 samples should be zero
        // for five samples we should have one
        // for ten samples we will have more than twice as many as for five since the counts fromt he samples start to influence
        // the estimated allele frequency
        final Genotype AB = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,1,1,1}, FAIRLY_CONFIDENT_PL);

        final List<VariantContext> vcsWithDifferentNumbersOfSamples = IntStream.range(1, 11)
                .mapToObj(n -> makeVC(alleles, Collections.nCopies(n, AB))).collect(Collectors.toList());
        final int[] counts = vcsWithDifferentNumbersOfSamples.stream().mapToInt(vc -> afCalc.getLog10PNonRef(vc).getAlleleCountAtMLE(B)).toArray();
        Assert.assertEquals(counts[0],0); // one sample
        Assert.assertEquals(counts[1],0); // two samples
        Assert.assertEquals(counts[4],2); // five samples
        Assert.assertTrue(counts[8] >= 3); // ten samples
    }

    @Test
    public void testApproximateMultiplicativeConfidence() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 1, 1, DEFAULT_PLOIDY);    //flat prior -- we will choose genotypes such that the posterior remains flat
        final List<Allele> alleles = Arrays.asList(A,B);

        final Genotype AA = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,2}, FAIRLY_CONFIDENT_PL);
        final Genotype BB = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {1,2}, FAIRLY_CONFIDENT_PL);

        final List<VariantContext> vcsWithDifferentNumbersOfSamples = new ArrayList<>();
        final List<Genotype> genotypeList = new ArrayList<>();

        for (int n = 0; n < 10; n++) {
            genotypeList.add(AA);
            genotypeList.add(BB);   //adding both keeps the flat prior.  Thus the posterior will equal the likelihood
            vcsWithDifferentNumbersOfSamples.add(makeVC(alleles, genotypeList));
        }

        // since we maintain a flat allele frequency distribution, the probability of being ref as each successive sample is added
        // is multiplied by the probability of any one.  Thus we get an arithmetic series in log space
        final double[] log10PRefs = vcsWithDifferentNumbersOfSamples.stream()
                .mapToDouble(vc -> afCalc.getLog10PNonRef(vc).getLog10LikelihoodOfAFEq0()).toArray();

        for (int n = 0; n < 9; n++) {
            Assert.assertEquals(log10PRefs[n+1] - log10PRefs[n], log10PRefs[0], 0.01);
        }
    }

    @Test
    public void testManyRefSamplesDontKillGoodVariant() {
        final AlleleFrequencyCalculator afCalc = new AlleleFrequencyCalculator(1, 0.1, 0.1, DEFAULT_PLOIDY);
        final List<Allele> alleles = Arrays.asList(A,B);
        final Genotype AA = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,2}, EXTREMELY_CONFIDENT_PL);
        final Genotype AB = genotypeWithObviousCall(DIPLOID, BIALLELIC, new int[] {0,1,1,1}, EXTREMELY_CONFIDENT_PL);
        for (final int numRef : new int[]{1, 10, 100, 1000, 10000, 100000}) {
            final List<Genotype> genotypeList = new ArrayList<>(Collections.nCopies(numRef, AA));
            genotypeList.add(AB);
            final VariantContext vc = makeVC(alleles, genotypeList);
            final double log10PRef = afCalc.getLog10PNonRef(vc).getLog10LikelihoodOfAFEq0();
            Assert.assertTrue(log10PRef < (-EXTREMELY_CONFIDENT_PL/10) + Math.log10(numRef) + 1);
        }
    }

    // make PLs that correspond to an obvious call i.e. one PL is relatively big and the rest are zero
    // alleleCounts is the GenotypeAlleleCounts format for the obvious genotype, with repeats but in no particular order
    private static int[] PLsForObviousCall(final int ploidy, final int numAlleles, final int[] alleleCounts, final int PL)   {
        final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, numAlleles);
        final int[] result = Collections.nCopies(glCalc.genotypeCount(), PL).stream().mapToInt(n->n).toArray();
        result[glCalc.alleleCountsToIndex(alleleCounts)] = 0;
        return result;
    }

    private static Genotype genotypeWithObviousCall(final int ploidy, final int numAlleles, final int[] alleles, final int PL) {
        return makeGenotype(ploidy, PLsForObviousCall(ploidy, numAlleles, alleles, PL));
    }
    //note the call is irrelevant to the AFCalculator, which only looks at PLs
    private static Genotype makeGenotype(final int ploidy, int ... pls) {
        return new GenotypeBuilder("sample" + sampleNameCounter++).alleles(Collections.nCopies(ploidy, Allele.NO_CALL)).PL(pls).make();
    }

    private static VariantContext makeVC(final List<Allele> alleles, final Genotype... genotypes) {
        return new VariantContextBuilder().chr("chr1").alleles(alleles).genotypes(genotypes).make();
    }

    private static VariantContext makeVC(final List<Allele> alleles, final Collection<Genotype> genotypes) {
        return new VariantContextBuilder().chr("chr1").alleles(alleles).genotypes(genotypes).make();
    }
}