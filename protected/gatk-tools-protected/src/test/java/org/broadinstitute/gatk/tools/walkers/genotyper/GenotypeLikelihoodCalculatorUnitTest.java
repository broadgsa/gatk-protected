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

package org.broadinstitute.gatk.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests {@link org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculators} and {@link org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculator}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeLikelihoodCalculatorUnitTest {

    @Test(dataProvider = "ploidyAndMaximumAlleleData")
    public void testPloidyAndMaximumAllele(final int ploidy, final int alleleCount) {
        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy, alleleCount);
        Assert.assertNotNull(calculator);
        Assert.assertEquals(calculator.ploidy(),ploidy);
        Assert.assertEquals(calculator.alleleCount(), alleleCount);
        Assert.assertEquals(calculator.genotypeCount(),calculateGenotypeCount(ploidy, alleleCount)," ploidy = " + ploidy + " alleleCount = " + alleleCount);
        final int genotypeCount = calculator.genotypeCount();
        final int testGenotypeCount = Math.min(30000,genotypeCount);
        for (int i = 0; i < testGenotypeCount; i++) {
            final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(i);
            Assert.assertNotNull(alleleCounts);
            if (i > 0)
                Assert.assertTrue(calculator.genotypeAlleleCountsAt(i - 1).compareTo(alleleCounts) < 0);
            final int[] alleleArray = new int[ploidy];
            int index = 0;
            for (int j = 0; j < alleleCounts.distinctAlleleCount(); j++)
                Arrays.fill(alleleArray, index, index += alleleCounts.alleleCountAt(j), alleleCounts.alleleIndexAt(j));
            final int[] alleleCountArray = new int[alleleCounts.distinctAlleleCount() << 1];
            alleleCounts.copyAlleleCounts(alleleCountArray,0);
            Assert.assertEquals(index,ploidy);
            Assert.assertEquals(calculator.allelesToIndex(alleleArray),i);
            Assert.assertEquals(calculator.alleleCountsToIndex(alleleCountArray),i);
        }
    }

    @Test(dataProvider = "ploidyAndMaximumAlleleAndReadCountsData", dependsOnMethods = "testPloidyAndMaximumAllele")
    public void testLikelihoodCalculation(final int ploidy, final int alleleCount, final int[] readCount) {
        final ReadLikelihoods<Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount,readCount);
        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy, alleleCount);
        final int genotypeCount = calculator.genotypeCount();
        final int testGenotypeCount = Math.min(30000,genotypeCount);
        final int sampleCount = readCount.length;
        for (int s = 0; s < sampleCount ; s++) {
            final ReadLikelihoods.Matrix<Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(s);
            final GenotypeLikelihoods genotypeLikelihoods = calculator.genotypeLikelihoods(sampleLikelihoods);
            final double[] genotypeLikelihoodsDoubles = genotypeLikelihoods.getAsVector();
            Assert.assertEquals(genotypeLikelihoodsDoubles.length,genotypeCount);
            for (int i = 0; i < testGenotypeCount; i++) {
                final GenotypeAlleleCounts genotypeAlleleCounts = calculator.genotypeAlleleCountsAt(i);
                Assert.assertNotNull(genotypeLikelihoods);
                final double[] readGenotypeLikelihoods = new double[sampleLikelihoods.readCount()];
                for (int r = 0; r < sampleLikelihoods.readCount(); r++) {
                    final double[] compoments = new double[genotypeAlleleCounts.distinctAlleleCount()];
                    for (int ar = 0; ar < genotypeAlleleCounts.distinctAlleleCount(); ar++) {
                        final int a = genotypeAlleleCounts.alleleIndexAt(ar);
                        final int aCount = genotypeAlleleCounts.alleleCountAt(ar);
                        final double readLk = sampleLikelihoods.get(a, r);
                        compoments[ar] = readLk + Math.log10(aCount);
                    }
                    readGenotypeLikelihoods[r] = MathUtils.approximateLog10SumLog10(compoments) - Math.log10(ploidy);
                }
                final double genotypeLikelihood = MathUtils.sum(readGenotypeLikelihoods);
                Assert.assertEquals(genotypeLikelihoodsDoubles[i], genotypeLikelihood, 0.0001);
            }
        }
    }

    @Test(dataProvider = "ploidyAndMaximumAlleleAndNewMaximumAlleleData")
    public void testGenotypeIndexMap(final int ploidy, final int oldAlleleCount, final int newAlleleCount) {
        final Random rnd = Utils.getRandomGenerator();
        final int maxAlleleCount = Math.max(oldAlleleCount,newAlleleCount);
        final int[] alleleMap = new int[newAlleleCount];
        final Map<Integer,Set<Integer>> reverseMap = new HashMap<>(oldAlleleCount);
        for (int i = 0; i < alleleMap.length; i++) {
            alleleMap[i] = rnd.nextInt(oldAlleleCount);
            if (reverseMap.get(alleleMap[i]) == null) reverseMap.put(alleleMap[i],new HashSet<Integer>(6));
            reverseMap.get(alleleMap[i]).add(i);
        }
        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy,maxAlleleCount);

        final int[] genotypeIndexMap = calculator.genotypeIndexMap(alleleMap);
        Assert.assertNotNull(genotypeIndexMap);
        Assert.assertEquals(genotypeIndexMap.length,GenotypeLikelihoodCalculators.genotypeCount(ploidy,newAlleleCount));

        final GenotypeLikelihoodCalculator oldCalculator = GenotypeLikelihoodCalculators.getInstance(ploidy,oldAlleleCount);
        final GenotypeLikelihoodCalculator newCalculator = GenotypeLikelihoodCalculators.getInstance(ploidy,newAlleleCount);

        for (int i = 0; i < genotypeIndexMap.length; i++) {
            final GenotypeAlleleCounts oldCounts = oldCalculator.genotypeAlleleCountsAt(genotypeIndexMap[i]);
            final GenotypeAlleleCounts newCounts = newCalculator.genotypeAlleleCountsAt(i);
            final int[] reverseCounts = new int[oldAlleleCount];
            for (int j = 0; j < newCounts.distinctAlleleCount(); j++) {
                final int newIndex = newCounts.alleleIndexAt(j);
                final int newRepeats = newCounts.alleleCountAt(j);
                final int expectedOldIndex = alleleMap[newIndex];
                final int oldIndexRank = oldCounts.alleleRankFor(expectedOldIndex);
                Assert.assertNotEquals(oldIndexRank, -1);
                final int oldIndex = oldCounts.alleleIndexAt(oldIndexRank);
                final int oldRepeats = oldCounts.alleleCountAt(oldIndexRank);
                Assert.assertEquals(oldIndex, expectedOldIndex);
                // not necessarily the same count if two or more new alleles map the same old allele.
                Assert.assertTrue(oldRepeats >= newRepeats);
                reverseCounts[oldIndex] += newRepeats;
            }
            for (int j = 0; j < oldAlleleCount; j++)
                Assert.assertEquals(oldCounts.alleleCountFor(j),reverseCounts[j]);
        }
    }


    // Simple inefficient calculation of the genotype count given the ploidy.
    private int calculateGenotypeCount(final int ploidy, final int alleleCount) {
        if (ploidy == 0)
            return 0;
        else if (ploidy == 1)
            return alleleCount;
        else if (ploidy == 2)
            return ((alleleCount) * (alleleCount + 1)) >> 1;
        else if (alleleCount == 0)
            return 0;
        else {
            return calculateGenotypeCount(ploidy - 1, alleleCount) +
                        calculateGenotypeCount(ploidy, alleleCount - 1);
        }
    }

    private static final int[] MAXIMUM_ALLELE = new int[] { 1, 2, 5, 6 };

    private static final int[] PLOIDY = new int[] { 1, 2, 3, 20 };

    private static final int[][] READ_COUNTS = new int[][] {
            { 10 , 100, 50 },
            { 0, 100, 10, 1 , 50 },
            { 1, 2, 3, 4, 20 },
            { 10, 0 },
    };

    @DataProvider(name="ploidyAndMaximumAlleleAndReadCountsData")
    public Object[][] ploidyAndMaximumAlleleAndReadCountsData() {
        final Object[][] result = new Object[PLOIDY.length * MAXIMUM_ALLELE.length * READ_COUNTS.length][];
        int index = 0;
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (final int[] k : READ_COUNTS)
                result[index++] = new Object[] { i, j, k };
        return result;
    }

    @DataProvider(name="ploidyAndMaximumAlleleData")
    public Object[][] ploidyAndMaximumAlleleData() {
        final Object[][] result = new Object[PLOIDY.length * MAXIMUM_ALLELE.length][];
        int index = 0;
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                    result[index++] = new Object[] { i, j };
        return result;
    }

    @DataProvider(name="ploidyAndMaximumAlleleAndNewMaximumAlleleData")
    public Object[][] ploidyAndMaximumAlleleAndNewMaximumAlleleData() {
        final List<Object[]> result = new ArrayList<>(PLOIDY.length * MAXIMUM_ALLELE.length * 20);
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (int k = 0; k < (i < 10? j * 2 : j + 1); k++)
                    result.add(new Object[] { i, j, k });
        return result.toArray(new Object[result.size()][]);
    }
}
