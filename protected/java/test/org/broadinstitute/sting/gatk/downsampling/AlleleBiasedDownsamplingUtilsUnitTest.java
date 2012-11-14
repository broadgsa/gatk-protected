/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.downsampling;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;


/**
 * Basic unit test for AlleleBiasedDownsamplingUtils
 */
public class AlleleBiasedDownsamplingUtilsUnitTest extends BaseTest {


    @Test
    public void testSmartDownsampling() {

        final int[] idealHetAlleleCounts = new int[]{0, 50, 0, 50};
        final int[] idealHomAlleleCounts = new int[]{0, 100, 0, 0};

        // no contamination, no removal
        testOneCase(0, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);

        // hom sample, het contaminant, different alleles
        testOneCase(5, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 5, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 0, 5, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);

        // hom sample, hom contaminant, different alleles
        testOneCase(10, 0, 0, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 10, 0, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);
        testOneCase(0, 0, 0, 10, 0.1, 100, idealHomAlleleCounts, idealHomAlleleCounts);

        // het sample, het contaminant, different alleles
        testOneCase(5, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 5, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);

        // het sample, hom contaminant, different alleles
        testOneCase(10, 0, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 10, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);

        // hom sample, het contaminant, overlapping alleles
        final int[] enhancedHomAlleleCounts = new int[]{0, 105, 0, 0};
        testOneCase(5, 5, 0, 0, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts);
        testOneCase(0, 5, 5, 0, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts);
        testOneCase(0, 5, 0, 5, 0.1, 100, idealHomAlleleCounts, enhancedHomAlleleCounts);

        // hom sample, hom contaminant, overlapping alleles
        testOneCase(0, 10, 0, 0, 0.1, 100, idealHomAlleleCounts, new int[]{0, 110, 0, 0});

        // het sample, het contaminant, overlapping alleles
        testOneCase(5, 5, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 5, 5, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 5, 0, 5, 0.1, 100, idealHetAlleleCounts, new int[]{0, 55, 0, 55});
        testOneCase(5, 0, 0, 5, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 5, 5, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);

        // het sample, hom contaminant, overlapping alleles
        testOneCase(0, 10, 0, 0, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
        testOneCase(0, 0, 0, 10, 0.1, 100, idealHetAlleleCounts, idealHetAlleleCounts);
    }

    private static void testOneCase(final int addA, final int addC, final int addG, final int addT, final double contaminationFraction,
                                    final int pileupSize, final int[] initialCounts, final int[] targetCounts) {

        final int[] actualCounts = initialCounts.clone();
        actualCounts[0] += addA;
        actualCounts[1] += addC;
        actualCounts[2] += addG;
        actualCounts[3] += addT;

        final int[] results = AlleleBiasedDownsamplingUtils.runSmartDownsampling(actualCounts, (int)(pileupSize * contaminationFraction));
        Assert.assertTrue(countsAreEqual(results, targetCounts));
    }

    private static boolean countsAreEqual(final int[] counts1, final int[] counts2) {
        for ( int i = 0; i < 4; i++ ) {
            if ( counts1[i] != counts2[i] )
                return false;
        }
        return true;
    }
}
