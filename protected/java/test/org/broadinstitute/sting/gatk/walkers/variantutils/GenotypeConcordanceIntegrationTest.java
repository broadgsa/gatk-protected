/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public class GenotypeConcordanceIntegrationTest extends WalkerTest {

    protected static final String emptyMd5 = "d41d8cd98f00b204e9800998ecf8427e";

    public static String baseTestString(String eval, String comp) {
        return "-T GenotypeConcordance -R " + b37KGReference + " --eval " + validationDataLocation + eval + " --comp " + validationDataLocation + comp + " -o %s";
    }

    @Test
    public void testIndelConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("NA12878.Jan2013.haplotypeCaller.subset.indels.vcf", "NA12878.Jan2013.bestPractices.subset.indels.vcf"),
                0,
                Arrays.asList("0f29a0c6dc44066228c8cb204fd53ec0")
        );

        executeTest("test indel concordance", spec);
    }
    
    @Test
    public void testNonoverlapingSamples() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("GenotypeConcordanceNonOverlapTest_Eval.vcf", "GenotypeConcordanceNonOverlapTest_Comp.vcf"),
                0,
                Arrays.asList("fc725022d47b4b5f8a6ef87f0f1ffe89")
        );

        executeTest("test non-overlapping samples", spec);
    }

    @Test
    public void testNonoverlappingSamplesMoltenized() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("GenotypeConcordanceNonOverlapTest_Eval.vcf", "GenotypeConcordanceNonOverlapTest_Comp.vcf"),
                0,
                Arrays.asList("")
        );

        executeTest("Test moltenized output",spec);
    }

    @Test
    public void testMultipleRecordsPerSite() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("GenotypeConcordance.multipleRecordsTest1.eval.vcf","GenotypeConcordance.multipleRecordsTest1.comp.vcf"),
                0,
                Arrays.asList("fdf2cac15775c613f596c27247a76570")
        );

        executeTest("test multiple records per site",spec);
    }

    @Test
    public void testGQFilteringEval() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("genotypeConcordanceFilterTest.vcf","genotypeConcordanceFilterTest.vcf") + " -gfe 'GQ<30'",
                0,
                Arrays.asList("b7b495ccfa6d50a6be3e095d3f6d3c52")
        );

        executeTest("Test filtering on the EVAL rod",spec);
    }

    @Test
    public void testFloatFilteringComp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("genotypeConcordanceFilterTest.vcf","genotypeConcordanceFilterTest.vcf") + " -gfc 'LX<0.50'",
                0,
                Arrays.asList("6406b16cde7960b8943edf594303afd6")
        );

        executeTest("Test filtering on the COMP rod", spec);
    }

    @Test
    public void testCombinedFilters() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString("genotypeConcordanceFilterTest.vcf","genotypeConcordanceFilterTest.vcf") + " -gfc 'LX<0.52' -gfe 'DP<5' -gfe 'GQ<37'",
                0,
                Arrays.asList("26ffd06215b6177acce0ea9f35d73d31")
        );

        executeTest("Test filtering on both rods",spec);
    }
}
