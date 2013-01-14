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
}
