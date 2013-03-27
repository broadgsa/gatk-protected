/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class KMerErrorCorrectorUnitTest extends BaseTest {
    @Test
    public void testMyData() {
        final KMerErrorCorrector corrector = new KMerErrorCorrector(3, 1, 2, 2);

        Assert.assertNotNull(corrector.toString());

        corrector.addKmers(
                "ATG", "ATG", "ATG", "ATG",
                "ACC", "ACC", "ACC",
                "AAA", "AAA",
                "CTG", // -> ATG
                "NNA", // -> AAA
                "CCC", // => ACC
                "NNN", // => null
                "NNC"  // => ACC [because of min count won't go to NNA]
        );

        testCorrection(corrector, "ATG", "ATG");
        testCorrection(corrector, "ACC", "ACC");
        testCorrection(corrector, "AAA", "AAA");
        testCorrection(corrector, "CTG", "ATG");
        testCorrection(corrector, "NNA", "AAA");
        testCorrection(corrector, "CCC", "ACC");
        testCorrection(corrector, "NNN", null);
        testCorrection(corrector, "NNC", "ACC");

        Assert.assertNotNull(corrector.toString());
    }

    private void testCorrection(final KMerErrorCorrector corrector, final String in, final String out) {
        Assert.assertEquals(corrector.getErrorCorrectedKmer(in), out);
        Assert.assertEquals(corrector.getErrorCorrectedKmer(in.getBytes()), out == null ? null : out.getBytes());
    }
}
