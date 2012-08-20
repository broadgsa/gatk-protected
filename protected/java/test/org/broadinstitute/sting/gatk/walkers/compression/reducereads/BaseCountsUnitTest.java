// our package
package org.broadinstitute.sting.gatk.walkers.compression.reducereads;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Basic unit test for BaseCounts in reduced reads
 */
public class BaseCountsUnitTest extends BaseTest {
    private class SingleTest {
        public String bases;
        public byte mostCountBase;
        public int mostCommonCount;

        private SingleTest(String bases, char mostCountBase, int mostCommonCount) {
            this.mostCommonCount = mostCommonCount;
            this.mostCountBase = (byte)mostCountBase;
            this.bases = bases;
        }
    }


    @DataProvider(name = "data")
    public Object[][] createData1() {
        List<SingleTest> params = new ArrayList<SingleTest>();

        params.add(new SingleTest("A", 'A', 1 ));
        params.add(new SingleTest("AA", 'A', 2 ));
        params.add(new SingleTest("AC", 'A', 1 ));
        params.add(new SingleTest("AAC", 'A', 2 ));
        params.add(new SingleTest("AAA", 'A', 3 ));
        params.add(new SingleTest("AAAN", 'A', 3 ));
        params.add(new SingleTest("AAANNNN", 'N', 4 ));
        params.add(new SingleTest("AACTG", 'A', 2 ));
        params.add(new SingleTest("D", 'D', 1 ));
        params.add(new SingleTest("DDAAD", 'D', 3));
        params.add(new SingleTest("", (char)BaseCounts.MAX_BASE_WITH_NO_COUNTS, 0 ));
        params.add(new SingleTest("AAIIIAI", 'I', 4 ));

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( SingleTest x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }



    @Test(dataProvider = "data", enabled = true)
    public void testCounting(SingleTest params) {
        BaseCounts counts = new BaseCounts();

        for ( byte base : params.bases.getBytes() )
            counts.incr(base);

        String name = String.format("Test-%s", params.bases);
        Assert.assertEquals(counts.totalCount(), params.bases.length(), name);
        Assert.assertEquals(counts.countOfMostCommonBase(), params.mostCommonCount, name);
        Assert.assertEquals((char)counts.baseWithMostCounts(), (char)params.mostCountBase, name);
    }
}