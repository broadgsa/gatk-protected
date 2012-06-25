package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.HashMap;
import java.util.Map;

/**
 * Given two BAMs with different read groups, it compares them based on ReduceReads metrics.
 * <p>
 * This is a test walker used for asserting that the ReduceReads procedure is not making blatant mistakes when compressing bam files.
 * </p>
 * <h2>Input</h2>
 * <p>
 * Two BAM files (using -I) with different read group IDs
 * </p>
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *  </pre>
 *
 * @author carneiro
 * @since 10/30/11
 */

@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckFilter.class})
public class CompareBAMWalker extends LocusWalker<Map<CompareBAMWalker.TestName, Boolean>, CompareBAMWalker.TestResults> {
    @Argument(required = true,  shortName = "rr",  fullName = "reduced_readgroup", doc = "The read group ID corresponding to the compressed BAM being tested") public String reducedReadGroupID;
    @Argument(required = false, shortName = "teq", fullName = "test_equal_bases",  doc = "Test if the bases marked as '=' are indeed ref bases.")              public boolean TEST_EQUAL_BASES = false;
    @Argument(required = false, shortName = "tbc", fullName = "test_base_counts",  doc = "Test if the base counts tag in consensus reads are accurate.")       public boolean TEST_BASE_COUNTS = false;
    @Argument(required = false, shortName = "mbq", fullName = "min_base_qual",     doc = "Minimum base quality to be considered.")                             public int MIN_BASE_QUAL    = 20;
    @Argument(required = false, shortName = "mmq", fullName = "min_mapping_qual",  doc = "Minimum mapping quality to be considered.")                          public int MIN_MAPPING_QUAL = 20;


    @Override
    public Map<TestName, Boolean> map (RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Map<TestName, Boolean> result = new HashMap<TestName, Boolean>();

        if (TEST_EQUAL_BASES) result.put(TestName.EQUAL_BASES, testEqualBases(ref, context));
        if (TEST_BASE_COUNTS) result.put(TestName.BASE_COUNTS, testBaseCounts(ref, context));

        return result;
    }

    @Override
    public TestResults reduceInit () {
        TestResults sum = new TestResults();      // a fresh new TestResults object to sum up the results of every object passed by MAP.

        if (TEST_EQUAL_BASES) sum.createTest(TestName.EQUAL_BASES);
        if (TEST_BASE_COUNTS) sum.createTest(TestName.BASE_COUNTS);

        return sum;
    }

    @Override
    public TestResults reduce (Map<TestName,Boolean> mapResult, TestResults sum) {
        for (TestName test : mapResult.keySet()) {
            if (mapResult.get(test))
                sum.reportSuccess(test);
            else
                sum.reportFailed(test);
        }

        return sum;
    }

    public void onTraversalDone (TestResults finalResults) {
        finalResults.report();
    }

    private boolean testEqualBases (ReferenceContext ref, AlignmentContext context) {
        return true;
    }

    private boolean testBaseCounts (ReferenceContext ref, AlignmentContext context) {

        return true;
    }

    public enum TestName {
        EQUAL_BASES ("testEqualBases"),
        BASE_COUNTS ("testBaseCounts");

        private String testName;

        TestName(String testName) {
            this.testName = testName;
        }

        public String getTestName() {
            return testName;
        }
    }

    public class TestResults {
        private Map<TestName, TestOutcome> testStats = new HashMap<TestName, TestOutcome>();

        public void createTest (TestName test) {
            testStats.put(test, new TestOutcome());
        }

        public void reportSuccess(TestName test) {
            if (testStats.containsKey(test))
                testStats.get(test).incPassed();
            else
                throw new ReviewedStingException("No such test: " + test);
        }

        public void reportFailed(TestName test) {
            if (testStats.containsKey(test))
                testStats.get(test).incFailed();
            else
                throw new ReviewedStingException("No such test: " + test);
        }

        public void report() {
            System.out.println();
            System.out.println(String.format("%20s\tPASS\tFAIL", ""));
            for (TestName test : testStats.keySet())
                System.out.println(String.format("%20s\t%d\t%d", test.getTestName(), testStats.get(test).getPassed(), testStats.get(test).getFailed()));
            System.out.println();
        }
    }

    private class TestOutcome {
        private long passed;
        private long failed;

        public long getPassed() {
            return passed;
        }

        public void incPassed() {
            this.passed++;
        }

        public long getFailed() {
            return failed;
        }

        public void incFailed() {
            this.failed++;
        }
    }

    private BaseCounts getFilteredBaseCounts(AlignmentContext context) {
        return getBaseCounts(context, MIN_BASE_QUAL, MIN_MAPPING_QUAL);
    }

    private BaseCounts getFullBaseCounts(AlignmentContext context) {
        return getBaseCounts(context, 3, 0);
    }

    private BaseCounts getBaseCounts(AlignmentContext context, int mbq, int mmq) {
        BaseCounts fullBaseCounts = new BaseCounts();
        for (String rg : context.getBasePileup().getReadGroups()) {
            if (!rg.equals(reducedReadGroupID)) {
                BaseCounts b = BaseCounts.createWithCounts(context.getBasePileup().getPileupForReadGroup(rg).getBaseAndMappingFilteredPileup(mbq, mmq).getBaseCounts());
                fullBaseCounts.add(b);
            }
        }
        return fullBaseCounts;
    }


}
