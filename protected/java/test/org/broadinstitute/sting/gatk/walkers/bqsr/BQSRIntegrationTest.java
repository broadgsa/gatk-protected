package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author ebanks
 * @since 7/16/12
 */
public class BQSRIntegrationTest extends WalkerTest {

    private static class BQSRTest {
        final String reference;
        final String interval;
        final String bam;
        final String args;
        final String md5;

        private BQSRTest(String reference, String bam, String interval, String args, String md5) {
            this.reference = reference;
            this.bam = bam;
            this.interval = interval;
            this.args = args;
            this.md5 = md5;
        }

        public String getCommandLine() {
            return " -T BaseRecalibrator" +
                    " -R " + reference +
                    " -I " + bam +
                    " -L " + interval +
                    args +
                    " -knownSites " + (reference.equals(b36KGReference) ? b36dbSNP129 : hg18dbSNP132) +
                    " --allow_potentially_misencoded_quality_scores" +  // TODO -- remove me when we get new SOLiD bams
                    " -o %s" +
                    " -sortAllCols";
        }

        @Override
        public String toString() {
            return String.format("BQSR(bam='%s', args='%s')", bam, args);
        }
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        String HiSeqBam = privateTestDir + "HiSeq.1mb.1RG.bam";
        String HiSeqInterval = "chr1:10,000,000-10,100,000";
        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, "", "2f250fecb930e0dfe0f63fe0fed3960b")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov ContextCovariate", "26c8d7226139a040557b1d3b1c8792f0")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov CycleCovariate", "9b43a1839cb6ea03aec1d96f15ca8efb")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --indels_context_size 4", "3159a9d136c45e4a65d46a23dc8fd3b5")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --low_quality_tail 5", "bb7262829effbbdbc8d88dd36f480368")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --quantizing_levels 6", "fbb002fa2b9197c4b555852dccc11562")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --mismatches_context_size 4", "7392acb71131a60a527ca32715fc59be")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", "", "49d4383896a90795d94138db1410a7df")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-10,200,000", "", "427448eff98cf194cc7217c0b1401e79")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.1RG.bam", "1:10,000,000-10,200,000", "", "50cd1a10b6ecb3d09f90f1e4a66da95d")},
                {new BQSRTest(b36KGReference, validationDataLocation + "originalQuals.1kg.chr1.1-1K.1RG.bam", "1:1-1,000", " -OQ", "1dc71561c9d0fb56f9876cb5043c5376")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-20,000,000", " --solid_recal_mode REMOVE_REF_BIAS", "13e8f032e76340b114847c90af0a1f8a")},
                {new BQSRTest(b36KGReference, privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam", "1:50,000-80,000", " --solid_nocall_strategy LEAVE_READ_UNRECALIBRATED", "03f58ae4f9d203034e895a3636fc108f")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:anyNameABCD,VCF " + privateTestDir + "vcfexample3.vcf", "49d4383896a90795d94138db1410a7df")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:bed " + validationDataLocation + "bqsrKnownTest.bed", "2db2ef8c2d63e167663d70340182f49a")},
        };
    }

    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.md5));
        executeTest("testBQSR-"+params.args, spec).getFirst();
    }

    @Test
    public void testBQSRFailWithoutDBSNP() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam" +
                        " -L 1:10,000,000-10,200,000" +
                        " -o %s",
                1, // just one output file
                UserException.CommandLineException.class);
        executeTest("testBQSRFailWithoutDBSNP", spec);
    }

    @Test
    public void testBQSRCSV() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam" +
                        " -knownSites " + b36dbSNP129 +
                        " -L 1:10,000,000-10,200,000" +
                        " -o /dev/null" +
                        " -sortAllCols" +
                        " --plot_pdf_file /dev/null" +
                        " --intermediate_csv_file %s",
                Arrays.asList("d1c38a3418979400630e2bca1140689c"));
        executeTest("testBQSR-CSVfile", spec);
    }

    @Test
    public void testBQSRFailWithSolidNoCall() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam" +
                        " -L 1:50,000-80,000" +
                        " --allow_potentially_misencoded_quality_scores" +  // TODO -- remove me when we get new SOLiD bams
                        " -o %s" +
                        " -sortAllCols",
                1, // just one output file
                UserException.class);
        executeTest("testBQSRFailWithSolidNoCall", spec);
    }

    private static class PRTest {
        final String args;
        final String md5;

        private PRTest(String args, String md5) {
            this.args = args;
            this.md5 = md5;
        }

        @Override
        public String toString() {
            return String.format("PrintReads(args='%s')", args);
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createPRTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{1, new PRTest(" -qq -1", "5226c06237b213b9e9b25a32ed92d09a")});
        tests.add(new Object[]{1, new PRTest(" -qq 6", "b592a5c62b952a012e18adb898ea9c33")});
        tests.add(new Object[]{1, new PRTest(" -DIQ", "8977bea0c57b808e65e9505eb648cdf7")});

        for ( final int nct : Arrays.asList(1, 2, 4) ) {
            tests.add(new Object[]{nct, new PRTest("", "ab2f209ab98ad3432e208cbd524a4c4a")});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PRTest")
    public void testPR(final int nct, PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + privateTestDir + "HiSeq.1mb.1RG.bam" +
                        " -nct " + nct +
                        " -BQSR " + privateTestDir + "HiSeq.20mb.1RG.table" +
                        params.args +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }
}
