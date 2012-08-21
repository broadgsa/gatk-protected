package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

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
                    " --no_plots" +
                    " -knownSites " + (reference.equals(b36KGReference) ? b36dbSNP129 : hg18dbSNP132) +
                    " -o %s";
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
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, "", "1cfc73371abb933ca26496745d105ff0")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov ContextCovariate", "ee5142776008741b1b2453b1258c6d99")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov CycleCovariate", "fbc520794f0f98d52159de956f7217f1")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --indels_context_size 4", "ab5b93794049c514bf8e407019d76b67")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --low_quality_tail 5", "81df636e3d0ed6f16113517e0169bc96")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --quantizing_levels 6", "ad3c47355448f8c45e172c6e1129c65d")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --mismatches_context_size 4", "fef7240140a9b6d6335ce009fa4edec5")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", "", "600652ee49b9ce1ca2d8ee2d8b7c8211")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-10,200,000", "", "769f95b9dcc78a405d3e6b191e5a19f5")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.1RG.bam", "1:10,000,000-10,200,000", "", "43fcba51264cc98bd8466d21e1b96766")},
                {new BQSRTest(b36KGReference, validationDataLocation + "originalQuals.1kg.chr1.1-1K.1RG.bam", "1:1-1,000", " -OQ", "48aaf9ac54b97eac3663882a59354ab2")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-20,000,000", " --solid_recal_mode REMOVE_REF_BIAS", "dac04b9e1e1c52af8d3a50c2e550fda9")},
                {new BQSRTest(b36KGReference, privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam", "1:50,000-80,000", " --solid_nocall_strategy LEAVE_READ_UNRECALIBRATED", "90d70542076715a8605a8d4002614b34")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:anyNameABCD,VCF " + privateTestDir + "vcfexample3.vcf", "600652ee49b9ce1ca2d8ee2d8b7c8211")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:bed " + validationDataLocation + "bqsrKnownTest.bed", "26a04f5a28c40750c603cbe8a926d7bd")},
        };
    }

    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.md5));
        executeTest("testBQSR-"+params.args, spec).getFirst();

        // TODO -- re-enable once parallelization is fixed in BaseRecalibrator
        //WalkerTestSpec specNT2 = new WalkerTestSpec(
        //        params.getCommandLine() + " -nt 2",
        //        Arrays.asList(params.md5));
        //executeTest("testBQSR-nt2-"+params.args, specNT2).getFirst();
    }

    @Test
    public void testBQSRFailWithoutDBSNP() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam" +
                        " -L 1:10,000,000-10,200,000" +
                        " --no_plots" +
                        " -o %s",
                1, // just one output file
                UserException.CommandLineException.class);
        executeTest("testBQSRFailWithoutDBSNP", spec);
    }

    @Test
    public void testBQSRFailWithSolidNoCall() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T BaseRecalibrator" +
                        " -R " + b36KGReference +
                        " -I " + privateTestDir + "NA19240.chr1.BFAST.SOLID.hasCSNoCall.bam" +
                        " -L 1:50,000-80,000" +
                        " --no_plots" +
                        " -o %s",
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
        return new Object[][]{
                {new PRTest("", "d2d6ed8667cdba7e56f5db97d6262676")},
                {new PRTest(" -qq -1", "b7053d3d67aba6d8892f0a60f0ded338")},
                {new PRTest(" -qq 6", "bfbf0855185b2b70aa35237fb71e4487")},
                {new PRTest(" -DIQ", "66aa65223f192ee39c1773aa187fd493")}
        };
    }

    @Test(dataProvider = "PRTest")
    public void testPR(PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + privateTestDir + "HiSeq.1mb.1RG.bam" +
                        " -BQSR " + privateTestDir + "HiSeq.1mb.1RG.table" +
                        params.args +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }
}
