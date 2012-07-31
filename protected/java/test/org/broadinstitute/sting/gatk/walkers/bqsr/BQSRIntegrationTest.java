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
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, "", "239ce3387b4540faf44ec000d844ccd1")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov ContextCovariate", "d69127341938910c38166dd18449598d")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov CycleCovariate", "b77e621bed1b0dc57970399a35efd0da")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --indels_context_size 4", "2697f38d467a7856c40abce0f778456a")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --low_quality_tail 5", "a55018b1643ca3964dbb50783db9f3e4")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --quantizing_levels 6", "54fe8d1f5573845e6a2aa9688f6dd950")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --mismatches_context_size 4", "6b518ad3c56d66c6f5ea812d058f5c4d")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", "", "3ddb9730f00ee3a612b42209ed9f7e03")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-10,200,000", "", "4cd4fb754e1ef142ad691cb35c74dc4c")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.1RG.bam", "1:10,000,000-10,200,000", "", "364eab693e5e4c7d18a77726b6460f3f")},
                {new BQSRTest(b36KGReference, validationDataLocation + "originalQuals.1kg.chr1.1-1K.1RG.bam", "1:1-1,000", " -OQ", "c449cfca61d605b534f0dce35581339d")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-20,000,000", " --solid_recal_mode REMOVE_REF_BIAS", "5268cb5a4b69335568751d5e5ab80d43")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:anyNameABCD,VCF " + privateTestDir + "vcfexample3.vcf", "3ddb9730f00ee3a612b42209ed9f7e03")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:bed " + validationDataLocation + "bqsrKnownTest.bed", "4a786ba42e38e7fd101947c34a6883ed")},
        };
    }

    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                params.getCommandLine(),
                Arrays.asList(params.md5));
        executeTest("testBQSR-"+params.args, spec).getFirst();

        WalkerTestSpec specNT2 = new WalkerTestSpec(
                params.getCommandLine() + " -nt 2",
                Arrays.asList(params.md5));
        executeTest("testBQSR-nt2-"+params.args, specNT2).getFirst();
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
