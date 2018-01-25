package org.broadinstitute.gatk.tools.walkers;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class RemoveAnnotationsIntegrationTest extends WalkerTest {
    @Test
    public void testBasicOperation() {
        String testDir = "/Users/bimber/Desktop/TestData/";

        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T RemoveAnnotations",
                                        "-R " + testDir + "18_MacaM.fasta",
                                        "-V " + testDir + "GBS_PstI.vcf.gz",
                                        "-L chr01",
                                        "-A AF",
                                        "-XA MAF",
                                        "-GA DP",
                                        "-XGA AD",
                                        "-ef",
                                        "-cgf",
                                        "--sitesOnly",
                                        "-o /Users/bimber/Desktop/removeAnn.vcf"
                                        //"-o %s"
                                ),
                                1,
                                Arrays.asList("7091cbeb47d041463806c8c8f98239a6")
                              );
        executeTest("testBasicOperation", spec);
    }
}
