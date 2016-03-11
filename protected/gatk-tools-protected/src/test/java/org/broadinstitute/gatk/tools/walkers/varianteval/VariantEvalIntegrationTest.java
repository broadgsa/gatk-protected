/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.varianteval;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VariantEvalIntegrationTest extends WalkerTest {
    private static String variantEvalTestDataRoot = privateTestDir + "VariantEval/";
    private static String fundamentalTestVCF = variantEvalTestDataRoot + "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";
    private static String fundamentalTestSNPsVCF = variantEvalTestDataRoot + "FundamentalsTest.annotated.db.subset.final.vcf";
    private static String fundamentalTestSNPsWithMLEVCF = variantEvalTestDataRoot + "FundamentalsTest.annotated.db.subset.final.withMLE.vcf";
    private static String fundamentalTestSNPsSplit1of2VCF = variantEvalTestDataRoot + "FundamentalsTest.annotated.db.subset.final.split_1_of_2.vcf";
    private static String fundamentalTestSNPsSplit2of2VCF = variantEvalTestDataRoot + "FundamentalsTest.annotated.db.subset.final.split_2_of_2.vcf";
    private static String fundamentalTestSNPsOneSampleVCF = variantEvalTestDataRoot + "FundamentalsTest.annotated.db.subset.final.NA12045.vcf";

    private static String cmdRoot = "-T VariantEval" +
            " -R " + b36KGReference;

    @Test
    public void testFunctionClassWithSnpeff() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + validationDataLocation + "snpEff2.0.5.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-noEV",
                                        "-EV TiTvVariantEvaluator",
                                        "-noST",
                                        "-ST FunctionalClass",
                                        "-L " + validationDataLocation + "snpEff2.0.5.AFR.unfiltered.VariantAnnotator.output.vcf",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("7091cbeb47d041463806c8c8f98239a6")
                              );
        executeTest("testFunctionClassWithSnpeff", spec);
    }

    @Test
    public void testStratifySamplesAndExcludeMonomorphicSites() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + variantEvalTestDataRoot + "/CEU.trio.callsForVE.vcf",
                                        "-noEV",
                                        "-EV TiTvVariantEvaluator",
                                        "-ST Sample",
                                        "-L " + variantEvalTestDataRoot + "/CEU.trio.callsForVE.vcf",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("7a09b8a6759ccee5da55f1f85a43fe9c")
                              );
        executeTest("testStratifySamplesAndExcludeMonomorphicSites", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndels() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "--dbsnp " + b37dbSNP132,
                                        "--eval " + fundamentalTestVCF,
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-L " + fundamentalTestVCF,
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("f70da7be5d4d8305f3e4433c9004aee4")
                              );
        executeTest("testFundamentalsCountVariantsSNPsandIndels", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNovelty() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("40abbc9be663aed8ee1158f832463ca8")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNovelty", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Novelty",
                        "-ST Filter",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("106a0e8753e839c0a2c030eb4b165fa9")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithNoveltyAndFilter", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithCpG() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST CpG",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("bca988c81a761f12627610e5a3bab5a0")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithCpG", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClasses() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST FunctionalClass",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("7ca5c0c5e79ba6cd1e5102ced851a1b4")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithFunctionalClass", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Degeneracy",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("a6a31f658ad1e76c79190ada758f157c")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithDegeneracy", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithSample() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST Sample",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("c1a3df6f89f5ddf7b7c296eb944f3fdd")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithSample", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthSelect",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("48652b360ce031aa2f9004c9bae6bda5")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithJexlExpression", spec);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST JexlExpression",
                        "-select 'DP < 20'",
                        "-selectName DepthLt20",
                        "-select 'DP > 20'",
                        "-selectName DepthGt20",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("c3521b18388aff7f53691a63619b3b07")
        );
        executeTest("testFundamentalsCountVariantsSNPsandIndelsWithMultipleJexlExpressions", spec);
    }

    @Test
    public void testFundamentalsCountVariantsNoCompRod() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--eval " + fundamentalTestVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-L " + fundamentalTestVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("90a46e045f3fe8b22f102acaaeec0201")
        );
        executeTest("testFundamentalsCountVariantsNoCompRod", spec);
    }

    @Test
    public void testSelect1() {
        String extraArgs = "-L 1:1-10,000,000";
        String tests = cmdRoot +
                " --dbsnp " + b36dbSNP129 +
                " --eval " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
                " --comp:comp_genotypes " + privateTestDir + "yri.trio.gatk.ug.head.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + extraArgs + " -ST CpG -o %s",
                1, Arrays.asList("eaa3708d9db22fca0844a652bb73b82f"));
        executeTestParallel("testSelect1", spec);
    }

    @Test
    public void testVEMendelianViolationEvaluator() {
        String vcfFile = "/MendelianViolationEval.vcf";
        String pedFile = "/MendelianViolationEval.ped";

        WalkerTestSpec spec = new WalkerTestSpec("-T VariantEval -R "+b37KGReference+" --eval " + variantEvalTestDataRoot + vcfFile + " -ped "+ variantEvalTestDataRoot + pedFile +" -noEV -EV MendelianViolationEvaluator -L 1:10109-10315 -o %s -mvq 0 -noST",
                1,
                Arrays.asList("c56e19d0647d826485d8a3b559d5c56d"));
        executeTestParallel("testVEMendelianViolationEvaluator" + vcfFile, spec);
    }

    @Test
    public void testMVEvalFamilyStrat() {
        String vcfFile = "/PhaseByTransmission/PhaseByTransmission.IntegrationTest.TP.vcf";
        String pedFile = "/PhaseByTransmission/PhaseByTransmission.IntegrationTest.goodFamilies.ped";

        WalkerTestSpec spec = new WalkerTestSpec("-R "+b37KGReference+ " -T VariantEval -ped " + privateTestDir + pedFile + " -eval " + privateTestDir + vcfFile + " -noEV -noST -ST Family -EV MendelianViolationEvaluator -o %s",
                1,
                Arrays.asList("d599d3e6b308ac06b2c2e003cf596328"));
        executeTestParallel("testMVEvalFamilyStrat", spec);
    }


    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -selectName %s", cmd, select, name);
    }

    @Test(enabled = false) // no longer supported in the GATK
    public void testTranches() {
        String extraArgs = "-T VariantEval -R "+ hg18Reference +" --eval " + validationDataLocation + "GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.vcf -o %s -EV TiTvVariantEvaluator -L chr1 -noEV -ST CpG -tf " + privateTestDir + "tranches.6.txt";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("6af2b9959aa1778a5b712536de453952"));
        executeTestParallel("testTranches",spec);
    }

    @Test
    public void testCompOverlap() {
        String extraArgs = "-T VariantEval -R " + b37KGReference + " -L " + variantEvalTestDataRoot + "pacbio.hg19.intervals --comp:comphapmap " + comparisonDataLocation + "Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf --eval " + variantEvalTestDataRoot + "pacbio.ts.recalibrated.vcf -noEV -EV CompOverlap -sn NA12878 -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("d0d9208060e69e157dac3bf01bdd83b0"));
        executeTestParallel("testCompOverlap",spec);
    }

    @Test
    public void testEvalTrackWithoutGenotypes() {
        String extraArgs = "-T VariantEval -R " +
                           b37KGReference +
                           " -L 20" +
                           " --dbsnp " + b37dbSNP132 +
                           " --eval:evalBI " + variantEvalTestDataRoot + "ALL.20100201.chr20.bi.sites.vcf" +
                           " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("fe9dcf4933a645f55be1cb0e33497e49"));
        executeTestParallel("testEvalTrackWithoutGenotypes",spec);
    }

    @Test
    public void testEvalTrackWithoutGenotypesWithSampleFields() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + variantEvalTestDataRoot + "noGenotypes.vcf",
                        "-o %s"
                ),
                1,
                Arrays.asList("")); //There is no md5 because we only care that this completes without an exception.
        executeTest("testEvalTrackWithoutGenotypesWithSampleFields", spec);

    }

    @Test
    public void testMultipleEvalTracksWithoutGenotypes() {
        String extraArgs = "-T VariantEval -R " + b37KGReference +
                " -L 20" +
                " --dbsnp " + b37dbSNP132 +
                " --eval:evalBI " + variantEvalTestDataRoot + "ALL.20100201.chr20.bi.sites.vcf" +
                " --eval:evalBC " + variantEvalTestDataRoot + "ALL.20100201.chr20.bc.sites.vcf" +
                " -noST -ST Novelty -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("8dfdec264fcff9472bdee7d223fdb3ca"));
        executeTestParallel("testMultipleEvalTracksWithoutGenotypes",spec);
    }

    @Test
    public void testMultipleCompTracks() {
        String dbsnp = GATKDataLocation + "dbsnp_132_b37.vcf";

        String extraArgs =  "-T VariantEval" +
                           " -R " + b37KGReference +
                           " --comp " + variantEvalTestDataRoot + "ALL.phase1.chr20.broad.snps.genotypes.subset.vcf" +
                           " --eval " + variantEvalTestDataRoot + "NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.subset.vcf" +
                           " --dbsnp " + dbsnp +
                           " -L 20:10000000-10100000" +
                           " -noST -noEV -ST Novelty -EV CompOverlap" +
                           " -o %s";

        WalkerTestSpec spec = new WalkerTestSpec(extraArgs,1,Arrays.asList("44146b8d4ddbaeb9409c9b88eefe7f40"));
        executeTestParallel("testMultipleCompTracks",spec);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults1() {
        String md5 = "5f894d726cfaa0b29d7c11ff5bb9b3fd";

        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-sn NA12045",
                        "-noST",
                        "-L " + fundamentalTestSNPsVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList(md5)
        );
        executeTestParallel("testPerSampleAndSubsettedSampleHaveSameResults-subset", spec);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsOneSampleVCF,
                        "-noEV",
                        "-EV CompOverlap",
                        "-noST",
                        "-L " + fundamentalTestSNPsOneSampleVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList(md5)
        );
        executeTestParallel("testPerSampleAndSubsettedSampleHaveSameResults-onesample", spec2);
    }


    @Test
    public void testAlleleCountStrat() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST AlleleCount",
                        "-L " + fundamentalTestSNPsVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("0d51321693d4afc262e4059353993d12")
        );
        executeTest("testAlleleCountStrat", spec);
    }

    @Test
    public void testAlleleCountStratWithMLE() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsWithMLEVCF,
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST AlleleCount",
                        "-L " + fundamentalTestSNPsWithMLEVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("cadd4e32ab82f7e43d0e510dfe5e6dda")
        );
        executeTest("testAlleleCountStratWithMLE", spec);
    }

    @Test
    public void testMultipleEvalTracksAlleleCountWithMerge() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsSplit1of2VCF,
                        "--eval " + fundamentalTestSNPsSplit2of2VCF,
                        "--mergeEvals",
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST AlleleCount",
                        "-L " + fundamentalTestSNPsVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("0d51321693d4afc262e4059353993d12")
        );
        executeTest("testMultipleEvalTracksAlleleCountWithMerge", spec);
    }

    @Test
    public void testMultipleEvalTracksAlleleCountWithoutMerge() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "--dbsnp " + b37dbSNP132,
                        "--eval " + fundamentalTestSNPsSplit1of2VCF,
                        "--eval " + fundamentalTestSNPsSplit2of2VCF,
                        //"--mergeEvals", No merge with AC strat ==> error
                        "-noEV",
                        "-EV CountVariants",
                        "-noST",
                        "-ST AlleleCount",
                        "-L " + fundamentalTestSNPsVCF
                ),
                0,
                UserException.class
        );
        executeTest("testMultipleEvalTracksAlleleCountWithoutMerge", spec);
    }

    @Test
    public void testIntervalStrat() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-eval " + privateTestDir + "/withSymbolic.b37.vcf",
                                        "-noEV",
                                        "-EV CountVariants",
                                        "-noST",
                                        "-stratIntervals " + privateTestDir + "/overlapTest.bed",
                                        "-ST IntervalStratification",
                                        "-L 20",
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("74fc726760c0fcfe50c44c853756f261")
                              );
        executeTest("testIntervalStrat", spec);
    }

    @Test
    public void testModernVCFWithLargeIndels() {
        WalkerTestSpec spec = new WalkerTestSpec(
                                buildCommandLine(
                                        "-T VariantEval",
                                        "-R " + b37KGReference,
                                        "-eval " + validationDataLocation + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.vcf",
                                        "-L 20",
                                        "-D " + b37dbSNP132,
                                        "-o %s"
                                ),
                                1,
                                Arrays.asList("af317f1ea1b80e5d4bc4f2d8523ef73d")
                              );
        executeTest("testModernVCFWithLargeIndels", spec);
    }

    @Test
    public void testStandardIndelEval() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + validationDataLocation + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.vcf",
                        "-L 20",
                        "-noST -ST Sample -ST OneBPIndel -ST TandemRepeat",
                        "-noEV -EV IndelSummary -EV IndelLengthHistogram",
                        "-gold " + validationDataLocation + "/Mills_and_1000G_gold_standard.indels.b37.sites.vcf",
                        "-D " + b37dbSNP132,
                        "-o %s"
                ),
                1,
                Arrays.asList("7dc2d8983cb7d98b291ca2f60a9151b2")
        );
        executeTest("testStandardIndelEval", spec);
    }

    @Test
    public void testBadACValue() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + privateTestDir + "vcfexample.withBadAC.vcf",
                        "-noST -ST AlleleCount",
                        "-noEV -EV VariantSummary"
                ),
                0,
                UserException.class);
        executeTest("testBadACValue", spec);
    }


    @Test()
    public void testIncompatibleEvalAndStrat() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + validationDataLocation + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.vcf",
                        "-L 20 -noST -ST AlleleCount -noEV -EV VariantSummary"
                ),
                0,
                UserException.class);
        executeTest("testIncompatibleEvalAndStrat", spec);
    }

    public void testIncludingAC0(boolean includeAC0, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + privateTestDir + "/ac0.vcf",
                        "-L 20:81006 -noST -noEV -EV VariantSummary -o %s" + (includeAC0 ? " -keepAC0" : "")
                ),
                1,
                Arrays.asList(md5));
        executeTest("testIncludingAC0 keep ac 0 = " + includeAC0, spec);
    }

    @Test public void testWithAC0() { testIncludingAC0(true, "c786128cfe4d3e28cdbc15c5c838ad20"); }
    @Test public void testWithoutAC0() { testIncludingAC0(false, "7bc505c07d9aee49571ad4b3fc9f7feb"); }

    //
    // Test validation report is doing the right thing with sites only and genotypes files
    // where the validation comp has more genotypes than eval
    //
    @Test(dataProvider = "testValidationReportData")
    public void testValidationReport(final String name, final String eval, final String comp, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + eval,
                        "-comp " + comp,
                        "-L 20:10,000,000-10,000,010 -noST -noEV -EV ValidationReport -o %s"
                ),
                1,
                Arrays.asList(md5));
        executeTest("testValidationReport with " + name, spec);
    }

    @DataProvider(name = "testValidationReportData")
    public Object[][] testValidationReportData() {
        final String compGenotypes = privateTestDir + "/validationReportComp.vcf";
        final String compSites = privateTestDir + "/validationReportComp.noGenotypes.vcf";
        final String evalGenotypes = privateTestDir + "/validationReportEval.vcf";
        final String evalSites = privateTestDir + "/validationReportEval.noGenotypes.vcf";

        List<Object[]> tests = new ArrayList<Object[]>();
        tests.add(new Object[]{"sites/sites", evalSites, compSites, "0b32e19efce28087cdc7b58e17ed633a"});
        tests.add(new Object[]{"sites/genotypes", evalSites, compGenotypes, "e2ffecee4a3acd0da7dd7fe10a59b2bc"});
        tests.add(new Object[]{"genotypes/sites", evalGenotypes, compSites, "f0dbb848a94b451e42765b0cb9d09ee2"});
        tests.add(new Object[]{"genotypes/genotypes", evalGenotypes, compGenotypes, "73790b530595fcbd467a88475ea9717f"});
        return tests.toArray(new Object[][]{});
    }


    @Test
    public void testPrintMissingComp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T VariantEval",
                        "-R " + b37KGReference,
                        "-eval " + privateTestDir + "validationReportEval.noGenotypes.vcf",
                        "--comp " + privateTestDir + "validationReportComp.noGenotypes.vcf",
                        "-L 20",
                        "-EV PrintMissingComp"
                ),
                0,
                Arrays.asList("d41d8cd98f00b204e9800998ecf8427e")); // sato: make sure it doesn't throw a null pointer exception.
        executeTest("testPrintMissingComp", spec);

    }
}
