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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.tools.walkers.annotator.AnnotationUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GenotypeGVCFsIntegrationTest extends WalkerTest {
    private static String baseTestString(String args, String ref) {
        return "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + ref + args;
    }

    private static String baseBPResolutionString(String args) {
        return baseTestString(" -V " + privateTestDir + "gvcf.basepairResolution.vcf " + args, b37KGReference);
    }

    @Test
    public void testUpdatePGT() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "testUpdatePGT.vcf", b37KGReference),
                1,
                Collections.singletonList("326ec5afa27ade4d0c562aa227997d88"));
        executeTest("testUpdatePGT", spec);
    }

    @Test
    public void testUpdatePGTStrandAlleleCountsBySample() throws IOException{
        final String logFileName = "testUpdatePGTStrandAlleleCountsBySample.log";
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "testUpdatePGT.vcf -A StrandAlleleCountsBySample -log " + logFileName, b37KGReference),
                1,
                Collections.singletonList("b995bf820c3edff5b721338ea9cf44aa"));
        executeTest("testUpdatePGTStrandAlleleCountsBySample", spec);

        final File file = new File(logFileName);
        Assert.assertTrue(FileUtils.readFileToString(file).contains(AnnotationUtils.ANNOTATION_HC_WARN_MSG));
    }

    @Test
    public void combineSingleSamplePipelineGVCF() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
                        " -L 20:10,000,000-11,000,000", b37KGReference),
                1,
                Collections.singletonList("7b2a135e694f9d1190e041e6fd420123"));
        executeTest("combineSingleSamplePipelineGVCF", spec);
    }

    @Test
    public void testTetraploidRun() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "tetraploid-gvcf-1.vcf" +
                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + privateTestDir + "tetraploid-gvcf-3.vcf" +
                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals", b37KGReference),
                1,
                Collections.singletonList("d6ef5e411ac5829a12d825a0fefac883"));
        executeTest("testTetraploidRun", spec);
    }

    @Test
    public void testMixedPloidyRun() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "haploid-gvcf-1.vcf" +
                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + privateTestDir + "diploid-gvcf-3.vcf" +
                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals", b37KGReference),
                1,
                Collections.singletonList("b497f16cd9eb99e353d9430fe7f34635"));
        executeTest("testMixedPloidyRun", spec);
    }

    @Test(enabled= true)
    public void testMixedPloidyMaxNumPLValuesRun() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "haploid-gvcf-1.vcf" +
                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + privateTestDir + "diploid-gvcf-3.vcf" +
                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals" +
                        " -maxNumPLValues 3", b37KGReference),
                1,
                Collections.singletonList("8c8ebe2069977ba13024a95827c6c50d"));
        executeTest("testMixedPloidyMaxNumPLValuesRun", spec);
    }

    @Test
    public void combineSingleSamplePipelineGVCF_includeNonVariants() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
                        " --includeNonVariantSites -L 20:10,030,000-10,033,000 -L 20:10,386,000-10,386,500", b37KGReference),
                1,
                Collections.singletonList("a9ecd152ec4b5b541887a0aed016f40d"));
        spec.disableShadowBCF();
        executeTest("combineSingleSamplePipelineGVCF_includeNonVariants", spec);
    }

    @Test
    public void combineSingleSamplePipelineGVCFHierarchical() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "combine.single.sample.pipeline.combined.vcf" +
                        " -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
                        " -L 20:10,000,000-20,000,000", b37KGReference),
                1,
                Collections.singletonList("8a37077d9a52b9d556cdd19403e27635"));
        executeTest("combineSingleSamplePipelineGVCFHierarchical", spec);
    }

    @Test
    public void combineSingleSamplePipelineGVCF_addDbsnp() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
                        " -L 20:10,000,000-11,000,000 --dbsnp " + b37dbSNP132, b37KGReference),
                1,
                Collections.singletonList("181fcb5d240b9bd92e3c793ca5aa7954"));
        executeTest("combineSingleSamplePipelineGVCF_addDbsnp", spec);
    }

    @Test
    public void combineBPresGVCFs() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "NA12891.BPres.g.vcf" +
                " -V " + privateTestDir + "NA12892.BPres.g.vcf" +
                " -L 20:10433000-10436909", b37KGReference),
                1,
                Collections.singletonList("f342872f485e6978501facc78c354078"));
        executeTest("combineBPresGVCFs", spec);
    }

    @Test
    public void testJustOneSample() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -L 1:69485-69791 -o %s -R " + b37KGReference +
                " -V " + privateTestDir + "gvcfExample1.vcf",
                1,
                Collections.singletonList("df88bbf2eea39a06f2bcc47d9379e5fa"));
        executeTest("testJustOneSample", spec);
    }

    @Test
    public void testSamplesWithDifferentLs() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -L 1:69485-69791 -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "gvcfExample1.vcf" +
                        " -V " + privateTestDir + "gvcfExample2.vcf",
                1,
                Collections.singletonList("933c3ec48870c54f7f74b259272d6645"));
        executeTest("testSamplesWithDifferentLs", spec);
    }

    @Test
    public void testNoPLsException() {
        // Test with input files with (1) 0/0 and (2) ./.
        final String md5 = "91038469a8133feb05038528f8565840";
        final WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -L 1:1115550-1115551 -o %s -R " + hg19Reference +
                        " --variant " + privateTestDir + "combined_genotype_gvcf_exception.vcf",
                1,
                Collections.singletonList(md5));
        final WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -L 1:1115550-1115551 -o %s -R " + hg19Reference +
                        " --variant " + privateTestDir + "combined_genotype_gvcf_exception.nocall.vcf",
                1,
                Collections.singletonList(md5));
        executeTest("testNoPLsException.1", spec1);
        executeTest("testNoPLsException.2", spec2);
    }

    @Test
    public void testNDA() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseBPResolutionString("-nda"),
                1,
                Collections.singletonList("b92742a67c5b1718fb5674393b2008fd"));
        executeTest("testNDA", spec);
    }

    @Test
    public void testAllSitesNonBiallelic() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseBPResolutionString("-allSites"),
                1,
                Collections.singletonList("77ae407c276aa7070fa22b110160c684"));
        spec.disableShadowBCF();
        executeTest("testAllSitesNonBiallelic", spec);
    }

    @Test
    public void testMaxAltAlleles() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseBPResolutionString("-maxAltAlleles 1"),
                1,
                Collections.singletonList("312f01d41dff851f449d9b4ecf0e78d2"));
        executeTest("testMaxAltAlleles", spec);
    }

    @Test
    public void testStandardConf() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseBPResolutionString("-stand_call_conf 300"),
                1,
                Collections.singletonList("89582b2a32e3f04d3daa565a0c185003"));
        executeTest("testStandardConf", spec);
    }

    @Test
    public void testStrandAlleleCountsBySample() throws IOException {
        //HaplotypeCaller creates gVCF
        final String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
        final WalkerTestSpec specHaplotypeCaller = new WalkerTestSpec(
                "-T HaplotypeCaller --disableDithering " +
                        String.format("-R %s -I %s ", b37KGReference, CEUTRIO_BAM) +
                        "--no_cmdline_in_header -o %s -L 20:10130000-10134800 " +
                        "-ERC GVCF --sample_name NA12878 -variant_index_type LINEAR " +
                        "-variant_index_parameter 128000 -A StrandAlleleCountsBySample",
                1, Collections.singletonList("")
        );
        specHaplotypeCaller.disableShadowBCF(); //TODO: Remove when BaseTest.assertAttributesEquals() works with SAC
        final File gVCF = executeTest("testStrandAlleleCountsBySampleHaplotypeCaller", specHaplotypeCaller).getFirst().get(0);

        //Use gVCF from HaplotypeCaller
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + gVCF.getAbsolutePath(), b37KGReference),
                1,
                Collections.singletonList("f613d0bc0d45aafe53227bc0d13712f1"));
        spec.disableShadowBCF();  //TODO: Remove when BaseTest.assertAttributesEquals() works with SAC
        executeTest("testStrandAlleleCountsBySample", spec);
    }

    @Test
    public void testUniquifiedSamples() throws IOException {
        //two copies of 5 samples; will also test InbreedingCoeff calculation for uniquified samples
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
                        " -V:sample1B " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
                        " -V:sample2B " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
                        " -V:combined1 " + privateTestDir + "combine.single.sample.pipeline.combined.vcf" +
                        " -V:combined2 " + privateTestDir + "combine.single.sample.pipeline.combined.vcf" +
                        " --uniquifySamples", b37KGReference),
                1,
                Collections.singletonList("0c99b1b20fb035a5dada036bd4cf39e5"));
        executeTest("testUniquifiedSamples", spec);

    }

    /**
     * Returns a list of attribute values from a VCF file
     *
     * @param vcfFile VCF file
     * @param attributeName attribute name
     *
     * @throws IOException if the file does not exist or can not be opened
     *
     * @return list of attribute values
     */
    private List<String> getAttributeValues(final File vcfFile, final String attributeName) throws IOException {
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(vcfFile);
        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIteratorVCF);

        final List<String> attributeValues = new ArrayList<String>();
        while (lineIteratorVCF.hasNext()) {
            final String line = lineIteratorVCF.next();
            Assert.assertFalse(line == null);
            final VariantContext vc = codec.decode(line);

            for (final Genotype g : vc.getGenotypes()) {
                if (g.hasExtendedAttribute(attributeName)) {
                    attributeValues.add((String) g.getExtendedAttribute(attributeName));
                }
            }
        }

        return attributeValues;
    }

    /**
     * Section to test spanning deletions
     */
    @Test
    public void testSpanningDeletions() throws IOException {
        final String gvcf1 = privateTestDir + "spanningDel.1.g.vcf";
        final String gvcf2 = privateTestDir + "spanningDel.2.g.vcf";
        final String gvcf3 = privateTestDir + "spanningDel.3.g.vcf";

        // create the genotyped VCF to use as a basis for comparison against all of the combined versions
        // case 0: GenotypeGVCFs(1.g.vcf, 2.g.vcf, 3.g.vcf)
        final WalkerTestSpec genotypeBase = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + gvcf1 + " -V " + gvcf2 + " -V " + gvcf3,
                1,
                Collections.singletonList(""));
        genotypeBase.disableShadowBCF();
        final File genotypeBaseVCF = executeTest("genotypeBase", genotypeBase).getFirst().get(0);
        final List<VariantContext> BASE_VARIANT_CONTEXTS = getVariantContexts(genotypeBaseVCF);

        // case 1: GenotypeGVCFs(CombineGVCFs(1.g.vcf, 2.g.vcf), 3.g.vcf)
        final WalkerTestSpec combine12 = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + gvcf1 + " -V " + gvcf2,
                1,
                Collections.singletonList(""));
        combine12.disableShadowBCF();
        final File combined_gVCF12 = executeTest("combine12", combine12).getFirst().get(0);
        final WalkerTestSpec genotype12_3 = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + combined_gVCF12.getAbsolutePath() + " -V " + gvcf3,
                1,
                Collections.singletonList(""));
        genotype12_3.disableShadowBCF();
        final File genotype12_3VCF = executeTest("genotype12_3", genotype12_3).getFirst().get(0);
        final List<VariantContext> VARIANT_CONTEXTS_12_3 = getVariantContexts(genotype12_3VCF);
        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_12_3);

        // case 2: GenotypeGVCFs(CombineGVCFs(CombineGVCFs(1.g.vcf, 2.g.vcf), 3.g.vcf))
        final WalkerTestSpec combine12then3 = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + combined_gVCF12 + " -V " + gvcf3,
                1,
                Collections.singletonList(""));
        combine12then3.disableShadowBCF();
        final File combined_gVCF12then3 = executeTest("combined_gVCF12then3", combine12then3).getFirst().get(0);
        final WalkerTestSpec genotype12then3 = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + combined_gVCF12then3.getAbsolutePath(),
                1,
                Collections.singletonList(""));
        genotype12then3.disableShadowBCF();
        final File genotype12then3VCF = executeTest("genotype12then3", genotype12then3).getFirst().get(0);
        final List<VariantContext> VARIANT_CONTEXTS_12then3 = getVariantContexts(genotype12then3VCF);
        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_12then3);

        // case 3: GenotypeGVCFs(CombineGVCFs(CombineGVCFs(1.g.vcf, 3.g.vcf), 2.g.vcf))
        final WalkerTestSpec combine13 = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + gvcf1 + " -V " + gvcf3,
                1,
                Collections.singletonList(""));
        combine13.disableShadowBCF();
        final File combined_gVCF13 = executeTest("combine13", combine13).getFirst().get(0);
        final WalkerTestSpec combine13then2 = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + combined_gVCF13 + " -V " + gvcf2,
                1,
                Collections.singletonList(""));
        combine13then2.disableShadowBCF();
        final File combined_gVCF13then2 = executeTest("combined_gVCF13then2", combine13then2).getFirst().get(0);
        final WalkerTestSpec genotype13then2 = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + combined_gVCF13then2.getAbsolutePath(),
                1,
                Collections.singletonList(""));
        genotype13then2.disableShadowBCF();
        final File genotype13then2VCF = executeTest("genotype13then2", genotype13then2).getFirst().get(0);
        final List<VariantContext> VARIANT_CONTEXTS_13then2 = getVariantContexts(genotype13then2VCF);
        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_13then2);

        // case 4: GenotypeGVCFs(CombineGVCFs(1.g.vcf, 2.g.vcf, 3.g.vcf))
        final WalkerTestSpec combine123 = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + gvcf1 + " -V " + gvcf2 + " -V " + gvcf3,
                1,
                Collections.singletonList(""));
        combine123.disableShadowBCF();
        final File combined_gVCF123 = executeTest("combine123", combine123).getFirst().get(0);
        final WalkerTestSpec genotype123 = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + combined_gVCF123.getAbsolutePath(),
                1,
                Collections.singletonList(""));
        genotype123.disableShadowBCF();
        final File genotype123VCF = executeTest("genotype123", genotype123).getFirst().get(0);
        final List<VariantContext> VARIANT_CONTEXTS_123 = getVariantContexts(genotype123VCF);
        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_123);
    }

    /**
     * Returns a list of VariantContext records from a VCF file
     *
     * @param vcfFile VCF file
     *
     * @throws IOException if the file does not exist or can not be opened
     *
     * @return list of VariantContext records
     */
    private static List<VariantContext> getVariantContexts(final File vcfFile) throws IOException {
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(vcfFile);
        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIteratorVCF);

        final List<VariantContext> VCs = new ArrayList<>();
        while ( lineIteratorVCF.hasNext() ) {
            final String line = lineIteratorVCF.next();
            Assert.assertFalse(line == null);
            VCs.add(codec.decode(line));
        }

        return VCs;
    }

    private static void testVCsAreEqual(final List<VariantContext> VCs1, final List<VariantContext> VCs2) {
        Assert.assertEquals(VCs1.size(), VCs2.size(), "number of Variant Contexts");
        for ( int i = 0; i < VCs1.size(); i++ ) {
            final VariantContext vc1 = VCs1.get(i);
            final VariantContext vc2 = VCs2.get(i);
            Assert.assertEquals(vc1.toStringDecodeGenotypes(), vc2.toStringDecodeGenotypes());
        }
    }


    private static final String simpleSpanningDeletionsMD5 = "53f2b8991e49a47efc44b8e02ebb8d91";

    @Test
    public void testSpanningDeletionsMD5() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.1.g.vcf -V " + privateTestDir + "spanningDel.2.g.vcf",
                1,
                Collections.singletonList(simpleSpanningDeletionsMD5));
        spec.disableShadowBCF();
        executeTest("testSpanningDeletionsMD5", spec);
    }

    @Test
    public void testSpanningDeletionsFromCombinedGVCF() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.combined.g.vcf",
                1,
                Collections.singletonList(simpleSpanningDeletionsMD5));
        spec.disableShadowBCF();
        executeTest("testSpanningDeletionsFromCombinedGVCFMD5", spec);
    }

    @Test
    public void testMultipleSpanningDeletionsMD5() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.1.g.vcf -V " + privateTestDir + "spanningDel.2.g.vcf -V " + privateTestDir + "spanningDel.3.g.vcf",
                1,
                Collections.singletonList("907dfaa4d31c22705eadd5890ae23929"));
        spec.disableShadowBCF();
        executeTest("testMultipleSpanningDeletionsMD5", spec);
    }

    @Test
    public void testSpanningDeletionDoesNotGetGenotypedWithNoOtherAlleles() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.delOnly.g.vcf",
                1,
                Collections.singletonList("b923e5c6d5dbce62034178bd5234b932"));
        spec.disableShadowBCF();
        executeTest("testSpanningDeletionDoesNotGetGenotypedWithNoOtherAlleles", spec);
    }

    @Test
    public void testDeprecatedSpanningDeletionDoesNotGetGenotypedWithNoOtherAlleles() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.depr.delOnly.g.vcf",
                1,
                Collections.singletonList("01ae75dfe5c0c2350fcef0f4cdca36b2"));
        spec.disableShadowBCF();
        executeTest("testSpanningDeletionDoesNotGetGenotypedWithNoOtherAlleles", spec);
    }

    @Test
    public void testGenotypingSpanningDeletionOverSpan() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.delOverSpan.1.g.vcf -V " +
                        privateTestDir + "spanningDel.delOverSpan.2.g.vcf",
                0,
                Collections.singletonList(""));     // we do not care about the md5; we just want to make sure it doesn't blow up with an error
        spec.disableShadowBCF();
        executeTest("testGenotypingSpanningDeletionOverSpan", spec);
    }

    @Test
    public void testBadADPropagationHaploidBugTest() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "ad-bug-input.vcf",
                1,
                Collections.singletonList("4d6cbd8d666a43fc136d73de2b217719"));
        spec.disableShadowBCF();
        executeTest("testBadADPropagationHaploidBugTest", spec);
    }

    @Test
    public void testSAC() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "261_S01_raw_variants_gvcf.vcf",
                1,
                Collections.singletonList("ea96440b537dd1b2b25ea565dfaa71fc"));
        spec.disableShadowBCF();
        executeTest("testSAC", spec);
    }

    @Test
    public void testSACMultisampleTetraploid() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "tetraploid-multisample-sac.g.vcf",
                1,
                Collections.singletonList("c21c847ef794c11e249985a16893b2fa"));
        spec.disableShadowBCF();
        executeTest("testSACMultisampleTetraploid", spec);
    }

    @Test
    public void testSetZeroRGQsToNoCall() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "set.zero.RGQs.no.call.sample1.g.vcf" +
                        " -V " + privateTestDir + "set.zero.RGQs.no.call.sample2.g.vcf" +
                        " -L chr16:1279274-1279874 -allSites", hg19ReferenceWithChrPrefixInChromosomeNames),
                Collections.singletonList("903047b6262fcb82070556ff74f26a75"));
        spec.disableShadowBCF();
        executeTest("testSetZeroRGQsToNoCall", spec);
    }

    @Test
    public void testAlleleSpecificAnnotations() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard --disableDithering -V "
                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf -V " + privateTestDir + "NA12891.AS.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("29d6db0a93abd72d64fb1e82da65c715"));
        spec.disableShadowBCF();
        executeTest("testAlleleSpecificAnnotations", spec);
    }

    @Test
    public void testASMateRankSumAnnotation() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -A AS_MQMateRankSumTest --disableDithering -V "
                + privateTestDir + "NA12878.AS.MateRankSum.chr20snippet.g.vcf -V " + privateTestDir + "NA12891.AS.MateRankSum.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("2a330015a7db9f9aee9bc5b776698f73"));
        spec.disableShadowBCF();
        executeTest("testASMateRankSumAnnotation", spec);
    }

    @Test
    public void testASInsertSizeRankSumAnnotation() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard --disableDithering -V "
                + privateTestDir + "NA12878.AS.InsertSizeRankSum.chr20snippet.g.vcf -V " + privateTestDir + "NA12891.AS.InsertSizeRankSum.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("75aee1e0c8c3528180e344ec6c0d8ffd"));
        spec.disableShadowBCF();
        executeTest("testASInsertSizeRankSumAnnotation", spec);
    }


    @Test
    //make sure none of the assumptions about things being merged as lists break the single-sample case
    //This test file also doesn't have raw data, so test to make sure that doesn't make GenotypeGVCFs crash and burn
    //Note that AS_InbreedingCoeff and InbreedingCoeff may still differ for bialleleic sites for low number of samples because allele frequencies are derived differently
    public void testAlleleSpecificAnnotations_oneSample() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard --disableDithering -V "
                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("f4fa3acec2b21037368898e913b7a3fa"));
        spec.disableShadowBCF();
        executeTest("testAlleleSpecificAnnotations_oneSample", spec);
    }

    @Test
    //do at least 10 samples so InbreedingCoeff and AS_InbreedingCoeff are output
    public void testAlleleSpecificAnnotations_elevenSamples() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard --disableDithering -V "
                + privateTestDir + "multiSamples.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("4e90f6908248fac9b3ce3e545180a8e5"));
        spec.disableShadowBCF();
        executeTest("testAlleleSpecificAnnotations_elevenSamples", spec);
    }

    @Test
    public void testMonomorphicVCwithAlt() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -G AS_Standard -o %s --no_cmdline_in_header --disableDithering -V "
                + privateTestDir + "monomorphicGVCwithAlt.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("43953b3e75a4d470b65773b1b5bea066"));
        spec.disableShadowBCF();
        executeTest("testAlleleSpecificAnnotations", spec);
    }

    @Test
    public void testFractionInformativeReads() {
        final String cmd = "-T GenotypeGVCFs -R " + b37KGReference + " -G AS_Standard -o %s --no_cmdline_in_header -A FractionInformativeReads --disableDithering -V "
                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf -V " + privateTestDir + "NA12891.AS.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Collections.singletonList("0b1bbcc7d24f8b0945c97907b1cdd974"));
        spec.disableShadowBCF();
        executeTest("testAlleleSpecificAnnotations", spec);
    }

    //Regression test for https://github.com/broadinstitute/gsa-unstable/issues/1281
    @Test
    public void testGenotypingSpanningDeletionWithAllSites() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "spanningDel.genotyping.g.vcf -allSites", b37KGReference),
                Collections.singletonList("04cfe93e92444cbde80e13ca8b8c3913"));
        spec.disableShadowBCF();
        executeTest("testGenotypingSpanningDeletionWithAllSites", spec);
    }

    @Test
    public void testGenotypingSpanningDeletionAcrossLines() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "input-1_2256566.vcf", b37KGReference),
                Collections.singletonList("152c8e07e35c592868f43626f27365de"));
        spec.disableShadowBCF();
        executeTest("testGenotypingSpanningDeletionAcrossLines", spec);
    }

    @Test
    public void testNewQualNaNBugFix() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -newQual -V " + privateTestDir + "input-newqual-nan-bug-fix.vcf", b37KGReferenceWithDecoy),
                Collections.singletonList("ab5994dcaf9b2d41269b4ff4729b5e81"));
        spec.disableShadowBCF();
        executeTest("testNewQualNaNBugFix", spec);
    }

    @Test
    public void testHomRefHighMQ() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "NA18503.22.vcf -V " + privateTestDir + "NA18504.22.vcf -V " +
                        privateTestDir + "NA18505.22.vcf -allSites", b37KGReference),
                Collections.singletonList("6d253024246e1024b9b6e8f885f53799"));
        spec.disableShadowBCF();
        executeTest("testHomRefHighMQ", spec);
    }

    @Test
    public void testInputPrior() {
        final WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -V " + privateTestDir + "gvcfExample1.vcf -V " + privateTestDir + "gvcfExample2.vcf " +
                        "-L 1:69485-69791 -inputPrior 0.2 -inputPrior 0.2 -inputPrior 0.2 -inputPrior 0.2", b37KGReference),
                Collections.singletonList("860d133262160bbc75ce5849e5fa491f"));
        spec.disableShadowBCF();
        executeTest("testInputPrior", spec);
    }
}