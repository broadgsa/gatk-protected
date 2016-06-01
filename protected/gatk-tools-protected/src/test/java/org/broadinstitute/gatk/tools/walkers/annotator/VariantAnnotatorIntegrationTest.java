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

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {

    final static String REF = b37KGReference;
    final static String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final static String STANDARD_ANNOTATIONS = " -G Standard -G StandardUG ";

    public static String baseTestString() {
        return "-T VariantAnnotator -R " + b36KGReference + " --no_cmdline_in_header -o %s";
    }

    @Test
    public void testHasAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b65bf866457f000926b76d0f9d40065e"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("8e830da0bf34f1dc91bbc2fa64b8a518"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b63baada372925a76c3f279e16eb631d"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("6f5856bc2d31f8aae4131717e5ab0b16"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("58793dec36f8aec2cd8894898ece7c4e"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        // the genotype annotations in this file are actually out of order.  If you don't parse the genotypes
        // they don't get reordered.  It's a good test of the genotype ordering system.
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("27745920cc780b04e8f5acba79f868ca"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7ec5470f742f80cdfbfb203213bea8cc"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("afc47e4f253d0999961f26920be8e834"));

        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testExcludeAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "-XA FisherStrand -XA ReadPosRankSumTest --variant " + privateTestDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("a8a87e2a67436e14ed32ce9d355a3440"));
        executeTest("test exclude annotations", spec);
    }

    @Test
    public void testAskingStrandAlleleCountsBySample() throws IOException{
        File logFile = createTempFile("testAskingStrandAlleleCountsBySample.log", ".tmp");
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000 "
                        + "-A StrandAlleleCountsBySample -log " + logFile.getAbsolutePath(), 1,
                Arrays.asList("5c0fe344544a887acbb5cd83083d303b"));
        executeTest("test file has annotations, adding StrandAlleleCountsBySample annotation", spec);

        Assert.assertTrue(FileUtils.readFileToString(logFile).contains(AnnotationUtils.ANNOTATION_HC_WARN_MSG));
    }

    @Test
    public void testAskingGCContent() throws IOException{
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000 -A GCContent", 1,
                Arrays.asList("82238b65dbd085baaec68be2975a9bf8"));
        final File outputVCF = executeTest("test file has annotations, adding GCContent annotation", spec).getFirst().get(0);
        final VCFCodec codec = new VCFCodec();
        final VCFHeader header = (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new FileInputStream(outputVCF)));
        final VCFHeaderLine infoLineGC = header.getInfoHeaderLine(GATKVCFConstants.GC_CONTENT_KEY);
        // GC content must be a Float type
        Assert.assertTrue(infoLineGC.toString().contains("Type=Float"));
    }

    @Test
    public void testOverwritingHeader() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample4.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,001,292", 1,
                Arrays.asList("fc7958261af93681fde73c1fc6b578a0"));
        executeTest("test overwriting header", spec);
    }

    @Test
    public void testNoReads() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("2d535a48ec1c66aa9c707f6d498fc81d"));
        executeTest("not passing it any reads", spec);
    }

    @Test
    public void testDBTagWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --dbsnp " + b36dbSNP129 + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("a2a2f9ce8e6f9c933ad46906719ce402"));
        executeTest("getting DB tag with dbSNP", spec);
    }

    @Test
    public void testMultipleIdsWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --alwaysAppendDbsnpId --dbsnp " + b36dbSNP129 + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3withIDs.vcf -L " + privateTestDir + "vcfexample3withIDs.vcf", 1,
                Arrays.asList("ed351765d63d92b1913784fa47b3d859"));
        executeTest("adding multiple IDs with dbSNP", spec);
    }

    @Test
    public void testDBTagWithHapMap() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --comp:H3 " + privateTestDir + "fakeHM3.vcf" + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("114d8300ec0fac613bb2e82f8951adc0"));
        executeTest("getting DB tag with HM3", spec);
    }

    @Test
    public void testDBTagWithTwoComps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --comp:H3 " + privateTestDir + "fakeHM3.vcf --comp:foo " + privateTestDir + "fakeHM3.vcf " + STANDARD_ANNOTATIONS + " --variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("af59185b6a03d4147d2755019dcc6bf9"));
        executeTest("getting DB tag with 2 comps", spec);
    }

    @Test
    public void testNoQuals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "noQual.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L " + privateTestDir + "noQual.vcf -A QualByDepth", 1,
                Arrays.asList("b6321f3ce7a60d083be64d5ec9a54c1b"));
        executeTest("test file doesn't have QUALs", spec);
    }

    @Test
    public void testUsingExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + privateTestDir + "targetAnnotations.vcf" + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty.vcf -E foo.AF -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("7bcbd8ad8388f371d1a990fde67d3273"));
        executeTest("using expression", spec);
    }

    @Test
    public void testUsingExpressionAlleleMisMatch() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resourceAlleleConcordance --resource:foo " + privateTestDir + "targetAnnotations.vcf" + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty-mod.vcf -E foo.AF -L " + privateTestDir + "vcfexample3empty-mod.vcf", 1,
                Arrays.asList("76f716569e33d88914e479b70e08ac88"));
        executeTest("using expression allele mismatch", spec);
    }

    @Test
    public void testUsingExpressionMultiAllele() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + privateTestDir + "targetAnnotations-multiAllele.vcf" + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty-multiAllele.vcf -E foo.AF -E foo.AC -L " + privateTestDir + "vcfexample3empty-multiAllele.vcf", 1,
                Arrays.asList("0e3fc86349f5fd28159d00d22d278e84"));
        executeTest("using expression with multi-alleles", spec);
    }

    @Test
    public void testFilterInExpression(){
        /* The order of filters in the output seems platform-dependent. May need to change htsjdk to make the order consistent across platforms. [Sato] */
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + privateTestDir + "annotationResourceWithFilter.vcf" + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty-multiAllele.vcf -E foo.FILTER -L " + privateTestDir + "vcfexample3empty-multiAllele.vcf", 1,
                Arrays.asList("6fe67e72232165a829fde7c3b12c2275"));
        executeTest("annotate a vcf with the FILTER field of another vcf", spec);
    }

    @Test
    public void testUsingExpressionWithID() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + privateTestDir + "targetAnnotations.vcf" + STANDARD_ANNOTATIONS + "--variant " + privateTestDir + "vcfexample3empty.vcf -E foo.ID -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("94b03ee63604ab8d61aacfd3297c5dca"));
        executeTest("using expression with ID", spec);
    }

    @Test
    public void testTabixAnnotationsAndParallelism() {
        final String MD5 = "c5beea399dadbba66a7fc46036eeafe5";
        for ( String file : Arrays.asList("CEU.exon.2010_03.sites.vcf", "CEU.exon.2010_03.sites.vcf.gz")) {
            WalkerTestSpec spec = new WalkerTestSpec(
                    baseTestString() + " -A HomopolymerRun --variant:vcf " + validationDataLocation + file + " -L " + validationDataLocation + "CEU.exon.2010_03.sites.vcf --no_cmdline_in_header", 1,
                    Arrays.asList(MD5));
            executeTest("Testing lookup vcf tabix vs. vcf tribble", spec);
        }

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -A HomopolymerRun -nt 2 --variant:vcf " + validationDataLocation + "CEU.exon.2010_03.sites.vcf -L " + validationDataLocation + "CEU.exon.2010_03.sites.vcf --no_cmdline_in_header", 1,
                Arrays.asList(MD5));

        executeTest("Testing lookup vcf tabix vs. vcf tribble plus parallelism", spec);
    }

    @Test
    public void testSnpEffAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + hg19Reference + " --no_cmdline_in_header -o %s -A SnpEff --variant " +
            validationDataLocation + "1kg_exomes_unfiltered.AFR.unfiltered.vcf --snpEffFile  " + validationDataLocation +
            "snpEff2.0.5.AFR.unfiltered.vcf -L 1:1-1,500,000 -L 2:232,325,429",
            1,
            Arrays.asList("db0c5f273583a54d0cefc4b3c01aae9a")
        );
        executeTest("Testing SnpEff annotations", spec);
    }

    @Test
    public void testSnpEffAnnotationsUnsupportedVersionGATKMode() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + b37KGReference + " --no_cmdline_in_header -o %s -A SnpEff " +
            "--variant " + privateTestDir + "vcf4.1.example.vcf " +
            "--snpEffFile  " + privateTestDir + "snpEff_unsupported_version_gatk_mode.vcf " +
            "-L 1:10001292-10012424",
            1,
            Arrays.asList("5f952cfcd25653edcffd3916d52e94ec")
        );
        executeTest("Testing SnpEff annotations (unsupported version, GATK mode)", spec);
    }

    @Test
    public void testSnpEffAnnotationsUnsupportedVersionNoGATKMode() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + b37KGReference + " --no_cmdline_in_header -o %s -A SnpEff " +
            "--variant " + privateTestDir + "vcf4.1.example.vcf " +
            "--snpEffFile  " + privateTestDir + "snpEff_unsupported_version_no_gatk_mode.vcf " +
            "-L 1:10001292-10012424",
            1,
            Arrays.asList("86fb28c2e3886eda194026b3b7d07c77")
        );
        executeTest("Testing SnpEff annotations (unsupported version, no GATK mode)", spec);
    }

    @Test(enabled = true)
    public void testTDTAnnotation() {
        final String MD5 = "cf65fdfbcd822279e84326989b2f6378";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A TransmissionDisequilibriumTest --variant:vcf " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + privateTestDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing TDT annotation ", spec);
    }


    @Test(enabled = true)
    public void testChromosomeCountsPed() {
        final String MD5 = "71045f20e4cd9f5fdaa367f6d7324e59";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A ChromosomeCounts --variant:vcf " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + privateTestDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing ChromosomeCounts annotation with PED file", spec);
    }

    @Test(enabled = true)
    public void testInbreedingCoeffPed() {
        final String MD5 = "0cf7c115316950abc0213935b20a653a";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A InbreedingCoeff --variant:vcf " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + privateTestDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing InbreedingCoeff annotation with PED file", spec);
    }

    @Test(enabled = true)
    public void testAlleleTrimming() {
        final String MD5 = "5db0cd72ad0ffc711afb95df58de31fa";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A InbreedingCoeff --variant:vcf " + privateTestDir + "alleleTrim.vcf.gz" +
                        " -L 1:26608870-26608875 -no_cmdline_in_header --resource:exac " + privateTestDir +  "exacAlleleTrim.vcf.gz  -E exac.AC_Adj" +
                        "  -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing allele trimming annotation", spec);
    }

    @Test
    public void testStrandBiasBySample() throws IOException {
        // pipeline 1: create variant via HalotypeCaller with no default annotations
        final String base = String.format("-T HaplotypeCaller --disableDithering -R %s -I %s", REF, CEUTRIO_BAM) + " --no_cmdline_in_header -o %s -L 20:10130000-10134800";
        final WalkerTestSpec spec = new WalkerTestSpec(base, 1, Arrays.asList(""));
        final File outputVCF = executeTest("testStrandBiasBySample", spec).getFirst().get(0);

        // pipeline 2: create variant via HalotypeCaller; include StrandBiasBySample, exclude FisherStrand annotation
        //             re-Annotate the variant with VariantAnnotator using FisherStrand annotation
        final String baseNoFS = String.format("-T HaplotypeCaller --disableDithering -R %s -I %s", REF, CEUTRIO_BAM) + " --no_cmdline_in_header -o %s -L 20:10130000-10134800 -XA FisherStrand -A StrandBiasBySample";
        final WalkerTestSpec specNoFS = new WalkerTestSpec(baseNoFS, 1, Arrays.asList(""));
        specNoFS.disableShadowBCF();
        final File outputVCFNoFS = executeTest("testStrandBiasBySample component stand bias annotation", specNoFS).getFirst().get(0);

        final String baseAnn = String.format("-T VariantAnnotator -R %s -V %s", REF, outputVCFNoFS.getAbsolutePath()) + " --no_cmdline_in_header -o %s -L 20:10130000-10134800 -A FisherStrand";
        final WalkerTestSpec specAnn = new WalkerTestSpec(baseAnn, 1, Arrays.asList(""));
        specAnn.disableShadowBCF();
        final File outputVCFAnn = executeTest("testStrandBiasBySample re-annotation of FisherStrand", specAnn).getFirst().get(0);

        // confirm that the FisherStrand values are identical for the two pipelines
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(outputVCF);
        final LineIterator lineIterator = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIterator);

        final VCFCodec codecAnn = new VCFCodec();
        final FileInputStream sAnn = new FileInputStream(outputVCFAnn);
        final LineIterator lineIteratorAnn = codecAnn.makeSourceFromStream(new PositionalBufferedStream(sAnn));
        codecAnn.readHeader(lineIteratorAnn);

        while( lineIterator.hasNext() && lineIteratorAnn.hasNext() ) {
            final String line = lineIterator.next();
            Assert.assertFalse(line == null);
            final VariantContext vc = codec.decode(line);

            final String lineAnn = lineIteratorAnn.next();
            Assert.assertFalse(lineAnn == null);
            final VariantContext vcAnn = codecAnn.decode(lineAnn);

            Assert.assertTrue(vc.hasAttribute("FS"));
            Assert.assertTrue(vcAnn.hasAttribute("FS"));
            Assert.assertEquals(vc.getAttributeAsDouble("FS", 0.0), vcAnn.getAttributeAsDouble("FS", -1.0));
        }

        Assert.assertFalse(lineIterator.hasNext());
        Assert.assertFalse(lineIteratorAnn.hasNext());
    }

    @Test
    public void testStrandAlleleCountsBySample() {
        final String MD5 = "564aeeefad92353d66dbb2a2222d5108";
        final WalkerTestSpec spec = new WalkerTestSpec(
                "-T HaplotypeCaller --disableDithering " +
                String.format("-R %s -I %s ", REF, CEUTRIO_BAM) +
                "--no_cmdline_in_header -o %s -L 20:10130000-10134800 " +
                "-A StrandBiasBySample -A StrandAlleleCountsBySample",
                1, Arrays.asList(MD5)
        );
        spec.disableShadowBCF();
        executeTest("testStrandAlleleCountsBySample", spec);
    }

    @Test(enabled = false)
    public void testQualByDepth() throws IOException {

        /*
        NOTE: Currently, no HC tests explicitly set their PairHMM implementation, such that outputs may vary depending
        on the underlying system. You may run into this test failing due to one of many PairHMM versions activating.

        Specifically for this test, the site chr 20:10130722 seems to vary in the HC calculated QUAL and QD. This has
        NOTHING to do with a problem in the VariantAnnotator, and is an issue with the way PairHMM is used in tests.

        According to Karthik, things to look for in the logs that may help indicate that the failure is NOT something
        you changed, but the PairHMM switching on you:

        - Using un-vectorized C++ implementation of PairHMM
        - Using AVX accelerated implementation of PairHMM
        - Using SSE4.1 accelerated implementation of PairHMM
        - libVectorLoglessPairHMM unpacked successfully from GATK jar file
        - Using vectorized implementation of PairHMM
        - (Nothing, which I think means another implementation is being used)

        Upon any failures in any HC related tests for now, your best bet is to contact Mauricio to make sure that
        everything still looks ok.

        TODO: Stabilize the PairHMM from arbitrarily running different versions.
         */

        final String base = String.format("-T HaplotypeCaller --disableDithering -R %s -I %s", REF, CEUTRIO_BAM) + " --no_cmdline_in_header -o %s -L 20:10130000-10134800";
        final WalkerTestSpec spec = new WalkerTestSpec(base, 1, Arrays.asList("707be4798b1e14e3d6827a49104be120"));
        final File outputVCF = executeTest("testQualByDepth", spec).getFirst().get(0);

        final String baseNoQD = String.format("-T HaplotypeCaller --disableDithering -R %s -I %s", REF, CEUTRIO_BAM) + " --no_cmdline_in_header -o %s -L 20:10130000-10134800 -XA QualByDepth";
        final WalkerTestSpec specNoQD = new WalkerTestSpec(baseNoQD, 1, Arrays.asList("7e582b422a5de47706daefaae17b8245"));
        specNoQD.disableShadowBCF();
        final File outputVCFNoQD = executeTest("testQualByDepth calling without QD", specNoQD).getFirst().get(0);

        final String baseAnn = String.format("-T VariantAnnotator -R %s -V %s", REF, outputVCFNoQD.getAbsolutePath()) + " --no_cmdline_in_header -o %s -L 20:10130000-10134800 -A QualByDepth";
        final WalkerTestSpec specAnn = new WalkerTestSpec(baseAnn, 1, Arrays.asList("cc4e2b12872b2b25ab1c106516b2ac6a"));
        specAnn.disableShadowBCF();
        final File outputVCFAnn = executeTest("testQualByDepth re-annotation of QD", specAnn).getFirst().get(0);

        // confirm that the QD values are present in the new file for all biallelic variants
        // QD values won't be identical because some filtered reads are missing during re-annotation

        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(outputVCF);
        final LineIterator lineIterator = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIterator);

        final VCFCodec codecAnn = new VCFCodec();
        final FileInputStream sAnn = new FileInputStream(outputVCFAnn);
        final LineIterator lineIteratorAnn = codecAnn.makeSourceFromStream(new PositionalBufferedStream(sAnn));
        codecAnn.readHeader(lineIteratorAnn);

        while( lineIterator.hasNext() && lineIteratorAnn.hasNext() ) {
            final String line = lineIterator.next();
            Assert.assertFalse(line == null);
            final VariantContext vc = codec.decode(line);

            final String lineAnn = lineIteratorAnn.next();
            Assert.assertFalse(lineAnn == null);
            final VariantContext vcAnn = codecAnn.decode(lineAnn);

	    Assert.assertTrue(vc.hasAttribute("QD"));
	    Assert.assertTrue(vcAnn.hasAttribute("QD"));
        }

        Assert.assertFalse(lineIterator.hasNext());
        Assert.assertFalse(lineIteratorAnn.hasNext());
    }

    @Test
    public void testHomopolymerRunWindow() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + hg19ReferenceWithChrPrefixInChromosomeNames + " -A HomopolymerRun --variant:vcf " + privateTestDir + "problem_del.vcf  " +
                        "-U ALLOW_SEQ_DICT_INCOMPATIBILITY -L chr18:44382010-44384010 --reference_window_stop 59 --no_cmdline_in_header -o %s", 1,
                Arrays.asList("f3166721ef0380636590b3e860aa06af"));
        executeTest("Testing testHomopolymerRunWindow", spec);
    }

    @Test
    public void testHomopolymerRunTooBig() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + hg19ReferenceWithChrPrefixInChromosomeNames + " -A HomopolymerRun --variant:vcf " + privateTestDir + "problem_del.vcf  " +
                        "-U ALLOW_SEQ_DICT_INCOMPATIBILITY -L chr18:44382010-44384010 --no_cmdline_in_header -o %s", 1,
                Arrays.asList("3dba997d03779781c82a25ace69c838e"));
        executeTest("Testing HomopolymerRunTooBig", spec);
    }
}
