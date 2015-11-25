/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class HaplotypeCallerIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final static String NA12878_BAM = privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final static String NA12878_CHR20_BAM = validationDataLocation + "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam";
    final static String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final static String NA12878_RECALIBRATED_BAM = privateTestDir + "NA12878.100kb.BQSRv2.example.bam";
    final static String NA12878_PCRFREE = privateTestDir + "PCRFree.2x250.Illumina.20_10_11.bam";
    final static String NA12878_PCRFREE250_ADAPTER_TRIMMED = privateTestDir + "PCRFree.2x250.b37_decoy.NA12878.adapter_trimmed-10000000-11000000.bam";
    final static String CEUTRIO_MT_TEST_BAM = privateTestDir + "CEUTrio.HiSeq.b37.MT.1_50.bam";
    final static String INTERVALS_FILE = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.test.intervals";
    final static String GGA_INTERVALS_FILE = privateTestDir + "haplotype-caller-reduced-test-interval.list";
    final static String HMM_SUB_IMPLEMENTATION = "UNVECTORIZED";
    final static String ALWAYS_LOAD_VECTOR_HMM = "-alwaysloadVectorHMM";

    private void HCTest(String bam, String args, String md5) throws IOException {
        final String base = String.format("-T HaplotypeCaller --contamination_fraction_to_filter 0.05 --disableDithering --pcr_indel_model NONE --maxReadsInRegionPerSample 1000 --minReadsPerAlignmentStart 5 --maxProbPropagationDistance 50 --activeProbabilityThreshold 0.002 -pairHMMSub %s %s -R %s -I %s -L %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, bam, INTERVALS_FILE) + " --no_cmdline_in_header -o %s -minPruning 3";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        final File outputVCF =  executeTest("testHaplotypeCaller: args=" + args, spec).getFirst().get(0);
        Assert.assertFalse(FileUtils.readFileToString(outputVCF).contains(VCFConstants.MAPPING_QUALITY_ZERO_KEY));
    }

    @Test
    public void testHaplotypeCallerMultiSample() throws IOException {
        HCTest(CEUTRIO_BAM, "", "bdc5a942a7df2d4c61388ecc74a163cf");
    }

    @Test
    public void testHaplotypeCallerSingleSample() throws IOException {

        HCTest(NA12878_BAM, "", "5941f88fa1372eb598880a4711ac08de");
    }

    @Test
    public void testHaplotypeCallerMultiSampleHaploid() throws IOException {
        HCTest(CEUTRIO_BAM, "-ploidy 1", "4c32640204a11fe46a23d97b1012bca2");
    }

    @Test
    public void testHaplotypeCallerSingleSampleHaploid() throws IOException {
        HCTest(NA12878_BAM, "-ploidy 1", "e09c5fd862f86f69289d847e3b293e19");
    }

    @Test
    public void testHaplotypeCallerSingleSampleTetraploid() throws IOException {
        HCTest(NA12878_BAM, "-ploidy 4", "bbae9f77683591200b2034d5461e6425");
    }

    @Test
    public void testHaplotypeCallerMinBaseQuality() throws IOException {
        HCTest(NA12878_BAM, "-mbq 15", "5941f88fa1372eb598880a4711ac08de");
    }

    @Test
    public void testHaplotypeCallerMinBaseQualityHaploid() throws IOException {
        HCTest(NA12878_BAM, "-mbq 15 -ploidy 1", "e09c5fd862f86f69289d847e3b293e19");
    }

    @Test
    public void testHaplotypeCallerMinBaseQualityTetraploid() throws IOException {
        HCTest(NA12878_BAM, "-mbq 15 -ploidy 4", "bbae9f77683591200b2034d5461e6425");
    }

    @Test
    public void testHaplotypeCallerGraphBasedSingleSample() throws IOException {
        HCTest(NA12878_BAM, "-likelihoodEngine GraphBased", "9a4ba54146449914b1fbaed5465b692e");
    }

    @Test
    public void testHaplotypeCallerGraphBasedMultiSampleHaploid() throws IOException {
        HCTest(CEUTRIO_BAM, "-likelihoodEngine GraphBased -ploidy 1", "5333e3ff4d7fc53624eab801250be4f0");
    }

    @Test
    public void testHaplotypeCallerGraphBasedMultiSample() throws IOException {
        HCTest(CEUTRIO_BAM, "-likelihoodEngine GraphBased", "e9759bf254b7b4d06e4a68c94a417cf1");
    }

    @Test
    public void testHaplotypeCallerSingleSampleWithDbsnp() throws IOException {
        HCTest(NA12878_BAM, "-D " + b37dbSNP132, "7466c4d9ce6a1d23bcb428fc9f446843");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGA() throws IOException {
        HCTest(CEUTRIO_BAM, "--max_alternate_alleles 3 -gt_mode GENOTYPE_GIVEN_ALLELES -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf" +
                " -isr INTERSECTION -L " + GGA_INTERVALS_FILE,
                "c2fe156912e59626c0393b9b47c9419e");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGAHaploid() throws IOException {
        HCTest(CEUTRIO_BAM, "--max_alternate_alleles 3 -gt_mode GENOTYPE_GIVEN_ALLELES -ploidy 1 -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf -isr INTERSECTION -L 20:10080000-10100000",
                "49f737983375c740361daa2d48ed0249");
    }

    @Test
    public void testHaplotypeCallerMultiSampleGGATetraploid() throws IOException {
        HCTest(CEUTRIO_BAM, "--max_alternate_alleles 3 -gt_mode GENOTYPE_GIVEN_ALLELES -ploidy 4 -alleles " + validationDataLocation + "combined.phase1.chr20.raw.indels.sites.vcf -isr INTERSECTION -L 20:10080000-10100000",
                "6ba0fa930899f70fee9a7b1161508f93");
    }

    @Test
    public void testHaplotypeCallerInsertionOnEdgeOfContig() throws IOException {
        HCTest(CEUTRIO_MT_TEST_BAM, "-L MT:1-10", "da1f6b9a7e5913910531b00f3b35ce06");
    }

    private void HCTestIndelQualityScores(String bam, String args, String md5) {
        final String base = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, bam) + " -L 20:10,005,000-10,025,000 --no_cmdline_in_header -o %s -minPruning 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
        executeTest("testHaplotypeCallerIndelQualityScores: args=" + args, spec);
    }

    @Test
    public void testHaplotypeCallerSingleSampleIndelQualityScores() {
        HCTestIndelQualityScores(NA12878_RECALIBRATED_BAM, "", "22028b104dd60b23a399e5e3a877a1fb");
    }

    private void HCTestNearbySmallIntervals(String bam, String args, String md5) {
        try {
            final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(b37KGReference));
            final GenomeLocParser parser = new GenomeLocParser(fasta.getSequenceDictionary());

            final String base = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, bam) + " -L 20:10,001,603-10,001,642 -L 20:10,001,653-10,001,742 --no_cmdline_in_header -o %s";
            final WalkerTestSpec spec = new WalkerTestSpec(base + " " + args, Arrays.asList(md5));
            for( final File vcf : executeTest("testHaplotypeCallerNearbySmallIntervals: args=" + args, spec).getFirst() ) {
                if( containsDuplicateRecord(vcf, parser) ) {
                    throw new IllegalStateException("Duplicate records detected but there should be none.");
                }
            }
        } catch( FileNotFoundException e ) {
            throw new IllegalStateException("Could not find the b37 reference file.");
        }
    }

    private boolean containsDuplicateRecord( final File vcf, final GenomeLocParser parser ) {
        final List<Pair<GenomeLoc, HaplotypeCallerGenotypingEngine.Event>> VCs = new ArrayList<>();
        try {
            for( final VariantContext vc :  GATKVCFUtils.readVCF(vcf).getSecond() ) {
                VCs.add(new Pair<>(parser.createGenomeLoc(vc), new HaplotypeCallerGenotypingEngine.Event(vc)));
            }
        } catch( IOException e ) {
            throw new IllegalStateException("Somehow the temporary VCF from the integration test could not be read.");
        }

        final Set<Pair<GenomeLoc, HaplotypeCallerGenotypingEngine.Event>> VCsAsSet = new HashSet<>(VCs);
        return VCsAsSet.size() != VCs.size(); // The se will remove duplicate Events.
    }


    @Test
    public void testHaplotypeCallerNearbySmallIntervals() {
        HCTestNearbySmallIntervals(NA12878_BAM, "", "6e55de7bf49ecdc71f6d4c9565a19853");
    }

    // This problem bam came from a user on the forum and it spotted a problem where the ReadClipper
    // was modifying the GATKSamRecord and that was screwing up the traversal engine from map call to
    // map call. So the test is there for consistency but not for correctness. I'm not sure we can trust
    // any of the calls in that region because it is so messy.
    @Test
    public void HCTestProblematicReadsModifiedInActiveRegions() {
        final String base = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, privateTestDir + "haplotype-problem-4.bam") + " --no_cmdline_in_header -o %s -minPruning 3 -L 4:49139026-49139965";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("f8774326a268f7d394b333818d15d05c"));
        executeTest("HCTestProblematicReadsModifiedInActiveRegions: ", spec);
    }

    @Test
    public void HCTestStructuralIndels() {
        final String base = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, privateTestDir + "AFR.structural.indels.bam") + " --no_cmdline_in_header -o %s -minPruning 6 -L 20:8187565-8187800 -L 20:18670537-18670730";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Arrays.asList("034c6b151dfde537d5843f70880bf8a4"));
        executeTest("HCTestStructuralIndels: ", spec);
    }

    @Test
    public void HCTestDoesNotFailOnBadRefBase() {
        // don't care about the output - just want to make sure it doesn't fail
        final String base = String.format("-T HaplotypeCaller --disableDithering -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, privateTestDir + "NA12878.readsOverBadBase.chr3.bam") + " --no_cmdline_in_header -o /dev/null -L 3:60830000-60840000 --minPruning 3 -stand_call_conf 2 -stand_emit_conf 2";
        final WalkerTestSpec spec = new WalkerTestSpec(base, Collections.<String>emptyList());
        executeTest("HCTestDoesNotFailOnBadRefBase: ", spec);
    }

    @Test
    public void HCTestDanglingTailMergingForDeletions() throws IOException {
        final String base = String.format("-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, NA12878_BAM) + " --no_cmdline_in_header -o %s -L 20:10130740-10130800 --allowNonUniqueKmersInRef";
        final WalkerTestSpec spec = new WalkerTestSpec(base, 1, Arrays.asList(""));
        final File outputVCF = executeTest("HCTestDanglingTailMergingForDeletions", spec).getFirst().get(0);

        // confirm that the call is the correct one
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(outputVCF);
        final LineIterator lineIterator = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIterator);
        final String line = lineIterator.next();
        Assert.assertFalse(line == null);
        final VariantContext vc = codec.decode(line);
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertTrue(vc.getReference().basesMatch("ATGTATG"));
        Assert.assertTrue(vc.getAlternateAllele(0).basesMatch("A"));
    }


    private static final String LEFT_ALIGNMENT_BAMOUT_TEST_INPUT = privateTestDir + "/bamout-indel-left-align-bugfix-input.bam";

    private static final String LEFT_ALIGNMENT_BAMOUT_TEST_OUTPUT = privateTestDir + "/bamout-indel-left-align-bugfix-expected-output.bam";

    @Test
    public void testLeftAlignmentBamOutBugFix() {
        final String base = String.format("-T HaplotypeCaller -pairHMMSub %s %s -R %s -I %s", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, REF, LEFT_ALIGNMENT_BAMOUT_TEST_INPUT)
                + " --no_cmdline_in_header -bamout %s -o /dev/null -L 1:11740000-11740700 --allowNonUniqueKmersInRef";
        final WalkerTestSpec spec = new WalkerTestSpec(base, 1, Arrays.asList("c1840293b4565ce1cef393c6a0d5fc9a"));
        executeTest("LeftAlignmentBamOutBugFix", spec);
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // test dbSNP annotation
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void HCTestDBSNPAnnotationWGS() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub " +  HMM_SUB_IMPLEMENTATION + " " + ALWAYS_LOAD_VECTOR_HMM + " -R " + b37KGReference + " --no_cmdline_in_header -I " + NA12878_PCRFREE + " -o %s -L 20:10,090,000-10,100,000 -D " + b37dbSNP132, 1,
                Arrays.asList("ab816ff6facc08acb19c55bd6e828f02"));
        executeTest("HC calling with dbSNP ID annotation on WGS intervals", spec);
    }

    @Test
    public void HCTestDBSNPAnnotationWEx() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller --disableDithering --pcr_indel_model NONE -pairHMMSub " +  HMM_SUB_IMPLEMENTATION + " " + ALWAYS_LOAD_VECTOR_HMM + " -R " + b37KGReference + " --no_cmdline_in_header -I " + NA12878_PCRFREE + " -o %s -L 20:10,100,000-11,000,000 -D " + b37dbSNP132
                        + " -L " + hg19Intervals + " -isr INTERSECTION", 1,
                Arrays.asList("fd1a539e14902f6957eb939aac1412f0"));
        executeTest("HC calling with dbSNP ID annotation on WEx intervals", spec);
    }

    @Test
    public void HCTestDBSNPAnnotationWGSGraphBased() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller -likelihoodEngine GraphBased --disableDithering --pcr_indel_model NONE -pairHMMSub " +  HMM_SUB_IMPLEMENTATION + " " + ALWAYS_LOAD_VECTOR_HMM + " -R " + b37KGReference + " --no_cmdline_in_header -I " + NA12878_PCRFREE + " -o %s -L 20:10,090,000-10,100,000 -D " + b37dbSNP132, 1,
                Arrays.asList("1810c35e8298dff4aa1b7b04fb5f4962"));
        executeTest("HC calling with dbSNP ID annotation on WGS intervals", spec);
    }

    @Test
    public void HCTestDBSNPAnnotationWExGraphBased() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller -likelihoodEngine GraphBased --disableDithering --pcr_indel_model NONE -pairHMMSub " +  HMM_SUB_IMPLEMENTATION + " " + ALWAYS_LOAD_VECTOR_HMM + " -R " + b37KGReference + " --no_cmdline_in_header -I " + NA12878_PCRFREE + " -o %s -L 20:10,000,000-11,000,000 -D " + b37dbSNP132
                        + " -L " + hg19Intervals + " -isr INTERSECTION", 1,
                Arrays.asList("8a06a53388d73ec667b353379f3b351e"));
        executeTest("HC calling with dbSNP ID annotation on WEx intervals", spec);
    }

    @Test
    public void HCTestGraphBasedPCRFreePositiveLogLkFix() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller -likelihoodEngine GraphBased --disableDithering --pcr_indel_model NONE -pairHMMSub " +  HMM_SUB_IMPLEMENTATION + " " + ALWAYS_LOAD_VECTOR_HMM + " -R " + hg19Reference + " --no_cmdline_in_header -I " + NA12878_PCRFREE250_ADAPTER_TRIMMED + " -o %s -L 20:10,024,000-10,024,500 "
                        , 1,
                Arrays.asList(""));
        executeTest("HCTestGraphBasedPCRFreePositiveLogLkFix", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // test PCR indel model
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void HCTestAggressivePcrIndelModelWGS() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller --disableDithering --pcr_indel_model AGGRESSIVE -pairHMMSub " + HMM_SUB_IMPLEMENTATION +  " " + ALWAYS_LOAD_VECTOR_HMM + " -R " + b37KGReference + " --no_cmdline_in_header -I " + NA12878_BAM + " -o %s -L 20:10,270,000-10,300,000", 1,
                Arrays.asList("be71ab16d55437cbe3005ea3b93cece6"));
        executeTest("HC calling with aggressive indel error modeling on WGS intervals", spec);
    }

    @Test
    public void HCTestConservativePcrIndelModelWGS() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T HaplotypeCaller --disableDithering --pcr_indel_model CONSERVATIVE -pairHMMSub " + HMM_SUB_IMPLEMENTATION +  " " + ALWAYS_LOAD_VECTOR_HMM  + " -R " + b37KGReference + " --no_cmdline_in_header -I " + NA12878_BAM + " -o %s -L 20:10,270,000-10,300,000", 1,
                Arrays.asList("8a5185d0a9400c8b3f4b12da65181c4b"));
        executeTest("HC calling with conservative indel error modeling on WGS intervals", spec);
    }

    @Test
    public void testNoSuchEdgeBugFix() {
        final String commandLine = String.format("-T HaplotypeCaller --pcr_indel_model NONE -pairHMMSub %s %s -R %s -I %s -L %s -dontTrimActiveRegions -ERC GVCF " +
                "-likelihoodEngine GraphBased -variant_index_type %s -variant_index_parameter %d",
                HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, b37KGReferenceWithDecoy, privateTestDir + "graphbased_no_such_edge_bug.bam", privateTestDir + "graphbased_no_such_edge_bug.intervals.bed",
                GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest("testGraphBasedNoSuchEdgeBugFix", spec);
    }

    // This test takes longer than 15 secs ... ~ 25-35 ,
    @Test
    public void testLackSensitivityDueToBadHaplotypeSelectionFix() {
        final String commandLine = String.format("-T HaplotypeCaller -pairHMMSub %s %s -R %s -I %s -L %s --no_cmdline_in_header --maxNumHaplotypesInPopulation 16",
                HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, b37KGReferenceWithDecoy, privateTestDir + "hc-lack-sensitivity.bam", privateTestDir + "hc-lack-sensitivity.interval_list");
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("02ec3090c4a6359fa10e6d8b30e1d5a2"));
        spec.disableShadowBCF();
        executeTest("testLackSensitivityDueToBadHaplotypeSelectionFix", spec);
    }

    @Test
    public void testMissingKeyAlternativeHaplotypesBugFix() {
        final String commandLine = String.format("-T HaplotypeCaller -pairHMMSub %s %s -R %s -I %s -L %s --no_cmdline_in_header ",
                HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, b37KGReferenceWithDecoy, privateTestDir + "lost-alt-key-hap.bam", privateTestDir + "lost-alt-key-hap.interval_list");
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine + " -o %s", Arrays.asList("13c9f094b9c54960dc2fd3a1815a2645"));
        spec.disableShadowBCF();
        executeTest("testMissingKeyAlternativeHaplotypesBugFix", spec);
    }

    @Test
    public void testDifferentIndelLocationsDueToSWExactDoubleComparisonsFix() {
        final String SHORT_INTERVAL = "12:7342264-7342464";
        final String LONG_INTERVAL = "12:7342270-7342824";
        final String TEST_BAM = privateTestDir + "sw_epsilon_test.bam";
        final String REFERENCE = b37KGReference;
        final String DBSNP = b37dbSNP138;
        final String commandLineWithoutInterval = String.format("-T HaplotypeCaller -pairHMMSub %s %s -I %s -R %s -D %s "
                + "-variant_index_type  LINEAR -variant_index_parameter 128000 --no_cmdline_in_header "
                + "-stand_call_conf 10.0 -stand_emit_conf 10.0", HMM_SUB_IMPLEMENTATION, ALWAYS_LOAD_VECTOR_HMM, TEST_BAM, REFERENCE, DBSNP);

        final String commandLineShortInterval = commandLineWithoutInterval + " -L " + SHORT_INTERVAL;
        final String commandLineLongInterval = commandLineWithoutInterval + " -L " + LONG_INTERVAL;

        //README: update MD5s accordingly when needed
        // but please make sure that both outputs get the same variant,
        // alleles all with DBSNP ids
        // We test here that change in active region size does not have an effect in placement of indels.
        final String md5 = "6121d05f96eca3b1dbe3a881d968b6c5";
        final WalkerTestSpec shortSpec = new WalkerTestSpec(commandLineShortInterval + " -o %s",Arrays.asList(md5));
        executeTest("testDifferentIndelLocationsDueToSWExactDoubleComparisonsFix::shortInterval",shortSpec);
        final WalkerTestSpec longSpec = new WalkerTestSpec(commandLineLongInterval + " -o %s",Arrays.asList(md5));
        executeTest("testDifferentIndelLocationsDueToSWExactDoubleComparisonsFix::longInterval",longSpec);
    }

    @Test
    public void testHaplotypeCallerPairHMMException(){
        executeTest("HaplotypeCallerPairHMMException",
                new WalkerTest.WalkerTestSpec(
                        " -T HaplotypeCaller" +
                                " --contamination_fraction_to_filter 0.05 --disableDithering --pcr_indel_model NONE --maxReadsInRegionPerSample 1000 " +
                                " --minReadsPerAlignmentStart 5 --maxProbPropagationDistance 50 --activeProbabilityThreshold 0.002 " +
                                " --no_cmdline_in_header -minPruning 3 -pairHMM VECTOR_LOGLESS_CACHING -pairHMMSub TEST_BEYOND_CAPABILITIES " +
                                ALWAYS_LOAD_VECTOR_HMM +
                                " -R " + REF +
                                " -I " + NA12878_BAM +
                                " -L " + INTERVALS_FILE +
                                " -o %s",
                        1, UserException.HardwareFeatureException.class));
    }

    @Test
    public void testHaplotypeCallerDcovException(){
        executeTest("HaplotypeCallerDcovException",
                new WalkerTest.WalkerTestSpec(
                        " -T HaplotypeCaller" +
                                " --contamination_fraction_to_filter 0.05 --disableDithering --pcr_indel_model NONE --maxReadsInRegionPerSample 1000 " +
                                " --minReadsPerAlignmentStart 5 --maxProbPropagationDistance 50 --activeProbabilityThreshold 0.002 " +
                                " --no_cmdline_in_header -minPruning 3 -pairHMM VECTOR_LOGLESS_CACHING -pairHMMSub " + HMM_SUB_IMPLEMENTATION +
                                " -dcov 50" +
                                " -R " + REF +
                                " -I " + NA12878_BAM +
                                " -L " + INTERVALS_FILE +
                                " -o %s",
                        1, UserException.CommandLineException.class));
    }

    @Test 
    public void testHaplotypeCallerMergeVariantsViaLDException(){
        executeTest("HaplotypeCallerMergeVariantsViaLDException",
                new WalkerTest.WalkerTestSpec(
                        " -T HaplotypeCaller" +
                                " -R " + REF +
                                " -I " + NA12878_BAM +
                                " -L " + INTERVALS_FILE +
				                " --mergeVariantsViaLD " +
                                " -o %s",
                        1, UserException.DeprecatedArgument.class));
    }

    @Test
    public void testHaplotypeCallerTandemRepeatAnnotator() throws IOException{
        HCTest(NA12878_BAM, " -L 20:10001000-10010000 -A TandemRepeatAnnotator -XA MappingQualityZero -XA SpanningDeletions", "03738462e7f0b6f149f40b790a3a7261");
    }

    @Test
    public void testHBaseCountsBySample() throws IOException{
        HCTest(NA12878_BAM, " -L 20:10001000-10010000 -A BaseCountsBySample", "8739d92898113436666f1cff3bc39bc5");
    }
}

