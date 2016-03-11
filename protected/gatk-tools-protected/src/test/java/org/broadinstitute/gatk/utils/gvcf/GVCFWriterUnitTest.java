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

package org.broadinstitute.gatk.utils.gvcf;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class GVCFWriterUnitTest extends BaseTest {
    private static class MockWriter implements VariantContextWriter {
        final List<VariantContext> emitted = new ArrayList<>();
        boolean headerWritten = false;
        boolean closed = false;

        @Override
        public void writeHeader(VCFHeader header) {
            headerWritten = true;
        }

        @Override
        public void close() {
            closed = true;
        }

        @Override
        public void add(VariantContext vc) {
            emitted.add(vc);
        }

        public boolean checkError(){
            return false;
        }
    }

    private MockWriter mockWriter;
    private List<Integer> standardPartition = Arrays.asList(1, 10, 20);
    private Allele REF = Allele.create("N", true);
    private Allele ALT = Allele.create("A");
    private List<Allele> ALLELES = Arrays.asList(REF, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
    private final String SAMPLE_NAME = "XXYYZZ";

    @BeforeMethod
    public void setUp() throws Exception {
        mockWriter = new MockWriter();
    }

    @Test
    public void testHeaderWriting() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);
        writer.writeHeader(new VCFHeader());
        Assert.assertTrue(mockWriter.headerWritten);
    }

    @Test
    public void testClose() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);
        writer.close();
        Assert.assertTrue(mockWriter.closed);
    }

    @Test
    public void testCloseWithoutClosingUnderlyingWriter() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);
        writer.close(false);
        Assert.assertFalse(mockWriter.closed);
    }

    private VariantContext makeHomRef(final String contig, final int start, final int GQ) {
        final VariantContextBuilder vcb = new VariantContextBuilder("test", contig, start, start, ALLELES);
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, REF));
        gb.GQ(GQ);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    private VariantContext makeHomRefAlt(final String contig, final int start, final int GQ) {
        final VariantContextBuilder vcb = new VariantContextBuilder("test", contig, start, start, Arrays.asList(REF, ALT));
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, REF));
        gb.GQ(GQ);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    private VariantContext makeNonRef(final String contig, final int start, final int GQ) {
        final VariantContextBuilder vcb = new VariantContextBuilder("test", contig, start, start, Arrays.asList(REF, ALT));
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, ALT));
        gb.GQ(GQ);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    private VariantContext makeDeletion(final String contig, final int start, final int size) {
        final String del = Utils.dupString("A", size);
        final String alt = del.substring(0, 1);
        final VariantContext vc = GATKVariantContextUtils.makeFromAlleles("test", contig, start, Arrays.asList(del, alt));
        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(vc.getReference(), vc.getAlternateAllele(0)));
        gb.GQ(50);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    @Test
    public void testCloseEmitsLastVariant() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 1);
    }

    @Test
    public void testCloseDoesntEmitsLastVariantWhenNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeNonRef("20", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 1);
    }

    @Test
    public void testCrossingContigBoundaryRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        writer.add(makeHomRef("20", 2, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef("21", 3, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);

        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(1), "21", 3, 3, false);
    }

    @Test
    public void testCrossingContigBoundaryToLowerPositionsRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 30, 30));
        writer.add(makeHomRef("20", 31, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef("21", 10, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 30, 31, false);
        writer.add(makeNonRef("21", 11, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(1), "21", 10, 10, false);
        assertGoodVC(mockWriter.emitted.get(2), "21", 11, 11, true);
    }

    @Test
    public void testCrossingContigBoundaryFromNonRefToLowerPositionsRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeNonRef("20", 20, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        writer.add(makeHomRef("21", 10, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 20, 20, true);
        writer.add(makeNonRef("21", 11, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(1), "21", 10, 10, false);
        assertGoodVC(mockWriter.emitted.get(2), "21", 11, 11, true);
    }

    @Test
    public void testCrossingContigBoundaryNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        writer.add(makeHomRef("20", 2, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeNonRef("21", 3, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "21", 3, 3, true);
    }

    @Test
    public void testCrossingContigBoundaryNonRefThenNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeNonRef("20", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        writer.add(makeNonRef("21", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 1, true);
        assertGoodVC(mockWriter.emitted.get(1), "21", 1, 1, true);
    }

    private void assertGoodVC(final VariantContext vc, final String contig, final int start, final int stop, final boolean nonRef) {
        Assert.assertEquals(vc.getChr(), contig);
        Assert.assertEquals(vc.getStart(), start);
        Assert.assertEquals(vc.getEnd(), stop);
        if ( nonRef ) {
            Assert.assertNotEquals(vc.getAlternateAllele(0), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        } else {
            Assert.assertEquals(vc.getNAlleles(), 2);
            Assert.assertEquals(vc.getAlternateAllele(0), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
            Assert.assertEquals(vc.getAttributeAsInt(VCFConstants.END_KEY, -1), stop);
            Assert.assertTrue(vc.hasGenotypes());
            Assert.assertTrue(vc.hasGenotype(SAMPLE_NAME));
            Assert.assertEquals(vc.getGenotypes().size(), 1);
            final Genotype g = vc.getGenotype(SAMPLE_NAME);
            Assert.assertEquals(g.hasAD(), false);
            Assert.assertEquals(g.hasLikelihoods(), true);
            Assert.assertEquals(g.hasPL(), true);
            Assert.assertEquals(g.getPL().length == 3, true);
            Assert.assertEquals(g.hasDP(), true);
            Assert.assertEquals(g.hasGQ(), true);
        }
    }

    @Test
    public void testVariantForcesNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        writer.add(makeHomRef("20", 2, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeNonRef("20", 3, 30));
        writer.add(makeHomRef("20", 4, 30));
        writer.add(makeHomRef("20", 5, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 3, 3, true);
        writer.close();
        assertGoodVC(mockWriter.emitted.get(2), "20", 4, 5, false);
    }

    @Test
    public void testEmittingTwoBands() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef("20", 3, 50));
        writer.add(makeHomRef("20", 4, 50));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 3, 4, false);
    }

    @Test
    public void testNonContiguousBlocks() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        writer.add(makeHomRef("20", 10, 0));
        writer.add(makeHomRef("20", 11, 0));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 10, 11, false);
    }

    @Test
    public void testDeletion() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        writer.add(makeDeletion("20", 3, 3));
        writer.add(makeHomRef("20", 4, 0));
        writer.add(makeHomRef("20", 5, 0));
        writer.add(makeHomRef("20", 6, 0));
        writer.add(makeHomRef("20", 7, 0));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 3, 5, true);
        assertGoodVC(mockWriter.emitted.get(2), "20", 6, 7, false);
    }

    @Test
    public void testHomRefAlt() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        writer.add(makeHomRefAlt("20", 3, 0));
        writer.add(makeHomRef("20", 4, 0));
        writer.add(makeHomRef("20", 5, 0));
        writer.add(makeHomRef("20", 6, 0));
        writer.add(makeHomRef("20", 7, 0));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        Assert.assertFalse(mockWriter.emitted.get(1).hasAttribute("END"));
        Assert.assertFalse(mockWriter.emitted.get(1).hasAttribute("BLOCK_SIZE"));
        assertGoodVC(mockWriter.emitted.get(2), "20", 4, 7, false);
    }

    @DataProvider(name = "BandPartitionData")
    public Object[][] makeBandPartitionData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{null, false});
        tests.add(new Object[]{Collections.emptyList(), false});
        tests.add(new Object[]{Arrays.asList(1), true});
        tests.add(new Object[]{Arrays.asList(1, 10), true});
        tests.add(new Object[]{Arrays.asList(1, 10, 30), true});
        tests.add(new Object[]{Arrays.asList(10, 1, 30), false});
        tests.add(new Object[]{Arrays.asList(-1, 1), false});
        tests.add(new Object[]{Arrays.asList(1, null, 10), false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BandPartitionData")
    public void testMyData(final List<Integer> partitions, final boolean expectedGood) {
        try {
            GVCFWriter.parsePartitions(partitions,2);
            Assert.assertTrue(expectedGood, "Expected to fail but didn't");
        } catch ( Exception e ) {
            Assert.assertTrue(! expectedGood, "Expected to succeed but failed with message " + e.getMessage());
        }
    }
}
