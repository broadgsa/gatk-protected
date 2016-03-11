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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.tools.walkers.genotyper.*;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.genotyper.SampleListUtils;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReferenceConfidenceModelUnitTest extends BaseTest {
    GenomeLocParser parser;
    final String RGID = "ID1";
    GATKSAMReadGroupRecord rg;
    final String sample = "NA12878";
    final SampleList samples = SampleListUtils.singletonList(sample);
    SAMFileHeader header;
    ReferenceConfidenceModel model;

    @BeforeClass
    public void setUp() throws Exception {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        rg = new GATKSAMReadGroupRecord(RGID);
        rg.setSample(sample);
        header.addReadGroup(rg);
        parser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @BeforeMethod
    public void setupModel() {
        model = new ReferenceConfidenceModel(parser, samples, header, 10);
    }

    @DataProvider(name = "CalcNIndelInformativeReadsData")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        { // very basic testing
            final String ref  = "ACGT";
            final String read = "ACGT";
            tests.add(new Object[]{read, ref, 1, Arrays.asList(1, 1, 1, 0)});
            tests.add(new Object[]{read, ref, 2, Arrays.asList(1, 1, 0, 0)});
            tests.add(new Object[]{read, ref, 3, Arrays.asList(1, 0, 0, 0)});
            tests.add(new Object[]{read, ref, 4, Arrays.asList(0, 0, 0, 0)});
        }

        { // actually interesting case where some sites aren't informative
            final String ref   = "NNAAAANN";
            final String read1 = "NNA";
            final String read2 = "NNAA";
            final String read3 = "NNAAA";
            final String read4 = "NNAAAA";
            final String read5 = "NNAAAAN";
            tests.add(new Object[]{read1, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read2, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read3, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read4, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read5, ref, 1, Arrays.asList(1, 1, 1, 1, 1, 1, 0, 0)});
        }

        {
            for ( final String repeatUnit : Arrays.asList("A", "CA", "TAG", "TAGC", "TCAGA")) {
                final String anchor = Utils.dupString("N", repeatUnit.length());
                for ( int nUnits = 1; nUnits < 10; nUnits++ ) {
                    final String repeat = Utils.dupString(repeatUnit, nUnits);
                    final String ref = anchor + repeat + anchor;
                    for ( int readLen = repeatUnit.length(); readLen < repeat.length(); readLen++ ) {
                        final String read = anchor + repeat.substring(0, readLen);
                        final List<Integer> expected = new LinkedList<>();
                        for ( int i = 0; i < anchor.length(); i++ ) expected.add(1);
                        for ( int i = 0; i < repeat.length(); i++ ) expected.add(readLen == repeat.length() ? 1 : 0);
                        for ( int i = 0; i < anchor.length(); i++ ) expected.add(0);
                        tests.add(new Object[]{read, ref, repeatUnit.length(), expected});

                        final List<Integer> result = new ArrayList<>(Collections.nCopies(ref.length() - anchor.length(), 1));
                        result.addAll(Collections.nCopies(anchor.length(), 0));
                        tests.add(new Object[]{ref, ref, repeatUnit.length(), result});
                    }
                }

            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CalcNIndelInformativeReadsData")
    public void testCalcNIndelInformativeReads(final String readBases, final String ref, final int maxIndelSize, final List<Integer> expected ) {
        final byte qual = (byte)30;
        final byte[] quals = Utils.dupBytes(qual, readBases.length());

        for ( int i = 0; i < readBases.getBytes().length; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(readBases.getBytes(), quals, readBases.length() + "M");
            final GenomeLoc loc = new UnvalidatingGenomeLoc("20", 0, i, i);
            final ReadBackedPileup pileup = new ReadBackedPileupImpl(loc, Collections.singletonList(read), i);
            final int actual = model.calcNIndelInformativeReads(pileup, i, ref.getBytes(), maxIndelSize);
            Assert.assertEquals(actual, (int)expected.get(i), "failed at position " + i);
        }
    }

    @Test
    public void testClose() {
        model.close();
    }

    @Test
    public void testWorstGL() {
        final GenotypeLikelihoods gq10 = GenotypeLikelihoods.fromPLField("0,10,100");
        final GenotypeLikelihoods gq20 = GenotypeLikelihoods.fromPLField("0,20,200");
        final GenotypeLikelihoods gq0 = GenotypeLikelihoods.fromPLField("20,0,200");

        Assert.assertSame(model.getGLwithWorstGQ(gq10, gq20), gq10);
        Assert.assertSame(model.getGLwithWorstGQ(gq20, gq10), gq10);
        Assert.assertSame(model.getGLwithWorstGQ(gq10, gq0), gq0);
        Assert.assertSame(model.getGLwithWorstGQ(gq0, gq10), gq0);
    }

    @Test
    public void testIndelLikelihoods() {
        GenotypeLikelihoods prev = model.getIndelPLs(HomoSapiensConstants.DEFAULT_PLOIDY,0);
        Assert.assertEquals(prev.getAsPLs(), new int[]{0, 0, 0});
        Assert.assertEquals(-10 * prev.getLog10GQ(GenotypeType.HOM_REF), 0.0);

        for ( int i = 1; i <= ReferenceConfidenceModel.MAX_N_INDEL_INFORMATIVE_READS; i++ ) {
            final GenotypeLikelihoods current = model.getIndelPLs(HomoSapiensConstants.DEFAULT_PLOIDY,i);
            final double prevGQ = -10 * prev.getLog10GQ(GenotypeType.HOM_REF);
            final double currGQ = -10 * current.getLog10GQ(GenotypeType.HOM_REF);
            Assert.assertTrue(prevGQ < currGQ, "GQ Failed with prev " + prev + " curr " + current + " at " + i);
            Assert.assertTrue(prev.getAsPLs()[1] < current.getAsPLs()[1], "het PL failed with prev " + prev + " curr " + current + " at " + i);
            Assert.assertTrue(prev.getAsPLs()[2] < current.getAsPLs()[2], "hom-var PL Failed with prev " + prev + " curr " + current + " at " + i);
//            logger.warn("result at " + i + " is " + current);
            prev = current;
        }
    }

    @Test
    public void testOverlappingVariantContext() {
        final VariantContext vc10 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 10, Arrays.asList("A", "C"));
        final VariantContext vc13 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 13, Arrays.asList("A", "C"));
        final VariantContext vc12_15 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 12, Arrays.asList("ACAT", "A"));
        final VariantContext vc18 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 18, Arrays.asList("A", "ACAT"));

        final List<VariantContext> calls = Arrays.asList(vc13, vc12_15, vc18, vc10);

        checkOverlapping(8, calls, null);
        checkOverlapping(9, calls, null);
        checkOverlapping(10, calls, vc10);
        checkOverlapping(11, calls, null);
        checkOverlapping(12, calls, vc12_15);
        checkOverlapping(13, calls, vc13);
        checkOverlapping(14, calls, vc12_15);
        checkOverlapping(15, calls, vc12_15);
        checkOverlapping(16, calls, null);
        checkOverlapping(17, calls, null);
        checkOverlapping(18, calls, vc18);
        checkOverlapping(19, calls, null);
        checkOverlapping(20, calls, null);
    }

    private void checkOverlapping(final int pos, Collection<VariantContext> calls, final VariantContext expected) {
        final GenomeLoc loc = parser.createGenomeLoc(parser.getContigs().getSequences().get(0).getSequenceName(), pos, pos);
        final VariantContext actual = model.getOverlappingVariantContext(loc, calls);
        Assert.assertEquals(actual, expected);
    }

    //
    // test reference calculation
    //
    private class RefConfData {
        final String ref;
        final int extension;
        final Haplotype refHap;
        final GenomeLoc refLoc, paddedRefLoc;
        final ActiveRegion region;
        int readCounter = 0;

        private RefConfData(String ref, int extension) {
            this.ref = ref;
            this.extension = extension;

            refLoc = parser.createGenomeLoc("chr1", getStart(), getEnd());
            paddedRefLoc = parser.createGenomeLoc("chr1", getStart() - extension, getEnd() + extension);
            region = new ActiveRegion(getRefLoc(), parser, extension);
            final String pad = Utils.dupString("N", extension);
            refHap = ReferenceConfidenceModel.createReferenceHaplotype(getActiveRegion(), (pad + ref + pad).getBytes(), getPaddedRefLoc());
        }

        public GenomeLoc getRefLoc() { return refLoc; }
        public GenomeLoc getPaddedRefLoc() { return paddedRefLoc; }
        public ActiveRegion getActiveRegion() { return region; }
        public Haplotype getRefHap() { return refHap; }
        public int getStart() { return 100; }
        public int getEnd() { return getStart() + getRefLength() - 1; }
        public byte[] getRefBases() { return ref.getBytes(); }
        public int getRefLength() { return ref.length(); }

        public GATKSAMRecord makeRead(final int start, final int length) {
            final byte[] quals = Utils.dupBytes((byte)30, length);
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read " + readCounter++, 0, start + getStart(), ref.substring(start, start + length).getBytes(), quals, length + "M");
            read.setReadGroup(rg);
            return read;
        }
    }


    @DataProvider(name = "RefConfidenceData")
    public Object[][] makeRefConfidenceData() {
        List<Object[]> tests = new ArrayList<>();

        for ( int i = 0; i < 10; i++ ) {
            for ( final int extension : Arrays.asList(0, 10) ) {
                tests.add(new Object[]{i, extension});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RefConfidenceData")
    public void testRefConfidenceBasic(final int nReads, final int extension) {
        final RefConfData data = new RefConfData("ACGTAACCGGTT", extension);
        final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
        final List<VariantContext> calls = Collections.emptyList();

        for ( int i = 0; i < nReads; i++ ) {
            data.getActiveRegion().add(data.makeRead(0, data.getRefLength()));
        }

        final ReadLikelihoods<Haplotype> likelihoods = HaplotypeCaller.createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();
        final List<Integer> expectedDPs = Collections.nCopies(data.getActiveRegion().getLocation().size(), nReads);
        final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, genotypingModel, calls);
        checkReferenceModelResult(data, contexts, expectedDPs, calls);
    }

    @Test
    public void testRefConfidencePartialReads() {

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();
        final String ref = "ACGTAACCGGTT";
        for ( int readLen = 3; readLen < ref.length(); readLen++ ) {
            for ( int start = 0; start < ref.length() - readLen; start++ ) {
                final RefConfData data = new RefConfData(ref, 0);
                final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
                final List<VariantContext> calls = Collections.emptyList();

                data.getActiveRegion().add(data.makeRead(start, readLen));
                final ReadLikelihoods<Haplotype> likelihoods = HaplotypeCaller.createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

                final List<Integer> expectedDPs = new ArrayList<>(Collections.nCopies(data.getActiveRegion().getLocation().size(), 0));
                for ( int i = start; i < readLen + start; i++ ) expectedDPs.set(i, 1);
                final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, genotypingModel, calls);
                checkReferenceModelResult(data, contexts, expectedDPs, calls);
            }
        }
    }

    @Test
    public void testRefConfidenceWithCalls() {
        final RefConfData xxxdata = new RefConfData("ACGTAACCGGTT", 0);
        final int start = xxxdata.getStart();
        final int stop = xxxdata.getEnd();

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();

        for ( int nReads = 0; nReads < 2; nReads++ ) {

            final VariantContext vcStart = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start, Arrays.asList("A", "C"));
            final VariantContext vcEnd = GATKVariantContextUtils.makeFromAlleles("test", "chr1", stop, Arrays.asList("A", "C"));
            final VariantContext vcMiddle = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 2, Arrays.asList("A", "C"));
            final VariantContext vcDel = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 4, Arrays.asList("AAC", "A"));
            final VariantContext vcIns = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 8, Arrays.asList("G", "GCG"));

            final List<VariantContext> allCalls = Arrays.asList(vcStart, vcEnd, vcMiddle, vcDel, vcIns);

            for ( int n = 1; n <= allCalls.size(); n++ ) {
                for ( final List<VariantContext> calls : Utils.makePermutations(allCalls, n, false) ) {
//                    logger.warn("Executing " + n + " " + calls.size());
                    final RefConfData data = new RefConfData("ACGTAACCGGTT", 0);
                    final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
                    for ( int i = 0; i < nReads; i++ ) {
                        data.getActiveRegion().add(data.makeRead(0, data.getRefLength()));
                    }

                    final ReadLikelihoods<Haplotype> likelihoods = HaplotypeCaller.createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

                    final List<Integer> expectedDPs = Collections.nCopies(data.getActiveRegion().getLocation().size(), nReads);
                    final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, genotypingModel, calls);
                    checkReferenceModelResult(data, contexts, expectedDPs, calls);
                }
            }
        }
    }

    private void checkReferenceModelResult(final RefConfData data, final List<VariantContext> contexts, final List<Integer> expectedDPs, final List<VariantContext> calls) {
        Assert.assertNotNull(contexts);

        final GenomeLoc loc = data.getActiveRegion().getExtendedLoc();
        final List<Boolean> seenBP = new ArrayList<>(Collections.nCopies(data.getActiveRegion().getLocation().size(), false));

        for ( int i = 0; i < loc.size(); i++ ) {
            final GenomeLoc curPos = parser.createGenomeLoc(loc.getContig(), loc.getStart() + i);
            final VariantContext call = model.getOverlappingVariantContext(curPos, calls);
            final VariantContext refModel = model.getOverlappingVariantContext(curPos, contexts);

            if ( ! data.getActiveRegion().getLocation().containsP(curPos) ) {
                // part of the extended interval, but not the full interval
                Assert.assertNull(refModel);
                continue;
            }

            if ( call != null ) {
                if (call.isVariant() && refModel.getType() ==  VariantContext.Type.SYMBOLIC ) {
                    //Assert.assertEquals(refModel, call, "Should have found call " + call + " but found " + refModel + " instead");
                    Assert.assertTrue(call.getReference().length() > 1); // must be a deletion.
                    Assert.assertTrue(call.getStart() < refModel.getStart()); // the deletion must not start at the same position
                    Assert.assertEquals(call.getReference().getBaseString().substring(refModel.getStart() - call.getStart(),
                                refModel.getStart() - call.getStart() + 1), refModel.getReference().getBaseString(), "" + data.getRefHap()); // the reference must be the same.
                    Assert.assertTrue(refModel.getGenotype(0).getGQ() <= 0); // No confidence in the reference hom-ref call across the deletion
                    Assert.assertEquals(refModel.getAlleles().size(),2); // the reference and the lonelly <NON_REF>
                    Assert.assertEquals(refModel.getAlleles().get(1), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                } else {
                    Assert.assertEquals(refModel, call, "Should have found call " + call + " but found " + refModel + " instead");
                }

            } else {
                final int expectedDP = expectedDPs.get(curPos.getStart() - data.getActiveRegion().getLocation().getStart());
                Assert.assertEquals(refModel.getStart(), loc.getStart() + i);
                Assert.assertEquals(refModel.getEnd(), loc.getStart() + i);
                Assert.assertFalse(refModel.hasLog10PError());
                Assert.assertEquals(refModel.getAlternateAlleles().size(), 1);
                Assert.assertEquals(refModel.getAlternateAllele(0), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                Assert.assertTrue(refModel.hasGenotype(sample));

                final Genotype g = refModel.getGenotype(sample);
                Assert.assertTrue(g.hasAD());
                Assert.assertTrue(g.hasDP());
                Assert.assertEquals(g.getDP(), expectedDP);
                Assert.assertTrue(g.hasGQ());
                Assert.assertTrue(g.hasPL());
            }

            final VariantContext vc = call == null ? refModel : call;
            if ( curPos.getStart() == vc.getStart() ) {
                for ( int pos = vc.getStart(); pos <= vc.getEnd(); pos++ ) {
                    final int j = pos - data.getActiveRegion().getLocation().getStart();
                    Assert.assertFalse(seenBP.get(j));
                    seenBP.set(j, true);
                }
            }
        }

        for ( int i = 0; i < seenBP.size(); i++ ) {
            Assert.assertEquals((boolean)seenBP.get(i), true);
        }
    }
}