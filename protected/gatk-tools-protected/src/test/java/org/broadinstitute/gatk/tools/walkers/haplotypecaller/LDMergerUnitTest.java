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

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.gatk.utils.haplotype.EventMap;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

public class LDMergerUnitTest extends BaseTest {
    LDMerger merger;
    GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        genomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(new File(b37KGReference)));
    }

    @BeforeMethod
    public void setUp() throws Exception {
        merger = new LDMerger();
    }

    @Test
    public void testCreateMergedVariantContext() {
        logger.warn("Executing testCreateMergedVariantContext");

        final byte[] ref = "AATTCCGGAATTCCGGAATT".getBytes();
        final GenomeLoc refLoc = genomeLocParser.createGenomeLoc("2", 1700, 1700 + ref.length);

        // SNP + SNP = simple MNP
        VariantContext thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        VariantContext nextVC = new VariantContextBuilder().loc("2", 1704, 1704).alleles("C","G").make();
        VariantContext truthVC = new VariantContextBuilder().loc("2", 1703, 1704).alleles("TC","GG").source("merged").make();
        VariantContext mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // SNP + ref + SNP = MNP with ref base gap
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","G").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","GCG").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + SNP
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TAAAAA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","G").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","TAAAAACG").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // SNP + insertion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","CAAAAA").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","GCCAAAAA").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // deletion + SNP
        thisVC = new VariantContextBuilder().loc("2", 1703, 1704).alleles("TC","T").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","G").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","TG").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // SNP + deletion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1706).alleles("TCCG","GCC").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + deletion = MNP
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1704, 1706).alleles("CCG","ACC").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + deletion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TAAAAA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1706).alleles("TCCG","TAAAAACC").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + insertion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","CA").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","TACCA").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // deletion + deletion
        thisVC = new VariantContextBuilder().loc("2", 1701, 1702).alleles("AT","A").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1701, 1706).alleles("ATTCCG","ATCC").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // deletion + insertion (abutting)
        thisVC = new VariantContextBuilder().loc("2", 1701, 1702).alleles("AT","A").make();
        nextVC = new VariantContextBuilder().loc("2", 1702, 1702).alleles("T","GCGCGC").make();
        truthVC = new VariantContextBuilder().loc("2", 1701, 1702).alleles("AT","AGCGCGC").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // complex + complex
        thisVC = new VariantContextBuilder().loc("2", 1703, 1704).alleles("TC","AAA").make();
        nextVC = new VariantContextBuilder().loc("2", 1706, 1707).alleles("GG","AC").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1707).alleles("TCCGG","AAACAC").source("merged").make();
        mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());
    }

    @Test
    public void testInsertionDeletionBecomingNullAllele() {
        final byte[] ref = "CAAA".getBytes();
        final GenomeLoc refLoc = genomeLocParser.createGenomeLoc("2", 1700, 1700 + ref.length);

        // insertion + deletion results in a null allele, should return false
        final VariantContext thisVC = new VariantContextBuilder().loc("2", 1700, 1701).alleles("CA","C").make();
        final VariantContext nextVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("A","AA").make();
        final VariantContext mergedVC = merger.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        Assert.assertNull(mergedVC,  "Insertion deletion becoming a null allele should return a null variant context");
    }

    /**
     * Just returns a given R2 value for testing
     */
    private static class MockLDCalculator extends HaplotypeLDCalculator {
        private final double R2;

        private MockLDCalculator(double r2) {
            R2 = r2;
        }

        @Override
        protected double computeProbOfBeingPhased(VariantContext first, VariantContext second) {
            return R2;
        }
    }

    @DataProvider(name = "R2MergerData")
    public Object[][] makeR2MergerData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        final double thres = LDMerger.MERGE_EVENTS_PROB_PHASED_THRESHOLD;
        for ( final double r2 : Arrays.asList(0.0, thres - 0.01, thres + 0.01, 1.0) ) {
            tests.add(new Object[]{"ACGT", "CCGC", 2, "4M", "ACGT", "CCGC", r2, r2 >= thres});
            tests.add(new Object[]{"ACGT", "AGGC", 2, "4M", "CGT", "GGC", r2, r2 >= thres});
            tests.add(new Object[]{"ACGT", "ACCC", 2, "4M", "GT", "CC", r2, r2 >= thres});
            tests.add(new Object[]{"ACGT", "ACCGTT", 2, "2M1I1M1I1M", "CG", "CCGT", r2, r2 >= thres});
            tests.add(new Object[]{"ACGT", "AGCT", 2, "4M", "CG", "GC", r2, r2 >= thres});
            tests.add(new Object[]{"ACAGT", "AAGC", 2, "1M1D3M", "ACAGT", "AAGC", r2, r2 >= thres});
            tests.add(new Object[]{"ACAGT", "AAT", 2, "1M1D1M1D1M", "ACAG", "AA", r2, r2 >= thres});

            // cannot be merged -- only 1 event
            tests.add(new Object[]{"AAA", "ACA", 1, "3M", null, null, r2, false});

            final int dist = LDMerger.MAX_DISTANCE_BETWEEN_SNPS_TO_MERGE + 2;
            tests.add(new Object[]{Utils.dupString("A", dist), "C" + Utils.dupString("A", dist - 2) + "C", 2, dist + "M", null, null, r2, false});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "R2MergerData")
    public void testR2Merger(final String refS, final String hapS, int nEvents, final String cigar, final String expectedMergedRef, final String expectedMergedAlt, final double r2, final boolean expectMerge) {
        final Haplotype ref = new Haplotype(refS.getBytes(), true, 0, TextCigarCodec.decode(refS.length() + "M"));
        final Haplotype hap = new Haplotype(hapS.getBytes(), false, 0, TextCigarCodec.decode(cigar));
        final GenomeLoc loc = new UnvalidatingGenomeLoc("1", 0, 1, ref.length());

        final List<Haplotype> haplotypes = Arrays.asList(ref, hap);
        final TreeSet<Integer> vcStarts = EventMap.buildEventMapsForHaplotypes(haplotypes, ref.getBases(), loc, false);
        final MockLDCalculator r2Calc = new MockLDCalculator(r2);

        Assert.assertEquals(vcStarts.size(), nEvents);
        final boolean merged = merger.mergeConsecutiveEventsBasedOnLDOnce(haplotypes, r2Calc, 1, vcStarts, ref.getBases(), loc);
        Assert.assertEquals(merged, expectMerge);
        Assert.assertEquals(vcStarts.size(), expectMerge ? 1 : nEvents);
        if ( expectMerge ) {
            final VariantContext vc = hap.getEventMap().getVariantContexts().iterator().next();
            Assert.assertTrue(vc.isBiallelic());
            Assert.assertEquals(vc.getReference().getDisplayString(), expectedMergedRef);
            Assert.assertEquals(vc.getAlternateAllele(0).getDisplayString(), expectedMergedAlt);
        }
    }

    @Test
    public void testR2MergerWithThirdHapWithoutEvent() {
        final String refS = "ACGT";
        final String hapS = "CCGA";
        final String cigar = "4M";
        final Haplotype ref = new Haplotype(refS.getBytes(), true, 0, TextCigarCodec.decode(refS.length() + "M"));
        final Haplotype hap1 = new Haplotype(hapS.getBytes(), false, 0, TextCigarCodec.decode(cigar));
        final Haplotype hap2 = new Haplotype("ACGA".getBytes(), false, 0, TextCigarCodec.decode(cigar));
        final GenomeLoc loc = new UnvalidatingGenomeLoc("1", 0, 1, ref.length());

        final List<Haplotype> haplotypes = Arrays.asList(ref, hap1, hap2);
        final TreeSet<Integer> vcStarts = EventMap.buildEventMapsForHaplotypes(haplotypes, ref.getBases(), loc, false);
        final MockLDCalculator r2Calc = new MockLDCalculator(1.0);

        Assert.assertEquals(vcStarts.size(), 2);
        final boolean merged = merger.mergeConsecutiveEventsBasedOnLDOnce(haplotypes, r2Calc, 1, vcStarts, ref.getBases(), loc);
        Assert.assertEquals(merged, true);
        Assert.assertEquals(vcStarts.size(), 1);

        final VariantContext vc = hap1.getEventMap().getVariantContexts().iterator().next();
        Assert.assertTrue(vc.isBiallelic());
        Assert.assertEquals(vc.getReference().getDisplayString(), "ACGT");
        Assert.assertEquals(vc.getAlternateAllele(0).getDisplayString(), "CCGA");

        Assert.assertEquals(hap2.getEventMap().size(), 0);
    }

    @Test
    public void testR2MergerWithMultipleAllelesAtSites() {
        final String refS = "ACGT";
        final String hapS = "TCGA";
        final String cigar = "4M";
        final Haplotype ref = new Haplotype(refS.getBytes(), true, 0, TextCigarCodec.decode(refS.length() + "M"));
        final Haplotype hap1 = new Haplotype(hapS.getBytes(), false, 0, TextCigarCodec.decode(cigar));

        final GenomeLoc loc = new UnvalidatingGenomeLoc("1", 0, 1, ref.length());
        for (final String hap2S : Arrays.asList("GCGA", "TCGG")) {
            final Haplotype hap2 = new Haplotype(hap2S.getBytes(), false, 0, TextCigarCodec.decode(cigar));

            final List<Haplotype> haplotypes = Arrays.asList(ref, hap1, hap2);
            final TreeSet<Integer> vcStarts = EventMap.buildEventMapsForHaplotypes(haplotypes, ref.getBases(), loc, false);
            final MockLDCalculator r2Calc = new MockLDCalculator(1.0);

            Assert.assertEquals(vcStarts.size(), 2);
            final boolean merged = merger.mergeConsecutiveEventsBasedOnLDOnce(haplotypes, r2Calc, 1, vcStarts, ref.getBases(), loc);
            Assert.assertEquals(merged, false);
            Assert.assertEquals(vcStarts.size(), 2);
        }
    }
}