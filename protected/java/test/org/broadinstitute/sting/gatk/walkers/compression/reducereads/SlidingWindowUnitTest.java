/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.UnvalidatingGenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

public class SlidingWindowUnitTest extends BaseTest {

    private static final int variantRegionLength = 1000;
    private static final int globalStartPosition = 1000000;

    private static boolean[] createBitset(final List<FinishedGenomeLoc> locs) {
        final boolean[] variantRegionBitset = new boolean[variantRegionLength];
        for ( FinishedGenomeLoc loc : locs ) {
            final int stop = loc.getStop() - globalStartPosition;
            for ( int i = loc.getStart() - globalStartPosition; i <= stop; i++ )
                variantRegionBitset[i] = true;
        }
        return variantRegionBitset;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    //// This section tests the findVariantRegions() method and related functionality ////
    //////////////////////////////////////////////////////////////////////////////////////

    private static final FinishedGenomeLoc loc90to95 = new FinishedGenomeLoc("1", 0, 1000090, 1000095, false);
    private static final FinishedGenomeLoc loc96to99 = new FinishedGenomeLoc("1", 0, 1000096, 1000099, false);
    private static final FinishedGenomeLoc loc100to110 = new FinishedGenomeLoc("1", 0, 1000100, 1000110, false);
    private static final FinishedGenomeLoc loc999 = new FinishedGenomeLoc("1", 0, 1000999, 1000999, false);

    private class FindVariantRegionsTest {
        public List<FinishedGenomeLoc> locs, expectedResult;
        public boolean[] variantRegionBitset;

        private FindVariantRegionsTest(final List<FinishedGenomeLoc> locs) {
            this.locs = locs;
            this.expectedResult = locs;
            variantRegionBitset = createBitset(locs);
        }

        private FindVariantRegionsTest(final List<FinishedGenomeLoc> locs, final List<FinishedGenomeLoc> expectedResult) {
            this.locs = locs;
            this.expectedResult = expectedResult;
            variantRegionBitset = createBitset(locs);
        }
    }

    @DataProvider(name = "findVariantRegions")
    public Object[][] createFindVariantRegionsData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<FinishedGenomeLoc>asList(loc90to95))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<FinishedGenomeLoc>asList(loc90to95, loc100to110))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<FinishedGenomeLoc>asList(loc90to95, loc96to99, loc100to110), Arrays.<FinishedGenomeLoc>asList(new FinishedGenomeLoc("1", 0, 1000090, 1000110, false)))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<FinishedGenomeLoc>asList(loc90to95, loc999))});
        tests.add(new Object[]{new FindVariantRegionsTest(Arrays.<FinishedGenomeLoc>asList(loc999))});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "findVariantRegions", enabled = true)
    public void testFindVariantRegions(FindVariantRegionsTest test) {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, globalStartPosition);
        final CompressionStash locs = slidingWindow.findVariantRegions(0, variantRegionLength, test.variantRegionBitset, true);
        int index = 0;
        for ( final FinishedGenomeLoc loc : locs ) {
            Assert.assertTrue(loc.equals(test.expectedResult.get(index++)));
        }
    }

    @Test(enabled = true)
    public void testNoClosingRegions() {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, globalStartPosition);
        final CompressionStash locs = slidingWindow.findVariantRegions(0, variantRegionLength, createBitset(Arrays.<FinishedGenomeLoc>asList(loc90to95, loc999)), false);
        Assert.assertEquals(locs.size(), 1);
        Assert.assertEquals(locs.iterator().next(), loc90to95);
    }


    /////////////////////////////////////////////////////////////////////////////
    //// This section tests the markSites() method and related functionality ////
    /////////////////////////////////////////////////////////////////////////////

    @Test(enabled = true)
    public void testMarkedSitesClass() {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, globalStartPosition);
        final SlidingWindow.MarkedSites markedSites = slidingWindow.new MarkedSites();

        markedSites.updateRegion(100, 100);
        Assert.assertEquals(markedSites.getStartLocation(), 100);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 100);

        markedSites.updateRegion(300, 100);
        Assert.assertEquals(markedSites.getStartLocation(), 300);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 100);

        markedSites.getVariantSiteBitSet()[10] = true;
        markedSites.updateRegion(290, 100);
        Assert.assertEquals(markedSites.getStartLocation(), 290);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 100);
        Assert.assertFalse(markedSites.getVariantSiteBitSet()[10]);

        markedSites.getVariantSiteBitSet()[20] = true;
        markedSites.updateRegion(290, 100);
        Assert.assertEquals(markedSites.getStartLocation(), 290);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 100);
        Assert.assertTrue(markedSites.getVariantSiteBitSet()[20]);

        markedSites.updateRegion(300, 100);
        Assert.assertEquals(markedSites.getStartLocation(), 300);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 100);

        markedSites.getVariantSiteBitSet()[95] = true;
        markedSites.updateRegion(390, 20);
        Assert.assertEquals(markedSites.getStartLocation(), 390);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 20);
        Assert.assertTrue(markedSites.getVariantSiteBitSet()[5]);

        markedSites.updateRegion(340, 60);
        Assert.assertEquals(markedSites.getStartLocation(), 340);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 60);

        markedSites.getVariantSiteBitSet()[20] = true;
        markedSites.updateRegion(350, 60);
        Assert.assertEquals(markedSites.getStartLocation(), 350);
        Assert.assertEquals(markedSites.getVariantSiteBitSet().length, 60);
        Assert.assertTrue(markedSites.getVariantSiteBitSet()[10]);
    }

    @Test(enabled = true)
    public void testMarkVariantRegion() {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, globalStartPosition);
        SlidingWindow.MarkedSites markedSites = slidingWindow.new MarkedSites();
        markedSites.updateRegion(100, 100);

        slidingWindow.markVariantRegion(markedSites, 40);
        Assert.assertEquals(countTrueBits(markedSites.getVariantSiteBitSet()), 21);

        slidingWindow.markVariantRegion(markedSites, 5);
        Assert.assertEquals(countTrueBits(markedSites.getVariantSiteBitSet()), 37);

        slidingWindow.markVariantRegion(markedSites, 95);
        Assert.assertEquals(countTrueBits(markedSites.getVariantSiteBitSet()), 52);
    }

    private static int countTrueBits(final boolean[] bitset) {
        int count = 0;
        for ( final boolean bit : bitset ) {
            if ( bit )
                count++;
        }
        return count;
    }


    /////////////////////////////////////////////////////////////////
    //// This section tests the consensus creation functionality ////
    /////////////////////////////////////////////////////////////////

    private static final int readLength = 100;
    private static final int testRegionSize = 1000;
    private final List<GATKSAMRecord> basicReads = new ArrayList<GATKSAMRecord>(20);
    private IndexedFastaSequenceFile seq;
    private SAMFileHeader header;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());

        final int readFrequency = 20;

        basicReads.clear();
        for ( int i = 0; i < testRegionSize; i += readFrequency ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "basicRead" + i, 0, globalStartPosition + i, readLength);
            read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            read.setMappingQuality(30);
            read.setReadNegativeStrandFlag(i % 40 == 20);
            basicReads.add(read);
        }
    }

    private class ConsensusCreationTest {
        public final int expectedNumberOfReads, expectedNumberOfReadsWithHetCompression;
        public final List<GATKSAMRecord> myReads = new ArrayList<GATKSAMRecord>(20);

        private ConsensusCreationTest(final List<GenomeLoc> locs, final boolean readsShouldBeLowQuality, final boolean variantBaseShouldBeLowQuality, final int expectedNumberOfReads, final int expectedNumberOfReadsWithHetCompression) {
            this.expectedNumberOfReads = expectedNumberOfReads;
            this.expectedNumberOfReadsWithHetCompression = expectedNumberOfReadsWithHetCompression;

            // first, add the basic reads to the collection
            myReads.addAll(basicReads);

            // then add the permuted reads
            for ( final GenomeLoc loc : locs )
                myReads.add(createVariantRead(loc, readsShouldBeLowQuality, variantBaseShouldBeLowQuality, CigarOperator.M));
        }

        private ConsensusCreationTest(final List<GenomeLoc> locs, final CigarOperator operator, final int expectedNumberOfReads, final int expectedNumberOfReadsWithHetCompression) {
            this.expectedNumberOfReads = expectedNumberOfReads;
            this.expectedNumberOfReadsWithHetCompression = expectedNumberOfReadsWithHetCompression;

            // first, add the basic reads to the collection
            myReads.addAll(basicReads);

            // then add the permuted reads
            for ( final GenomeLoc loc : locs )
                myReads.add(createVariantRead(loc, false, false, operator));
        }

        private GATKSAMRecord createVariantRead(final GenomeLoc loc, final boolean readShouldBeLowQuality,
                                                final boolean variantBaseShouldBeLowQuality, final CigarOperator operator) {

            final int startPos = loc.getStart() - 50;

            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead" + startPos, 0, startPos, readLength);

            final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
            // create a mismatch if requested
            if ( operator == CigarOperator.M )
                bases[50] = 'C';
            read.setReadBases(bases);

            final byte[] baseQuals = Utils.dupBytes((byte) 30, readLength);
            if ( variantBaseShouldBeLowQuality )
                baseQuals[50] = (byte)10;
            read.setBaseQualities(baseQuals);
            final byte mappingQual = readShouldBeLowQuality ? (byte)10 : (byte)30;
            read.setMappingQuality(mappingQual);

            if ( operator != CigarOperator.M ) {
                final List<CigarElement> elements = new ArrayList<CigarElement>(3);
                elements.add(new CigarElement(operator == CigarOperator.D ? 50 : 51, CigarOperator.M));
                elements.add(new CigarElement(1, operator));
                elements.add(new CigarElement(operator == CigarOperator.D ? 50 : 48, CigarOperator.M));
                read.setCigar(new Cigar(elements));
            }

            return read;
        }
    }

    private static final GenomeLoc loc290 = new UnvalidatingGenomeLoc("1", 0, 1000290, 1000290);
    private static final GenomeLoc loc295 = new UnvalidatingGenomeLoc("1", 0, 1000295, 1000295);
    private static final GenomeLoc loc309 = new UnvalidatingGenomeLoc("1", 0, 1000309, 1000309);
    private static final GenomeLoc loc310 = new UnvalidatingGenomeLoc("1", 0, 1000310, 1000310);
    private static final GenomeLoc loc1100 = new UnvalidatingGenomeLoc("1", 0, 1001100, 1001100);

    @DataProvider(name = "ConsensusCreation")
    public Object[][] createConsensusCreationTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // test high quality reads and bases
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(), false, false, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290), false, false, 9, 5)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc295), false, false, 10, 10)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc309), false, false, 10, 10)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc310), false, false, 11, 11)});

        // test low quality reads
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(), true, false, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290), true, false, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc295), true, false, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc309), true, false, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc310), true, false, 1, 1)});

        // test low quality bases
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(), false, true, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290), false, true, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc295), false, true, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc309), false, true, 1, 1)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc310), false, true, 1, 1)});

        // test mixture
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc1100), true, false, 2, 2)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc1100), false, true, 3, 3)});

        // test I/D operators
        // TODO -- uncomment this test when the deletion bug is fixed!
        // tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290), CigarOperator.D, 9, 5)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc295), CigarOperator.D, 10, 10)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc309), CigarOperator.D, 10, 10)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc310), CigarOperator.D, 11, 11)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290), CigarOperator.I, 9, 9)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc295), CigarOperator.I, 10, 10)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc309), CigarOperator.I, 10, 10)});
        tests.add(new Object[]{new ConsensusCreationTest(Arrays.<GenomeLoc>asList(loc290, loc310), CigarOperator.I, 11, 11)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ConsensusCreation", enabled = true)
    public void testConsensusCreationTest(ConsensusCreationTest test) {
        // test WITHOUT het compression allowed
        SlidingWindow slidingWindow = new SlidingWindow("1", 0, 10, header, new GATKSAMReadGroupRecord("test"), 0, 0.05, 0.05, 20, 20, 100, ReduceReads.DownsampleStrategy.Normal, false, false);
        for ( final GATKSAMRecord read : test.myReads )
            slidingWindow.addRead(read);
        Pair<Set<GATKSAMRecord>, CompressionStash> result = slidingWindow.close();

        Assert.assertEquals(result.getFirst().size(), test.expectedNumberOfReads);

        // test WITH het compression allowed
        slidingWindow = new SlidingWindow("1", 0, 10, header, new GATKSAMReadGroupRecord("test"), 0, 0.05, 0.05, 20, 20, 100, ReduceReads.DownsampleStrategy.Normal, false, true);
        for ( final GATKSAMRecord read : test.myReads )
            slidingWindow.addRead(read);
        result = slidingWindow.close();

        Assert.assertEquals(result.getFirst().size(), test.expectedNumberOfReadsWithHetCompression);
    }


    ///////////////////////////////////////////////////////////
    //// This section tests the downsampling functionality ////
    ///////////////////////////////////////////////////////////

    private class DSTest {
        public final int dcov;

        private DSTest(final int dcov) {
            this.dcov = dcov;
        }
    }

    @DataProvider(name = "Downsampling")
    public Object[][] createDownsamplingTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( int i = 1; i < basicReads.size() + 10; i++ )
            tests.add(new Object[]{new DSTest(i)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "Downsampling", enabled = true)
    public void testDownsamplingTest(DSTest test) {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, 10, header, new GATKSAMReadGroupRecord("test"), 0, 0.05, 0.05, 20, 20, test.dcov, ReduceReads.DownsampleStrategy.Normal, false, false);
        final List<GATKSAMRecord> result = slidingWindow.downsampleVariantRegion(basicReads);

        Assert.assertEquals(result.size(), Math.min(test.dcov, basicReads.size()));
    }


    //////////////////////////////////////////////////////////////
    //// This section tests the consensus base quals accuracy ////
    //////////////////////////////////////////////////////////////

    private class QualsTest {
        public final List<Integer> quals;
        public final List<GATKSAMRecord> myReads = new ArrayList<GATKSAMRecord>(5);

        private QualsTest(final List<Integer> quals) {
            this.quals = quals;
            for ( int i = 0; i < quals.size(); i++ ) {
                final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "basicRead" + i, 0, globalStartPosition, 1);
                read.setReadBases(new byte[]{(byte)'A'});
                read.setBaseQualities(new byte[]{quals.get(i).byteValue()});
                read.setMappingQuality(30);
                myReads.add(read);
            }
        }
    }

    @DataProvider(name = "ConsensusQuals")
    public Object[][] createConsensusQualsData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int[] quals = new int[]{ 0, 5, 10, 15, 20, 30, 40, 50 };

        for ( final int qual1 : quals ) {
            for ( final int qual2 : quals ) {
                for ( final int qual3 : quals ) {
                    tests.add(new Object[]{new QualsTest(Arrays.asList(qual1, qual2, qual3))});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private static final byte minUsableConsensusQual = 10;

    @Test(dataProvider = "ConsensusQuals", enabled = true)
    public void testConsensusQualsTest(QualsTest test) {
        final SlidingWindow slidingWindow = new SlidingWindow("1", 0, 10, header, new GATKSAMReadGroupRecord("test"), 0, 0.05, 0.05, minUsableConsensusQual, 20, 100, ReduceReads.DownsampleStrategy.Normal, false, false);
        for ( final GATKSAMRecord read : test.myReads )
            slidingWindow.addRead(read);
        final Pair<Set<GATKSAMRecord>, CompressionStash> result = slidingWindow.close();

        Assert.assertEquals(result.getFirst().size(), 1);
        final GATKSAMRecord read = result.getFirst().iterator().next();
        final int actualBaseQual = read.getReducedCount(0) * read.getBaseQualities()[0];
        final int expectedBaseQual = qualSum(test.quals);
        Assert.assertEquals(actualBaseQual, expectedBaseQual);
    }

    private static int qualSum(final List<Integer> quals) {
        int goodBases = 0;
        int sum = 0;
        for ( final int qual : quals ) {
            if ( qual >= minUsableConsensusQual ) {
                goodBases++;
                sum += qual;
            }
        }

        // handle a low quality consensus
        if ( sum == 0 ) {
            for ( final int qual : quals ) {
                goodBases++;
                sum += qual;
            }
        }

        return sum - (sum % goodBases);
    }



}
