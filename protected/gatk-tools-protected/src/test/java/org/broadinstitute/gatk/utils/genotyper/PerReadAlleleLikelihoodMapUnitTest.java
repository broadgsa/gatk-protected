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

package org.broadinstitute.gatk.utils.genotyper;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.BaseTest;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.Utils;
import java.util.Map;
import java.util.List;
import org.testng.Assert;
import org.testng.annotations.Test;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;

public class PerReadAlleleLikelihoodMapUnitTest extends BaseTest {

    // example genome loc parser for this test, can be deleted if you don't use the reference
    private GenomeLocParser genomeLocParser;

    // example fasta index file, can be deleted if you don't use the reference
    private IndexedFastaSequenceFile seq;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    @Test()
    public void testMultiAlleleWithHomLiks() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GenomeLoc myLocation = genomeLocParser.createGenomeLoc("1", 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte)'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other read the first base is a C

            // set the read's bases and quals
            read.setReadBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(myLocation, reads, 0);
        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        Allele base_C = Allele.create(BaseUtils.Base.C.base);
        Allele base_T = Allele.create(BaseUtils.Base.T.base);

        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        for ( final PileupElement e : pileup ) {
            for ( final Allele allele : Arrays.asList(base_A,base_C,base_T) ) {
                Double likelihood = allele == base_A ? -0.04 : -3.0;
                perReadAlleleLikelihoodMap.add(e,allele,likelihood);
            }
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.size(),pileup.depthOfCoverage());
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().keySet().size(),3);
        Map<Allele,List<GATKSAMRecord>> shouldBeAllA = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(shouldBeAllA.get(base_A).size(),pileup.depthOfCoverage());
        Assert.assertEquals(shouldBeAllA.get(base_C).size(),0);
        Assert.assertEquals(shouldBeAllA.get(base_T).size(),0);
    }


    @Test()
    public void testMultiAlleleWithHetLiks() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GenomeLoc myLocation = genomeLocParser.createGenomeLoc("1", 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte)'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other base is a C

            // set the read's bases and quals
            read.setReadBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(myLocation, reads, 0);
        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        Allele base_C = Allele.create(BaseUtils.Base.C.base);
        Allele base_T = Allele.create(BaseUtils.Base.T.base);

        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        int idx = 0;
        for ( final PileupElement e : pileup ) {
            for ( final Allele allele : Arrays.asList(base_A,base_C,base_T) ) {
                Double likelihood;
                if ( idx % 2 == 0 )
                    likelihood = allele == base_A ? -0.04 : -3.0;
                else
                    likelihood = allele == base_C ? -0.04 : -3.0;
                perReadAlleleLikelihoodMap.add(e,allele,likelihood);
            }
            idx++;
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.size(),pileup.depthOfCoverage());
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().keySet().size(),3);
        Map<Allele,List<GATKSAMRecord>> halfAhalfC = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(halfAhalfC.get(base_A).size(),pileup.depthOfCoverage()/2);
        Assert.assertEquals(halfAhalfC.get(base_C).size(),pileup.depthOfCoverage()/2);
        Assert.assertEquals(halfAhalfC.get(base_T).size(),0);

        // make sure the likelihoods are retrievable

        idx = 0;
        for ( final PileupElement e : pileup ) {
            Assert.assertTrue(perReadAlleleLikelihoodMap.containsPileupElement(e));
            Map<Allele,Double> likelihoods = perReadAlleleLikelihoodMap.getLikelihoodsAssociatedWithPileupElement(e);
            for ( final Allele allele : Arrays.asList(base_A,base_C,base_T) ) {
                Double expLik;
                if ( idx % 2 == 0 )
                    expLik = allele == base_A ? -0.04 : -3.0;
                else
                    expLik = allele == base_C ? -0.04 : -3.0;
                Assert.assertEquals(likelihoods.get(allele),expLik);
            }
            idx++;
        }

        // and test downsampling for good measure

        final List<GATKSAMRecord> excessReads = new LinkedList<GATKSAMRecord>();
        int prevSize = perReadAlleleLikelihoodMap.size();
        for ( int i = 0; i < 10 ; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myExcessRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte)'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other base is a C

            // set the read's bases and quals
            read.setReadBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            for ( final Allele allele : Arrays.asList(base_A,base_C,base_T) ) {
                perReadAlleleLikelihoodMap.add(read,allele,allele==base_A ? -0.04 : -3.0);
            }
            Assert.assertEquals(perReadAlleleLikelihoodMap.size(),1+prevSize);
            prevSize = perReadAlleleLikelihoodMap.size();
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.size(),pileup.depthOfCoverage()+10);
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().get(base_A).size(),60);
        perReadAlleleLikelihoodMap.performPerAlleleDownsampling(0.1);
        Assert.assertEquals(perReadAlleleLikelihoodMap.size(),(int) (0.9*(pileup.depthOfCoverage()+10)));

        Map<Allele,List<GATKSAMRecord>> downsampledStrat = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(downsampledStrat.get(base_A).size(),(int) (pileup.depthOfCoverage()/2) - 1);
        Assert.assertEquals(downsampledStrat.get(base_C).size(),(int) (pileup.depthOfCoverage()/2));
        Assert.assertEquals(downsampledStrat.get(base_T).size(),0);
    }

    @DataProvider(name = "PoorlyModelledReadData")
    public Object[][] makePoorlyModelledReadData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{10, 0.1, false, Arrays.asList(0.0)});
        tests.add(new Object[]{10, 0.1, true, Arrays.asList(-10.0)});
        tests.add(new Object[]{10, 0.1, false, Arrays.asList(0.0, -10.0)});
        tests.add(new Object[]{10, 0.1, true, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.1, false, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.01, true, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.01, false, Arrays.asList(-5.0, -10.0, -3.0)});
        tests.add(new Object[]{100, 0.01, false, Arrays.asList(-5.0, -10.0, -2.0)});
        tests.add(new Object[]{100, 0.01, true, Arrays.asList(-5.0, -10.0, -4.2)});
        tests.add(new Object[]{100, 0.001, true, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.001, false, Arrays.asList(-5.0, -10.0, 0.0)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PoorlyModelledReadData")
    public void testPoorlyModelledRead(final int readLen, final double maxErrorRatePerBase, final boolean expected, final List<Double> log10likelihoods) {
        final byte[] bases = Utils.dupBytes((byte)'A', readLen);
        final byte[] quals = Utils.dupBytes((byte) 40, readLen);

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, readLen + "M");

        final PerReadAlleleLikelihoodMap map = new PerReadAlleleLikelihoodMap();
        final boolean actual = map.readIsPoorlyModelled(read, log10likelihoods, maxErrorRatePerBase);
        Assert.assertEquals(actual, expected);
    }


    @DataProvider(name = "RemovingPoorlyModelledReadData")
    public Object[][] makeRemovingPoorlyModelledReadData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        final int readLen = 10;
        for ( int nReads = 0; nReads < 4; nReads++ ) {
            for ( int nBad = 0; nBad <= nReads; nBad++ ) {
                final int nGood = nReads - nBad;
                tests.add(new Object[]{readLen, nReads, nBad, nGood});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RemovingPoorlyModelledReadData")
    public void testRemovingPoorlyModelledReads(final int readLen, final int nReads, final int nBad, final int nGood) {
        final PerReadAlleleLikelihoodMap map = new PerReadAlleleLikelihoodMap();
        final Set<GATKSAMRecord> goodReads = new HashSet<GATKSAMRecord>();
        final Set<GATKSAMRecord> badReads = new HashSet<GATKSAMRecord>();
        for ( int readI = 0; readI < nReads; readI++ ) {
            final boolean bad = readI < nBad;
            final double likelihood = bad ? -100.0 : 0.0;

            final byte[] bases = Utils.dupBytes((byte)'A', readLen);
            final byte[] quals = Utils.dupBytes((byte) 40, readLen);

            final Allele allele = Allele.create(Utils.dupString("A", readI+1));

            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, readLen + "M");
            read.setReadName("readName" + readI);
            map.add(read, allele, likelihood);
            (bad ? badReads : goodReads).add(read);
        }

        final List<GATKSAMRecord> removedReads = map.filterPoorlyModelledReads(0.01);
        Assert.assertEquals(removedReads.size(), nBad, "nBad " + nBad + " nGood " + nGood);
        Assert.assertEquals(new HashSet<GATKSAMRecord>(removedReads), badReads, "nBad " + nBad + " nGood " + nGood);
        Assert.assertEquals(map.size(), nGood, "nBad " + nBad + " nGood " + nGood);
        Assert.assertTrue(map.getStoredElements().containsAll(goodReads), "nBad " + nBad + " nGood " + nGood);
        Assert.assertEquals(map.getStoredElements().size(), nGood, "nBad " + nBad + " nGood " + nGood);
    }

    @DataProvider(name = "MostLikelyAlleleData")
    public Object[][] makeMostLikelyAlleleData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Allele a = Allele.create("A");
        final Allele c = Allele.create("C");
        final Allele g = Allele.create("G");

        tests.add(new Object[]{Arrays.asList(a), Arrays.asList(Arrays.asList(0.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c), Arrays.asList(Arrays.asList(0.0, -1.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c), Arrays.asList(Arrays.asList(-1.0, 0.0)), c, c});
        tests.add(new Object[]{Arrays.asList(a, c, g), Arrays.asList(Arrays.asList(0.0, 0.0, -10.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c, g), Arrays.asList(Arrays.asList(0.0, 0.0, -10.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c, g),
                Arrays.asList(
                        Arrays.asList(0.0, -10.0, -10.0),
                        Arrays.asList(-100.0, 0.0, -10.0)),
                c, a});
        tests.add(new Object[]{Arrays.asList(a, c, g),
                Arrays.asList(
                        Arrays.asList(0.0, -10.0, -10.0),
                        Arrays.asList(-20.0, 0.0, -100.0)),
                c, a});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MostLikelyAlleleData")
    public void testMostLikelyAllele(final List<Allele> alleles, final List<List<Double>> perReadlikelihoods, final Allele best, final Allele second) {
        final PerReadAlleleLikelihoodMap map = new PerReadAlleleLikelihoodMap();

        for ( int readI = 0; readI < perReadlikelihoods.size(); readI++ ) {
            final List<Double> likelihoods = perReadlikelihoods.get(readI);

            final byte[] bases = Utils.dupBytes((byte)'A', 10);
            final byte[] quals = Utils.dupBytes((byte) 30, 10);
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, quals, "10M");
            read.setReadName("readName" + readI);

            for ( int i = 0; i < alleles.size(); i++ ) {
                final Allele allele = alleles.get(i);
                final double likelihood = likelihoods.get(i);
                map.add(read, allele, likelihood);
            }
        }

        final MostLikelyAllele mla = map.getMostLikelyDiploidAlleles();
        Assert.assertEquals(mla.getMostLikelyAllele(), best);
        Assert.assertEquals(mla.getSecondMostLikelyAllele(), second);
    }
}