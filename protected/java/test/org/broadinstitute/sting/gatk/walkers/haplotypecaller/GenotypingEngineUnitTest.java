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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/15/12
 */

import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Unit tests for GenotypingEngine
 */
public class GenotypingEngineUnitTest extends BaseTest {

    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    @Test
    public void testFindHomVarEventAllelesInSample() {
        final List<Allele> eventAlleles = new ArrayList<Allele>();
        eventAlleles.add( Allele.create("A", true) );
        eventAlleles.add( Allele.create("C", false) );
        final List<Allele> haplotypeAlleles = new ArrayList<Allele>();
        haplotypeAlleles.add( Allele.create("AATA", true) );
        haplotypeAlleles.add( Allele.create("AACA", false) );
        haplotypeAlleles.add( Allele.create("CATA", false) );
        haplotypeAlleles.add( Allele.create("CACA", false) );
        final List<Haplotype> haplotypes = new ArrayList<Haplotype>();
        haplotypes.add(new Haplotype("AATA".getBytes()));
        haplotypes.add(new Haplotype("AACA".getBytes()));
        haplotypes.add(new Haplotype("CATA".getBytes()));
        haplotypes.add(new Haplotype("CACA".getBytes()));
        final List<Allele> haplotypeAllelesForSample = new ArrayList<Allele>();
        haplotypeAllelesForSample.add( Allele.create("CATA", false) );
        haplotypeAllelesForSample.add( Allele.create("CACA", false) );
        final List<List<Haplotype>> alleleMapper = new ArrayList<List<Haplotype>>();
        List<Haplotype> Aallele = new ArrayList<Haplotype>();
        Aallele.add(haplotypes.get(0));
        Aallele.add(haplotypes.get(1));
        List<Haplotype> Callele = new ArrayList<Haplotype>();
        Callele.add(haplotypes.get(2));
        Callele.add(haplotypes.get(3));
        alleleMapper.add(Aallele);
        alleleMapper.add(Callele);
        final List<Allele> eventAllelesForSample = new ArrayList<Allele>();
        eventAllelesForSample.add( Allele.create("C", false) );
        eventAllelesForSample.add( Allele.create("C", false) );

        if(!compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper, haplotypes))) {
            logger.warn("calc alleles = " + GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper, haplotypes));
            logger.warn("expected alleles = " + eventAllelesForSample);
        }
        Assert.assertTrue(compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper, haplotypes)));
    }

    @Test
    public void testFindHetEventAllelesInSample() {
        final List<Allele> eventAlleles = new ArrayList<Allele>();
        eventAlleles.add( Allele.create("A", true) );
        eventAlleles.add( Allele.create("C", false) );
        eventAlleles.add( Allele.create("T", false) );
        final List<Allele> haplotypeAlleles = new ArrayList<Allele>();
        haplotypeAlleles.add( Allele.create("AATA", true) );
        haplotypeAlleles.add( Allele.create("AACA", false) );
        haplotypeAlleles.add( Allele.create("CATA", false) );
        haplotypeAlleles.add( Allele.create("CACA", false) );
        haplotypeAlleles.add( Allele.create("TACA", false) );
        haplotypeAlleles.add( Allele.create("TTCA", false) );
        haplotypeAlleles.add( Allele.create("TTTA", false) );
        final List<Haplotype> haplotypes = new ArrayList<Haplotype>();
        haplotypes.add(new Haplotype("AATA".getBytes()));
        haplotypes.add(new Haplotype("AACA".getBytes()));
        haplotypes.add(new Haplotype("CATA".getBytes()));
        haplotypes.add(new Haplotype("CACA".getBytes()));
        haplotypes.add(new Haplotype("TACA".getBytes()));
        haplotypes.add(new Haplotype("TTCA".getBytes()));
        haplotypes.add(new Haplotype("TTTA".getBytes()));
        final List<Allele> haplotypeAllelesForSample = new ArrayList<Allele>();
        haplotypeAllelesForSample.add( Allele.create("TTTA", false) );
        haplotypeAllelesForSample.add( Allele.create("AATA", true) );
        final List<List<Haplotype>> alleleMapper = new ArrayList<List<Haplotype>>();
        List<Haplotype> Aallele = new ArrayList<Haplotype>();
        Aallele.add(haplotypes.get(0));
        Aallele.add(haplotypes.get(1));
        List<Haplotype> Callele = new ArrayList<Haplotype>();
        Callele.add(haplotypes.get(2));
        Callele.add(haplotypes.get(3));
        List<Haplotype> Tallele = new ArrayList<Haplotype>();
        Tallele.add(haplotypes.get(4));
        Tallele.add(haplotypes.get(5));
        Tallele.add(haplotypes.get(6));
        alleleMapper.add(Aallele);
        alleleMapper.add(Callele);
        alleleMapper.add(Tallele);
        final List<Allele> eventAllelesForSample = new ArrayList<Allele>();
        eventAllelesForSample.add( Allele.create("A", true) );
        eventAllelesForSample.add( Allele.create("T", false) );

        if(!compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper, haplotypes))) {
            logger.warn("calc alleles = " + GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper, haplotypes));
            logger.warn("expected alleles = " + eventAllelesForSample);
        }
        Assert.assertTrue(compareAlleleLists(eventAllelesForSample, GenotypingEngine.findEventAllelesInSample(eventAlleles, haplotypeAlleles, haplotypeAllelesForSample, alleleMapper, haplotypes)));
    }

    private boolean compareAlleleLists(List<Allele> l1, List<Allele> l2) {
        if( l1.size() != l2.size() ) {
            return false; // sanity check
        }

        for( int i=0; i < l1.size(); i++ ){
            if ( !l2.contains(l1.get(i)) )
                return false;
        }
        return true;
    }

    
    private class BasicGenotypingTestProvider extends TestDataProvider {
        byte[] ref;
        byte[] hap;
        Map<Integer,Byte> expected;

        public BasicGenotypingTestProvider(String refString, String hapString, Map<Integer, Byte> expected) {
            super(BasicGenotypingTestProvider.class, String.format("Haplotype to VCF test: ref = %s, alignment = %s", refString,hapString));
            ref = refString.getBytes();
            hap = hapString.getBytes();
            this.expected = expected;
        }
        
        public Map<Integer,VariantContext> calcAlignment() {
            final SWPairwiseAlignment alignment = new SWPairwiseAlignment(ref, hap);
            return GenotypingEngine.generateVCsFromAlignment( new Haplotype(hap), alignment.getAlignmentStart2wrt1(), alignment.getCigar(), ref, hap, genomeLocParser.createGenomeLoc("4",1,1+ref.length), "name");
        }
    }

    @DataProvider(name = "BasicGenotypingTestProvider")
    public Object[][] makeBasicGenotypingTests() {

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(21 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG", "ATCTCGCATCGCGAGCATCGCCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1, (byte)'M');
            map.put(20, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider("AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            map.put(30 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "ACCTCGCATCGCGAGCATCGTTACTAGCCGATG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            map.put(28 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCCATAG", map);
        }

        return BasicGenotypingTestProvider.getTests(BasicGenotypingTestProvider.class);
    }
    
    @Test(dataProvider = "BasicGenotypingTestProvider", enabled = true)
    public void testHaplotypeToVCF(BasicGenotypingTestProvider cfg) {
        Map<Integer,VariantContext> calculatedMap = cfg.calcAlignment();
        Map<Integer,Byte> expectedMap = cfg.expected;
        logger.warn(String.format("Test: %s", cfg.toString()));
        if(!compareVCMaps(calculatedMap, expectedMap)) {
            logger.warn("calc map = " + calculatedMap);
            logger.warn("expected map = " + expectedMap);
        }
        Assert.assertTrue(compareVCMaps(calculatedMap, expectedMap));
    }

    /**
     * Tests that we get the right values from the R^2 calculation
     */
    @Test
    public void testCalculateR2LD() {
        logger.warn("Executing testCalculateR2LD");

        Assert.assertEquals(GenotypingEngine.calculateR2LD(1,1,1,1), 0.0, 0.00001);
        Assert.assertEquals(GenotypingEngine.calculateR2LD(100,100,100,100), 0.0, 0.00001);
        Assert.assertEquals(GenotypingEngine.calculateR2LD(1,0,0,1), 1.0, 0.00001);
        Assert.assertEquals(GenotypingEngine.calculateR2LD(100,0,0,100), 1.0, 0.00001);
        Assert.assertEquals(GenotypingEngine.calculateR2LD(1,2,3,4), (0.1 - 0.12) * (0.1 - 0.12) / (0.3 * 0.7 * 0.4 * 0.6), 0.00001);
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
        VariantContext mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // SNP + ref + SNP = MNP with ref base gap
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","G").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","GCG").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + SNP
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TAAAAA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","G").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","TAAAAACG").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // SNP + insertion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","CAAAAA").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","GCCAAAAA").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // deletion + SNP
        thisVC = new VariantContextBuilder().loc("2", 1703, 1704).alleles("TC","T").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","G").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","TG").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // SNP + deletion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","G").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1706).alleles("TCCG","GCC").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + deletion = MNP
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1704, 1706).alleles("CCG","ACC").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + deletion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TAAAAA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1706).alleles("TCCG","TAAAAACC").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // insertion + insertion
        thisVC = new VariantContextBuilder().loc("2", 1703, 1703).alleles("T","TA").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1705).alleles("C","CA").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1705).alleles("TCC","TACCA").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // deletion + deletion
        thisVC = new VariantContextBuilder().loc("2", 1701, 1702).alleles("AT","A").make();
        nextVC = new VariantContextBuilder().loc("2", 1705, 1706).alleles("CG","C").make();
        truthVC = new VariantContextBuilder().loc("2", 1701, 1706).alleles("ATTCCG","ATCC").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // deletion + insertion (abutting)
        thisVC = new VariantContextBuilder().loc("2", 1701, 1702).alleles("AT","A").make();
        nextVC = new VariantContextBuilder().loc("2", 1702, 1702).alleles("T","GCGCGC").make();
        truthVC = new VariantContextBuilder().loc("2", 1701, 1702).alleles("AT","AGCGCGC").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());

        // complex + complex
        thisVC = new VariantContextBuilder().loc("2", 1703, 1704).alleles("TC","AAA").make();
        nextVC = new VariantContextBuilder().loc("2", 1706, 1707).alleles("GG","AC").make();
        truthVC = new VariantContextBuilder().loc("2", 1703, 1707).alleles("TCCGG","AAACAC").source("merged").make();
        mergedVC = GenotypingEngine.createMergedVariantContext(thisVC, nextVC, ref, refLoc);
        logger.warn(truthVC + " == " + mergedVC);
        Assert.assertTrue(truthVC.hasSameAllelesAs(mergedVC));
        Assert.assertEquals(truthVC.getStart(), mergedVC.getStart());
        Assert.assertEquals(truthVC.getEnd(), mergedVC.getEnd());
    }
    
    /**
     * Private function to compare Map of VCs, it only checks the types and start locations of the VariantContext
     */
    private boolean compareVCMaps(Map<Integer, VariantContext> calc, Map<Integer, Byte> expected) {
        if( !calc.keySet().equals(expected.keySet()) ) { return false; } // sanity check
        for( Integer loc : expected.keySet() ) {
            Byte type = expected.get(loc);
            switch( type ) {
                case 'I':
                    if( !calc.get(loc).isSimpleInsertion() ) { return false; }
                    break;
                case 'D':
                    if( !calc.get(loc).isSimpleDeletion() ) { return false; }
                    break;
                case 'M':
                    if( !(calc.get(loc).isMNP() || calc.get(loc).isSNP()) ) { return false; }
                    break;
                default:
                    return false;
            }
        }
        return true;
    }
}