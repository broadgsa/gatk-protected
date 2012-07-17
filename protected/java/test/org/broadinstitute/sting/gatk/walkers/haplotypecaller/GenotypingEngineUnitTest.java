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
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
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
        final ArrayList<Haplotype> haplotypes = new ArrayList<Haplotype>();
        haplotypes.add(new Haplotype("AATA".getBytes()));
        haplotypes.add(new Haplotype("AACA".getBytes()));
        haplotypes.add(new Haplotype("CATA".getBytes()));
        haplotypes.add(new Haplotype("CACA".getBytes()));
        final List<Allele> haplotypeAllelesForSample = new ArrayList<Allele>();
        haplotypeAllelesForSample.add( Allele.create("CATA", false) );
        haplotypeAllelesForSample.add( Allele.create("CACA", false) );
        final ArrayList<ArrayList<Haplotype>> alleleMapper = new ArrayList<ArrayList<Haplotype>>();
        ArrayList<Haplotype> Aallele = new ArrayList<Haplotype>();
        Aallele.add(haplotypes.get(0));
        Aallele.add(haplotypes.get(1));
        ArrayList<Haplotype> Callele = new ArrayList<Haplotype>();
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
        final ArrayList<Haplotype> haplotypes = new ArrayList<Haplotype>();
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
        final ArrayList<ArrayList<Haplotype>> alleleMapper = new ArrayList<ArrayList<Haplotype>>();
        ArrayList<Haplotype> Aallele = new ArrayList<Haplotype>();
        Aallele.add(haplotypes.get(0));
        Aallele.add(haplotypes.get(1));
        ArrayList<Haplotype> Callele = new ArrayList<Haplotype>();
        Callele.add(haplotypes.get(2));
        Callele.add(haplotypes.get(3));
        ArrayList<Haplotype> Tallele = new ArrayList<Haplotype>();
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
        HashMap<Integer,Byte> expected;
        GenotypingEngine ge = new GenotypingEngine(false, 0, false);

        public BasicGenotypingTestProvider(String refString, String hapString, HashMap<Integer, Byte> expected) {
            super(BasicGenotypingTestProvider.class, String.format("Haplotype to VCF test: ref = %s, alignment = %s", refString,hapString));
            ref = refString.getBytes();
            hap = hapString.getBytes();
            this.expected = expected;
        }
        
        public HashMap<Integer,VariantContext> calcAlignment() {
            final SWPairwiseAlignment alignment = new SWPairwiseAlignment(ref, hap);
            return ge.generateVCsFromAlignment( alignment.getAlignmentStart2wrt1(), alignment.getCigar(), ref, hap, genomeLocParser.createGenomeLoc("4",1,1+ref.length), "name", 0);
        }
    }

    @DataProvider(name = "BasicGenotypingTestProvider")
    public Object[][] makeBasicGenotypingTests() {

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(21 + contextSize, (byte)'M');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG", "ATCTCGCATCGCGAGCATCGCCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(1, (byte)'M');
            map.put(20, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider("AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
            map.put(2 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            map.put(30 + contextSize, (byte)'D');
            final String context = Utils.dupString('G', contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "ACCTCGCATCGCGAGCATCGTTACTAGCCGATG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            HashMap<Integer, Byte> map = new HashMap<Integer, Byte>();
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
        HashMap<Integer,VariantContext> calculatedMap = cfg.calcAlignment();
        HashMap<Integer,Byte> expectedMap = cfg.expected;
        logger.warn(String.format("Test: %s", cfg.toString()));
        if(!compareVCMaps(calculatedMap, expectedMap)) {
            logger.warn("calc map = " + calculatedMap);
            logger.warn("expected map = " + expectedMap);
        }
        Assert.assertTrue(compareVCMaps(calculatedMap, expectedMap));
    }

    /**
     * Tests that we get the right values from the binomial distribution
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
    
    /**
     * Private function to compare HashMap of VCs, it only checks the types and start locations of the VariantContext
     */
    private boolean compareVCMaps(HashMap<Integer, VariantContext> calc, HashMap<Integer, Byte> expected) {
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