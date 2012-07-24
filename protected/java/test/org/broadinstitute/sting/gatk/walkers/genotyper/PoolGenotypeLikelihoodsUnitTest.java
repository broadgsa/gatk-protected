/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintStream;
import java.util.*;


public class PoolGenotypeLikelihoodsUnitTest {

    final UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();
    final Logger logger = Logger.getLogger(Walker.class);
    private static final boolean VERBOSE = false;
    private static final boolean SIMULATE_NOISY_PILEUP = false;
    private static final int NUM_SIMULATED_OBS = 10;

    @Test
    public void testStoringLikelihoodElements() {


        // basic test storing a given PL vector in a PoolGenotypeLikelihoods object and then retrieving it back

        int ploidy = 20;
        int numAlleles = 4;
        int res = GenotypeLikelihoods.numLikelihoods(numAlleles, ploidy);
        //       System.out.format("Alt Alleles: %d, Ploidy: %d, #Likelihoods: %d\n", numAltAlleles, ploidy, res);

        List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create("T",true));
        alleles.add(Allele.create("C",false));
        alleles.add(Allele.create("A",false));
        alleles.add(Allele.create("G",false));

        double[] gls = new double[res];

        for (int k=0; k < gls.length; k++)
            gls[k]= (double)k;

        PoolGenotypeLikelihoods gl = new PoolSNPGenotypeLikelihoods(alleles, gls,ploidy, null, false,true);
        double[] glnew = gl.getLikelihoods();

        Assert.assertEquals(gls, glnew);
    }

    @Test
    public void testElementStorageCache() {
        // compare cached element storage with compuationally hard-coded iterative computation

        for (int ploidy = 2; ploidy < 10; ploidy++) {
            for (int nAlleles = 2; nAlleles < 10; nAlleles++)
                Assert.assertEquals(PoolGenotypeLikelihoods.getNumLikelihoodElements(nAlleles,ploidy),
                        GenotypeLikelihoods.numLikelihoods(nAlleles, ploidy));
        }

    }

    @Test
    public void testVectorToLinearIndex() {

        // create iterator, compare linear index given by iterator with closed form function
        int numAlleles = 4;
        int ploidy = 2;
        PoolGenotypeLikelihoods.SumIterator iterator = new PoolGenotypeLikelihoods.SumIterator(numAlleles, ploidy);

        while(iterator.hasNext()) {
            System.out.format("\n%d:",iterator.getLinearIndex());
            int[] a =  iterator.getCurrentVector();
            for (int aa: a)
                System.out.format("%d ",aa);


            int computedIdx = PoolGenotypeLikelihoods.getLinearIndex(a, numAlleles, ploidy);
            System.out.format("Computed idx = %d\n",computedIdx);
            iterator.next();
        }

    }
    @Test
    public void testSubsetToAlleles() {

        int ploidy = 2;
        int numAlleles = 4;
        int res = GenotypeLikelihoods.numLikelihoods(numAlleles, ploidy);
        //       System.out.format("Alt Alleles: %d, Ploidy: %d, #Likelihoods: %d\n", numAltAlleles, ploidy, res);

        List<Allele> originalAlleles = new ArrayList<Allele>();
        originalAlleles.add(Allele.create("T",true));
        originalAlleles.add(Allele.create("C",false));
        originalAlleles.add(Allele.create("A",false));
        originalAlleles.add(Allele.create("G",false));

        double[] oldLikelihoods = new double[res];

        for (int k=0; k < oldLikelihoods.length; k++)
            oldLikelihoods[k]= (double)k;

        List<Allele> allelesToSubset = new ArrayList<Allele>();
        allelesToSubset.add(Allele.create("A",false));
        allelesToSubset.add(Allele.create("C",false));

        double[] newGLs = PoolGenotypeLikelihoods.subsetToAlleles(oldLikelihoods, ploidy,
                originalAlleles, allelesToSubset);


        /*
            For P=2, N=4, default iteration order:
                0:2 0 0 0
                1:1 1 0 0
                2:0 2 0 0
                3:1 0 1 0
                4:0 1 1 0
                5:0 0 2 0
                6:1 0 0 1
                7:0 1 0 1
                8:0 0 1 1
                9:0 0 0 2

            For P=2,N=2, iteration order is:
                0:2 0
                1:1 1
                2:0 2

            From first list, if we're extracting alleles 2 and 1, we need all elements that have zero at positions 0 and 3.
            These are only elements {2,4,5}. Since test is flipping alleles 2 and 1, order is reversed.
  */
        Assert.assertEquals(newGLs,new double[]{5.0,4.0,2.0});
    }
    @Test
    public void testIndexIterator() {
        int[] seed = new int[]{1,2,3,4};
        PoolGenotypeLikelihoods.SumIterator iterator = runIterator(seed,-1);
        // Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),prod(seed)-1);

        seed = new int[]{1,0,1,1};
        iterator = runIterator(seed,-1);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),prod(seed)-1);

        seed = new int[]{5};
        iterator = runIterator(seed,-1);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),prod(seed)-1);

        // Diploid, # alleles = 4
        seed = new int[]{2,2,2,2};
        iterator = runIterator(seed,2);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),9);

        // Diploid, # alleles = 2
        seed = new int[]{2,2};
        iterator = runIterator(seed,2);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),2);

        // Diploid, # alleles = 3
        seed = new int[]{2,2,2};
        iterator = runIterator(seed,2);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),5);

        // Triploid, # alleles = 2
        seed = new int[]{3,3};
        iterator = runIterator(seed,3);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),3);
        // Triploid, # alleles = 3
        seed = new int[]{3,3,3};
        iterator = runIterator(seed,3);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),9);

        // Triploid, # alleles = 4
        seed = new int[]{3,3,3,3};
        iterator = runIterator(seed,3);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),19);

        // 8-ploid, # alleles = 6
        seed = new int[]{8,8,8,8,8,8};
        iterator = runIterator(seed,8);
        //  Assert.assertTrue(compareIntArrays(iterator.getCurrentVector(), seed));
        Assert.assertEquals(iterator.getLinearIndex(),1286);


    }

    private PoolGenotypeLikelihoods.SumIterator runIterator(int[] seed, int restrictSumTo) {
        PoolGenotypeLikelihoods.SumIterator iterator = new PoolGenotypeLikelihoods.SumIterator(seed, restrictSumTo);

        while(iterator.hasNext()) {
            int[] a =  iterator.getCurrentVector();
            int idx = PoolGenotypeLikelihoods.getLinearIndex(a, a.length, restrictSumTo);
            if (VERBOSE)   {
                System.out.format("%d:",iterator.getLinearIndex());
                for (int i=0; i < seed.length; i++)
                    System.out.format("%d ",a[i]);
                System.out.format(" LI:%d\n", idx);
            }
            iterator.next();
        }

        return iterator;

    }

    private static int prod(int[] x) {
        int prod = 1;
        for (int xx : x) {
            prod *= (1+xx);
        }
        return prod;
    }

    @Test
    public void testErrorModel() {
        final ArtificialReadPileupTestProvider refPileupTestProvider = new ArtificialReadPileupTestProvider(1,"ref");
        final byte minQ = 5;
        final byte maxQ = 40;
        final byte refByte = refPileupTestProvider.getRefByte();
        final byte altByte = refByte == (byte)'T'? (byte) 'C': (byte)'T';
        final String refSampleName = refPileupTestProvider.getSampleNames().get(0);
        final List<Allele> trueAlleles = new ArrayList<Allele>();
        trueAlleles.add(Allele.create(refByte, true));

        final VariantContext refVC = new VariantContextBuilder("test","chr1",5, 5,
                trueAlleles).genotypes(GenotypeBuilder.create(refSampleName, trueAlleles)).make();
        final int[] matchArray = {95, 995, 9995, 10000};
        final int[] mismatchArray = {1,5,10,20};
        if (VERBOSE) System.out.println("Running SNP error model test");

        for (int matches: matchArray) {
            for (int mismatches: mismatchArray) {
                // get artificial alignment context for ref sample - no noise
                Map<String,AlignmentContext> refContext = refPileupTestProvider.getAlignmentContextFromAlleles(0, new String(new byte[]{altByte}), new int[]{matches, mismatches}, false, 30);
                final ReadBackedPileup refPileup = refContext.get(refSampleName).getBasePileup();
                final ErrorModel emodel = new ErrorModel(minQ,maxQ, (byte)20, refPileup, refVC, 0.0);
                final double[] errorVec = emodel.getErrorModelVector().getProbabilityVector();

                final double mlEst = -10.0*Math.log10((double)mismatches/(double)(matches+mismatches));
                final int peakIdx = (int)Math.round(mlEst);
                if (VERBOSE) System.out.format("Matches:%d Mismatches:%d maxV:%d peakIdx:%d\n",matches, mismatches, MathUtils.maxElementIndex(errorVec),peakIdx);
                Assert.assertEquals(MathUtils.maxElementIndex(errorVec),peakIdx);

            }
        }


    }

    @Test
    public void testIndelErrorModel() {
        final ArtificialReadPileupTestProvider refPileupTestProvider = new ArtificialReadPileupTestProvider(1,"ref");
        final byte minQ = 5;
        final byte maxQ = 40;
        final byte refByte = refPileupTestProvider.getRefByte();
        final String altBases = "TCA";
        final String refSampleName = refPileupTestProvider.getSampleNames().get(0);
        final List<Allele> trueAlleles = new ArrayList<Allele>();
        trueAlleles.add(Allele.create(Allele.NULL_ALLELE_STRING, true));
        trueAlleles.add(Allele.create("TC", false));

        final String fw = new String(refPileupTestProvider.getReferenceContext().getForwardBases());
        final VariantContext refInsertionVC = new VariantContextBuilder("test","chr1",refPileupTestProvider.getReferenceContext().getLocus().getStart(),
                refPileupTestProvider.getReferenceContext().getLocus().getStart(), trueAlleles).
                genotypes(GenotypeBuilder.create(refSampleName, trueAlleles)).referenceBaseForIndel(refByte).make();


        final int[] matchArray = {95, 995, 9995, 10000};
        final int[] mismatchArray = {1,5,10,20};

        if (VERBOSE) System.out.println("Running indel error model test");
        for (int matches: matchArray) {
            for (int mismatches: mismatchArray) {
                // get artificial alignment context for ref sample - no noise
                // CASE 1: Test HET insertion
                // Ref sample has TC insertion but pileup will have TCA inserted instead to test mismatches
                Map<String,AlignmentContext> refContext = refPileupTestProvider.getAlignmentContextFromAlleles(altBases.length(), altBases, new int[]{matches, mismatches}, false, 30);
                final ReadBackedPileup refPileup = refContext.get(refSampleName).getBasePileup();
                final ErrorModel emodel = new ErrorModel(minQ,maxQ, (byte)20, refPileup, refInsertionVC, 0.0);
                final double[] errorVec = emodel.getErrorModelVector().getProbabilityVector();

                final double mlEst = -10.0*Math.log10((double)mismatches/(double)(matches+mismatches));
                final int peakIdx = (int)Math.round(mlEst);
                if (VERBOSE) System.out.format("Matches:%d Mismatches:%d peakIdx:%d\n",matches, mismatches, peakIdx);
                Assert.assertEquals(MathUtils.maxElementIndex(errorVec),peakIdx);

                // CASE 2: Test HET deletion

            }
        }

        // create deletion VC
        final int delLength = 4;
        final List<Allele> delAlleles = new ArrayList<Allele>();
        delAlleles.add(Allele.create(fw.substring(1,delLength+1), true));
        delAlleles.add(Allele.create(Allele.NULL_ALLELE_STRING, false));

        final VariantContext refDeletionVC =  new VariantContextBuilder("test","chr1",refPileupTestProvider.getReferenceContext().getLocus().getStart(),
                refPileupTestProvider.getReferenceContext().getLocus().getStart()+delLength, delAlleles).
                genotypes(GenotypeBuilder.create(refSampleName, delAlleles)).referenceBaseForIndel(refByte).make();

        for (int matches: matchArray) {
            for (int mismatches: mismatchArray) {
                // get artificial alignment context for ref sample - no noise
                // CASE 1: Test HET deletion
                // Ref sample has 4bp deletion but pileup will have 3 bp deletion instead to test mismatches
                Map<String,AlignmentContext> refContext = refPileupTestProvider.getAlignmentContextFromAlleles(-delLength+1, altBases, new int[]{matches, mismatches}, false, 30);
                final ReadBackedPileup refPileup = refContext.get(refSampleName).getBasePileup();
                final ErrorModel emodel = new ErrorModel(minQ,maxQ, (byte)20, refPileup, refDeletionVC, 0.0);
                final double[] errorVec = emodel.getErrorModelVector().getProbabilityVector();

                final double mlEst = -10.0*Math.log10((double)mismatches/(double)(matches+mismatches));
                final int peakIdx = (int)Math.round(mlEst);
                if (VERBOSE) System.out.format("Matches:%d Mismatches:%d peakIdx:%d\n",matches, mismatches, peakIdx);
                Assert.assertEquals(MathUtils.maxElementIndex(errorVec),peakIdx);

                // CASE 2: Test HET deletion

            }
        }

    }

    @Test
    public void testAddPileupToPoolGL() {

        // dummy error model - Q=infinity FAPP so that there's no source of uncertainty
        final double[] emv = new double[SAMUtils.MAX_PHRED_SCORE+1];
        
        // error rate for noisy tests
        final int PHRED_SITE_ERROR_RATE = 20;

        Arrays.fill(emv, Double.NEGATIVE_INFINITY);
        emv[SAMUtils.MAX_PHRED_SCORE] = 0;

        final int numSamples = 1;

        // have a high quality site say Q40 site, and create artificial pileups for one single sample, at coverage N, with given
        // true pool AC = x.

        final ArtificialReadPileupTestProvider readPileupTestProvider = new ArtificialReadPileupTestProvider(numSamples,"sample", (byte)SAMUtils.MAX_PHRED_SCORE);
        final ErrorModel noiselessErrorModel = new ErrorModel(emv);

        final double[] emverr = new double[SAMUtils.MAX_PHRED_SCORE+1];
        Arrays.fill(emverr, Double.NEGATIVE_INFINITY);
        emverr[PHRED_SITE_ERROR_RATE] = 0;
        final ErrorModel Q30ErrorModel = new ErrorModel(emverr);


        final int eventLength = 0; // test snp only
        final byte refByte = readPileupTestProvider.getRefByte();
        final byte altByte = refByte == (byte)'T'? (byte) 'C': (byte)'T';

        final int refIdx = BaseUtils.simpleBaseToBaseIndex(refByte);
        final int altIdx = BaseUtils.simpleBaseToBaseIndex(altByte);

        final List<Allele> allAlleles = new ArrayList<Allele>();  // this contains only ref Allele up to now
        final Set<String> laneIDs = new TreeSet<String>();
        laneIDs.add(GenotypeLikelihoodsCalculationModel.DUMMY_LANE);

        final HashMap<String, ErrorModel> noiselessErrorModels = new HashMap<String, ErrorModel>();

        // build per-lane error model for all lanes present in ref sample
        for (String laneID : laneIDs)
            noiselessErrorModels.put(laneID, noiselessErrorModel);

        final HashMap<String, ErrorModel> noisyErrorModels = new HashMap<String, ErrorModel>();

        // build per-lane error model for all lanes present in ref sample
        for (String laneID : laneIDs)
            noisyErrorModels.put(laneID, Q30ErrorModel);

        for (byte b: BaseUtils.BASES) {
            if (refByte == b)
                allAlleles.add(Allele.create(b,true));
            else
                allAlleles.add(Allele.create(b, false));
        }

        PrintStream out = null;
        if (SIMULATE_NOISY_PILEUP) {
            try {
                out = new PrintStream(new File("/humgen/gsa-scr1/delangel/GATK/Sting_unstable_mac/GLUnitTest.table"));
    //                            out = new PrintStream(new File("/Users/delangel/GATK/Sting_unstable/GLUnitTest.table"));
            }
            catch (Exception e) {}
            // write header
            out.format("Depth\tPoolPloidy\tACTrue\tACEst\tREF\tALTTrue\tALTEst\n");
        }
        final int[] depthVector = {1000,10000};
        //final double[] alleleFrequencyVector = {0.01,0.1,0.5,1.0};
        final int[] spVector = {10,100};
        //final int[] spVector = {1};
        for (int depth : depthVector) {
            for (int nSamplesPerPool : spVector) {
                final int ploidy = 2*nSamplesPerPool;
                for (int ac =2; ac <=ploidy; ac++) {

                    // simulate pileup with given AC and depth
                    int altDepth = (int)Math.round( (double)ac/(double)ploidy * (double)depth);
                    final int[] numReadsPerAllele = {depth-altDepth,altDepth};
                    final Map<String,AlignmentContext> alignmentContextMap =
                            readPileupTestProvider.getAlignmentContextFromAlleles(eventLength, new String(new byte[]{altByte}), numReadsPerAllele);

                    // get now likelihoods for this

                    final PoolSNPGenotypeLikelihoods GL = new PoolSNPGenotypeLikelihoods(allAlleles, null, nSamplesPerPool*2, noiselessErrorModels, false, true);
                    final int nGoodBases = GL.add(alignmentContextMap.get("sample0000").getBasePileup(), true, false, UAC.MIN_BASE_QUALTY_SCORE);
                    if (VERBOSE) {
                        System.out.format("Depth:%d, AC:%d, altDepth:%d, samplesPerPool:%d\nGLs:", depth,ac,altDepth, nSamplesPerPool);
                       System.out.println(GL.toString());
                    }
                    Assert.assertEquals(nGoodBases, depth);
                    Pair<int[],Double> mlPair = GL.getMostLikelyACCount();

                    // Most likely element has to be conformation REF = nSamples-AC,ALT = AC
                    if (ac == 0) {
                        Assert.assertEquals(mlPair.first[refIdx],ploidy);
                    } else {
                        Assert.assertEquals(mlPair.first[altIdx],ac);
                        Assert.assertEquals(mlPair.first[refIdx],ploidy-ac);
                    }


                    // simulate now pileup with base error rate
                    if (SIMULATE_NOISY_PILEUP) {
                        System.out.format("Depth:%d, AC:%d, altDepth:%d, samplesPerPool:%d\n", depth,ac,altDepth, nSamplesPerPool);

                         for (int k=0; k < NUM_SIMULATED_OBS; k++) {
                            final Map<String,AlignmentContext> noisyAlignmentContextMap =
                                    readPileupTestProvider.getAlignmentContextFromAlleles(eventLength, new String(new byte[]{altByte}), numReadsPerAllele,
                                            true, PHRED_SITE_ERROR_RATE);

                            // get now likelihoods for this

                            final PoolSNPGenotypeLikelihoods noisyGL = new PoolSNPGenotypeLikelihoods(allAlleles, null, nSamplesPerPool*2, noisyErrorModels, false,true);
                            noisyGL.add(noisyAlignmentContextMap.get("sample0000").getBasePileup(), true, false, UAC.MIN_BASE_QUALTY_SCORE);
                            mlPair = noisyGL.getMostLikelyACCount();

                            // Most likely element has to be conformation REF = nSamples-AC,ALT = AC
                            int acEst;
                            if (ac == 0) {
                                acEst =  mlPair.first[refIdx];
                            } else {
                                acEst = mlPair.first[altIdx];
                            }
                            byte altEst = BaseUtils.baseIndexToSimpleBase(MathUtils.maxElementIndex(mlPair.first));
                            out.format("%d\t%d\t%d\t%d\t%c\t%c\t%c\n",depth, ploidy, ac, acEst, refByte, altByte, altEst);

                        }
                     }
                }
            }


        }
        if (SIMULATE_NOISY_PILEUP)
            out.close();


    }



}
