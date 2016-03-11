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

package org.broadinstitute.gatk.tools.walkers.genotyper;

import htsjdk.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintStream;
import java.util.*;


public class GeneralPloidyGenotypeLikelihoodsUnitTest {

    final UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();
    final Logger logger = Logger.getLogger(Walker.class);
    private static final boolean VERBOSE = false;
    private static final boolean SIMULATE_NOISY_PILEUP = false;
    private static final int NUM_SIMULATED_OBS = 10;

    void PoolGenotypeLikelihoodsUnitTest() {
        UAC.minQualityScore = 5;
        UAC.maxQualityScore = 40;
        UAC.phredScaledPrior = (byte)20;
        UAC.minPower = 0.0;

    }
    @Test
    public void testStoringLikelihoodElements() {


        // basic test storing a given PL vector in a GeneralPloidyGenotypeLikelihoods object and then retrieving it back

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

        GeneralPloidyGenotypeLikelihoods gl = new GeneralPloidySNPGenotypeLikelihoods(alleles, gls,ploidy, null, false,true);
        double[] glnew = gl.getLikelihoods();

        Assert.assertEquals(gls, glnew);
    }

    @Test
    public void testElementStorageCache() {
        // compare cached element storage with compuationally hard-coded iterative computation

        for (int ploidy = 2; ploidy < 10; ploidy++) {
            for (int nAlleles = 2; nAlleles < 10; nAlleles++)
                Assert.assertEquals(GeneralPloidyGenotypeLikelihoods.getNumLikelihoodElements(nAlleles, ploidy),
                        GenotypeLikelihoods.numLikelihoods(nAlleles, ploidy));
        }

    }

    @Test
    public void testVectorToLinearIndex() {

        // create iterator, compare linear index given by iterator with closed form function
        int numAlleles = 4;
        int ploidy = 2;
        GeneralPloidyGenotypeLikelihoods.SumIterator iterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(numAlleles, ploidy);

        while(iterator.hasNext()) {
            System.out.format("\n%d:",iterator.getLinearIndex());
            int[] a =  iterator.getCurrentVector();
            for (int aa: a)
                System.out.format("%d ",aa);


            int computedIdx = GeneralPloidyGenotypeLikelihoods.getLinearIndex(a, numAlleles, ploidy);
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

        double[] newGLs = GeneralPloidyGenotypeLikelihoods.subsetToAlleles(oldLikelihoods, ploidy,
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
        GeneralPloidyGenotypeLikelihoods.SumIterator iterator = runIterator(seed,-1);
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

    private GeneralPloidyGenotypeLikelihoods.SumIterator runIterator(int[] seed, int restrictSumTo) {
        GeneralPloidyGenotypeLikelihoods.SumIterator iterator = new GeneralPloidyGenotypeLikelihoods.SumIterator(seed, restrictSumTo);

        while(iterator.hasNext()) {
            int[] a =  iterator.getCurrentVector();
            int idx = GeneralPloidyGenotypeLikelihoods.getLinearIndex(a, a.length, restrictSumTo);
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
                final ErrorModel emodel = new ErrorModel(UAC, refPileup, refVC, refPileupTestProvider.getReferenceContext());
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
        final byte refByte = refPileupTestProvider.getRefByte();
        final String altBases = "TCA";
        final String refSampleName = refPileupTestProvider.getSampleNames().get(0);
        final List<Allele> trueAlleles = new ArrayList<Allele>();
        trueAlleles.add(Allele.create(refByte, true));
        trueAlleles.add(Allele.create((char)refByte + "TC", false));

        final String fw = new String(refPileupTestProvider.getReferenceContext().getForwardBases());
        final VariantContext refInsertionVC = new VariantContextBuilder("test","chr1",refPileupTestProvider.getReferenceContext().getLocus().getStart(),
                refPileupTestProvider.getReferenceContext().getLocus().getStart(), trueAlleles).
                genotypes(GenotypeBuilder.create(refSampleName, trueAlleles)).make();


        final int[] matchArray = {95, 995, 9995, 10000};
        final int[] mismatchArray = {1,5,10,20};

        if (VERBOSE) System.out.println("Running indel error model test");
        for (int matches: matchArray) {
            for (int mismatches: mismatchArray) {
                // get artificial alignment context for ref sample - no noise
                // CASE 1: Test HET insertion
                // Ref sample has TC insertion but pileup will have TCA inserted instead to test mismatches
                Map<String,AlignmentContext> refContext = refPileupTestProvider.getAlignmentContextFromAlleles(1+altBases.length(), altBases, new int[]{matches, mismatches}, false, 30);
                final ReadBackedPileup refPileup = refContext.get(refSampleName).getBasePileup();
                final ErrorModel emodel = new ErrorModel(UAC, refPileup, refInsertionVC, refPileupTestProvider.getReferenceContext());
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
        delAlleles.add(Allele.create(fw.substring(0,delLength+1), true));
        delAlleles.add(Allele.create(refByte, false));

        final VariantContext refDeletionVC =  new VariantContextBuilder("test","chr1",refPileupTestProvider.getReferenceContext().getLocus().getStart(),
                refPileupTestProvider.getReferenceContext().getLocus().getStart()+delLength, delAlleles).
                genotypes(GenotypeBuilder.create(refSampleName, delAlleles)).make();

        for (int matches: matchArray) {
            for (int mismatches: mismatchArray) {
                // get artificial alignment context for ref sample - no noise
                // CASE 1: Test HET deletion
                // Ref sample has 4bp deletion but pileup will have 3 bp deletion instead to test mismatches
                Map<String,AlignmentContext> refContext = refPileupTestProvider.getAlignmentContextFromAlleles(-delLength+1, altBases, new int[]{matches, mismatches}, false, 30);
                final ReadBackedPileup refPileup = refContext.get(refSampleName).getBasePileup();
                final ErrorModel emodel = new ErrorModel(UAC, refPileup, refDeletionVC, refPileupTestProvider.getReferenceContext());
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

         // all first ref allele
        allAlleles.add(Allele.create(refByte,true));
        for (byte b: BaseUtils.BASES) {
            if (refByte != b)
                allAlleles.add(Allele.create(b, false));
        }

        final int refIdx = 0;
        int altIdx = -1;

        for (int k=0; k < allAlleles.size(); k++)
            if (altByte == allAlleles.get(k).getBases()[0]) {
                altIdx = k;
                break;
            }



        PrintStream out = null;
        if (SIMULATE_NOISY_PILEUP) {
            try {
                out = new PrintStream(new File("GLUnitTest.table"));
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

                    final GeneralPloidySNPGenotypeLikelihoods GL = new GeneralPloidySNPGenotypeLikelihoods(allAlleles, null, nSamplesPerPool*2, noiselessErrorModels, false, true);
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

                            final GeneralPloidySNPGenotypeLikelihoods noisyGL = new GeneralPloidySNPGenotypeLikelihoods(allAlleles, null, nSamplesPerPool*2, noisyErrorModels, false,true);
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
