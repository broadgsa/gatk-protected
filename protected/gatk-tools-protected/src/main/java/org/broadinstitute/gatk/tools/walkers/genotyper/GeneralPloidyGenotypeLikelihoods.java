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
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.ExactACcounts;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.ExactACset;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.*;

public abstract class GeneralPloidyGenotypeLikelihoods {
    protected final int numChromosomes;
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    protected static final boolean VERBOSE = false;
    protected static final double qualVec[] = new double[SAMUtils.MAX_PHRED_SCORE+1];

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods;
    protected double[][] logMismatchProbabilityArray;

    protected final int nSamplesPerPool;
    protected final HashMap<String, ErrorModel> perLaneErrorModels;
    protected final int likelihoodDim;
    protected final boolean ignoreLaneInformation;
    protected final double LOG10_PLOIDY;
    protected boolean hasReferenceSampleData;

    protected final int nAlleles;
    protected final List<Allele> alleles;

    private static final double MIN_LIKELIHOOD = Double.NEGATIVE_INFINITY;
            
    private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
    private static final int MAX_NUM_SAMPLES_PER_POOL = 1000;

    private static final boolean FAST_GL_COMPUTATION = true;
    // constructor with given logPL elements
    public GeneralPloidyGenotypeLikelihoods(final List<Allele> alleles, final double[] logLikelihoods, final int ploidy,
                                            final HashMap<String, ErrorModel> perLaneErrorModels, final boolean ignoreLaneInformation) {
        this.alleles = alleles;
        this.nAlleles = alleles.size();
        numChromosomes = ploidy;
        nSamplesPerPool = numChromosomes/2;
        this.perLaneErrorModels = perLaneErrorModels;
        this.ignoreLaneInformation = ignoreLaneInformation;

        // check if at least one lane has actual data
        if (perLaneErrorModels == null || perLaneErrorModels.isEmpty())
            hasReferenceSampleData = false;
        else {
            for (Map.Entry<String,ErrorModel> elt : perLaneErrorModels.entrySet()) {
                if (elt.getValue().hasData()) {
                    hasReferenceSampleData = true;
                    break;
                }
            }
        }
        // check sizes
        if (nAlleles > MAX_NUM_ALLELES_TO_CACHE)
            throw new UserException("No support for this number of alleles");

        if (nSamplesPerPool > MAX_NUM_SAMPLES_PER_POOL)
            throw new UserException("No support for such large number of samples per pool");

        likelihoodDim = GenotypeLikelihoods.numLikelihoods(nAlleles, numChromosomes);

        if (logLikelihoods == null){
            log10Likelihoods = new double[likelihoodDim]; 
            Arrays.fill(log10Likelihoods, MIN_LIKELIHOOD);
        } else {
            if (logLikelihoods.length != likelihoodDim)
                throw new ReviewedGATKException("BUG: inconsistent parameters when creating GeneralPloidyGenotypeLikelihoods object");

            log10Likelihoods = logLikelihoods; //.clone(); // is clone needed?
        }
        fillCache();
        LOG10_PLOIDY = Math.log10((double)numChromosomes);
   }


    /**
     * Crucial inner class that handles addressing elements of pool likelihoods. We store likelihoods as a map
     * of form int[] -> double (to be more precise, IntArrayWrapper -> Double).
     * For a given ploidy (chromosome count) and number of alleles, we need a form to iterate deterministically
     * across all possible allele conformations.
     * Problem equivalent to listing in deterministic order all possible ways in which N integers will sum to P,
     * where N is number of alleles and P is number of chromosomes.
     * There's an option to list all integers so that sum will be UP to P.
     * For example, with P=2,N=2, restrictSumTo = 2 iterator will produce
     * [2 0] [1 1] [0 2]
     *
     *
     */
    public static class SumIterator {
        private int[] currentState;
        private final int[] finalState;
        private final int restrictSumTo;
        private final int dim;
        private boolean hasNext;
        private int linearIndex;
        private int currentSum;

        /**
         * Default constructor. Typical use case: restrictSumTo = -1 if there's no sum restriction, or will generate int[]
         * vectors so that all add to this value.
         *
         * @param finalState                    End state - typically we should set value to (P,P,P,...)
         * @param restrictSumTo                 See above
         */
        public SumIterator(final int[] finalState,final int restrictSumTo) {
            this.finalState = finalState;
            this.dim = finalState.length;
            this.restrictSumTo = restrictSumTo;
            currentState = new int[dim];
            reset();

        }

        /**
         * Shortcut constructor for common use case: iterator will produce 
         * all vectors of length numAlleles whose sum = numChromosomes
         * @param numAlleles              Number of alleles
         * @param numChromosomes          Ploidy
         */
        public SumIterator(final int numAlleles, final int numChromosomes) {
            this(getInitialStateVector(numAlleles,numChromosomes), numChromosomes);            
        }


        private static int[] getInitialStateVector(final int nAlleles, final int numChromosomes) {
            int[] initialState = new int[nAlleles];
            Arrays.fill(initialState,numChromosomes);
            return initialState;
        }
        
        public void setInitialStateVector(final int[] stateVector) {
            if (restrictSumTo > 0) {
                // check that desired vector is valid
                if (MathUtils.sum(stateVector) != restrictSumTo)
                    throw new ReviewedGATKException("BUG: initial state vector nor compatible with sum iterator");

                final int numAlleles = currentState.length;
                final int ploidy = restrictSumTo;

                linearIndex = GeneralPloidyGenotypeLikelihoods.getLinearIndex(stateVector, numAlleles, ploidy);
            }
            else
                throw new ReviewedGATKException("BUG: Not supported");

        }
        public void next() {
            int initialDim = (restrictSumTo > 0)?1:0;
            hasNext = next(finalState, initialDim);
            if (hasNext)
                linearIndex++;
        }

        private boolean next(final int[] finalState, int initialDim) {
            boolean hasNextState = false;
            for (int currentDim=initialDim; currentDim < finalState.length; currentDim++) {
                final int x = currentState[currentDim]+1;

                if (x > finalState[currentDim] || (currentSum >= restrictSumTo && initialDim > 0)) {
                    // update vector sum, and reset position
                    currentSum -= currentState[currentDim];
                    currentState[currentDim] = 0;
                    if (currentDim >= dim-1) {
                        hasNextState = false;
                        break;
                    }
                }
                else {
                    currentState[currentDim] = x;
                    hasNextState = true;
                    currentSum++;
                    break;
                }
            }
            if (initialDim > 0) {
                currentState[0] = restrictSumTo - currentSum;
            }
            return hasNextState;
        }

        public void reset() {
            Arrays.fill(currentState, 0);
            if (restrictSumTo > 0)
                currentState[0] = restrictSumTo;
            hasNext = true;
            linearIndex = 0;
            currentSum = 0;
        }
        public int[] getCurrentVector() {
            return currentState;
        }
        
        public int[] getCurrentAltVector() {
            return Arrays.copyOfRange(currentState,1,currentState.length);
        }
      /*  public int getCurrentSum() {
            return currentSum;
        }
        */
        public int getLinearIndex() {
            return linearIndex;
        }

        public boolean hasNext() {
            return hasNext;
        }
    }

    public List<Allele> getAlleles() { return alleles;}

    /**
     * Returns an array of log10 likelihoods for each genotype conformation, with ordering determined by SumIterator class.
     *
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }





    /**
     * Set particular element of logPL vector
     * @param idx          index of allele count conformation to modify
     * @param pl                    Likelihood to associate with map
     */
    public void setLogPLs(final int idx, final double pl) {
            log10Likelihoods[idx] = pl;
    }

    public void renormalize() {
        log10Likelihoods = MathUtils.normalizeFromLog10(log10Likelihoods,false,true);
    }
    /** Compute most likely AC conformation based on currently stored PL's - just loop through log PL map and output max value
     *
     * @return vector with most likely allele count, ordered according to this object's alleles
     */
    public Pair<int[],Double> getMostLikelyACCount() {

        int[] mlInd = null;
        double maxVal = Double.NEGATIVE_INFINITY;

        final SumIterator iterator = new SumIterator(alleles.size(),numChromosomes);

        int idx = 0;
        while (iterator.hasNext()) {
            double pl = log10Likelihoods[idx++];
            if (pl > maxVal) {
                maxVal = pl;
                mlInd = iterator.getCurrentVector().clone();
                
            }
            iterator.next();
        }
        if (VERBOSE) {
            System.out.println(GATKVCFConstants.MLE_ALLELE_COUNT_KEY + ": " + Arrays.toString(mlInd));
        }
        return new Pair<int[], Double>(mlInd,maxVal);
    }

    /**
     * Given set of alleles with corresponding vector of likelihoods, subset to a new set of alleles
     *
     * @param oldLikelihoods        Vector of PL's corresponding to original alleles
     * @param numChromosomes        Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param allelesToSubset       Alleles to subset
     * @return                      Vector of new PL's, ordered according to SumIterator's ordering
     */
    public static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
                                                   final List<Allele> originalAlleles, final List<Allele> allelesToSubset) {

        int newPLSize = GeneralPloidyGenotypeLikelihoods.getNumLikelihoodElements(allelesToSubset.size(), numChromosomes);
        double[] newPLs = new double[newPLSize];

        int idx = 0;
        // First fill boolean array stating whether each original allele is present in new mapping
        final boolean [] allelePresent = new boolean[originalAlleles.size()];
        for ( Allele allele : originalAlleles )
            allelePresent[idx++] = allelesToSubset.contains(allele);

        // compute mapping from old idx to new idx
        // This might be needed in case new allele set is not ordered in the same way as old set
        // Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}. Permutation key = [2,1]

        int[] permutationKey = new int[allelesToSubset.size()];
        for (int k=0; k < allelesToSubset.size(); k++)
            // for each allele to subset, find corresponding index in original allele list
            permutationKey[k] = originalAlleles.indexOf(allelesToSubset.get(k));


        if (VERBOSE) {
            System.out.println("permutationKey:"+Arrays.toString(permutationKey));
        }

        final SumIterator iterator = new SumIterator(originalAlleles.size(),numChromosomes);

        while (iterator.hasNext()) {
            // for each entry in logPL table, associated originally with allele count stored in vec[],
            // see if this allele count conformation will be present in new logPL table.
            // For entry to be present, elements in dimensions not present in requested allele list have to have count = 0
            int[] pVec = iterator.getCurrentVector();
            double pl = oldLikelihoods[iterator.getLinearIndex()];
            
            boolean keyPresent = true;
            for (int k=0; k < allelePresent.length; k++)
                if ( pVec[k]>0 && !allelePresent[k] )
                    keyPresent = false;

            if (keyPresent) {// skip to next entry in logPLs if this conformation is not present in subset

                final int[] newCount = new int[allelesToSubset.size()];
    
                // map from old allele mapping count to new allele mapping
                // In pseudo-Matlab notation: newCount = vec[permutationKey] for permutationKey vector
                for (idx = 0; idx < newCount.length; idx++)
                    newCount[idx] =  pVec[permutationKey[idx]];
    
                // get corresponding index from new count
                int outputIdx = GeneralPloidyGenotypeLikelihoods.getLinearIndex(newCount, allelesToSubset.size(), numChromosomes);
                newPLs[outputIdx] = pl;
                if (VERBOSE) {
                    System.out.println("Old Key:"+Arrays.toString(pVec));
                    System.out.println("New Key:"+Arrays.toString(newCount));
                }
            }
            iterator.next();
        }

        return  newPLs;
    }

    public static int getLinearIndex(int[] vectorIdx, int numAlleles, int ploidy) {

        if (ploidy <= 0)
            return 0;

        int linearIdx = 0;
        int cumSum = ploidy;
        for (int k=numAlleles-1;k>=1; k--) {
            int idx = vectorIdx[k];
            // how many blocks are before current position
            if (idx == 0)
                continue;
            for (int p=0; p < idx; p++)
                linearIdx += getNumLikelihoodElements( k, cumSum-p);
            
            cumSum -= idx;
        }

        return linearIdx;
        
    }

    /**
     * Given a scalar index, what's the alelle count conformation corresponding to it?
     * @param nAlleles                    Number of alleles
     * @param numChromosomes              Ploidy
     * @param PLindex                     Index to query
     * @return                            Allele count conformation, according to iteration order from SumIterator
     */
    public static int[] getAlleleCountFromPLIndex(final int nAlleles, final int numChromosomes, final int PLindex) {

        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(numChromosomes, nAlleles);
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);
        return alleleCounts.alleleCountsByIndex(nAlleles - 1);
    }

    /*
    * a cache of the PL ivector sizes as a function of # of alleles and pool sizes
    */
    
    public static int getNumLikelihoodElements(int numAlleles, int ploidy) {
        return GenotypeLikelihoodVectorSizes[numAlleles][ploidy];
    }

    private final static int[][] GenotypeLikelihoodVectorSizes = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE, 2*MAX_NUM_SAMPLES_PER_POOL);

    private static int[][] fillGLVectorSizeCache(int maxAlleles, int maxPloidy) {
        
        int[][] cache = new int[maxAlleles][maxPloidy];
        for (int numAlleles=1; numAlleles < maxAlleles; numAlleles++) {
            for (int ploidy=0; ploidy < maxPloidy; ploidy++) {

                if (numAlleles == 1)
                    cache[numAlleles][ploidy] = 1;
                else if (ploidy == 1)
                    cache[numAlleles][ploidy] = numAlleles;
                else {
                    int acc =0;
                    for (int k=0; k <= ploidy; k++ )
                        acc += cache[numAlleles-1][ploidy-k];

                    cache[numAlleles][ploidy] = acc;
                }
            }
        }
        return cache;
    }

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        StringBuilder s = new StringBuilder(1000);

        s.append("Alleles:");
        for (Allele a: this.alleles){
            s.append(a.getDisplayString());
            s.append(",");
        }
        s.append("\nGLs:\n");
        SumIterator iterator = new SumIterator(nAlleles,numChromosomes);
        while (iterator.hasNext()) {
            if (!Double.isInfinite(getLikelihoods()[iterator.getLinearIndex()])) {

                s.append("Count [");
                StringBuilder b = new StringBuilder(iterator.getCurrentVector().length*2);
                for (int it:iterator.getCurrentVector()) {
                    b.append(it);
                    b.append(",");
                }
                s.append(b.toString());
                s.append(String.format("] GL=%4.3f\n",this.getLikelihoods()[iterator.getLinearIndex()]) );
            }
            iterator.next();
        }
        return s.toString();
    }


    public void computeLikelihoods(ErrorModel errorModel,
        List<Allele> alleleList, List<Integer> numObservations, ReadBackedPileup pileup) {

        if (FAST_GL_COMPUTATION) {
            //  queue up elements to be computed. Assumptions:
            // GLs distributions are unimodal
            // GLs are continuous
            // Hence, once an AC conformation is computed, we queue up its immediate topological neighbors.
            // If neighbors fall below maximum - threshold, we don't queue up THEIR own neighbors
            // and we repeat until queue is empty
            // queue of AC conformations to process
            final LinkedList<ExactACset> ACqueue = new LinkedList<ExactACset>();
            // mapping of ExactACset indexes to the objects
            final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<ExactACcounts, ExactACset>(likelihoodDim);
            // add AC=0 to the queue
            final int[] zeroCounts = new int[nAlleles];
            zeroCounts[0] = numChromosomes;

            ExactACset zeroSet =
                    new ExactACset(1, new ExactACcounts(zeroCounts));

            ACqueue.add(zeroSet);
            indexesToACset.put(zeroSet.getACcounts(), zeroSet);

            // keep processing while we have AC conformations that need to be calculated
            double maxLog10L = Double.NEGATIVE_INFINITY;
            while ( !ACqueue.isEmpty() ) {
                // compute log10Likelihoods
                final ExactACset ACset = ACqueue.remove();
                final double log10LofKs = calculateACConformationAndUpdateQueue(ACset, errorModel, alleleList, numObservations, maxLog10L, ACqueue, indexesToACset, pileup);

                // adjust max likelihood seen if needed
                maxLog10L = Math.max(maxLog10L, log10LofKs);
                // clean up memory
                indexesToACset.remove(ACset.getACcounts());
                if ( VERBOSE )
                    System.out.printf(" *** removing used set=%s%n", ACset.getACcounts());

             }


        }   else {
            int plIdx = 0;
            SumIterator iterator = new SumIterator(nAlleles, numChromosomes);
            while (iterator.hasNext()) {
                ExactACset ACset =
                       new ExactACset(1, new ExactACcounts(iterator.getCurrentVector()));
                // for observed base X, add Q(jX,k) to likelihood vector for all k in error model
                //likelihood(jA,jC,jG,jT) = logsum(logPr (errorModel[k],nA*Q(jA,k) +  nC*Q(jC,k) + nG*Q(jG,k) + nT*Q(jT,k))
                getLikelihoodOfConformation(ACset, errorModel, alleleList, numObservations, pileup);

                setLogPLs(plIdx++, ACset.getLog10Likelihoods()[0]);
                iterator.next();
            }
        }
        // normalize PL's
        renormalize();

    }

    private double calculateACConformationAndUpdateQueue(final ExactACset set,
                                                         final ErrorModel errorModel,
                                                         final List<Allele> alleleList,
                                                         final List<Integer> numObservations,
                                                         final double  maxLog10L,
                                                         final LinkedList<ExactACset> ACqueue,
                                                         final HashMap<ExactACcounts,
                                                                 ExactACset> indexesToACset,
                                                         final ReadBackedPileup pileup) {
        // compute likelihood of set
        getLikelihoodOfConformation(set, errorModel, alleleList, numObservations, pileup);
        final double log10LofK = set.getLog10Likelihoods()[0];
        
        // log result in PL vector
        int idx = getLinearIndex(set.getACcounts().getCounts(), nAlleles, numChromosomes);
        setLogPLs(idx, log10LofK);

        // can we abort early because the log10Likelihoods are so small?
        if ( log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY ) {
            if ( VERBOSE )
                System.out.printf(" *** breaking early set=%s log10L=%.2f maxLog10L=%.2f%n", set.getACcounts(), log10LofK, maxLog10L);
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        // by convention, ACcounts contained in set have full vector of possible pool ac counts including ref count.
        final int ACwiggle = numChromosomes - set.getACsum() + set.getACcounts().getCounts()[0];
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;


        // add conformations for other cases
        for ( int allele = 1; allele < nAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // is this a valid conformation?
            int altSum = (int)MathUtils.sum(ACcountsClone) - ACcountsClone[0];
            ACcountsClone[0] = numChromosomes - altSum;
            if (ACcountsClone[0] < 0)
                continue;


            updateACset(ACcountsClone, ACqueue, indexesToACset);
        }
        return log10LofK;

    }

    /**
     * Abstract methods, must be implemented in subclasses
     *
     * @param ACset       Count to compute
     * @param errorModel    Site-specific error model object
     * @param alleleList    List of alleles
     * @param numObservations Number of observations for each allele
     * @param pileup        Read backed pileup in case it's necessary
     */
    public abstract void getLikelihoodOfConformation(final ExactACset ACset,
                                                     final ErrorModel errorModel,
                                                     final List<Allele> alleleList,
                                                     final List<Integer> numObservations,
                                                     final ReadBackedPileup pileup);


    public abstract int add(ReadBackedPileup pileup, UnifiedArgumentCollection UAC);

    // Static methods
    public static void updateACset(final int[] newSetCounts,
                                    final LinkedList<ExactACset> ACqueue,
                                    final HashMap<ExactACcounts, ExactACset> indexesToACset) {

        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            ExactACset newSet = new ExactACset(1, index);
            indexesToACset.put(index, newSet);
            ACqueue.add(newSet);     
            if (VERBOSE)
                System.out.println(" *** Adding set to queue:" + index.toString());
        }

    }
    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------


    //
    // Constant static data
    //

    static {
        // cache 10^(-k/10)
        for (int j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++)
            qualVec[j] = Math.pow(10.0,-(double)j/10.0);
    }

    private void fillCache() {
        // cache Q(j,k) = log10(j/2N*(1-ek) + (2N-j)/2N*ek) for j = 0:2N

        logMismatchProbabilityArray = new double[1+numChromosomes][1+SAMUtils.MAX_PHRED_SCORE];
        for (int i=0; i <= numChromosomes; i++) {
            for (int j=0; j <= SAMUtils.MAX_PHRED_SCORE; j++) {
                double phi = (double)i/numChromosomes;
                logMismatchProbabilityArray[i][j] = Math.log10(phi * (1.0-qualVec[j]) + qualVec[j]/3.0 * (1.0-phi));
            }
        }
    }

}

