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

package org.broadinstitute.gatk.utils.pairhmm;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.QualityUtils;

import java.util.Arrays;

import static org.broadinstitute.gatk.utils.pairhmm.PairHMMModel.*;

/**
 * Created with IntelliJ IDEA.
 * User: bradt
 * Date: 6/11/13
 */
public class ArrayLoglessPairHMM extends PairHMM {
    private static final double INITIAL_CONDITION = Math.pow(2, 1020);
    private static final double INITIAL_CONDITION_LOG10 = Math.log10(INITIAL_CONDITION);

    // we divide e by 3 because the observed base could have come from any of the non-observed alleles
    protected static final double TRISTATE_CORRECTION = 3.0;

    protected double[][] transition = null; // The transition probabilities cache
    protected double[][] prior = null;      // The prior probabilities cache

    // Array declarations for arrays implementation
    private double[] currentMatchArray = null;
    private double[] currentDeleteArray = null;
    private double[] currentInsertArray = null;
    private double[] parentMatchArray = null;
    private double[] parentDeleteArray = null;
    private double[] parentInsertArray = null;
    private double[] grandparentMatchArray = null;
    private double[] grandparentDeleteArray = null;
    private double[] grandparentInsertArray = null;

    // When successive haplotypes have a common prefix, these arrays store cached info from the previous haplotype; for reading
    private double[] matchCacheArray = null;
    private double[] deleteCacheArray = null;
    private double[] insertCacheArray = null;

    // These arrays store cache info for use with the next haplotype; for writing
    private double[] nextMatchCacheArray = null;
    private double[] nextDeleteCacheArray = null;
    private double[] nextInsertCacheArray = null;

    // Used when caching to store our intermediate sum at point of first difference bw successive haplotypes
    private double partialSum;


    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize(final int readMaxLength, final int haplotypeMaxLength ) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        transition = PairHMMModel.createTransitionMatrix(maxReadLength);
        prior = new double[paddedMaxReadLength][paddedMaxHaplotypeLength];

        // Initialize all arrays
        // Final Cell of array is a padding cell, initialized to zero.
        currentMatchArray = new double[paddedMaxReadLength];
        currentDeleteArray = new double[paddedMaxReadLength];
        currentInsertArray = new double[paddedMaxReadLength];

        parentMatchArray = new double[paddedMaxReadLength];
        parentDeleteArray = new double[paddedMaxReadLength];
        parentInsertArray = new double[paddedMaxReadLength];

        grandparentMatchArray = new double[paddedMaxReadLength];
        grandparentDeleteArray = new double[paddedMaxReadLength];
        grandparentInsertArray = new double[paddedMaxReadLength];

        // Initialize the special arrays used for caching when successive haplotypes have a common prefix
        matchCacheArray = new double[paddedMaxReadLength];
        deleteCacheArray = new double[paddedMaxReadLength];
        insertCacheArray = new double[paddedMaxReadLength];

        nextMatchCacheArray = new double[paddedMaxReadLength];
        nextDeleteCacheArray = new double[paddedMaxReadLength];
        nextInsertCacheArray = new double [paddedMaxReadLength];
    }


    /**
     * {@inheritDoc}
     */
    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               int hapStartIndex,
                                                               final boolean recacheReadValues,
                                                               final int nextHapStartIndex) {

        if ( ! constantsAreInitialized) {
            initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);

            // note that we initialized the constants
            constantsAreInitialized = true;
        }
        initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);

        // Some housekeeping to be done if we are starting a new read
        if (recacheReadValues) {
            hapStartIndex = 0;

            initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP);
            // note that we initialized the constants
            constantsAreInitialized = true;

            // Read length may have changed, so we need to set zero-value padding at the appropriate position.
            padMatchAndInsertArrays(readBases.length);
        }

        // if we have not cached from a previous haplotype, clear any info we may have accumulated in a previous HMM iteration
        if (hapStartIndex == 0) {
            clearPreviouslyCachedInfo(readBases.length);

            // Haplotype length may have changed, so we need to set initial-value padding at the appropriate position.
            padDeleteArrays(haplotypeBases.length, readBases.length);
        }

        // We build up our solution by looking at position [0] in the match, insert arrays. Need to set these to 0 before we start.
        clearArraySolutionPosition();

        // Some parameters to control behavior during the dynamic programming loop
        final int maxDiagonals = readBases.length + haplotypeBases.length - hapStartIndex - 1;   // Number of diagonals for a matrix  = rows + cols - 1;
        int startFill;                                                                           // The lower bound of the array indices we want to over-write
        int endFill;                                                                             // The upper bound of the array indices we want to over-write
        final int cacheSumIndex = nextHapStartIndex - hapStartIndex + readBases.length - 1;      // This array will contain the partial sum to cache for the next haplotype
        double finalArraySumProbabilities = partialSum;                                          // The final answer prior to log10 correction

        // Perform dynamic programming using arrays, as if over diagonals of a hypothetical read/haplotype alignment matrix
        for (int i = 1; i <= maxDiagonals; i++) {
            // set the bounds for cells we wish to fill in the arrays
            startFill = Math.max(readBases.length - i, 0);
            endFill = Math.min(maxDiagonals - i + 1, readBases.length);

            // apply any previously cached array information
            if (i <= readBases.length)
                applyPreviouslyCachedInfo(startFill);

            // fill in the cells for our current arrays
            updateArrays(readBases.length, hapStartIndex, nextHapStartIndex, startFill, endFill, i);

            // final probability is the log10 sum of the last element in the Match and Insertion state arrays
            // this way we ignore all paths that ended in deletions! (huge)
            // but we have to sum all the paths ending in the M and I arrays, because they're no longer extended.
            // Where i > readBases.length, array[0] corresponds to bottom row of a [read] x [haplotype] matrix. Before this, they carries the 0's we set above.
            finalArraySumProbabilities += currentInsertArray[0] + currentMatchArray[0];

            // Partial sum for caching the next haplotype:
            // At the position of the last similar base between this haplotype and the next one...
            // ...remember the partial sum, so that we can start here on the next hap.
            if (i == cacheSumIndex)
                partialSum = finalArraySumProbabilities;

            rotateArrayReferences();
        }
        // The cache arrays we wrote for this haplotype will be read for the next haplotype.
        rotateCacheArrays();

        //return result
        return Math.log10(finalArraySumProbabilities) - INITIAL_CONDITION_LOG10;
    }

    /**
     * Initializes the matrix that holds all the constants related to the editing
     * distance between the read and the haplotype.
     *
     * @param haplotypeBases the bases of the haplotype
     * @param readBases      the bases of the read
     * @param readQuals      the base quality scores of the read
     * @param startIndex     where to start updating the distanceMatrix (in case this read is similar to the previous read)
     */
    public void initializePriors(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {

        // initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProb(qual) : (QualityUtils.qualToErrorProb(qual) / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION)) );
            }
        }
    }

    /**
     * Initializes the matrix that holds all the constants related to quality scores.
     *
     * @param insertionGOP   insertion quality scores of the read
     * @param deletionGOP    deletion quality scores of the read
     * @param overallGCP     overall gap continuation penalty
     */
    @Requires({
            "insertionGOP != null",
            "deletionGOP != null",
            "overallGCP != null"
    })
    @Ensures("constantsAreInitialized")
    protected static void initializeProbabilities(final double[][] transition, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP) {
        PairHMMModel.qualToTransProbs(transition,insertionGOP,deletionGOP,overallGCP);
    }

    /**
     * Pad the ends of the Match and Insert arrays with 0.
     * Analogous to setting zeros in the first row in the Match, Insert matrices of N2MemoryPairHMM.
     *
     * @param padPosition Which index in the arrays we wish to pad
     */
    private void padMatchAndInsertArrays(final int padPosition) {
        grandparentMatchArray[padPosition] = 0;
        grandparentInsertArray[padPosition] = 0;
        parentMatchArray[padPosition] = 0;
        parentInsertArray[padPosition] = 0;
        currentMatchArray[padPosition] = 0;
        currentInsertArray[padPosition] = 0;
        matchCacheArray[padPosition] = 0;
        insertCacheArray[padPosition] = 0;
        nextMatchCacheArray[padPosition] = 0;
        nextInsertCacheArray[padPosition] = 0;
    }

    /**
     * Pad the Delete arrays with an intial value. Let's us have free deletions at the beginning of the alignment.
     * Analogous to padding the first row of the Delete matrix of N2MemoryPairHMM.
     *
     * @param haplotypeLength   The length of the present haplotype. Necessary for calculating initial padding value
     * @param padPosition       Which index in the arrays we wish to pad
     */
    private void padDeleteArrays(final int haplotypeLength, final int padPosition) {
        final double initialValue = INITIAL_CONDITION / haplotypeLength;

        // Pad the deletion arrays. Akin to padding the first row in the deletion matrix
        parentDeleteArray[padPosition] = initialValue;
        grandparentDeleteArray[padPosition] = initialValue;
        currentDeleteArray[padPosition] = initialValue;
        deleteCacheArray[padPosition] = initialValue;
        nextDeleteCacheArray[padPosition] = initialValue;
    }

    /**
     * We build up our solution by looking at position [0] in the match, insert arrays. Need to set these to 0 before we start.
     *
     */
    private void clearArraySolutionPosition() {
        grandparentMatchArray[0] = 0;
        grandparentInsertArray[0] = 0;
        parentMatchArray[0] = 0;
        parentInsertArray[0] = 0;
        currentMatchArray[0] = 0;
        currentInsertArray[0] = 0;
    }

    /**
     * Clears cached information saved from the last haplotype,
     * allowing us to start at the beginning of the present haplotype with intitial values of 0.
     *
     * @param fillLength How much of the cache arrays do we need to zero
     */
    private void clearPreviouslyCachedInfo(final int fillLength) {
        Arrays.fill(matchCacheArray, 0, fillLength, 0);
        Arrays.fill(deleteCacheArray, 0, fillLength, 0);
        Arrays.fill(insertCacheArray, 0, fillLength, 0);

        partialSum = 0;
    }

    /**
     * Applies cached information saved from the last haplotype,
     * allowing us to start in the middle of the present haplotype.
     *
     * @param indK the index in the arrays we wish to update with cached info
     */
    private void applyPreviouslyCachedInfo(int indK) {
        // apply caching info necessary for calculating current DELETE array values
        parentMatchArray[indK] = matchCacheArray[indK];
        parentDeleteArray[indK] = deleteCacheArray[indK];

        // apply caching info necessary for calculating current MATCH array values
        grandparentMatchArray[indK + 1] = matchCacheArray[indK + 1];
        grandparentDeleteArray[indK + 1] = deleteCacheArray[indK + 1];
        grandparentInsertArray[indK + 1] = insertCacheArray[indK + 1];
    }

    /**
     * Records the mid-process state of one location in the read/haplotype alignment.
     * Writes new cache information for use with the next haplotype we see.
     *
     * @param indK the index in the cache arrays we wish to store information in
     */
    private void recordNewCacheInfo(int indK) {
        nextMatchCacheArray[indK] = currentMatchArray[indK];
        nextDeleteCacheArray[indK] = currentDeleteArray[indK];
        nextInsertCacheArray[indK] = currentInsertArray[indK];
    }

    /**
     * Update the HMM arrays for the current diagonal.
     *
     * @param readLength        The length of the read
     * @param hapStartIndex     An offset that tells us if we are starting in the middle of the present haplotype
     * @param nextHapStartIndex An offset that tells us which base in the NEXT haplotype we need to look at to record new caching info
     * @param startFill         The lower bound of the array indices we want to over-write
     * @param endFill           The upper bound of the array indices we want to over-write
     * @param iii               The index indicating which diagonal of the read/haplotype alignment we are working on
     */
    private void updateArrays(final int readLength,
                              final int hapStartIndex,
                              final int nextHapStartIndex,
                              final int startFill,
                              final int endFill,
                              final int iii) {

        // The coordinate in our priors and transition matrices corresponding to a given position in the read/haplotype alignment
        int matrixRow;
        int matrixCol;

        int arrayIndex;
        for (arrayIndex = startFill; arrayIndex < endFill; arrayIndex++) {
            // translate the array position into a row, column in the priors and transition matrices
            matrixRow = readLength - arrayIndex - 1;
            matrixCol = iii - matrixRow - 1 + hapStartIndex;

            // update cell for each of our current arrays. Prior, transition matrices are padded +1 row,col
            updateArrayCell(arrayIndex, prior[matrixRow+1][matrixCol+1], transition[matrixRow+1]);

            // Set up caching for the next haplotype
            // At the position of the final similar base between this haplotype and the next one, remember the mid-array values
            if (matrixCol == nextHapStartIndex - 1)
                recordNewCacheInfo(arrayIndex);
        }
    }

    /**
     * Updates a cell in the HMM arrays
     *
     * @param indK             index in the arrays to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param transition       an array with the six transition relevant to this location
     */
    private void updateArrayCell( final int indK, final double prior, final double[] transition) {
        currentMatchArray[indK] = prior * ( grandparentMatchArray[indK + 1] * transition[matchToMatch] +
                grandparentInsertArray[indK + 1] * transition[indelToMatch] +
                grandparentDeleteArray[indK + 1] * transition[indelToMatch] );
        currentInsertArray[indK] = parentMatchArray[indK + 1] * transition[matchToInsertion] + parentInsertArray[indK + 1] * transition[insertionToInsertion];
        currentDeleteArray[indK] = parentMatchArray[indK] * transition[matchToDeletion] + parentDeleteArray[indK] * transition[deletionToDeletion];
    }

    /**
     * To prepare for the next diagonal in our loop, each array must be bumped to an older generation
     *
     */
    private void rotateArrayReferences() {
        double[] tempMatchArray = grandparentMatchArray;
        double[] tempDeleteArray = grandparentDeleteArray;
        double[] tempInsertArray = grandparentInsertArray;

        grandparentMatchArray = parentMatchArray;
        grandparentDeleteArray = parentDeleteArray;
        grandparentInsertArray = parentInsertArray;

        parentMatchArray = currentMatchArray;
        parentDeleteArray = currentDeleteArray;
        parentInsertArray = currentInsertArray;

        currentMatchArray = tempMatchArray;
        currentDeleteArray = tempDeleteArray;
        currentInsertArray = tempInsertArray;
    }

    /**
     * To prepare for the next haplotype, the caching info we wrote is copied into the cach-read arrays
     *
     */
    private void rotateCacheArrays() {
        matchCacheArray = nextMatchCacheArray.clone();
        deleteCacheArray = nextDeleteCacheArray.clone();
        insertCacheArray = nextInsertCacheArray.clone();
    }
}
