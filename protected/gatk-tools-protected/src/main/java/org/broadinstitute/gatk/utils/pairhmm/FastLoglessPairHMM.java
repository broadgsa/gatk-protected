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

import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.broadinstitute.gatk.utils.pairhmm.PairHMMModel.*;

/**
 * Fast partial PairHMM backed on the standard Logless PairHMM
 *
 */
public class FastLoglessPairHMM extends LoglessPairHMM  implements FlexibleHMM {

    /**
     * Initial read length capacity.
     */
    private static final int INITIAL_READ_LENGTH_CAPACITY = 200;

    /**
     * Initial haplotype length capacity.
     */
    private static final int INITIAL_HAPLOTYPE_LENGTH_CAPACITY = 400;


    /**
     * Holds the current read capacity.
     * <p>It can only go up overtime.</p>
     */
    private int readCapacity = INITIAL_READ_LENGTH_CAPACITY;

    /**
     * Holds the current haplotype length capacity.
     * <p>It can only go up overtime.</p>
     */
    private int haplotypeCapacity = INITIAL_HAPLOTYPE_LENGTH_CAPACITY;

    private int maxToCol;
    private int haplotypeLength;

    private double[][] logTransition;

    private boolean indelToIndelIsConstant;

    /**
     * Returns the currently loaded read base qualities.
     *
     * @throws IllegalStateException if no read was previously loaded using {@link #loadRead}.
     * @return never {@code null}.
     */
    public byte[] getReadQuals() {
        if (readQuals == null)
            throw new IllegalStateException("no read was loaded onto the pairhmm calculator");
        return readQuals;
    }

    /**
     * Returns the currently loaded read insertion qualities.
     *
     * @throws IllegalStateException if no read was previously loaded using {@link #loadRead}.
     * @return never {@code null}.
     */
    @SuppressWarnings("unused")
    public byte[] getReadInsQuals() {
        if (readQuals == null)
            throw new IllegalStateException("no read was loaded onto the pairhmm calculator");
        return readInsQuals;
    }

    /**
     * Returns the currently loaded read deletion qualities.
     *
     * @throws IllegalStateException if no read was previously loaded using {@link #loadRead}.
     * @return never {@code null}.
     */
    @SuppressWarnings("unused")
    public byte[] getReadDelQuals() {
        if (readQuals == null)
            throw new IllegalStateException("no read was loaded onto the pairhmm calculator");
        return readDelQuals;
    }

    /**
     * Returns the currently loaded read gap extension penalty..
     *
     * @throws IllegalStateException if no read was previously loaded using {@link #loadRead}.
     * @return never {@code null}.
     */
    @SuppressWarnings("unused")
    public byte[] getReadGepQuals() {
        if (readQuals == null)
            throw new IllegalStateException("no read was loaded onto the pairhmm calculator");
        return readGepQuals;
    }

    /**
     * Creates a new pair-hmm calculator instance give the gap continuation penalty.
     *
     * @param gcp the gap-continuation penalty.
     */
    public FastLoglessPairHMM(final byte gcp) {
        constantGCP = gcp;
        initialize(readCapacity,haplotypeCapacity);
    }

    @Override
    public byte getGapExtensionPenalty() {
        return constantGCP;
    }


    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10(final byte[] haplotypeBases,
                   final byte[] readBases,
                   final byte[] readQuals,
                   final byte[] insertionGOP,
                   final byte[] deletionGOP,
                   final byte[] overallGCP,
                   final int hapStartIndex,
                   final boolean recacheReadValues, final int nextHapStartIndex) {
        this.readBases = readBases;
        this.haplotypeBases = haplotypeBases;
        this.haplotypeLength = haplotypeBases.length;
        return super.subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases,readBases,readQuals,
                insertionGOP,deletionGOP,overallGCP,hapStartIndex,recacheReadValues,nextHapStartIndex);
    }

    /**
     * Implement the last step summation to calculate the total likelihood.
     *
     * @param row number of the last row of the pair-hmm where the likelihood values are present.
     * @param fromCol inclusive first column to include in the summation.
     * @param toCol exclusive last column to include in the summation.
     * @return 0 or less.
     */
    protected double finalLikelihoodCalculation(final int row,
                                                final int fromCol, final int toCol) {

        final double divider = Math.max(1,2 *(toCol - fromCol));
        final double dividerInverse = 1.0 / divider;
        double finalLikelihood = 0;

        for (int j = fromCol; j < toCol; j++) {
            finalLikelihood += matchMatrix[row][j] * dividerInverse;
            finalLikelihood += insertionMatrix[row][j] * dividerInverse;
        }
        return StrictMath.log10(finalLikelihood) - INITIAL_CONDITION_LOG10 + StrictMath.log10(divider);
    }

    /**
     * Initialize the matrix values for a problem including the trailing end of the read.
     *
     * <p>
     *     Notice that you can improve performance by omitting filling reusable values from
     *     previous haplotype calculations. You can set {@code haplotypeStartOffset} to skill
     *     those columns.
     * </p>
     *
     * @param readStart inclusive first position of the read used in the calculations.
     * @param readEnd exclusive last position of the read considered in the calculations.
     * @param haplotypeStartOffset offset of the haplotype right after the reusable prefix
     *                             from previous calls.
     *
     *
     */
    protected void initializeMatrixValuesForTrailingProblem(final int readStart, final int readEnd,
                                                            final int haplotypeStartOffset) {

        @SuppressWarnings("all")
        final int zeroRow = readStart;
        final int toRow = readEnd + 1;
        final int toCol = haplotypeLength + 1;

        // fill first row with -Inf fot M and I but not for Deletion if leading
        // to allow for free deletions at the beginning.
        if (readStart == 0) {
            // First row initialization:
            Arrays.fill(matchMatrix[zeroRow],haplotypeStartOffset,toCol,0);
            Arrays.fill(deletionMatrix[zeroRow],haplotypeStartOffset,toCol,INITIAL_CONDITION);

            if (haplotypeStartOffset == 0)
                for (int i = zeroRow + 1; i < toRow; i++)
                    insertionMatrix[i][0] = matchMatrix[i][0] = deletionMatrix[i][0] = 0;

        } else {
            Arrays.fill(matchMatrix[zeroRow], Math.max(1,haplotypeStartOffset), toCol,0);
            Arrays.fill(insertionMatrix[zeroRow], haplotypeStartOffset, toCol,0);
            if (haplotypeStartOffset == 0) {
                matchMatrix[zeroRow][0] = INITIAL_CONDITION;
                deletionMatrix[zeroRow][0] = 0;
            }
            if (haplotypeStartOffset <= 1) deletionMatrix[zeroRow][1] = matchMatrix[zeroRow][1] * transition[zeroRow][matchToDeletion];
            for (int i = Math.max(haplotypeStartOffset,2); i < toCol; i++) {
                deletionMatrix[zeroRow][i] = deletionMatrix[zeroRow][i - 1]
                        * transition[zeroRow][deletionToDeletion];
            }

            if (haplotypeStartOffset == 0) {
                matchMatrix[zeroRow + 1][0] = deletionMatrix[zeroRow + 1][0] = 0;
                insertionMatrix[zeroRow + 1][0] = matchMatrix[zeroRow][0] * transition[zeroRow + 1][matchToInsertion];


                for (int i = zeroRow + 2; i < toRow; i++) {
                  matchMatrix[i][0] = deletionMatrix[i][0] = 0;
                  insertionMatrix[i][0] = insertionMatrix[i - 1][0]
                          * transition[i][insertionToInsertion];
                }
            }
        }
    }

    /**
     * Initializes calculation matrices give the characteristics of the next and previous problems.
     * @param currentProblem reference to the Lk calculation problem we are dealing currently.
     * @param previousProblem reference to the Lk calculation problem that has been solved just before.
     *
     */
    protected void initializeMatrixValues(final Problem currentProblem, final Problem previousProblem) {
        if (previousProblem != null &&
                previousProblem.readStart == currentProblem.readStart &&
                previousProblem.hapStart == currentProblem.hapStart &&
                maxToCol >= currentProblem.hapEnd + 1)
            return;

        final int zeroRow = currentProblem.readStart;
        final int zeroCol = currentProblem.hapStart;
        final int toRow = currentProblem.readEnd + 1;
        final int toCol = currentProblem.hapEnd + 1;
        maxToCol = toCol;

        // fill first row with -Inf fot M and I but not for Deletion if leading
        // to allow for free deletions at the beginning.
        if (currentProblem.leading) {
            // First row initialization:
            Arrays.fill(matchMatrix[zeroRow],zeroCol,toCol,0);
            Arrays.fill(deletionMatrix[zeroRow],zeroCol,toCol,INITIAL_CONDITION);

            for (int i = zeroRow + 1; i < toRow; i++)
                insertionMatrix[i][zeroCol] = matchMatrix[i][zeroCol] = deletionMatrix[i][zeroCol] = 0;

        } else { // If not leading set the appropriate matching 1.0 prob and
            // deletion + extension.

            Arrays.fill(matchMatrix[zeroRow], zeroCol + 1, toCol,0);
            Arrays.fill(insertionMatrix[zeroRow], zeroCol, toCol,0);
            matchMatrix[zeroRow][zeroCol] = INITIAL_CONDITION;
            deletionMatrix[zeroRow][zeroCol] = 0;
            deletionMatrix[zeroRow][zeroCol + 1] = matchMatrix[zeroRow][zeroCol] * transition[zeroRow][matchToDeletion];
            for (int i = zeroCol + 2; i < toCol; i++) {
                deletionMatrix[zeroRow][i] = deletionMatrix[zeroRow][i - 1]
                        * transition[zeroRow][deletionToDeletion];
            }

            matchMatrix[zeroRow + 1][zeroCol] = deletionMatrix[zeroRow + 1][zeroCol] = 0;
            insertionMatrix[zeroRow + 1][zeroCol] = matchMatrix[zeroRow][zeroCol] * transition[zeroRow + 1][matchToInsertion];

            for (int i = zeroRow + 2; i < toRow; i++) {
                matchMatrix[i][zeroCol] = deletionMatrix[i][zeroCol] = 0;
                insertionMatrix[i][zeroCol] = insertionMatrix[i - 1][zeroCol]
                        * transition[i][insertionToInsertion];
            }
        }
    }

    /**
     * Constant gap-continuation-penalty.
     */
    private final byte constantGCP;

    /**
     * Currently loaded haplotype base sequence.
     */
    private byte[] haplotypeBases;

    /**
     * Currently loaded read base sequence.
     */
    private byte[] readBases;

    /**
     * Read qualities.
     */
    private byte[] readQuals;

    /**
     * Read insertion qualities.
     */
    private byte[] readInsQuals;

    /**
     * Read deletion qualities.
     */
    private byte[] readDelQuals;

    /**
     * Read gap-extension-penalties.
     */
    private byte[] readGepQuals;

    /**
     * Cached results.
     */
    private Map<Problem, Double> cachedResults = new HashMap<>();

    /**
     * Loads the read that is going to be evaluated in following calls to {@link #calculateLocalLikelihoods}.
     *
     * @param read the target read.
     * @throws NullPointerException if {@code read} is null.
     */
    @Override
    public void loadRead(final GATKSAMRecord read) {
        loadRead(read.getReadBases(),read.getBaseQualities(),read.getBaseInsertionQualities(),read.getBaseDeletionQualities(),read.getMappingQuality());
    }

    /**
     * Loads the read that is going to be evaluated in following calls to {@link #calculateLocalLikelihoods}.
     *
     * @param readBases the read bases.
     * @param readQuals the read base call quality scores.
     * @param readInsQuals the read insertion quality scores.
     * @param readDelQuals  the read deletion quality scores.
     * @param mq the read mapping quality score.
     * @throws NullPointerException if any of the arrays passed is {@code null}.
     * @throws IllegalArgumentException if the arrays passed have incompatible sizes.
     */
    public void loadRead(final byte[] readBases, final byte[] readQuals, final byte[] readInsQuals, final byte[] readDelQuals, int mq) {
        // TODO This is a copy&paste from PairHMM*Engine read data preparation code.
        // TODO It is simply to difficult to share the code without changing that class
        if (readBases.length != readQuals.length) throw new IllegalArgumentException("the read quality array length does not match the read base array length");
        if (readBases.length != readInsQuals.length) throw new IllegalArgumentException("the read insert quality array length does not match the read base array length");
        if (readBases.length != readDelQuals.length) throw new IllegalArgumentException("the read deletion quality length does not match the read base array length");
        maxToCol = 0;

        if (readBases.length > readCapacity) {
            readCapacity = readBases.length;
            initialize(readCapacity,haplotypeCapacity);
        }
        paddedReadLength = readBases.length + 1;
        final byte[] overallGCP = new byte[readBases.length];
        Arrays.fill(overallGCP, constantGCP); // Is there a way to derive

        for (int kkk = 0; kkk < readQuals.length; kkk++) {
            readQuals[kkk] = (byte) Math.min(0xff & readQuals[kkk],
                    mq); // cap base quality by mapping
            readQuals[kkk] = (byte) (readQuals[kkk] < BASE_QUALITY_SCORE_THRESHOLD ? QualityUtils.MIN_USABLE_Q_SCORE
                    : Math.max(QualityUtils.MIN_USABLE_Q_SCORE,readQuals[kkk]));
            readInsQuals[kkk] = (byte) Math.max(QualityUtils.MIN_USABLE_Q_SCORE,readInsQuals[kkk]);
            readDelQuals[kkk] = (byte) Math.max(QualityUtils.MIN_USABLE_Q_SCORE,readDelQuals[kkk]);
        }
        this.readBases = readBases;
        this.readQuals = readQuals;
        this.readInsQuals = readInsQuals;
        this.readDelQuals = readDelQuals;
        this.readGepQuals = overallGCP;

        initializeProbabilities(transition,readInsQuals, readDelQuals, overallGCP);
        PairHMMModel.qualToTransProbsLog10(logTransition,readInsQuals,readDelQuals,overallGCP);

        indelToIndelIsConstant = true;
        for (final double d : overallGCP)
            if (d != overallGCP[0]) {
                indelToIndelIsConstant = false;
                break;
            }

        if (haplotypeBases != null)
            fillPriorsTable(0);
        cachedResults.clear();
    }

    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength )  {
        super.initialize(readMaxLength,haplotypeMaxLength);
        logTransition = PairHMMModel.createTransitionMatrix(readMaxLength);
    }

    @Override
    public void loadHaplotypeBases(final byte[] haplotypeBases) {
        if (readBases == null)
            throw new IllegalStateException(
                    "no read was loaded before the haplotype");
        this.haplotypeBases = haplotypeBases.clone();
        haplotypeLength = haplotypeBases.length;
        paddedHaplotypeLength = haplotypeLength;
        if (haplotypeCapacity < haplotypeLength) {
            haplotypeCapacity = haplotypeLength;
            initialize(readCapacity,haplotypeCapacity);
            initializeProbabilities(transition, readInsQuals, readDelQuals, readGepQuals);
        }
        initializePriors(this.haplotypeBases, readBases, readQuals, 0);
    }


    /**
     * Changes only the suffix of the currently loaded haplotype.
     *
     * <p>
     *     If from is 0, this is equivalent to call {@link #loadHaplotypeBases(byte[])} directly.
     * </p>
     * @param from first position on the current haplotype to substitute with the new suffix.
     *             It can be up to the length of the haplotype in such case this operation is in
     *             effect just extending that haplotype.
     * @param suffix the new bases for the end part of the current haplotype.
     * @param suffixFrom inclusive first position of the actual suffix within the {@code suffix} array.
     * @param suffixTo exclusive last position of the actual suffix within the {@code suffix} array.
     *
     * @throws IllegalStateException if no read was loaded with {@link #loadRead}.
     * @throws IllegalArgumentException if from is more than 0 but no haplotype was loaded previously or if indices passed are inconsistent.
     * @throws ArrayIndexOutOfBoundsException if indices passed are outside valid ranges.
     */
    public void changeHaplotypeSuffix(final int from, final byte[] suffix, final int suffixFrom, final int suffixTo) {
        if (readBases == null)
            throw new IllegalStateException(
                    "no read was loaded before the haplotype");
        if (haplotypeBases == null && from > 0)
            throw new IllegalArgumentException("from cannot be larger than 0 if no haplotype bases was previously loaded");
        if (suffixFrom < 0)
            throw new ArrayIndexOutOfBoundsException("the suffix from index cannot be negative");
        if (suffixTo > suffix.length)
            throw new ArrayIndexOutOfBoundsException("the suffix to index cannot be larger than the suffix array length");
        if (suffixFrom > suffixTo)
            throw new IllegalArgumentException("the suffix to index cannot be smaller than the suffix from index");
        if (from > haplotypeLength)
            throw new IllegalArgumentException("the from index cannot be greater than the current haplotype length");
        if (from < 0)
            throw new IllegalArgumentException("the from index cannot be negative");

        int startIndex = from;
        if (haplotypeBases == null) {
            haplotypeBases = Arrays.copyOfRange(suffix,suffixFrom,suffixTo);
            haplotypeLength = suffixTo - suffixFrom;
        } else {
            final int newLength = from + suffixTo - suffixFrom;
            if (haplotypeBases.length < newLength)
                haplotypeBases = Arrays.copyOf(haplotypeBases,newLength);
            System.arraycopy(suffix,suffixFrom,haplotypeBases,from,newLength - from);
            haplotypeLength = newLength;
        }
        paddedHaplotypeLength = haplotypeLength + 1;
        if (haplotypeCapacity < haplotypeLength) {
            haplotypeCapacity = haplotypeLength;
            initialize(readCapacity,haplotypeCapacity);
            initializeProbabilities(transition, readInsQuals, readDelQuals, readGepQuals);
            startIndex = 0;
        }
        //startIndex = 0;
        fillPriorsTable(startIndex);
    }

    /**
     * Returns the bases of the current haplotype.
     *
     * @throws IllegalStateException if no haplotype was loaded previously
     * @return never {@code null}
     */
    public byte[] getHaplotypeBases() {
        if (haplotypeBases == null)
            throw new IllegalStateException();
        return Arrays.copyOfRange(haplotypeBases,0,haplotypeLength);
    }

    /**
     * Returns a debug representation of the pair-hmm.
     * @return never {@code null}.
     */
    public String toString() {
        return "" + haplotypeLength + ":" + new String(Arrays.copyOfRange(haplotypeBases,0,haplotypeLength));
    }

    @Override
    protected void initializePriors(final byte[] hapBases, final byte[] readBases, final byte[] baseQuals, final int idx) {
        haplotypeBases = hapBases;
        haplotypeLength = haplotypeBases.length;
        this.readBases = readBases;
        this.readQuals = baseQuals;
        fillPriorsTable(idx);
    }

    /**
     * Fills the prior table up.
     *
     * <p>
     *     It accepts an argument to save unnecessary prefix filling up.
     * </p>
     *
     * @param idx first position in the haplotype to start filling from.
     */
    protected void fillPriorsTable(final int idx) {
        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = idx; j < haplotypeLength; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProb(qual) : (QualityUtils.qualToErrorProb(qual) / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION)) );
            }
        }
    }


    /**
     * Decorates haplotype set with their likelihoods as compared with the currently loaded read.
     *
     *
     * @param readStart inclusive start position of the targeted section of the read.
     * @param readEnd exclusive end position just beyond the targeted section of the read.
     * @param haplotypes in/out set of haplotypes.
     */
    public void calculateLocalLikelihoods(final int readStart, final int readEnd, final PairHMMReadyHaplotypes haplotypes) {
        final PairHMMReadyHaplotypes.Iterator entryIterator = haplotypes.iterator();
        boolean isFirst = true;
        while (entryIterator.hasNext()) {
            entryIterator.next();
            final int startIndex = entryIterator.startIndex();
            final byte[] bases = entryIterator.bases();
            changeHaplotypeSuffix(startIndex,bases,startIndex,bases.length);
            final double likelihood = calculateLikelihood(readStart, readEnd, startIndex, isFirst);
            isFirst = false;
            entryIterator.setLikelihood(likelihood);
        }
    }



    @Override
    public double calculateLocalLikelihood(final int readStart, final int readEnd,
                                           final int hapStart, final int hapEnd, final boolean kmerMatch) {
        if (readBases == null || haplotypeBases == null)
            throw new IllegalStateException("read or haplotype was not loaded");
        final int hapSegmentLength = hapEnd - hapStart;
        final int readSegmentLength = readEnd - readStart;
        // trivial case when there is a single base match.
        if (kmerMatch) {
            return calculateLocalLikelihoodsExactMatch(readStart, hapStart, hapSegmentLength, readSegmentLength);
        } else if (hapSegmentLength == readSegmentLength) {
            if (hapSegmentLength == 0) {
                return calculateLocalLikelihoodEmptySquare(readStart, readEnd);
            } else if (hapSegmentLength == 1) {
                return calculateLocalLikelihoodSingleBase(readStart, readEnd, hapStart);
            } else { // general (slower) solution.
                return calculateLocalLikelihoodsGeneral(readStart, readEnd, hapStart, hapEnd);
            }
        } else if (hapSegmentLength == 0) { // must be full insertion we
            return calculateLocalLikelihoodInsertion(readStart, readEnd);
        } else if (readSegmentLength == 0) { // full deletion.
            return calculateLocalLikelihoodDeletion(readStart, hapStart, hapEnd);
        } else { // general (slower) solution.
            return calculateLocalLikelihoodsGeneral(readStart, readEnd, hapStart, hapEnd);
        }
    }

    /**
     * Fast likelihood when the pair-hmm represents a deletion in the read.
     */
    private double calculateLocalLikelihoodDeletion(final int readStart, final int hapStart, final int hapEnd) {
        if (readStart == 0 || readStart >= readBases.length) // Deletion at the beginning or end not have a cost.
            return 0;
        return logTransition[readStart + 1][matchToDeletion]
                    + logTransition[readStart + 1][deletionToDeletion] * (hapEnd - hapStart - 1);
    }

    /**
     * Fast likelihood when the pair-hmm represents a insertion in the read.
     */
    private double calculateLocalLikelihoodInsertion(final int readStart, final int readEnd) {
        double result = logTransition[readStart + 1][matchToInsertion];

        if (indelToIndelIsConstant)
            result += logTransition[readStart + 1][insertionToInsertion] * (readEnd - readStart - 1);
        else
            for (int i = readStart + 1; i < readEnd; i++)
                result += logTransition[i + 1][insertionToInsertion];

        if (readEnd < readBases.length)
            result += logTransition[readEnd + 1][indelToMatch];
        return result;
    }

    /**
     * Single base mismatch fast likelihood calculation.
     */
    private double calculateLocalLikelihoodSingleBase(final int readStart, final int readEnd, final int hapStart) {
        double result = INITIAL_CONDITION;
        result *= prior[readStart + 1][hapStart + 1];
        if (readStart > 0) {
            result *= transition[readStart + 1][matchToMatch];
        }
        if (readEnd < readBases.length) {
            result *= transition[readEnd + 1][matchToMatch];
        }
        return StrictMath.log10(result) - INITIAL_CONDITION_LOG10;
    }

    /**
     * Empty square Pair-hmm.
     */
    private double calculateLocalLikelihoodEmptySquare(final int readStart, final int readEnd) {
        double result = INITIAL_CONDITION;
        if (readStart > 0 && readEnd < readBases.length) {
            result *= transition[readStart + 1][matchToMatch];
        }
        return StrictMath.log10(result) - INITIAL_CONDITION_LOG10;
    }

    /**
     * Likelihood assuming that there is a exact match between both sequences: read and haplotype
     */
    private double calculateLocalLikelihoodsExactMatch(final int readStart, final int hapStart, final int hapSegmentLength, final int readSegmentLength) {
        double result = INITIAL_CONDITION;
        if (hapSegmentLength == 1) {
            result *= prior[readStart + 1][hapStart + 1];
        } else {
            for (int i = 0; i < readSegmentLength; i++) {
                result *= prior[readStart + i + 1][hapStart + i + 1];
                if (i > 0) {
                    result *= transition[readStart + i + 1][matchToMatch];
                }
            }
        }
        return StrictMath.log10(result) - INITIAL_CONDITION_LOG10;
    }

    /**
     * Revert to a general pair-hmm solution.
     */
    private double calculateLocalLikelihoodsGeneral(final int readStart, final int readEnd, final int hapStart, final int hapEnd) {
        final Problem p = new Problem(readStart, readEnd, hapStart, hapEnd);
        final Double cachedCost = cachedResults.get(p);
        if (cachedCost != null) {
            return cachedCost;
        }
        double cost = calculateLocalLikelihoodGeneral(p);
        cachedResults.put(p, cost);
        return cost;
    }

    /**
     * Resolve the regular full pair-hmm.
     *
     * <p>
     * With the possibility of reuse the previous haplotype common prefix by using
     * a startIndex which is greater than 0.
     */
    private double calculateLikelihood(final int readStart, final int readEnd, final int startIndex, final boolean initializeEdges) {
        final int edgeStart = initializeEdges ? 0 : startIndex + 1;
        initializeMatrixValuesForTrailingProblem(readStart, readEnd, edgeStart);
        updateTable(readStart + 1, readEnd + 1, startIndex + 1, haplotypeLength + 1);
        if (readEnd == readBases.length)
            return finalLikelihoodCalculation(readEnd,0,haplotypeLength + 1) - (readStart == 0 ? StrictMath.log10(haplotypeLength) : 0);
        else {
            final double divider = 3.0;
            final double dividerInverted = 1.0 / divider;
            return StrictMath.log10(matchMatrix[readEnd][haplotypeLength]
                * transition[readEnd][matchToMatch] * dividerInverted +
                insertionMatrix[readEnd][haplotypeLength]
                        * transition[readEnd][indelToMatch] * dividerInverted +
                deletionMatrix[readEnd][haplotypeLength]
                        * transition[readEnd][indelToMatch] * dividerInverted) - INITIAL_CONDITION_LOG10 + StrictMath.log10(divider);
        }
    }


    private double calculateLocalLikelihoodGeneral(final Problem p) {

        initializeMatrixValues(p,null);
     // int fromCol = p.hapStart + 1;
     //   if (previousProblem == null) {
     //       fromCol = p.hapStart + 1;
     //   } else {
     //       final int sharedPrefix = previousProblem.followerStartIndex(p);
     //       if (sharedPrefix >= 0)
     //           fromCol = sharedPrefix + 1;
     //       else
     //           fromCol = p.hapStart + 1;
     //   }
     //   previousProblem = p;

        updateTable(p.readStart + 1, p.readEnd + 1,
                p.hapStart + 1, p.hapEnd + 1);

        if (p.trailing) {
            return finalLikelihoodCalculation(p.readEnd,p.hapStart,p.hapEnd + 1)
                    - (p.leading ? StrictMath.log10(p.hapEnd - p.hapStart) : 0);
        } else {
            final double divider = 3.0;
            final double dividerInverted = 1.0 / divider;
            return StrictMath.log10(matchMatrix[p.readEnd][p.hapEnd]
                    * transition[p.readEnd][matchToMatch] * dividerInverted +
                    insertionMatrix[p.readEnd][p.hapEnd]
                            * transition[p.readEnd][indelToMatch] * dividerInverted +
                    deletionMatrix[p.readEnd][p.hapEnd]
                            * transition[p.readEnd][indelToMatch] * dividerInverted) - INITIAL_CONDITION_LOG10 + StrictMath.log10(divider);
        }
    }

    private void updateTable(final int rowFrom, final int rowTo,
                             final int colFrom, final int colTo) {

        for (int i = rowFrom; i < rowTo; i++) {
            for (int j = colFrom; j < colTo; j++) {
                updateCell(i, j, prior[i][j], transition[i]);
            }
        }

    }

    /**
     * Holds the properties of a pair-hmm computational problem.
     */
    public class Problem {
        private final byte[] haplotypeSegment;
        private final int readStart;
        private final int readEnd;
        private final int hapStart;
        private final int hapEnd;
        private final int hashCode;
        private final boolean trailing;
        private final boolean leading;

        /**
         * Construct a new project object.
         * @param start inclusive start position on the read to consider.
         * @param end exclusive after last position on the read to consider.
         * @param hapStart inclusive start position on the haplotype to consider.
         * @param hapEnd exclusive after last position on the haplotype to consider.
         */
        public Problem(final int start, final int end, final int hapStart,
                       final int hapEnd) {
            if (start < 0 || start > readBases.length)
                throw new IllegalArgumentException("bad start index " + start);
            if (end < start || end > readBases.length)
                throw new IllegalArgumentException("bad end index " + end + " < " + start + " or " + end + " > " + readBases.length);
            if (hapStart < 0 || hapStart > haplotypeLength)
                throw new IllegalArgumentException("bad hap start index "
                        + hapStart + " is larger than the haplotypeLength " + haplotypeLength);
            if (hapEnd < hapStart || hapEnd > haplotypeLength)
                throw new IllegalArgumentException("bad hap end index "
                        + hapEnd + " outside [" + hapStart + ","
                        + haplotypeLength + "]");

            haplotypeSegment = Arrays.copyOfRange(haplotypeBases, hapStart, hapEnd);
            readStart = start;
            readEnd = end;
            this.hapStart = hapStart;
            this.hapEnd = hapEnd;
            trailing = readEnd == readBases.length;
            leading = readStart == 0;

            hashCode = ((start * 31 + end) * 31 + Arrays.hashCode(haplotypeSegment) * 31);
        }

        @Override
        public int hashCode() {
            return hashCode;
        }

        @Override
        public boolean equals(Object o) {
            if (o == this)
                return true;
            else if (o == null)
                return false;
            else if (o.getClass() != this.getClass())
                return false;
            else {
                final Problem p = (Problem) o;
                return (p.hashCode == this.hashCode) && (p.readStart == this.readStart) && (p.readEnd == this.readEnd) && Arrays.equals(haplotypeSegment, p.haplotypeSegment);
            }
        }


    }

    /**
     * Returns the currently loaded read base calls.
     * @return {@code never null}.
     */
    public byte[] getReadBases() {
        if (readBases == null)
            throw new IllegalStateException("no read was previously loaded.");
        return readBases;
    }


}
