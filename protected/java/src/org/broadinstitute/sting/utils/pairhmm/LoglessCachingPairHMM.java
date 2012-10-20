/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.pairhmm;

import org.broadinstitute.sting.utils.QualityUtils;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin, carneiro
 * Date: 10/16/12
 */

public class LoglessCachingPairHMM extends CachingPairHMM {

    protected static final double SCALE_FACTOR_LOG10 = 300.0;

    protected static final double [] firstRowConstantMatrix = {
            QualityUtils.qualToProb((byte) (DEFAULT_GOP + DEFAULT_GOP)),
            QualityUtils.qualToProb(DEFAULT_GCP),
            QualityUtils.qualToErrorProb(DEFAULT_GOP),
            QualityUtils.qualToErrorProb(DEFAULT_GCP),
            1.0,
            1.0
    };

    @Override
    public void initialize( final int READ_MAX_LENGTH, final int HAPLOTYPE_MAX_LENGTH ) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = READ_MAX_LENGTH + 2;
        final int Y_METRIC_LENGTH = HAPLOTYPE_MAX_LENGTH + 2;

        matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        for( int iii=0; iii < X_METRIC_LENGTH; iii++ ) {
            Arrays.fill(matchMetricArray[iii], 0.0);
            Arrays.fill(XMetricArray[iii], 0.0);
            Arrays.fill(YMetricArray[iii], 0.0);
        }

        // the initial condition
        matchMetricArray[1][1] = Math.pow(10.0, SCALE_FACTOR_LOG10); // Math.log10(1.0);

        constantMatrix = new double[X_METRIC_LENGTH][6];
        distanceMatrix = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        // fill in the first row
        for( int jjj = 2; jjj < Y_METRIC_LENGTH; jjj++ ) {
            updateCell(1, jjj, 1.0, firstRowConstantMatrix, matchMetricArray, XMetricArray, YMetricArray);
        }
    }

    @Override
    public double computeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                            final byte[] readBases,
                                                            final byte[] readQuals,
                                                            final byte[] insertionGOP,
                                                            final byte[] deletionGOP,
                                                            final byte[] overallGCP,
                                                            final int hapStartIndex,
                                                            final boolean recacheReadValues ) {

        if( recacheReadValues ) {
            initializeConstants( insertionGOP, deletionGOP, overallGCP );
        }
        initializeDistanceMatrix( haplotypeBases, readBases, readQuals, hapStartIndex );

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = readBases.length + 2;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 2;

        for (int i = 2; i < X_METRIC_LENGTH; i++) {
            for (int j = hapStartIndex+1; j < Y_METRIC_LENGTH; j++) {
                updateCell(i, j, distanceMatrix[i][j], constantMatrix[i], matchMetricArray, XMetricArray, YMetricArray);
            }
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = X_METRIC_LENGTH - 1;
        final int endJ = Y_METRIC_LENGTH - 1;
        return Math.log10( matchMetricArray[endI][endJ] + XMetricArray[endI][endJ] + YMetricArray[endI][endJ] ) - SCALE_FACTOR_LOG10;
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
    public void initializeDistanceMatrix( final byte[] haplotypeBases,
                                          final byte[] readBases,
                                          final byte[] readQuals,
                                          final int startIndex ) {

        // initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                distanceMatrix[i+2][j+2] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProb(qual) : QualityUtils.qualToErrorProb(qual) );
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
    public void initializeConstants( final byte[] insertionGOP,
                                     final byte[] deletionGOP,
                                     final byte[] overallGCP ) {

        final int l = insertionGOP.length;
        constantMatrix[1] = firstRowConstantMatrix;
        for (int i = 0; i < l; i++) {
            final int qualIndexGOP = Math.min(insertionGOP[i] + deletionGOP[i], Byte.MAX_VALUE);
            constantMatrix[i+2][0] = QualityUtils.qualToProb((byte) qualIndexGOP);
            constantMatrix[i+2][1] = QualityUtils.qualToProb(overallGCP[i]);
            constantMatrix[i+2][2] = QualityUtils.qualToErrorProb(insertionGOP[i]);
            constantMatrix[i+2][3] = QualityUtils.qualToErrorProb(overallGCP[i]);
            constantMatrix[i+2][4] = QualityUtils.qualToErrorProb(deletionGOP[i]);
            constantMatrix[i+2][5] = QualityUtils.qualToErrorProb(overallGCP[i]);
        }
        constantMatrix[l+1][4] = 1.0;
        constantMatrix[l+1][5] = 1.0;
    }

    /**
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions

     * @param indI             row index in the matrices to update
     * @param indJ             column index in the matrices to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param constants        an array with the six constants relevant to this location
     * @param matchMetricArray the matches likelihood matrix
     * @param XMetricArray     the insertions likelihood matrix
     * @param YMetricArray     the deletions likelihood matrix
     */
    private void updateCell( final int indI, final int indJ, final double prior, final double[] constants,
                             final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        matchMetricArray[indI][indJ] = prior * ( matchMetricArray[indI - 1][indJ - 1] * constants[0] +
                                                 XMetricArray[indI - 1][indJ - 1] * constants[1] +
                                                 YMetricArray[indI - 1][indJ - 1] * constants[1] );
        XMetricArray[indI][indJ] = matchMetricArray[indI - 1][indJ] * constants[2] + XMetricArray[indI - 1][indJ] * constants[3];
        YMetricArray[indI][indJ] = matchMetricArray[indI][indJ - 1] * constants[4] + YMetricArray[indI][indJ - 1] * constants[5];
    }
}
