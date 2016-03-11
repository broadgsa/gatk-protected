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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import htsjdk.variant.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class HaplotypeIndelErrorModel {

    private final int maxReadDeletionLength; // maximum length of deletion on a read
    private final double noDeletionProbability; // alpha for geometric distribution for deletion length
    private final int haplotypeSize;
    
    private final int BASE_QUAL_THRESHOLD = 6;

    private final int PATH_METRIC_TABLE_LENGTH;
    private final int RIGHT_ALIGN_INDEX;
    private final int LEFT_ALIGN_INDEX;


    private double deletionErrorProbabilities[];

    private double pathMetricArray[][];
    private int bestStateIndexArray[][];

    private final double logOneMinusInsertionStartProbability;
    private final double logInsertionStartProbability;
    private final double logInsertionEndProbability;
    private final double logOneMinusInsertionEndProbability;

    private boolean DEBUG = false;
    private boolean doSimpleCalculationModel = false;

    private static final double QUAL_ONE_HALF = -10*Math.log10(0.5);

    private static final int MAX_CACHED_QUAL = 60;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    static {
        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = QualityUtils.qualToProb((byte)k);


            baseMatchArray[k] =  probToQual(baseProb);
            baseMismatchArray[k] = (double)(k);
        }
    }

    public  HaplotypeIndelErrorModel(int mrdl, double insStart, double insEnd, double alpha, int haplotypeSize,
                                     boolean dosimple, boolean deb) {
        this(mrdl, insStart, insEnd, alpha, haplotypeSize);
        this.DEBUG = deb;
        this.doSimpleCalculationModel = dosimple;
    }
    public  HaplotypeIndelErrorModel(int mrdl, double insStart, double insEnd, double alpha, int haplotypeSize) {
        this.maxReadDeletionLength = mrdl;
        this.noDeletionProbability = 1-alpha;
        this.haplotypeSize = haplotypeSize;

        PATH_METRIC_TABLE_LENGTH = haplotypeSize+2;
        RIGHT_ALIGN_INDEX = PATH_METRIC_TABLE_LENGTH-1;
        LEFT_ALIGN_INDEX = 0;

        logOneMinusInsertionStartProbability = probToQual(1-insStart);
        logInsertionStartProbability = probToQual(insStart);
        logInsertionEndProbability = probToQual(insEnd);
        logOneMinusInsertionEndProbability = probToQual(1-insEnd);


        // fill in probability for read deletion of length k = C*exp(k-1)
        double prob = 1.0;
        double sumProb = 0.0;

        deletionErrorProbabilities = new double[maxReadDeletionLength+1];

        deletionErrorProbabilities[1] = noDeletionProbability;
        for (int k=2; k <= maxReadDeletionLength; k++) {
            deletionErrorProbabilities[k] = prob;
            sumProb = sumProb + prob;
            prob = prob*Math.exp(-1);
        }

        // now get everything in log domain, set normalizing constant so that probabilities sum to one
        deletionErrorProbabilities[1] = probToQual(deletionErrorProbabilities[1]);
        for (int k=2; k <= maxReadDeletionLength; k++) {
            deletionErrorProbabilities[k] = probToQual((1-noDeletionProbability)*deletionErrorProbabilities[k]/sumProb);
        }



    }

    public static double probToQual(double prob) {
        // TODO: see if I can use QualityUtils version, right now I don't want rounding or byte conversion
        return -10.0*Math.log10(prob);
    }

    public double computeReadLikelihoodGivenHaplotype(Haplotype haplotype, SAMRecord read) {

        long numStartClippedBases = 0;
        long numEndClippedBases = 0;


        byte[] unclippedReadQuals = read.getBaseQualities();
        byte[] unclippedReadBases = read.getReadBases();


        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i=0; i < read.getReadLength(); i++) {
            if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                numStartClippedBases++;
            else
                break;

        }
        for (int i=read.getReadLength()-1; i >= 0; i-- ){
            if (unclippedReadQuals[i] < BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }
        //System.out.format("numstart: %d numend: %d\n", numStartClippedBases, numEndClippedBases);
        if (numStartClippedBases + numEndClippedBases >= read.getReadLength()) {
            return 0;///Double.POSITIVE_INFINITY;
        }
        byte[] readBases = Arrays.copyOfRange(unclippedReadBases,(int)numStartClippedBases,
                (int)(read.getReadBases().length-numEndClippedBases));

        byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,(int)numStartClippedBases,
                (int)(read.getReadBases().length-numEndClippedBases));


        int readLength = readBases.length;

        // initialize path metric and traceback memories for Viterbi computation
        pathMetricArray = new double[readLength+1][PATH_METRIC_TABLE_LENGTH];
        bestStateIndexArray = new int[readLength+1][PATH_METRIC_TABLE_LENGTH];

        for (int k=1; k < PATH_METRIC_TABLE_LENGTH; k++)
            pathMetricArray[0][k] = 0;

 /*

         if (doSimpleCalculationModel) {

            // No Viterbi algorithm - assume no sequencing indel artifacts,

            // so we can collapse computations and pr(read | haplotype) is just probability of observing overlap
            // of read with haplotype.
            int haplotypeIndex = initialIndexInHaplotype;
            double c =  0.0;//deletionErrorProbabilities[1] +logOneMinusInsertionStartProbability;
            // compute likelihood of portion of base to the left of the haplotype
            for (int indR=readStartIdx-1; indR >= 0; indR--) {
                byte readBase = readBases[indR];
                byte readQual = readQuals[indR];
                if (readQual <= 2)
                    continue;
                double pBaseRead = getProbabilityOfReadBaseGivenXandI((byte)0, readBase, readQual, LEFT_ALIGN_INDEX, 0);

                // pBaseRead has -10*log10(Prob(base[i]|haplotype[i])
                pRead += pBaseRead;

            }
            //System.out.format("\nSt: %d Pre-Likelihood:%f\n",readStartIdx, pRead);

            for (int indR=readStartIdx; indR < readBases.length; indR++) {
                byte readBase = readBases[indR];
                byte readQual = readQuals[indR];

                byte haplotypeBase;
                if (haplotypeIndex < RIGHT_ALIGN_INDEX)
                    haplotypeBase = haplotype.getBases()[haplotypeIndex];
                else
                    haplotypeBase = (byte)0; // dummy

                double pBaseRead = getProbabilityOfReadBaseGivenXandI(haplotypeBase, readBase, readQual, haplotypeIndex, 0);
                if (haplotypeBase != 0)
                    pBaseRead += c;

                // pBaseRead has -10*log10(Prob(base[i]|haplotype[i])
                if (readQual > 3)
                    pRead += pBaseRead;
                haplotypeIndex++;
                if (haplotypeIndex >= haplotype.getBases().length)
                    haplotypeIndex = RIGHT_ALIGN_INDEX;
                //System.out.format("H:%c R:%c RQ:%d HI:%d %4.5f %4.5f\n", haplotypeBase, readBase, (int)readQual, haplotypeIndex, pBaseRead, pRead);
             }
            //System.out.format("\nSt: %d Post-Likelihood:%f\n",readStartIdx, pRead);

            if (DEBUG) {
                System.out.println(read.getReadName());
                System.out.print("Haplotype:");

                for (int k=0; k <haplotype.getBases().length; k++) {
                    System.out.format("%c ", haplotype.getBases()[k]);
                }
                System.out.println();

                System.out.print("Read bases: ");
                for (int k=0; k <readBases.length; k++) {
                    System.out.format("%c ", readBases[k]);
                }
                System.out.format("\nLikelihood:%f\n",pRead);

            }

            if (read.getReadName().contains("106880")) {

                System.out.println("aca");

                System.out.println("Haplotype:");

                for (int k=initialIndexInHaplotype; k <haplotype.getBases().length; k++) {
                    System.out.format("%c ", haplotype.getBases()[k]);
                }
                System.out.println();

                System.out.println("Read bases: ");
                for (int k=readStartIdx; k <readBases.length; k++) {
                    System.out.format("%c ", readBases[k]);
                }

            }
            return pRead;

        }
                */


        // Update path metric computations based on branch metric (Add/Compare/Select operations)
        // do forward direction first, ie from anchor to end of read
        // outer loop
        for (int indR=0; indR < readLength; indR++) {
            byte readBase = readBases[indR];
            byte readQual = readQuals[indR];

             for (int indX=LEFT_ALIGN_INDEX; indX <= RIGHT_ALIGN_INDEX; indX++) {


                byte haplotypeBase;
                if (indX > LEFT_ALIGN_INDEX && indX < RIGHT_ALIGN_INDEX)
                    haplotypeBase = haplotype.getBases()[indX-1];
                else
                    haplotypeBase = readBase;

                updatePathMetrics(haplotypeBase, indX, indR, readBase, readQual);
            }
        }



        // for debugging only: compute backtracking to find optimal route through trellis. Since I'm only interested
        // in log-likelihood of best state, this isn't really necessary.
        double bestMetric = MathUtils.arrayMin(pathMetricArray[readLength]);

        if (DEBUG) {



            System.out.println(read.getReadName());
            System.out.print("Haplotype:");

            for (int k=0; k <haplotype.getBases().length; k++) {
                System.out.format("%c ", haplotype.getBases()[k]);
            }
            System.out.println();

            System.out.print("Read bases: ");
            for (int k=0; k <readBases.length; k++) {
                System.out.format("%c ", readBases[k]);
            }
            System.out.println();

            System.out.print("Read quals: ");
            for (int k=0; k <readQuals.length; k++) {
                System.out.format("%d ", (int)readQuals[k]);
            }
            System.out.println();

            // start from last position of read, go backwards to find optimal alignment
            int[] bestIndexArray = new int[readLength];
            int bestIndex = MathUtils.minElementIndex(pathMetricArray[readLength]);
            bestIndexArray[readLength-1] = bestIndex;

            for (int k=readLength-2; k>=0; k--) {
                bestIndex = bestStateIndexArray[k][bestIndex];
                bestIndexArray[k] = bestIndex;
            }
            
            System.out.print("Alignment: ");
            for (int k=0; k <readBases.length; k++) {
                System.out.format("%d ", bestIndexArray[k]);
            }
            System.out.println();
        }
        // now just take optimum along all path metrics: that's the log likelihood of best alignment
        if (DEBUG)
            System.out.format("Likelihood: %5.4f\n", bestMetric);
        return bestMetric;

    }

    private void updatePathMetrics(byte haplotypeBase, int indX, int indR, byte readBase, byte readQual) {

        double bmetric;

        double bestMetric = Double.POSITIVE_INFINITY;
        int bestMetricIndex = -1;

        // compute metric for match/mismatch
        double pBaseRead;

        // workaround for reads whose bases quality = 0,
        if (readQual < 1)
            readQual = 1;

        if (readQual > MAX_CACHED_QUAL)
            readQual = MAX_CACHED_QUAL;

        double pBaseMatch =  baseMatchArray[(int)readQual];
        double pBaseMismatch = baseMismatchArray[(int)readQual];

        if (haplotypeBase == readBase)
            pBaseRead =  pBaseMatch;
        else
            pBaseRead = pBaseMismatch;





        if (indX == LEFT_ALIGN_INDEX) {
            // special treatment for L->L case  (read aligns to the right/left of haplotype) when Xold = X = R
            bestMetric = pathMetricArray[indR][indX] + QUAL_ONE_HALF; //1/2 in log scale
            bestMetricIndex = indX;
        }
        else {
            for (int indXold = indX-1; indXold >= indX-this.maxReadDeletionLength; indXold--) {

                if (indXold < 0)
                    break;

                // fetch path metric and add branch metric
                bmetric = pathMetricArray[indR][indXold] + deletionErrorProbabilities[indX-indXold] + pBaseRead;

                if (indXold == indX-1) {
                    // special case for exact walk down diagonal: need to consider that an insertion may have ended
               //     bmetric += logInsertionEndProbability;
                } else {
                    bmetric += logOneMinusInsertionStartProbability;
                }

                if (bmetric < bestMetric) {
                    bestMetric = bmetric;
                    bestMetricIndex = indXold;
                }
            }
            // additional case: a walk right (ie indXold = indX). Can be because of an insertion in the middle of reads,
            // or we're aligned to right of read
            bmetric = pathMetricArray[indR][indX]+pBaseRead;
            if (indX < RIGHT_ALIGN_INDEX) {
                bmetric += logInsertionStartProbability + logOneMinusInsertionEndProbability;
            }
            else {
                // anything extra to do for R->R case?
                bmetric = pathMetricArray[indR][indX] + QUAL_ONE_HALF;
            }

            if (bmetric < bestMetric) {
                bestMetric = bmetric;
                bestMetricIndex = indX;
            }


        }

         // record best path metric
        pathMetricArray[indR+1][indX] = bestMetric;
        bestStateIndexArray[indR+1][indX] = bestMetricIndex;

    }

    public double[] computeReadHaplotypeLikelihoods(ReadBackedPileup pileup, HashMap<Allele,Haplotype> haplotypesInVC){
        double[][] haplotypeLikehoodMatrix = new double[haplotypesInVC.size()][haplotypesInVC.size()];
        double readLikelihoods[][] = new double[pileup.getReads().size()][haplotypesInVC.size()];
        int i=0;
        for (GATKSAMRecord read : pileup.getReads()) {
            if(ReadUtils.is454Read(read)) {
                continue;
            }
            // for each read/haplotype combination, compute likelihoods, ie -10*log10(Pr(R | Hi))
            // = sum_j(-10*log10(Pr(R_j | Hi) since reads are assumed to be independent
            int j=0;
            for (Map.Entry<Allele,Haplotype> a: haplotypesInVC.entrySet()) {
                readLikelihoods[i][j]= computeReadLikelihoodGivenHaplotype(a.getValue(), read);
                if (DEBUG) {
                    System.out.print(read.getReadName()+" ");

                    System.out.format("%d %d S:%d US:%d E:%d UE:%d C:%s %3.4f\n",i, j, read.getAlignmentStart(),
                            read.getUnclippedStart(), read.getAlignmentEnd(), read.getUnclippedEnd(),
                            read.getCigarString(), readLikelihoods[i][j]);
                }
                j++;
            }
            i++;
        }

        for (i=0; i < haplotypesInVC.size(); i++) {
            for (int j=i; j < haplotypesInVC.size(); j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                double[] readLikelihood = new double[2]; // diploid sample
                for (int readIdx=0; readIdx < pileup.getReads().size(); readIdx++) {
                    readLikelihood[0] = -readLikelihoods[readIdx][i]/10;
                    readLikelihood[1] = -readLikelihoods[readIdx][j]/10;

                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+x0^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    // Second term is a constant added to both likelihoods so will be ignored
                    haplotypeLikehoodMatrix[i][j] += MathUtils.approximateLog10SumLog10(readLikelihood[0], readLikelihood[1]);

                }


            }
        }

        return getHaplotypeLikelihoods(haplotypeLikehoodMatrix);

    }

    private double[] getHaplotypeLikelihoods(double[][] haplotypeLikehoodMatrix) {
        int hSize = haplotypeLikehoodMatrix.length;
        double[] genotypeLikelihoods = new double[hSize*(hSize+1)/2];

        int k=0;
        double maxElement = Double.NEGATIVE_INFINITY;
        for (int i=0; i < hSize; i++) {
            for (int j=i; j < hSize; j++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
                if (haplotypeLikehoodMatrix[i][j] > maxElement)
                    maxElement = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize
        for (int i=0; i < genotypeLikelihoods.length; i++)
            genotypeLikelihoods[i] -= maxElement;

        return genotypeLikelihoods;

    }


}
