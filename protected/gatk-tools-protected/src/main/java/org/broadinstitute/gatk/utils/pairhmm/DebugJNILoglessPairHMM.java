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
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.broadinstitute.gatk.utils.pairhmm.PairHMMModel.*;


/**
 * Created with IntelliJ IDEA.
 * User: rpoplin, carneiro
 * Date: 10/16/12
 */
public class DebugJNILoglessPairHMM extends LoglessPairHMM {

    private static final boolean dumpSandboxOnly = false;        //simulates ifdef
    private static final boolean debug = false;  //simulates ifdef
    private static final boolean verify = !dumpSandboxOnly && (debug || true);        //simulates ifdef
    private static final boolean debug0_1 = false;       //simulates ifdef
    private static final boolean debug1 = false; //simulates ifdef
    private static final boolean debug2 = false;
    private static final boolean debug3 = false;
   
    //Debugging stats 
    private int numCalls = 0;
    private int numComputeLikelihoodCalls = 0;
    protected HashMap<String, BufferedWriter> filenameToWriter = new HashMap<String, BufferedWriter>();

    private JNILoglessPairHMM jniPairHMM = null;
    public DebugJNILoglessPairHMM(final PairHMM.HMM_IMPLEMENTATION hmmType, PairHMM.HMM_SUB_IMPLEMENTATION pairHMMSub, final boolean alwaysLoadVectorLoglessPairHMMLib) {
        super();
        switch(hmmType) {
            case VECTOR_LOGLESS_CACHING:
                jniPairHMM = new VectorLoglessPairHMM(pairHMMSub, alwaysLoadVectorLoglessPairHMMLib);
                break;
            default:
                throw new UserException.BadArgumentValue("pairHMM","Specified JNIPairHMM implementation is unrecognized or incompatible with the HaplotypeCaller. Acceptable options are VECTOR_LOGLESS_CACHING");
        }
    }

    @Override
    public void close()
    {
        jniPairHMM.close();
        debugClose();
    }

    //Used only when testing parts of the compute kernel
    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) {
        if(verify)
            super.initialize(readMaxLength, haplotypeMaxLength);
        if(debug3)
        {
            System.out.println("Java: alloc initialized readMaxLength : "+readMaxLength+" haplotypeMaxLength : "+haplotypeMaxLength);
            debugDump("lengths_java.txt", String.format("%d %d\n",readMaxLength, haplotypeMaxLength),
                    true);
        }
        if(debug2)
            jniInitialize(readMaxLength, haplotypeMaxLength);
    }

    private HashMap<Haplotype,Integer> haplotypeToHaplotypeListIdxMap = null;
    //Used to transfer data to JNI
    //Since the haplotypes are the same for all calls to computeLikelihoods within a region, transfer the haplotypes only once to the JNI per region
    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize( final List<Haplotype> haplotypes, final Map<String, List<GATKSAMRecord>> perSampleReadList,
            final int readMaxLength, final int haplotypeMaxLength ) {
        if(verify)
        {
            super.initialize(haplotypes, perSampleReadList, readMaxLength, haplotypeMaxLength);
            jniPairHMM.initialize(haplotypes, perSampleReadList, readMaxLength, haplotypeMaxLength);
            haplotypeToHaplotypeListIdxMap = jniPairHMM.getHaplotypeToHaplotypeListIdxMap();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void finalizeRegion()
    {
        if(!dumpSandboxOnly)
            jniPairHMM.finalizeRegion();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void computeLikelihoods( final ReadLikelihoods.Matrix<Haplotype> likelihoods, final List<GATKSAMRecord> processedReads, final Map<GATKSAMRecord, byte[]> GCPArrayMap ) {
        // (re)initialize the pairHMM only if necessary
        final int readMaxLength = verify ? findMaxReadLength(processedReads) : 0;
        final int haplotypeMaxLength = verify ? findMaxHaplotypeLength(likelihoods.alleles()) : 0;
        if(verify)
        {
            if (!initialized || readMaxLength > maxReadLength || haplotypeMaxLength > maxHaplotypeLength)
            { initialize(readMaxLength, haplotypeMaxLength); }
            if ( ! initialized )
                throw new IllegalStateException("Must call initialize before calling jniComputeLikelihoods in debug/verify mode");
        }
        int readListSize = processedReads.size();
        int numHaplotypes = likelihoods.alleleCount();
        int numTestcases = readListSize*numHaplotypes;
        if(debug0_1)
            System.out.println("Java numReads "+readListSize+" numHaplotypes "+numHaplotypes);
        int idx = 0;
        for(final GATKSAMRecord read : processedReads)
        {
            final byte [] overallGCP = GCPArrayMap.get(read);
            if(debug0_1)
                System.out.println("Java read length "+read.getReadBases().length);
            if(debug3)
            {
                for(int i=0;i<read.getReadBases().length;++i)
                {
                    debugDump("reads_java.txt",String.format("%d\n",(int)read.getReadBases()[i]),true);
                    debugDump("reads_java.txt",String.format("%d\n",(int)read.getBaseQualities()[i]),true);
                    debugDump("reads_java.txt",String.format("%d\n",(int)read.getBaseInsertionQualities()[i]),true);
                    debugDump("reads_java.txt",String.format("%d\n",(int)read.getBaseDeletionQualities()[i]),true);
                    debugDump("reads_java.txt",String.format("%d\n",(int)overallGCP[i]),true);
                }
            }
            ++idx;
        }

        if(verify)
        {
            idx = 0;
            for (final Haplotype h : likelihoods.alleles())       //order is important - access in same order always
            {
                byte[] haplotypeBases = h.getBases();
                if(debug0_1)
                    System.out.println("Java haplotype length "+haplotypeBases.length);
                if(debug3)
                {
                    for(int i=0;i<haplotypeBases.length;++i)
                        debugDump("haplotype_bases_java.txt",String.format("%d\n",(int)haplotypeBases[i]),true);
                }
                ++idx;
            }
        }
        double[] likelihoodArray = null;
        PerReadAlleleLikelihoodMap likelihoodMap  = null;
        if(verify)
        {
            jniPairHMM.computeLikelihoods(likelihoods, processedReads, GCPArrayMap);
            likelihoodArray = jniPairHMM.getLikelihoodArray();
            //to compare values
            super.computeLikelihoods(likelihoods, processedReads, GCPArrayMap);
        }
        else
        {
            likelihoodMap = new PerReadAlleleLikelihoodMap();
            likelihoodArray = new double[numTestcases];
            for(int i=0;i<numTestcases;++i)
                likelihoodArray[i] = -0.5;
        }
        if(verify || dumpSandboxOnly)
        {
            boolean toDump = dumpSandboxOnly ? true : false;
            if(verify)
            {
                //re-order values in likelihoodArray
                double[] tmpArray = new double[numHaplotypes];
                idx = 0;
                int idxInsideHaplotypeList = 0;
                int readIdx = 0;
                for(final GATKSAMRecord read : processedReads)
                {
                    for(int j=0;j<numHaplotypes;++j)
                        tmpArray[j] = likelihoodArray[readIdx+j];
                    for (final Haplotype haplotype : likelihoods.alleles())//order is important - access in same order always
                    {
                        idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.get(haplotype);
                        likelihoodArray[idx] = tmpArray[idxInsideHaplotypeList];
                        ++idx;
                    }
                    readIdx += numHaplotypes;
                }
                //for floating point values, no exact equality
                //check whether numbers are close in terms of abs_error or relative_error
                //For very large values, relative_error is relevant
                //For very small values, abs_error is relevant
                for(int i=0;i<likelihoodArray.length;++i)
                {
                    double abs_error = Math.abs(likelihoodArray[i] - mLikelihoodArray[i]);
                    double relative_error = 0;
                    if(mLikelihoodArray[i] == 0)
                        relative_error = 0;
                    else
                        relative_error = Math.abs(abs_error/mLikelihoodArray[i]);
                    if(abs_error > 1e-5 && relative_error > 1e-5)
                    {
                        toDump = true;
                        break;
                    }
                }
            }
            //if numbers are not close, then dump out the data that produced the inconsistency
            if(toDump)
            {
                idx = 0;
                System.out.println("Dump : Java numReads "+readListSize+" numHaplotypes "+numHaplotypes);
                boolean firstLine = true;
                for(final GATKSAMRecord read : processedReads)
                {
                    byte [] overallGCP = GCPArrayMap.get(read);
                    byte[] tmpByteArray = new byte[read.getReadBases().length];
                    for (final Haplotype haplotype : likelihoods.alleles())       //order is important - access in same order always
                    {
                        byte[] haplotypeBases = haplotype.getBases();
                        debugDump("debug_dump.txt",new String(haplotypeBases)+" ",true);
                        debugDump("debug_dump.txt",new String(read.getReadBases())+" ",true);
                        for(int k=0;k<read.getReadBases().length;++k)
                            tmpByteArray[k] = (byte)((int)((read.getBaseQualities())[k]) + 33);
                        debugDump("debug_dump.txt",new String(tmpByteArray)+" ",true);
                        for(int k=0;k<read.getReadBases().length;++k)
                            tmpByteArray[k] = (byte)((int)((read.getBaseInsertionQualities())[k]) + 33);
                        debugDump("debug_dump.txt",new String(tmpByteArray)+" ",true);
                        for(int k=0;k<read.getReadBases().length;++k)
                            tmpByteArray[k] = (byte)((int)((read.getBaseDeletionQualities())[k]) + 33);
                        debugDump("debug_dump.txt",new String(tmpByteArray)+" ",true);
                        for(int k=0;k<read.getReadBases().length;++k)
                            tmpByteArray[k] = (byte)((int)(overallGCP[k]) + 33);
                        debugDump("debug_dump.txt",new String(tmpByteArray),true);
                        if(firstLine)
                        {
                            debugDump("debug_dump.txt",String.format(" %d %d\n",readListSize, numHaplotypes), true);
                            firstLine = false;
                        }
                        else
                            debugDump("debug_dump.txt","\n",true);
                        if(verify)
                            debugDump("debug_results.txt",String.format("%e %e\n",mLikelihoodArray[idx],likelihoodArray[idx]),true);
                        else
                            if(dumpSandboxOnly)
                                likelihoods.set(likelihoods.alleleIndex(haplotype),likelihoods.readIndex(read), likelihoodArray[idx]);
                        ++idx;
                    }
                }
            }   
            debugClose();
        }
        ++numComputeLikelihoodCalls;
        //if(numComputeLikelihoodCalls == 5)
        //jniPairHMM.close();
        //System.exit(0);
    }

    //Used to test parts of the compute kernel separately    
    private native void jniInitialize( final int readMaxLength, final int haplotypeMaxLength);
    private native static void jniInitializeProbabilities( final double[][] transition, final byte[] insertionGOP,
            final byte[] deletionGOP, final byte[] overallGCP);
    private native double jniInitializePriorsAndUpdateCells( boolean doInitialization, final int paddedReadLength,
            final int paddedHaplotypeLength, final byte[] readBases, final byte[] haplotypeBases, final byte[] readQuals,
            final int hapStartIndex);
    private native double jniSubComputeReadLikelihoodGivenHaplotypeLog10( final int readLength, final int haplotypeLength,
            final byte[] readBases, final byte[] haplotypeBases, final byte[] readQuals, final byte[] insertionGOP,
            final byte[] deletionGOP, final byte[] overallGCP, final int hapStartIndex);
    /**
     * {@inheritDoc}
     */
    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,  final byte[] readBases,
            final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP,
            final int hapStartIndex, final boolean recacheReadValues, final int nextHapStartIndex) {
        //System.out.println("#### START STACK TRACE ####");
        //for (StackTraceElement ste : Thread.currentThread().getStackTrace()) {
        //System.out.println(ste);
        //}
        //System.out.println("#### END STACK TRACE ####");
        //
        if(debug1)
            jniSubComputeReadLikelihoodGivenHaplotypeLog10(readBases.length, haplotypeBases.length,
                    readBases, haplotypeBases, readQuals,
                    insertionGOP, deletionGOP, overallGCP,
                    hapStartIndex);

        boolean doInitialization = (previousHaplotypeBases == null || previousHaplotypeBases.length != haplotypeBases.length);
        if (doInitialization) {
            final double initialValue = INITIAL_CONDITION / haplotypeBases.length;
            // set the initial value (free deletions in the beginning) for the first row in the deletion matrix
            for( int j = 0; j < paddedHaplotypeLength; j++ ) {
                deletionMatrix[0][j] = initialValue;
            }
        }

        if ( ! constantsAreInitialized || recacheReadValues ) {
            initializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP); 
            if(debug3)
            {
                System.out.println("Java: initializeProbabilities lengths : "+insertionGOP.length+" padded "+paddedReadLength+" "+paddedHaplotypeLength);
                for(int i=0;i<insertionGOP.length;++i)
                    for(int j=0;j<6;++j)
                        debugDump("transitions_java.txt",String.format("%e\n",transition[i+1][j]),true);
            }
            if(debug2)
                jniInitializeProbabilities(transition, insertionGOP, deletionGOP, overallGCP); 

            // note that we initialized the constants
            constantsAreInitialized = true;
        }

        if(debug3)
            System.out.println("Java: initializePriors : lengths "+readBases.length+" "+haplotypeBases.length+" padded "+paddedReadLength+" "+paddedHaplotypeLength + " doNotUseTristateCorrection "+doNotUseTristateCorrection);
        initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex);

        for (int i = 1; i < paddedReadLength; i++) {
            // +1 here is because hapStartIndex is 0-based, but our matrices are 1 based
            for (int j = hapStartIndex+1; j < paddedHaplotypeLength; j++) {
                updateCell(i, j, prior[i][j], transition[i]);
            }
        }

        // final probability is the log10 sum of the last element in the Match and Insertion state arrays
        // this way we ignore all paths that ended in deletions! (huge)
        // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
        final int endI = paddedReadLength - 1;
        double finalSumProbabilities = 0.0;
        for (int j = 1; j < paddedHaplotypeLength; j++) {
            finalSumProbabilities += matchMatrix[endI][j] + insertionMatrix[endI][j];
        }
        if(debug2)
            jniInitializePriorsAndUpdateCells(doInitialization, paddedReadLength, paddedHaplotypeLength,
                    readBases, haplotypeBases, readQuals,
                    hapStartIndex);
        if(debug)
            debugDump("return_values_java.txt",String.format("%e\n",Math.log10(finalSumProbabilities) - INITIAL_CONDITION_LOG10),true);
        ++numCalls;
        //if(numCalls == 100)
        //{
        //debugClose();
        //System.exit(0);
        //}
        return Math.log10(finalSumProbabilities) - INITIAL_CONDITION_LOG10;
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
    protected void initializePriors(final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals, final int startIndex) {

        // initialize the pBaseReadLog10 matrix for all combinations of read x haplotype bases
        // Abusing the fact that java initializes arrays with 0.0, so no need to fill in rows and columns below 2.

        if(debug3)
            System.out.println("hapStartIndex "+startIndex);

        for (int i = 0; i < readBases.length; i++) {
            final byte x = readBases[i];
            final byte qual = readQuals[i];
            for (int j = startIndex; j < haplotypeBases.length; j++) {
                final byte y = haplotypeBases[j];
                prior[i+1][j+1] = ( x == y || x == (byte) 'N' || y == (byte) 'N' ?
                        QualityUtils.qualToProb(qual) : (QualityUtils.qualToErrorProb(qual) / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION)) );
                if(debug3)
                    debugDump("priors_java.txt",String.format("%e\n",prior[i+1][j+1]),true);
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
     * Updates a cell in the HMM matrix
     *
     * The read and haplotype indices are offset by one because the state arrays have an extra column to hold the
     * initial conditions

     * @param indI             row index in the matrices to update
     * @param indJ             column index in the matrices to update
     * @param prior            the likelihood editing distance matrix for the read x haplotype
     * @param transition        an array with the six transition relevant to this location
     */
    protected void updateCell( final int indI, final int indJ, final double prior, final double[] transition) {

        matchMatrix[indI][indJ] = prior * ( matchMatrix[indI - 1][indJ - 1] * transition[matchToMatch] +
                insertionMatrix[indI - 1][indJ - 1] * transition[indelToMatch] +
                deletionMatrix[indI - 1][indJ - 1] * transition[indelToMatch] );
        insertionMatrix[indI][indJ] = matchMatrix[indI - 1][indJ] * transition[matchToInsertion] + insertionMatrix[indI - 1][indJ] * transition[insertionToInsertion];
        deletionMatrix[indI][indJ] = matchMatrix[indI][indJ - 1] * transition[matchToDeletion] + deletionMatrix[indI][indJ - 1] * transition[deletionToDeletion];
        if(debug3)
        {
            debugDump("matrices_java.txt",String.format("%e\n",matchMatrix[indI][indJ]),true);
            debugDump("matrices_java.txt",String.format("%e\n",insertionMatrix[indI][indJ]),true);
            debugDump("matrices_java.txt",String.format("%e\n",deletionMatrix[indI][indJ]),true);
        }
    }

    protected void debugDump( String filename, String s, boolean toAppend ) {
        try {
            File file = new File(filename);
            if (!file.exists())
                file.createNewFile();
            BufferedWriter currWriter = filenameToWriter.get(filename);
            if(currWriter == null)
            {
                FileWriter fw = new FileWriter(file, toAppend);
                currWriter = new BufferedWriter(fw);
                filenameToWriter.put(filename, currWriter);
            }
            currWriter.write(s);
        }
        catch(IOException e)
        {
            e.printStackTrace();
        }
    }

    protected void debugClose()  {
        for(Map.Entry<String, BufferedWriter> currEntry : filenameToWriter.entrySet()) {
            BufferedWriter currWriter = currEntry.getValue();
            try
            {
                currWriter.flush();
                currWriter.close();
            }
            catch(IOException e)
            {
                e.printStackTrace();

            }
        }
        filenameToWriter.clear();
    }
}
