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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.samtools.SAMUtils;
import htsjdk.variant.variantcontext.Allele;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.AlleleList;
import org.broadinstitute.gatk.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pairhmm.*;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.TandemRepeatFinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

public class PairHMMLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {
    private final static Logger logger = Logger.getLogger(PairHMMLikelihoodCalculationEngine.class);

    private final byte constantGCP;

    private final double log10globalReadMismappingRate;

    private final PairHMM.HMM_IMPLEMENTATION hmmType;
    private final PairHMM.HMM_SUB_IMPLEMENTATION hmmSubType;
    private final boolean alwaysLoadVectorLoglessPairHMMLib;
    private final boolean noFpga;

    private final ThreadLocal<PairHMM> pairHMMThreadLocal = new ThreadLocal<PairHMM>() {
        @Override
        protected PairHMM initialValue() {
            switch (hmmType) {
                case EXACT: return new Log10PairHMM(true);
                case ORIGINAL: return new Log10PairHMM(false);
                case LOGLESS_CACHING:
                    if (noFpga || !CnyPairHMM.isAvailable())
                        return new LoglessPairHMM();
                    else
                        return new CnyPairHMM();
                case VECTOR_LOGLESS_CACHING:
                    try
                    {
                        return new VectorLoglessPairHMM(hmmSubType, alwaysLoadVectorLoglessPairHMMLib);
                    }
                    catch(UnsatisfiedLinkError ule)
                    {
                        logger.warn("Failed to load native library for VectorLoglessPairHMM - using Java implementation of LOGLESS_CACHING");
                        return new LoglessPairHMM();
                    }
                case DEBUG_VECTOR_LOGLESS_CACHING:
                    return new DebugJNILoglessPairHMM(PairHMM.HMM_IMPLEMENTATION.VECTOR_LOGLESS_CACHING, hmmSubType, alwaysLoadVectorLoglessPairHMMLib);
                case ARRAY_LOGLESS:
                    if (noFpga || !CnyPairHMM.isAvailable())
                        return new ArrayLoglessPairHMM();
                    else
                        return new CnyPairHMM();
                default:
                    throw new UserException.BadArgumentValue("pairHMM", "Specified pairHMM implementation is unrecognized or incompatible with the HaplotypeCaller. Acceptable options are ORIGINAL, EXACT, CACHING, LOGLESS_CACHING, and ARRAY_LOGLESS.");
            }
        }
    };
//    Attempted to do as below, to avoid calling pairHMMThreadLocal.get() later on, but it resulted in a NullPointerException
//    private final PairHMM pairHMM = pairHMMThreadLocal.get();

    private final static boolean WRITE_LIKELIHOODS_TO_FILE = false;
    private final static String LIKELIHOODS_FILENAME = "likelihoods.txt";
    private final PrintStream likelihoodsStream;

    public enum PCR_ERROR_MODEL {
        /** no specialized PCR error model will be applied; if base insertion/deletion qualities are present they will be used */
        NONE(null),
        /** a most aggressive model will be applied that sacrifices true positives in order to remove more false positives */
        HOSTILE(1.0),
        /** a more aggressive model will be applied that sacrifices true positives in order to remove more false positives */
        AGGRESSIVE(2.0),
        /** a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives */
        CONSERVATIVE(3.0);

        private final Double rateFactor;

        /** rate factor is applied to the PCR error model.  Can be null to imply no correction */
        PCR_ERROR_MODEL(Double rateFactor) {
            this.rateFactor = rateFactor;
        }
        private Double getRateFactor() { return rateFactor; }
        private boolean hasRateFactor() { return rateFactor != null; }
    }

    private final PCR_ERROR_MODEL pcrErrorModel;

    /**
     * The expected rate of random sequencing errors for a read originating from its true haplotype.
     *
     * For example, if this is 0.01, then we'd expect 1 error per 100 bp.
     */
    private final static double EXPECTED_ERROR_RATE_PER_BASE = 0.02;

    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP the gap continuation penalty to use with the PairHMM
     * @param hmmType the type of the HMM to use
     * @param hmmSubType the type of the machine dependent sub-implementation of HMM to use
     * @param alwaysLoadVectorLoglessPairHMMLib always load the vector logless HMM library instead of once
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param noFpga disable FPGA acceleration
     * @param pcrErrorModel model to correct for PCR indel artifacts
     */
    public PairHMMLikelihoodCalculationEngine( final byte constantGCP, final PairHMM.HMM_IMPLEMENTATION hmmType, final PairHMM.HMM_SUB_IMPLEMENTATION hmmSubType,
                                               final boolean alwaysLoadVectorLoglessPairHMMLib, final double log10globalReadMismappingRate, final boolean noFpga, final PCR_ERROR_MODEL pcrErrorModel ) {
        this.hmmType = hmmType;
        this.hmmSubType = hmmSubType;
        this.alwaysLoadVectorLoglessPairHMMLib = alwaysLoadVectorLoglessPairHMMLib;
        this.constantGCP = constantGCP;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.noFpga = noFpga;
        this.pcrErrorModel = pcrErrorModel;

        initializePCRErrorModel();

        if ( WRITE_LIKELIHOODS_TO_FILE ) {
            try {
                likelihoodsStream = new PrintStream(new FileOutputStream(new File(LIKELIHOODS_FILENAME)));
            } catch ( FileNotFoundException e ) {
                throw new RuntimeException(e);
            }
        } else {
            likelihoodsStream = null;
        }
    }

    @Override
    public void close() {
        if ( likelihoodsStream != null ) likelihoodsStream.close();
        pairHMMThreadLocal.get().close();
    }

    private void capMinimumReadQualities(GATKSAMRecord read, byte[] readQuals, byte[] readInsQuals, byte[] readDelQuals) {
        for( int kkk = 0; kkk < readQuals.length; kkk++ ) {
            readQuals[kkk] = (byte) Math.min( 0xff & readQuals[kkk], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
            readQuals[kkk] = ( readQuals[kkk] < PairHMM.BASE_QUALITY_SCORE_THRESHOLD ? QualityUtils.MIN_USABLE_Q_SCORE : readQuals[kkk] );
            readInsQuals[kkk] = ( readInsQuals[kkk] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : readInsQuals[kkk] );
            readDelQuals[kkk] = ( readDelQuals[kkk] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : readDelQuals[kkk] );
        }
    }

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */
    private List<GATKSAMRecord> modifyReadQualities(final List<GATKSAMRecord> reads) {
        final List<GATKSAMRecord> result = new ArrayList<>(reads.size());

        for (final GATKSAMRecord read : reads) {
            final byte[] readBases = read.getReadBases();

            // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
            final byte[] readQuals = read.getBaseQualities().clone();
            final byte[] readInsQuals = read.getBaseInsertionQualities().clone();
            final byte[] readDelQuals = read.getBaseDeletionQualities().clone();

            applyPCRErrorModel(readBases, readInsQuals, readDelQuals);
            capMinimumReadQualities(read, readQuals, readInsQuals, readDelQuals);

            // Create a new copy of the read and sets its base qualities to the modified versions.
            // Pack this into a new list for return
            result.add(GATKSAMRecord.createQualityModifiedRead(read, readBases, readQuals, readInsQuals, readDelQuals));
        }
        return result;
    }

    /**
     * Initialize our pairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     *
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *
     * @param haplotypes a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    private void initializePairHMM(final List<Haplotype> haplotypes, final Map<String, List<GATKSAMRecord>> perSampleReadList) {
        int X_METRIC_LENGTH = 0;
        for( final Map.Entry<String, List<GATKSAMRecord>> sample : perSampleReadList.entrySet() ) {
            for( final GATKSAMRecord read : sample.getValue() ) {
                final int readLength = read.getReadLength();
                if( readLength > X_METRIC_LENGTH ) { X_METRIC_LENGTH = readLength; }
            }
        }
        int Y_METRIC_LENGTH = 0;
        for( final Haplotype h : haplotypes ) {
            final int haplotypeLength = h.getBases().length;
            if( haplotypeLength > Y_METRIC_LENGTH ) { Y_METRIC_LENGTH = haplotypeLength; }
        }

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        pairHMMThreadLocal.get().initialize(haplotypes, perSampleReadList, X_METRIC_LENGTH, Y_METRIC_LENGTH);
    }

    private void finalizePairHMM()
    {
        pairHMMThreadLocal.get().finalizeRegion();
    }

    @Override
    public ReadLikelihoods<Haplotype> computeReadLikelihoods( final AssemblyResultSet assemblyResultSet, final SampleList samples, final Map<String, List<GATKSAMRecord>> perSampleReadList ) {

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // configure the HMM
        initializePairHMM(haplotypeList, perSampleReadList);

        // Add likelihoods for each sample's reads to our result
        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.sampleCount();
        for (int s = 0; s < sampleCount; s++) {
            final ReadLikelihoods.Matrix<Haplotype> sampleLikelihoods = result.sampleMatrix(s);
            computeReadLikelihoods(sampleLikelihoods);
        }

        result.normalizeLikelihoods(false, log10globalReadMismappingRate);
        result.filterPoorlyModeledReads(EXPECTED_ERROR_RATE_PER_BASE);
        finalizePairHMM();
        return result;
    }

    private void computeReadLikelihoods( final ReadLikelihoods.Matrix<Haplotype> likelihoods) {

        // Modify the read qualities by applying the PCR error model and capping the minimum base,insertion,deletion qualities
        final List<GATKSAMRecord> processedReads = modifyReadQualities(likelihoods.reads());

        final Map<GATKSAMRecord,byte[]> gapContinuationPenalties = buildGapContinuationPenalties(processedReads,constantGCP);
        // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
        pairHMMThreadLocal.get().computeLikelihoods(likelihoods,processedReads,gapContinuationPenalties);

        if (WRITE_LIKELIHOODS_TO_FILE)
            writeDebugLikelihoods(likelihoods);
    }

    private Map<GATKSAMRecord, byte[]> buildGapContinuationPenalties(final List<GATKSAMRecord> processedReads, final byte gcp) {
        final Map<GATKSAMRecord,byte[]> result = new HashMap<>(processedReads.size());
        for (final GATKSAMRecord read : processedReads) {
            final byte[] readGcpArray = new byte[read.getReadLength()];
            Arrays.fill(readGcpArray,gcp);
            result.put(read,readGcpArray);
        }
        return result;
    }

    private void writeDebugLikelihoods(final ReadLikelihoods.Matrix<Haplotype> likelihoods) {
        final List<GATKSAMRecord> reads = likelihoods.reads();
        final List<Haplotype> haplotypes = likelihoods.alleles();
        final int haplotypeCount = haplotypes.size();
        final int readCount = reads.size();
        for (int r = 0; r < readCount; r++)
            for (int a = 0; a < haplotypeCount; a++)
                writeDebugLikelihoods(reads.get(r),haplotypes.get(a),likelihoods.get(a,r));
        likelihoodsStream.flush();
    }

    private void writeDebugLikelihoods(final GATKSAMRecord processedRead, final Haplotype haplotype, final double log10l){
        likelihoodsStream.printf("%s %s %s %s %s %s %f%n",
                    haplotype.getBaseString(),
                    new String(processedRead.getReadBases() ),
                    SAMUtils.phredToFastq(processedRead.getBaseQualities()),
                    SAMUtils.phredToFastq(processedRead.getBaseInsertionQualities() ),
                    SAMUtils.phredToFastq(processedRead.getBaseDeletionQualities() ),
                    SAMUtils.phredToFastq(constantGCP),
                    log10l);
    }

    @Requires({"alleleOrdering.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == alleleOrdering.size()"})
    @Deprecated
    public static double[][] computeDiploidHaplotypeLikelihoods( final String sample,
                                                                 final ReadLikelihoods readLikelihoods,
                                                                 final List alleleOrdering,
                                                                 final boolean normalize ) {
        return computeDiploidHaplotypeLikelihoods(Collections.singleton(sample), readLikelihoods, alleleOrdering, normalize);
    }

    @Requires({"alleleOrdering.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == alleleOrdering.size()"})
    @Deprecated
    private static double[][] computeDiploidHaplotypeLikelihoods( final Set<String> samples,
                                                                 final ReadLikelihoods readLikelihoods,
                                                                 final List alleleOrdering,
                                                                 final boolean normalize) {

        final int numHaplotypes = alleleOrdering.size();
        final int[] alleleIndices = new int[alleleOrdering.size()];
        final ListIterator alleleIterator = alleleOrdering.listIterator();
        int nextAlleleIndex = 0;
        while (alleleIterator.hasNext())
            if ((alleleIndices[nextAlleleIndex++] = readLikelihoods.alleleIndex((Allele) alleleIterator.next())) == -1)
                throw new IllegalArgumentException("allele " + alleleIterator.previous() + " not found in likelihood collection ");

        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];

        // compute the diploid haplotype likelihoods
        for(final String sample : samples) {
            final int sampleIndex = readLikelihoods.sampleIndex(sample);
            if (sampleIndex == -1)
                throw new IllegalArgumentException("the sample provided is not in the likelihood collection");
            final ReadLikelihoods.Matrix sampleLikelihoods = readLikelihoods.sampleMatrix(sampleIndex);
            final int sampleReadCount = readLikelihoods.sampleReadCount(sampleIndex);
            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                final int iii_allele = alleleIndices[iii];
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    final int jjj_allele = alleleIndices[jjj];
                    double haplotypeLikelihood = 0.0;
                    for (int r = 0; r < sampleReadCount; r++) {
                        final double value = MathUtils.approximateLog10SumLog10(sampleLikelihoods.get(iii_allele,r),
                                sampleLikelihoods.get(jjj_allele,r)) + MathUtils.LOG_ONE_HALF;
                        haplotypeLikelihood += value;
                    }
                    haplotypeLikelihoodMatrix[iii][jjj] += haplotypeLikelihood;
                }
            }
        }

        // normalize the diploid likelihoods matrix
        return normalize ? normalizeDiploidLikelihoodMatrixFromLog10( haplotypeLikelihoodMatrix ) : haplotypeLikelihoodMatrix;
    }

    @Requires({"likelihoodMatrix.length == likelihoodMatrix[0].length"})
    @Ensures({"result.length == result[0].length", "result.length == likelihoodMatrix.length"})
    @Deprecated
    protected static double[][] normalizeDiploidLikelihoodMatrixFromLog10( final double[][] likelihoodMatrix ) {
        final int numHaplotypes = likelihoodMatrix.length;
        double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int index = 0;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ){
                genotypeLikelihoods[index++] = likelihoodMatrix[iii][jjj];
            }
        }
        genotypeLikelihoods = MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
        index = 0;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ){
                likelihoodMatrix[iii][jjj] = genotypeLikelihoods[index++];
            }
        }
        return likelihoodMatrix;
    }

    // --------------------------------------------------------------------------------
    //
    // Experimental attempts at PCR error rate modeling
    //
    // --------------------------------------------------------------------------------

    protected static final int MAX_STR_UNIT_LENGTH = 8;
    protected static final int MAX_REPEAT_LENGTH = 20;
    protected static final int MIN_ADJUSTED_QSCORE = 10;
    protected static final double INITIAL_QSCORE = 40.0;

    private byte[] pcrIndelErrorModelCache = null;

    private void initializePCRErrorModel() {
        if ( pcrErrorModel == PCR_ERROR_MODEL.NONE || !pcrErrorModel.hasRateFactor() )
            return;

        pcrIndelErrorModelCache = new byte[MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel.getRateFactor();

        for( int iii = 0; iii <= MAX_REPEAT_LENGTH; iii++ )
            pcrIndelErrorModelCache[iii] = getErrorModelAdjustedQual(iii, rateFactor);
    }

    protected static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor) {
        return (byte) Math.max(MIN_ADJUSTED_QSCORE, MathUtils.fastRound( INITIAL_QSCORE - Math.exp(((double) repeatLength) / (rateFactor * Math.PI)) + 1.0 ));
    }



    protected void applyPCRErrorModel(final byte[] readBases, final byte[] readInsQuals, final byte[] readDelQuals ) {
        if ( pcrErrorModel == PCR_ERROR_MODEL.NONE )
            return;

        final TandemRepeatFinder repeatFinder = new TandemRepeatFinder(readBases, MAX_STR_UNIT_LENGTH, MAX_REPEAT_LENGTH);
        for ( int iii = 1; iii < readBases.length; iii++ ) {
            final int repeatLength = repeatFinder.findMostRelevantTandemRepeatUnitAt(iii - 1).getRepeatCount();
            readInsQuals[iii-1] = (byte) Math.min(0xff & readInsQuals[iii-1], 0xff & pcrIndelErrorModelCache[repeatLength]);
            readDelQuals[iii-1] = (byte) Math.min(0xff & readDelQuals[iii-1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        }
    }
}
