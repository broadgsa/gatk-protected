/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.pairhmm.*;
import org.broadinstitute.sting.utils.recalibration.covariates.RepeatCovariate;
import org.broadinstitute.sting.utils.recalibration.covariates.RepeatLengthCovariate;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.variant.variantcontext.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

public class PairHMMLikelihoodCalculationEngine implements LikelihoodCalculationEngine {
    private final static Logger logger = Logger.getLogger(PairHMMLikelihoodCalculationEngine.class);

    public static final byte BASE_QUALITY_SCORE_THRESHOLD = (byte) 18; // Base quals less than this value are squashed down to min possible qual

    private final byte constantGCP;
    private final double log10globalReadMismappingRate;
    private final boolean DEBUG;

    private final PairHMM.HMM_IMPLEMENTATION hmmType;
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
		case JNI_LOGLESS_CACHING:
		    return new JNILoglessPairHMM();
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
        NONE,
        /** a more aggressive model will be applied that sacrifices true positives in order to remove more false positives */
        AGGRESSIVE,
        /** a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives */
        CONSERVATIVE
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
     * @param debug should we emit debugging information during the calculation?
     * @param hmmType the type of the HMM to use
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param noFpga disable FPGA acceleration
     */
    public PairHMMLikelihoodCalculationEngine( final byte constantGCP, final boolean debug, final PairHMM.HMM_IMPLEMENTATION hmmType, final double log10globalReadMismappingRate, final boolean noFpga, final PCR_ERROR_MODEL pcrErrorModel ) {
        this.hmmType = hmmType;
        this.constantGCP = constantGCP;
        this.DEBUG = debug;
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

    private void writeDebugLikelihoods(final GATKSAMRecord processedRead, final Haplotype haplotype, final double log10l){
        if ( WRITE_LIKELIHOODS_TO_FILE ) {
            likelihoodsStream.printf("%s %s %s %s %s %s %f%n",
                    haplotype.getBaseString(),
                    new String(processedRead.getReadBases() ),
                    SAMUtils.phredToFastq(processedRead.getBaseQualities() ),
                    SAMUtils.phredToFastq(processedRead.getBaseInsertionQualities() ),
                    SAMUtils.phredToFastq(processedRead.getBaseDeletionQualities() ),
                    SAMUtils.phredToFastq(constantGCP),
                    log10l);
        }
    }

    private Map<Allele, Haplotype> createAlleleMap(List<Haplotype> haplotypes){
        final int numHaplotypes = haplotypes.size();
        final Map<Allele, Haplotype> alleleMap = new LinkedHashMap<>(numHaplotypes);
        for ( final Haplotype haplotype : haplotypes ) {
            final Allele allele = Allele.create(haplotype, true);
            alleleMap.put(allele, haplotype);
        }
        return alleleMap;
    }

    private Map<GATKSAMRecord, byte[]> fillGCPArrays(List<GATKSAMRecord> reads){
        final Map<GATKSAMRecord, byte []> GCPArrayMap = new LinkedHashMap<>();
        for (GATKSAMRecord read: reads){
            byte [] GCPArray = new byte[read.getReadBases().length];
            Arrays.fill( GCPArray, constantGCP ); // Is there a way to derive empirical estimates for this from the data?
            GCPArrayMap.put(read, GCPArray);
        }
        return GCPArrayMap;
    }

    private void capMinimumReadQualities(GATKSAMRecord read, byte[] readQuals, byte[] readInsQuals, byte[] readDelQuals) {
        for( int kkk = 0; kkk < readQuals.length; kkk++ ) {
            readQuals[kkk] = (byte) Math.min( 0xff & readQuals[kkk], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
            readQuals[kkk] = ( readQuals[kkk] < BASE_QUALITY_SCORE_THRESHOLD ? QualityUtils.MIN_USABLE_Q_SCORE : readQuals[kkk] );
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
        List<GATKSAMRecord> processedReads = new LinkedList<>();
        for ( GATKSAMRecord read : reads ) {

            final byte[] readBases = read.getReadBases();

            // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
            final byte[] readQuals = read.getBaseQualities().clone();
            final byte[] readInsQuals = read.getBaseInsertionQualities().clone();
            final byte[] readDelQuals = read.getBaseDeletionQualities().clone();

            applyPCRErrorModel(readBases, readInsQuals, readDelQuals);
            capMinimumReadQualities(read, readQuals, readInsQuals, readDelQuals);

            // Create a new copy of the read and sets its base qualities to the modified versions.
            // Pack this into a new list for return
            final GATKSAMRecord processedRead = GATKSAMRecord.createQualityModifiedRead(read, readBases, readQuals, readInsQuals, readDelQuals);
            processedReads.add(processedRead);
        }
        return processedReads;
    }

    /**
     * Post-processing of the read/allele likelihoods.
     *
     * We send quality-capped reads to the pairHMM for evaluation, and it returns a map containing these capped reads.
     * We wish to return a map containing the original, unmodified reads.
     *
     * At the same time, we want to effectively set a lower cap on the reference score, based on the global mis-mapping rate.
     * This protects us from the case where the assembly has produced haplotypes
     * that are very divergent from reference, but are supported by only one read.  In effect
     * we capping how badly scoring the reference can be for any read by the chance that the read
     * itself just doesn't belong here
     *
     * @param perReadAlleleLikelihoodMap the original map returned by the PairHMM. Contains the processed reads, the haplotype Alleles, and their log10ls
     * @param reads Our original, unmodified reads
     * @param processedReads Reads whose minimum base,insertion,deletion qualities have been capped; these were actually used to derive log10ls
     * @param alleleHaplotypeMap The map associating the Allele and Haplotype versions of each haplotype
     *
     * @return processedReadAlleleLikelihoodMap; a new PRALM containing the original reads, and their haplotype log10ls including capped reference log10ls
     */
    private PerReadAlleleLikelihoodMap capReferenceHaplotypeLikelihoods(PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap, List<GATKSAMRecord> reads, List<GATKSAMRecord> processedReads, Map<Allele, Haplotype> alleleHaplotypeMap){

        // a new read/allele map, to contain the uncapped reads, haplotypes, and potentially the capped reference log10ls
        final PerReadAlleleLikelihoodMap processedReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();

        Allele refAllele = null;
        final int numReads = reads.size();
        for (int readIndex = 0; readIndex < numReads; readIndex++) {

            // Get the original and quality-modified read from their respective lists
            // Note that this requires both lists to have reads in the same order
            final GATKSAMRecord originalRead = reads.get(readIndex);
            final GATKSAMRecord processedRead = processedReads.get(readIndex);

            // keep track of the reference likelihood and the best non-ref likelihood
            double refLog10l = Double.NEGATIVE_INFINITY;
            double bestNonReflog10L = Double.NEGATIVE_INFINITY;

            for ( Allele allele : alleleHaplotypeMap.keySet() ) {
                final double log10l = perReadAlleleLikelihoodMap.getLikelihoodAssociatedWithReadAndAllele(processedRead, allele);
                final Haplotype haplotype = alleleHaplotypeMap.get(allele);
                if ( haplotype.isNonReference() )
                    bestNonReflog10L = Math.max(bestNonReflog10L, log10l);
                else {
                    refAllele = allele;
                    refLog10l = log10l;
                }
                writeDebugLikelihoods(processedRead, haplotype, log10l);

                // add the ORIGINAL (non-capped) read to the final map, along with the current haplotype and associated log10l
                processedReadAlleleLikelihoodMap.add(originalRead, allele, log10l);
            }

            // ensure that the reference haplotype is no worse than the best non-ref haplotype minus the global
            // mismapping rate.  This protects us from the case where the assembly has produced haplotypes
            // that are very divergent from reference, but are supported by only one read.  In effect
            // we capping how badly scoring the reference can be for any read by the chance that the read
            // itself just doesn't belong here
            final double worstRefLog10Allowed = bestNonReflog10L + log10globalReadMismappingRate;
            if ( refLog10l < (worstRefLog10Allowed) ) {
                processedReadAlleleLikelihoodMap.add(originalRead, refAllele, worstRefLog10Allowed);
            }
        }
        return processedReadAlleleLikelihoodMap;
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
    public Map<String, PerReadAlleleLikelihoodMap> computeReadLikelihoods( final AssemblyResultSet assemblyResultSet, final Map<String, List<GATKSAMRecord>> perSampleReadList ) {

        final List<Haplotype> haplotypes = assemblyResultSet.getHaplotypeList();
        // configure the HMM
        initializePairHMM(haplotypes, perSampleReadList);

        // Add likelihoods for each sample's reads to our stratifiedReadMap
        final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = new LinkedHashMap<>();
        for( final Map.Entry<String, List<GATKSAMRecord>> sampleEntry : perSampleReadList.entrySet() ) {
             // evaluate the likelihood of the reads given those haplotypes
            final PerReadAlleleLikelihoodMap map = computeReadLikelihoods(haplotypes, sampleEntry.getValue());

            map.filterPoorlyModelledReads(EXPECTED_ERROR_RATE_PER_BASE);
            stratifiedReadMap.put(sampleEntry.getKey(), map);
        }
	
	//Used mostly by the JNI implementation(s) to free arrays
	finalizePairHMM();

        return stratifiedReadMap;
    }


    public Map<String, PerReadAlleleLikelihoodMap> computeReadLikelihoods( final List<Haplotype> haplotypes, final Map<String, List<GATKSAMRecord>> perSampleReadList ) {

        // Add likelihoods for each sample's reads to our stratifiedReadMap
        final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = new LinkedHashMap<>();
        for( final Map.Entry<String, List<GATKSAMRecord>> sampleEntry : perSampleReadList.entrySet() ) {
            // evaluate the likelihood of the reads given those haplotypes
            final PerReadAlleleLikelihoodMap map = computeReadLikelihoods(haplotypes, sampleEntry.getValue());

            map.filterPoorlyModelledReads(EXPECTED_ERROR_RATE_PER_BASE);
            stratifiedReadMap.put(sampleEntry.getKey(), map);
        }

        return stratifiedReadMap;
    }

    private PerReadAlleleLikelihoodMap computeReadLikelihoods( final List<Haplotype> haplotypes, final List<GATKSAMRecord> reads) {

        // Modify the read qualities by applying the PCR error model and capping the minimum base,insertion,deletion qualities
        List<GATKSAMRecord> processedReads = modifyReadQualities(reads);

        // Get alleles corresponding to our haplotypees
        Map<Allele, Haplotype> alleleHaplotypeMap = createAlleleMap(haplotypes);

        // Get an array containing the constantGCP for each read in our modified read list
        Map<GATKSAMRecord,byte[]> GCPArrayMap = fillGCPArrays(processedReads);

        // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = pairHMMThreadLocal.get().computeLikelihoods(processedReads, alleleHaplotypeMap, GCPArrayMap);

        // Generate a new map containing the original, unmodified reads, and with minimal reference haplotype log10ls determined from the global mis-mapping rate

        return capReferenceHaplotypeLikelihoods(perReadAlleleLikelihoodMap, reads, processedReads, alleleHaplotypeMap);
    }

    @Requires({"alleleOrdering.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == alleleOrdering.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final String sample,
                                                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap,
                                                                 final List<Allele> alleleOrdering,
                                                                 final boolean normalize ) {
        return computeDiploidHaplotypeLikelihoods(Collections.singleton(sample), stratifiedReadMap, alleleOrdering, normalize);
    }

    @Requires({"alleleOrdering.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == alleleOrdering.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final Set<String> samples,
                                                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap,
                                                                 final List<Allele> alleleOrdering,
                                                                 final boolean normalize) {

        final int numHaplotypes = alleleOrdering.size();
        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            Arrays.fill(haplotypeLikelihoodMatrix[iii], Double.NEGATIVE_INFINITY);
        }

        // compute the diploid haplotype likelihoods
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            final Allele iii_allele = alleleOrdering.get(iii);
            for( int jjj = 0; jjj <= iii; jjj++ ) {
                final Allele jjj_allele = alleleOrdering.get(jjj);
                double haplotypeLikelihood = 0.0;
                for( final String sample : samples ) {
                    for( final Map.Entry<GATKSAMRecord, Map<Allele,Double>> entry : stratifiedReadMap.get(sample).getLikelihoodReadMap().entrySet() ) {
                        // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                        // First term is approximated by Jacobian log with table lookup.
                        haplotypeLikelihood += ReadUtils.getMeanRepresentativeReadCount( entry.getKey() ) *
                                ( MathUtils.approximateLog10SumLog10(entry.getValue().get(iii_allele), entry.getValue().get(jjj_allele)) + MathUtils.LOG_ONE_HALF );
                    }
                }
                haplotypeLikelihoodMatrix[iii][jjj] = haplotypeLikelihood;
            }
        }

        // normalize the diploid likelihoods matrix
        return normalize ? normalizeDiploidLikelihoodMatrixFromLog10( haplotypeLikelihoodMatrix ) : haplotypeLikelihoodMatrix;
    }

    @Requires({"likelihoodMatrix.length == likelihoodMatrix[0].length"})
    @Ensures({"result.length == result[0].length", "result.length == likelihoodMatrix.length"})
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
    // System to compute the best N haplotypes for genotyping
    //
    // --------------------------------------------------------------------------------
//
//    /**
//     * Helper function for selectBestHaplotypesFromEachSample that updates the score of haplotype haplotypeAsAllele
//     * @param map an annoying map object that moves us between the allele and haplotype representation
//     * @param haplotypeAsAllele the allele version of the haplotype
//     * @return the haplotype version, with its score incremented by 1 if its non-reference
//     */
//    private Haplotype updateSelectHaplotype(final Map<Allele, Haplotype> map, final Allele haplotypeAsAllele) {
//        final Haplotype h = map.get(haplotypeAsAllele); // TODO -- fixme when haplotypes are properly generic
//        if ( h.isNonReference() ) h.setScore(h.getScore() + 1); // ref is already at max value
//        return h;
//    }
//
//    /**
//     * Take the best N haplotypes and return them as a list
//     *
//     * Only considers the haplotypes selectedHaplotypes that were actually selected by at least one sample
//     * as it's preferred haplotype.  Takes the best N haplotypes from selectedHaplotypes in decreasing
//     * order of score (so higher score haplotypes are preferred).  The N we take is determined by
//     *
//     * N = min(2 * nSamples + 1, maxNumHaplotypesInPopulation)
//     *
//     * where 2 * nSamples is the number of chromosomes in 2 samples including the reference, and our workload is
//     * bounded by maxNumHaplotypesInPopulation as that number can grow without bound
//     *
//     * @param selectedHaplotypes a non-null set of haplotypes with scores >= 1
//     * @param nSamples the number of samples used to select the haplotypes
//     * @param maxNumHaplotypesInPopulation the maximum number of haplotypes we're allowed to take, regardless of nSamples
//     * @return a list of N or fewer haplotypes, with the reference haplotype first
//     */
//    private List<Haplotype> selectBestHaplotypesAccordingToScore(final Set<Haplotype> selectedHaplotypes, final int nSamples, final int maxNumHaplotypesInPopulation) {
//        final List<Haplotype> selectedHaplotypesList = new ArrayList<>(selectedHaplotypes);
//        Collections.sort(selectedHaplotypesList, new HaplotypeScoreComparator());
//        final int numChromosomesInSamplesPlusRef = 2 * nSamples + 1;
//        final int haplotypesToKeep = Math.min(numChromosomesInSamplesPlusRef, maxNumHaplotypesInPopulation);
//        final List<Haplotype> bestHaplotypes = selectedHaplotypesList.size() <= haplotypesToKeep ? selectedHaplotypesList : selectedHaplotypesList.subList(0, haplotypesToKeep);
//        if ( bestHaplotypes.get(0).isNonReference()) throw new IllegalStateException("BUG: reference haplotype should be first in list");
//        return bestHaplotypes;
//    }
//
//    /**
//     * Select the best haplotypes for genotyping the samples in stratifiedReadMap
//     *
//     * Selects these haplotypes by counting up how often each haplotype is selected as one of the most likely
//     * haplotypes per sample.  What this means is that each sample computes the diploid genotype likelihoods for
//     * all possible pairs of haplotypes, and the pair with the highest likelihood has each haplotype each get
//     * one extra count for each haplotype (so hom-var haplotypes get two counts).  After performing this calculation
//     * the best N haplotypes are selected (@see #selectBestHaplotypesAccordingToScore) and a list of the
//     * haplotypes in order of score are returned, ensuring that at least one of the haplotypes is reference.
//     *
//     * @param haplotypes a list of all haplotypes we're considering
//     * @param stratifiedReadMap a map from sample -> read likelihoods per haplotype
//     * @param maxNumHaplotypesInPopulation the max. number of haplotypes we can select from haplotypes
//     * @return a list of selected haplotypes with size <= maxNumHaplotypesInPopulation
//     */
//    public List<Haplotype> selectBestHaplotypesFromEachSample(final List<Haplotype> haplotypes, final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap, final int maxNumHaplotypesInPopulation) {
//        if ( haplotypes.size() < 2 ) throw new IllegalArgumentException("Must have at least 2 haplotypes to consider but only have " + haplotypes);
//
//        if ( haplotypes.size() == 2 ) return haplotypes; // fast path -- we'll always want to use 2 haplotypes
//
//        // all of the haplotypes that at least one sample called as one of the most likely
//        final Set<Haplotype> selectedHaplotypes = new HashSet<>();
//        selectedHaplotypes.add(findReferenceHaplotype(haplotypes)); // ref is always one of the selected
//
//        // our annoying map from allele -> haplotype
//        final Map<Allele, Haplotype> allele2Haplotype = new HashMap<>();
//        for ( final Haplotype h : haplotypes ) {
//            h.setScore(h.isReference() ? Double.MAX_VALUE : 0.0); // set all of the scores to 0 (lowest value) for all non-ref haplotypes
//            allele2Haplotype.put(Allele.create(h, h.isReference()), h);
//        }
//
//        // for each sample, compute the most likely pair of haplotypes
//        for ( final Map.Entry<String, PerReadAlleleLikelihoodMap> entry : stratifiedReadMap.entrySet() ) {
//            // get the two most likely haplotypes under a diploid model for this sample
//            final MostLikelyAllele mla = entry.getValue().getMostLikelyDiploidAlleles();
//
//            if ( mla != null ) { // there was something to evaluate in this sample
//                // note that there must be at least 2 haplotypes
//                final Haplotype best = updateSelectHaplotype(allele2Haplotype, mla.getMostLikelyAllele());
//                final Haplotype second = updateSelectHaplotype(allele2Haplotype, mla.getSecondMostLikelyAllele());
//
////            if ( DEBUG ) {
////                logger.info("Chose haplotypes " + best + " " + best.getCigar() + " and " + second + " " + second.getCigar() + " for sample " + entry.getKey());
////            }
//
//                // add these two haplotypes to the set of haplotypes that have been selected
//                selectedHaplotypes.add(best);
//                selectedHaplotypes.add(second);
//
//                // we've already selected all of our haplotypes, and we don't need to prune them down
//                if ( selectedHaplotypes.size() == haplotypes.size() && haplotypes.size() < maxNumHaplotypesInPopulation )
//                    break;
//            }
//        }
//
//        // take the best N haplotypes forward, in order of the number of samples that choose them
//        final int nSamples = stratifiedReadMap.size();
//        final List<Haplotype> bestHaplotypes = selectBestHaplotypesAccordingToScore(selectedHaplotypes, nSamples, maxNumHaplotypesInPopulation);
//
//        if ( DEBUG ) {
//            logger.info("Chose " + (bestHaplotypes.size() - 1) + " alternate haplotypes to genotype in all samples.");
//            for ( final Haplotype h : bestHaplotypes ) {
//                logger.info("\tHaplotype " + h.getCigar() + " selected for further genotyping" + (h.isNonReference() ? " found " + (int)h.getScore() + " haplotypes" : " as ref haplotype"));
//            }
//        }
//        return bestHaplotypes;
//    }
//
//    /**
//     * Find the haplotype that isRef(), or @throw ReviewedStingException if one isn't found
//     * @param haplotypes non-null list of haplotypes
//     * @return the reference haplotype
//     */
//    private static Haplotype findReferenceHaplotype( final List<Haplotype> haplotypes ) {
//        for( final Haplotype h : haplotypes ) {
//            if( h.isReference() ) return h;
//        }
//        throw new ReviewedStingException( "No reference haplotype found in the list of haplotypes!" );
//    }

    // --------------------------------------------------------------------------------
    //
    // Experimental attempts at PCR error rate modeling
    //
    // --------------------------------------------------------------------------------

    protected static final int MAX_STR_UNIT_LENGTH = 8;
    protected static final int MAX_REPEAT_LENGTH = 20;
    protected static final int MIN_ADJUSTED_QSCORE = 10;
    protected static final double INITIAL_QSCORE = 40.0;

    private byte[] pcrIndelErrorModelCache = new byte[MAX_REPEAT_LENGTH * MAX_STR_UNIT_LENGTH + 1];
    private final RepeatCovariate repeatCovariate = new RepeatLengthCovariate();

    private void initializePCRErrorModel() {
        if ( pcrErrorModel == PCR_ERROR_MODEL.NONE )
            return;

        repeatCovariate.initialize(MAX_STR_UNIT_LENGTH, MAX_REPEAT_LENGTH);

        pcrIndelErrorModelCache = new byte[MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel == PCR_ERROR_MODEL.AGGRESSIVE ? 2.0 : 3.0;

        for( int iii = 0; iii <= MAX_REPEAT_LENGTH; iii++ )
            pcrIndelErrorModelCache[iii] = getErrorModelAdjustedQual(iii, rateFactor);
    }

    protected static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor) {
        return (byte) Math.max(MIN_ADJUSTED_QSCORE, MathUtils.fastRound( INITIAL_QSCORE - Math.exp(((double) repeatLength) / (rateFactor * Math.PI)) + 1.0 ));
    }

    protected void applyPCRErrorModel( final byte[] readBases, final byte[] readInsQuals, final byte[] readDelQuals ) {
        if ( pcrErrorModel == PCR_ERROR_MODEL.NONE )
            return;

        for ( int iii = 1; iii < readBases.length; iii++ ) {
            final int repeatLength = repeatCovariate.findTandemRepeatUnits(readBases, iii-1).getSecond();
            readInsQuals[iii-1] = (byte) Math.min(0xff & readInsQuals[iii-1], 0xff & pcrIndelErrorModelCache[repeatLength]);
            readDelQuals[iii-1] = (byte) Math.min(0xff & readDelQuals[iii-1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        }
    }
}
