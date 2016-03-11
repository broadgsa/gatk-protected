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

import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Utility class that error-corrects reads.
 * Main idea: An error in a read will appear as a bubble in a k-mer (de Bruijn) graph and such bubble will have very low multiplicity.
 * Hence, read errors will appear as "sparse" kmers with very little support.
 * Historically, the most common approach to error-correct reads before assembly has been to first compute the kmer spectrum of the reads,
 * defined as the kmer composition of a set of reads along with the multiplicity of each kmer.
 * First-generation correctors like the Euler corrector (Pevzner 2001) mapped low frequency kmers (kmers appearing say below N times)
 * into high frequency ones that lied within a certain Hamming or edit distance.
 * This is doable, but has some drawbacks:
 * - Kmers used for error correction become tied to kmers used for graph building.
 * - Hence, large kmers (desireable for graph building because they can resolve repeats better) are a hindrance for error correction,
 * because they are seen less often.
 * - After error correction, there is no guarantee that a sequence of kmers corresponds to an "actual" read.
 *
 * An error-corrected set of reads also makes a much smoother graph without the need to resolving so many bubbles.
 *
 * Idea hence is to correct reads based on their kmer content, but in a context independent from graph building.
 * In order to do this, the following steps are taken:
 * - The k-mer spectrum of a set of reads in computed. However, we are at freedom to choose the most convenient k-mer size (typicially around
 * read length /2).
 * - We partition the set of observed k-mers into "solid" kmers which have multiplicity > M, and "insolid" ones otherwise (Pevzner 2001).
 *
 * - Main idea of the algorithm is to try to substitute a sequence of bases in a read by a sequence better supported by kmers.
 * - For each "unsolid" kmer observed in reads, we try to find a "solid" kmer within a maximum Hamming distance.
 * - If such solid kmer exists, then this unsolid kmer is "correctable", otherwise, uncorrectable.
 * - For each read, then:
 * -- Walk through  read and visit all kmers.
 * -- If kmer is solid, continue to next kmer.
 * -- If not, and if it's correctable (i.e. there exists a mapping from an unsolid kmer to a solid kmer within a given Hamming distance),
 *    add the bases and offsets corresponding to differing positions between unsolid and solid kmer to correction list.
 * -- At the end, each base in read will have a list of corrections associated with it. We can then choose to correct or not.
 *    If read has only consistent corrections, then we can correct base to common base in corrections.
 *
 *    TODO:
 *    todo Q: WHAT QUALITY TO USE??
 *    todo how do we deal with mate pairs?
 *
 *


 */
public class ReadErrorCorrector {
    private final static Logger logger = Logger.getLogger(ReadErrorCorrector.class);
    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    KMerCounter countsByKMer;

    Map<Kmer,Kmer> kmerCorrectionMap = new HashMap<>();
    Map<Kmer,Pair<int[],byte[]>> kmerDifferingBases = new HashMap<>();
    private final int kmerLength;
    private final boolean debug;
    private final boolean trimLowQualityBases;
    private final byte minTailQuality;
    private final int maxMismatchesToCorrect;
    private final byte qualityOfCorrectedBases;
    private final int maxObservationsForKmerToBeCorrectable;
    private final int maxHomopolymerLengthInRegion;
    private final int minObservationsForKmerToBeSolid;

    // default values, for debugging
    private final static boolean doInplaceErrorCorrection = false;    // currently not used, since we want corrected reads to be used only for assembly
    private final static int MAX_MISMATCHES_TO_CORRECT = 2;
    private final static byte QUALITY_OF_CORRECTED_BASES = 30; // what's a reasonable value here?
    private final static int MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE = 1;
    private final static boolean TRIM_LOW_QUAL_TAILS = false;
    private final static boolean DONT_CORRECT_IN_LONG_HOMOPOLYMERS = false;
    private final static int MAX_HOMOPOLYMER_THRESHOLD = 12;

    // debug counter structure
    private final ReadErrorCorrectionStats readErrorCorrectionStats = new ReadErrorCorrectionStats();

    /**
     * Create a new kmer corrector
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     * @param maxMismatchesToCorrect e >= 0
     * @param qualityOfCorrectedBases  Bases to be corrected will be assigned this quality
     */
    public ReadErrorCorrector(final int kmerLength,
                              final int maxMismatchesToCorrect,
                              final int maxObservationsForKmerToBeCorrectable,
                              final byte qualityOfCorrectedBases,
                              final int minObservationsForKmerToBeSolid,
                              final boolean trimLowQualityBases,
                              final byte minTailQuality,
                              final boolean debug,
                              final byte[] fullReferenceWithPadding) {
        if ( kmerLength < 1 ) throw new IllegalArgumentException("kmerLength must be > 0 but got " + kmerLength);
        if ( maxMismatchesToCorrect < 1 )
            throw new IllegalArgumentException("maxMismatchesToCorrect must be >= 1 but got " + maxMismatchesToCorrect);
        if ( qualityOfCorrectedBases < 2 || qualityOfCorrectedBases > QualityUtils.MAX_REASONABLE_Q_SCORE)
            throw new IllegalArgumentException("qualityOfCorrectedBases must be >= 2 and <= MAX_REASONABLE_Q_SCORE but got " + qualityOfCorrectedBases);

        countsByKMer = new KMerCounter(kmerLength);
        this.kmerLength = kmerLength;
        this.maxMismatchesToCorrect = maxMismatchesToCorrect;
        this.qualityOfCorrectedBases = qualityOfCorrectedBases;
        this.minObservationsForKmerToBeSolid = minObservationsForKmerToBeSolid;
        this.trimLowQualityBases = trimLowQualityBases;
        this.minTailQuality = minTailQuality;
        this.debug = debug;
        this.maxObservationsForKmerToBeCorrectable = maxObservationsForKmerToBeCorrectable;

        // when region has long homopolymers, we may want not to correct reads, since assessment is complicated,
        // so we may decide to skip error correction in these regions
        maxHomopolymerLengthInRegion = computeMaxHLen(fullReferenceWithPadding);
    }

    /**
     * Simple constructor with sensible defaults
     * @param kmerLength            K-mer length for error correction (not necessarily the same as for assembly graph)
     * @param minTailQuality        Minimum tail quality: remaining bases with Q's below this value are hard-clipped after correction
     * @param debug                 Output debug information
     */
    public ReadErrorCorrector(final int kmerLength, final byte minTailQuality, final int minObservationsForKmerToBeSolid, final boolean debug,final byte[] fullReferenceWithPadding) {
        this(kmerLength, MAX_MISMATCHES_TO_CORRECT, MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE, QUALITY_OF_CORRECTED_BASES, minObservationsForKmerToBeSolid, TRIM_LOW_QUAL_TAILS, minTailQuality, debug,fullReferenceWithPadding);
    }

    /**
    * Main entry routine to add all kmers in a read to the read map counter
    * @param read                        Read to add bases
    */
    @Requires("read != null")
    protected void addReadKmers(final GATKSAMRecord read) {
        if (DONT_CORRECT_IN_LONG_HOMOPOLYMERS && maxHomopolymerLengthInRegion > MAX_HOMOPOLYMER_THRESHOLD)
            return;

        final byte[] readBases = read.getReadBases();
        for (int offset = 0; offset <= readBases.length-kmerLength; offset++ )  {
            countsByKMer.addKmer(new Kmer(readBases,offset,kmerLength),1);

        }
    }

    /**
     * Correct a collection of reads based on stored k-mer counts
     * @param reads
     */
    public final List<GATKSAMRecord> correctReads(final Collection<GATKSAMRecord> reads) {

        final List<GATKSAMRecord> correctedReads = new ArrayList<>(reads.size());
        if (DONT_CORRECT_IN_LONG_HOMOPOLYMERS && maxHomopolymerLengthInRegion > MAX_HOMOPOLYMER_THRESHOLD) {
            // just copy reads into output and exit
            correctedReads.addAll(reads);
        }
        else {
            computeKmerCorrectionMap();
            for (final GATKSAMRecord read: reads) {
                final GATKSAMRecord correctedRead = correctRead(read);
                if (trimLowQualityBases)
                    correctedReads.add(ReadClipper.hardClipLowQualEnds(correctedRead, minTailQuality));
                else
                    correctedReads.add(correctedRead);
            }
            if (debug) {
                logger.info("Number of corrected bases:"+readErrorCorrectionStats.numBasesCorrected);
                logger.info("Number of corrected reads:"+readErrorCorrectionStats.numReadsCorrected);
                logger.info("Number of skipped reads:"+readErrorCorrectionStats.numReadsUncorrected);
                logger.info("Number of solid kmers:"+readErrorCorrectionStats.numSolidKmers);
                logger.info("Number of corrected kmers:"+readErrorCorrectionStats.numCorrectedKmers);
                logger.info("Number of uncorrectable kmers:"+readErrorCorrectionStats.numUncorrectableKmers);
            }
        }
        return correctedReads;
    }


    /**
     * Do actual read correction based on k-mer map. First, loop through stored k-mers to get a list of possible corrections
     * for each position in the read. Then correct read based on all possible consistent corrections.
     * @param inputRead                               Read to correct
     * @return                                        Corrected read (can be same reference as input if doInplaceErrorCorrection is set)
     */
    @Requires("inputRead != null")
    private GATKSAMRecord correctRead(final GATKSAMRecord inputRead) {
        // do actual correction
        boolean corrected = false;
        final byte[] correctedBases = inputRead.getReadBases();
        final byte[] correctedQuals = inputRead.getBaseQualities();

        // array to store list of possible corrections for read
        final CorrectionSet correctionSet = buildCorrectionMap(correctedBases);

        for (int offset = 0; offset < correctedBases.length; offset++) {
            final Byte b = correctionSet.getConsensusCorrection(offset);
            if (b != null && b != correctedBases[offset]) {
                correctedBases[offset] = b;
                correctedQuals[offset] = qualityOfCorrectedBases;
                corrected = true;
            }
            readErrorCorrectionStats.numBasesCorrected++;
        }

        if (corrected) {
            readErrorCorrectionStats.numReadsCorrected++;
            if (doInplaceErrorCorrection) {
                inputRead.setReadBases(correctedBases);
                inputRead.setBaseQualities(correctedQuals);
                return inputRead;
            }
            else {
                GATKSAMRecord correctedRead = new GATKSAMRecord(inputRead);

                //  do the actual correction
                // todo - do we need to clone anything else from read?
                correctedRead.setBaseQualities(inputRead.getBaseQualities());
                correctedRead.setIsStrandless(inputRead.isStrandless());
                correctedRead.setReadBases(inputRead.getReadBases());
                correctedRead.setReadString(inputRead.getReadString());
                correctedRead.setReadGroup(inputRead.getReadGroup());
                return correctedRead;
            }
        }
        else {
            readErrorCorrectionStats.numReadsUncorrected++;
            return inputRead;
        }
    }

    /**
     * Build correction map for each of the bases in read.
     * For each of the constituent kmers in read:
     * a) See whether the kmer has been mapped to a corrected kmer.
     * b) If so, get list of differing positions and corresponding bases.
     * c) Add then list of new bases to index in correction list.
     * Correction list is of read size, and holds a list of bases to correct.
     * @param correctedBases                        Bases to attempt to correct
     * @return                                      CorrectionSet object.
     */
    @Requires("correctedBases != null")
    private CorrectionSet buildCorrectionMap(final byte[] correctedBases) {
        // array to store list of possible corrections for read
        final CorrectionSet correctionSet = new CorrectionSet(correctedBases.length);

        for (int offset = 0; offset <= correctedBases.length-kmerLength; offset++ )  {
            final Kmer kmer = new Kmer(correctedBases,offset,kmerLength);
            final Kmer newKmer = kmerCorrectionMap.get(kmer);
            if (newKmer != null && !newKmer.equals(kmer)){
                final Pair<int[],byte[]> differingPositions = kmerDifferingBases.get(kmer);
                final int[] differingIndeces = differingPositions.first;
                final byte[] differingBases = differingPositions.second;

                for (int k=0; k < differingIndeces.length; k++) {
                    // get list of differing positions for corrected kmer
                    // for each of these, add correction candidate to correction set
                    correctionSet.add(offset + differingIndeces[k],differingBases[k]);
                }
            }
        }
        return correctionSet;
    }


    /**
     * Top-level entry point that adds a collection of reads to our kmer list.
     * For each read in list, its constituent kmers will be logged in our kmer table.
     * @param reads
     */
    @Requires("reads != null")
    public void addReadsToKmers(final Collection<GATKSAMRecord> reads) {
        for (final GATKSAMRecord read: reads)
            addReadKmers(read);

        if (debug)
            for ( final KMerCounter.CountedKmer countedKmer: countsByKMer.getCountedKmers() )
                logger.info(String.format("%s\t%d\n", countedKmer.kmer, countedKmer.count));
    }


    /**
     * For each kmer we've seen, do the following:
     * a) If kmer count > threshold1, this kmer is good, so correction map will be to itself.
     * b) If kmer count <= threshold2, this kmer is bad.
     *    In that case, loop through all other kmers. If kmer is good, compute distance, and get minimal distance.
     *    If such distance is < some threshold, map to this kmer, and record differing positions and bases.
     *
     */
    private void computeKmerCorrectionMap() {
        for (final KMerCounter.CountedKmer storedKmer : countsByKMer.getCountedKmers()) {
            if (storedKmer.getCount() >= minObservationsForKmerToBeSolid) {
                // this kmer is good: map to itself
                kmerCorrectionMap.put(storedKmer.getKmer(),storedKmer.getKmer());
                kmerDifferingBases.put(storedKmer.getKmer(),new Pair<>(new int[0],new byte[0])); // dummy empty array
                readErrorCorrectionStats.numSolidKmers++;
            }
            else if (storedKmer.getCount() <= maxObservationsForKmerToBeCorrectable) {
                // loop now thru all other kmers to find nearest neighbor
                final Pair<Kmer,Pair<int[],byte[]>> nearestNeighbor = findNearestNeighbor(storedKmer.getKmer(),countsByKMer,maxMismatchesToCorrect);

                // check if nearest neighbor lies in a close vicinity. If so, log the new bases and the correction map
                if (nearestNeighbor != null) { // ok, found close neighbor
                    kmerCorrectionMap.put(storedKmer.getKmer(), nearestNeighbor.first);
                    kmerDifferingBases.put(storedKmer.getKmer(), nearestNeighbor.second);
                    readErrorCorrectionStats.numCorrectedKmers++;
//                    if (debug)
//                        logger.info("Original kmer:"+storedKmer + "\tCorrected kmer:"+nearestNeighbor.first+"\tDistance:"+dist);
                }
                else
                    readErrorCorrectionStats.numUncorrectableKmers++;

            }
         }
    }

    /**
     * Finds nearest neighbor of a given k-mer, among a list of counted K-mers, up to a given distance.
     * If many k-mers share same closest distance, an arbitrary k-mer is picked
     * @param kmer                        K-mer of interest
     * @param countsByKMer                KMerCounter storing set of counted k-mers (may include kmer of interest)
     * @param maxDistance                 Maximum distance to search
     * @return                            Pair of values: closest K-mer in Hamming distance and list of differing bases.
     *                                      If no neighbor can be found up to given distance, returns null
     */
    @Requires({"kmer != null", "countsByKMer != null","maxDistance >= 1"})
    private Pair<Kmer,Pair<int[],byte[]>> findNearestNeighbor(final Kmer kmer,
                                                             final KMerCounter  countsByKMer,
                                                             final int maxDistance) {
        int minimumDistance = Integer.MAX_VALUE;
        Kmer closestKmer = null;

        final int[] differingIndeces = new int[maxDistance+1];
        final byte[] differingBases = new byte[maxDistance+1];

        final int[] closestDifferingIndices = new int[maxDistance+1];
        final byte[] closestDifferingBases = new byte[maxDistance+1];

        for (final KMerCounter.CountedKmer candidateKmer : countsByKMer.getCountedKmers()) {
            // skip if candidate set includes test kmer
            if (candidateKmer.getKmer().equals(kmer))
                continue;

            final int hammingDistance  = kmer.getDifferingPositions(candidateKmer.getKmer(), maxDistance, differingIndeces, differingBases);
            if (hammingDistance < 0) // can't compare kmer? skip
                continue;

            if (hammingDistance < minimumDistance)  {
                minimumDistance = hammingDistance;
                closestKmer = candidateKmer.getKmer();
                System.arraycopy(differingBases,0,closestDifferingBases,0,differingBases.length);
                System.arraycopy(differingIndeces,0,closestDifferingIndices,0,differingIndeces.length);
            }
        }
        return new Pair<>(closestKmer, new Pair<>(closestDifferingIndices,closestDifferingBases));
    }


    /**
     * experimental function to compute max homopolymer length in a given reference context
     * @param fullReferenceWithPadding                Reference context of interest
     * @return                                        Max homopolymer length in region
     */
    @Requires("fullReferenceWithPadding != null")
    private static int computeMaxHLen(final byte[] fullReferenceWithPadding) {

        int leftRun = 1;
        int maxRun = 1;
        for ( int i = 1; i < fullReferenceWithPadding.length; i++) {
            if ( fullReferenceWithPadding[i] == fullReferenceWithPadding[i-1] )
                leftRun++;
            else
                leftRun = 1;
            }
            if (leftRun > maxRun)
                maxRun = leftRun;


        return maxRun;
    }

    private static final class ReadErrorCorrectionStats {
        public int numReadsCorrected;
        public int numReadsUncorrected;
        public int numBasesCorrected;
        public int numSolidKmers;
        public int numUncorrectableKmers;
        public int numCorrectedKmers;
    }

    /**
     * Wrapper utility class that holds, for each position in read, a list of bytes representing candidate corrections.
     * So, a read ACAGT where the middle A has found to be errorful might look like:
     * 0: {}
     * 1: {}
     * 2: {'C','C','C'}
     * 3: {}
     * 4: {}
     *
     * It's up to the method getConsensusCorrection()  to decide how to use the correction sets for each position.
     * By default, only strict consensus is allowed right now.
     *
     */
    protected static class CorrectionSet {
        private final int size;
        private ArrayList<List<Byte>> corrections;

        /**
         * Main class constructor.
         * @param size      Size of correction set, needs to be set equal to the read being corrected
         */
        public CorrectionSet(final int size) {
            this.size = size;
            corrections = new ArrayList<>(size);
            for (int k=0; k < size; k++)
                corrections.add(k,new ArrayList<Byte>());
        }

        /**
         * Add a base to this correction set at a particular offset, measured from the start of the read
         * @param offset                          Offset from start of read
         * @param base                            base to be added to list of corrections at this offset
         */
        public void add(final int offset, final byte base) {
            if (offset >= size || offset < 0)
                throw new IllegalStateException("Bad entry into CorrectionSet: offset > size");
            if (!BaseUtils.isRegularBase(base))
                return; // no irregular base correction

            final List<Byte> storedBytes = corrections.get(offset);
            storedBytes.add(base);
        }

        /**
         * Get list of corrections for a particular offset
         * @param offset                            Offset of interest
         * @return                                  List of bases representing possible corrections at this offset
         */
        public List<Byte> get(final int offset) {
            if (offset >= size || offset < 0)
                throw new IllegalArgumentException("Illegal call of CorrectionSet.get(): offset must be < size");
            return corrections.get(offset);
        }

        /**
         * Get consensus correction for a particular offset. In this implementation, it just boils down to seeing if
         * byte list associated with offset has identical values. If so, return this base, otherwise return null.
         * @param offset
         * @return                                 Consensus base, or null if no consensus possible.
         */
        public Byte getConsensusCorrection(final int offset) {
            if (offset >= size || offset < 0)
                throw new IllegalArgumentException("Illegal call of CorrectionSet.getConsensusCorrection(): offset must be < size");
            final List<Byte> storedBytes = corrections.get(offset);
            if (storedBytes.isEmpty())
                return null;

            // todo - is there a cheaper/nicer way to compare if all elements in list are identical??
            final byte lastBase = storedBytes.remove(storedBytes.size()-1);
            for (final Byte b: storedBytes) {
                // strict correction rule: all bases must match
                if (b != lastBase)
                    return null;
            }

            // all bytes then are equal:
            return lastBase;

        }



    }
}
