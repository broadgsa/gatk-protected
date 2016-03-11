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

import com.google.java.contract.Ensures;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.AlleleList;
import org.broadinstitute.gatk.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.gatk.utils.genotyper.IndexedSampleList;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.pairhmm.ArrayLoglessPairHMM;
import org.broadinstitute.gatk.utils.pairhmm.Log10PairHMM;
import org.broadinstitute.gatk.utils.pairhmm.LoglessPairHMM;
import org.broadinstitute.gatk.utils.pairhmm.PairHMM;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.util.*;


public class PairHMMIndelErrorModel {
    public static final int BASE_QUAL_THRESHOLD = 20;

    private boolean DEBUG = false;

    private static final int MAX_CACHED_QUAL = 127;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    private static final int START_HRUN_GAP_IDX = 4;
    private static final int MAX_HRUN_GAP_IDX = 20;

    private static final byte MIN_GAP_OPEN_PENALTY = 30;
    private static final byte MIN_GAP_CONT_PENALTY = 10;
    private static final byte GAP_PENALTY_HRUN_STEP = 1; // each increase in hrun decreases gap penalty by this.

    private final byte[] GAP_OPEN_PROB_TABLE;
    private final byte[] GAP_CONT_PROB_TABLE;

    private final PairHMM pairHMM;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////

    static {
        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10, -k/10.);


            baseMatchArray[k] =  Math.log10(1-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public PairHMMIndelErrorModel(byte indelGOP, byte indelGCP, boolean deb, final PairHMM.HMM_IMPLEMENTATION hmmType ) {
        this.DEBUG = deb;

        switch (hmmType) {
            case EXACT:
                pairHMM = new Log10PairHMM(true);
                break;
            case ORIGINAL:
                pairHMM = new Log10PairHMM(false);
                break;
            case LOGLESS_CACHING:
                pairHMM = new LoglessPairHMM();
                break;
            case ARRAY_LOGLESS:
                pairHMM = new ArrayLoglessPairHMM();
                break;
            default:
                throw new UserException.BadArgumentValue("pairHMM", "Specified pairHMM implementation is unrecognized or incompatible with the UnifiedGenotyper. Acceptable options are ORIGINAL, EXACT, LOGLESS_CACHING, or ARRAY_LOGLESS.");
        }

        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new byte[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new byte[MAX_HRUN_GAP_IDX];

        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = indelGOP;
            GAP_CONT_PROB_TABLE[i] = indelGCP;
        }

        double step = GAP_PENALTY_HRUN_STEP/10.0;

        // initialize gop and gcp to their default values
        byte gop = indelGOP;
        byte gcp = indelGCP;

        // all of the following is computed in QUal-space
        for (int i=START_HRUN_GAP_IDX; i < MAX_HRUN_GAP_IDX; i++) {
            gop -= GAP_PENALTY_HRUN_STEP;
            if (gop < MIN_GAP_OPEN_PENALTY)
                gop = MIN_GAP_OPEN_PENALTY;

            gcp -= step;
            if(gcp < MIN_GAP_CONT_PENALTY)
                gcp = MIN_GAP_CONT_PENALTY;
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

    }

    static private void getContextHomopolymerLength(final byte[] refBytes, final int[] hrunArray) {
        // compute forward hrun length, example:
        // AGGTGACCCCCCTGAGAG
        // 001000012345000000
        hrunArray[0] = 0;
        int[] hforward = new int[hrunArray.length];
        int[] hreverse = new int[hrunArray.length];

        for (int i = 1; i < refBytes.length; i++) {
            if (refBytes[i] == refBytes[i-1])
                hforward[i] = hforward[i-1]+1;
            else
                hforward[i] = 0;
        }

        // do similar thing for reverse length, example:
        // AGGTGACCCCCCTGAGAG
        // 021000543210000000
        // and then accumulate with forward values.
        // Total:
        // AGGTGACCCCCCTGAGAG
        // 022000555555000000
        for (int i=refBytes.length-1; i > 0; i--) {
            if (refBytes[i-1] == refBytes[i])
                hreverse[i-1] += hreverse[i]+1;
        }

        for (int i = 1; i < refBytes.length; i++)
            hrunArray[i] = hforward[i]+hreverse[i];
    }


    private void fillGapProbabilities(final int[] hrunProfile,
                                      final byte[] contextLogGapOpenProbabilities,
                                      final byte[] contextLogGapContinuationProbabilities) {
        // fill based on lookup table
        for (int i = 0; i < hrunProfile.length; i++) {
            if (hrunProfile[i] >= MAX_HRUN_GAP_IDX) {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
            }
            else {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[hrunProfile[i]];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[hrunProfile[i]];
            }
        }
    }

    /**
     * Trims the haplotypes in the given map to the provided start/stop.
     *
     * @param haplotypeMap                      the input map
     * @param startLocationInRefForHaplotypes   the start location of the trim
     * @param stopLocationInRefForHaplotypes    the stop location of the trim
     * @param ref                               the reference context (used for debugging only, so can be null)
     * @return a non-null mapping corresponding to the trimmed version of the original;
     *   some elements may be lost if trimming cannot be performed on them (e.g. they fall outside of the region to keep)
     */
    protected static Map<Allele, Haplotype> trimHaplotypes(final Map<Allele, Haplotype> haplotypeMap,
                                                           long startLocationInRefForHaplotypes,
                                                           long stopLocationInRefForHaplotypes,
                                                           final ReferenceContext ref) {
        if ( haplotypeMap == null ) throw new IllegalArgumentException("The input allele to haplotype map cannot be null");

        final LinkedHashMap<Allele, Haplotype> trimmedHaplotypeMap = new LinkedHashMap<>();
        for (final Allele a: haplotypeMap.keySet()) {

            final Haplotype haplotype = haplotypeMap.get(a);

            if (stopLocationInRefForHaplotypes > haplotype.getStopPosition())
                stopLocationInRefForHaplotypes = haplotype.getStopPosition();

            if (startLocationInRefForHaplotypes < haplotype.getStartPosition())
                startLocationInRefForHaplotypes = haplotype.getStartPosition();
            else if (startLocationInRefForHaplotypes > haplotype.getStopPosition())
                startLocationInRefForHaplotypes = haplotype.getStopPosition();

            final long indStart = startLocationInRefForHaplotypes - haplotype.getStartPosition();
            final long indStop =  stopLocationInRefForHaplotypes - haplotype.getStartPosition();
            if ( indStart >= indStop )
                continue;

            // commented out here because we need to make this method static for unit testing
            //if (DEBUG)
            //    System.out.format("indStart: %d indStop: %d WinStart:%d WinStop:%d start: %d stop: %d\n",
            //            indStart, indStop, ref.getWindow().getStart(), ref.getWindow().getStop(), startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes);

            // get the trimmed haplotype-bases array and create a new haplotype based on it. Pack this into the new map
            final byte[] trimmedHaplotypeBases = Arrays.copyOfRange(haplotype.getBases(), (int)indStart, (int)indStop);
            final Haplotype trimmedHaplotype = new Haplotype(trimmedHaplotypeBases, haplotype.isReference());
            trimmedHaplotypeMap.put(a, trimmedHaplotype);
        }
        return trimmedHaplotypeMap;
    }


    public synchronized double[] computeDiploidReadHaplotypeLikelihoods(final ReadBackedPileup pileup,
                                                                        final LinkedHashMap<Allele, Haplotype> haplotypeMap,
                                                                        final ReferenceContext ref,
                                                                        final int eventLength,
                                                                        final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap,
                                                                        final double downsamplingFraction) {
        final int numHaplotypes = haplotypeMap.size();

        final double[][] readLikelihoods = computeGeneralReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, perReadAlleleLikelihoodMap);
        perReadAlleleLikelihoodMap.performPerAlleleDownsampling(downsamplingFraction);
        return getDiploidHaplotypeLikelihoods(numHaplotypes, readLikelihoods);
        
    }

    /**
     * Should we clip a downstream portion of a read because it spans off the end of a haplotype?
     *
     * @param read               the read in question
     * @param refWindowStop      the end of the reference window
     * @return true if the read needs to be clipped, false otherwise
     */
    protected static boolean mustClipDownstream(final GATKSAMRecord read, final int refWindowStop) {
        return ( !read.isEmpty() && read.getSoftStart() < refWindowStop && read.getSoftStart() + read.getReadLength() - 1 > refWindowStop );
    }

    /**
     * Should we clip a upstream portion of a read because it spans off the end of a haplotype?
     *
     * @param read               the read in question
     * @param refWindowStart     the start of the reference window
     * @return true if the read needs to be clipped, false otherwise
     */
    protected static boolean mustClipUpstream(final GATKSAMRecord read, final int refWindowStart) {
        return ( !read.isEmpty() && read.getSoftStart() < refWindowStart && read.getSoftEnd() > refWindowStart );
    }

    @Ensures("result != null && result.length == pileup.getNumberOfElements()")
    public synchronized double[][] computeGeneralReadHaplotypeLikelihoods(final ReadBackedPileup pileup,
                                                                          final LinkedHashMap<Allele, Haplotype> haplotypeMap, 
                                                                          final ReferenceContext ref,
                                                                          final int eventLength, 
                                                                          final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap) {
        final double readLikelihoods[][] = new double[pileup.getNumberOfElements()][haplotypeMap.size()];

        int readIdx=0;
        for (final PileupElement p: pileup) {

            // check if we've already computed likelihoods for this pileup element (i.e. for this read at this location)
            if (perReadAlleleLikelihoodMap.containsPileupElement(p)) {
                Map<Allele,Double> el = perReadAlleleLikelihoodMap.getLikelihoodsAssociatedWithPileupElement(p);
                int j=0;
                for (Allele a: haplotypeMap.keySet()) {
                    readLikelihoods[readIdx][j++] = el.get(a);
                }
            }
            else {
                // extra padding on candidate haplotypes to make sure reads are always strictly contained
                // in them - a value of 1 will in theory do but we use a slightly higher one just for safety sake, mostly
                // in case bases at edge of reads have lower quality.
                final int trailingBases = 3;
                final int refWindowStart = ref.getWindow().getStart() + trailingBases;
                final int refWindowStop  = ref.getWindow().getStop() - trailingBases;

                if (DEBUG) {
                    System.out.format("Read Name:%s, aln start:%d aln stop:%d orig cigar:%s\n",p.getRead().getReadName(), p.getRead().getAlignmentStart(), p.getRead().getAlignmentEnd(), p.getRead().getCigarString());
                }

                GATKSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());

                // if the read extends beyond the downstream (right) end of the reference window, clip it
                if ( mustClipDownstream(read, refWindowStop) )
                    read = ReadClipper.hardClipByReadCoordinates(read, refWindowStop - read.getSoftStart() + 1, read.getReadLength() - 1);

                // if the read extends beyond the upstream (left) end of the reference window, clip it
                if ( mustClipUpstream(read, refWindowStart) )
                    read = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, refWindowStart);

                if (read.isEmpty())
                    continue;

                // hard-clip low quality ends - this may introduce extra H elements in CIGAR string
                read = ReadClipper.hardClipLowQualEnds(read, (byte) BASE_QUAL_THRESHOLD );

                if (read.isEmpty())
                    continue;

                // get bases of candidate haplotypes that overlap with reads
                final long readStart = read.getSoftStart();
                final long readEnd = read.getSoftEnd();

                // see if we want to use soft clipped bases. Aligners may soft clip all bases at insertions because they don't match,
                // but they're actually consistent with the insertion!
                // Rule: if a read starts in interval [eventStart-eventLength,eventStart+1] and we are at an insertion, we'll use all soft clipped bases at the beginning.
                // Conversely, if a read ends at [eventStart,eventStart+eventLength] we'll use all soft clipped bases in the end of the read.
                final long eventStartPos = ref.getLocus().getStart();

                // compute total number of clipped bases (soft or hard clipped) and only use them if necessary
                final boolean softClips = useSoftClippedBases(read, eventStartPos, eventLength);
                final int numStartSoftClippedBases = softClips ? read.getAlignmentStart()- read.getSoftStart() : 0;
                final int numEndSoftClippedBases = softClips ? read.getSoftEnd()- read.getAlignmentEnd() : 0 ;
                final byte [] unclippedReadBases = read.getReadBases();
                final byte [] unclippedReadQuals = read.getBaseQualities();

                /**
                 * Compute genomic locations that candidate haplotypes will span.
                 * Read start and stop locations (variables readStart and readEnd) are the original unclipped positions from SAMRecord,
                 * adjusted by hard clips from Cigar string and by qual-based soft-clipping performed above.
                 * We will propose haplotypes that overlap the read with some padding.
                 * True read start = readStart + numStartSoftClippedBases - ReadUtils.getFirstInsertionOffset(read)
                 * Last term is because if a read starts with an insertion then these bases are not accounted for in readStart.
                 * trailingBases is a padding constant(=3) and we additionally add abs(eventLength) to both sides of read to be able to
                 * differentiate context between two haplotypes
                 */
                final int absEventLength = Math.abs(eventLength);
                long startLocationInRefForHaplotypes = Math.max(readStart + numStartSoftClippedBases - trailingBases - ReadUtils.getFirstInsertionOffset(read) - absEventLength, 0);
                long stopLocationInRefForHaplotypes = readEnd - numEndSoftClippedBases + trailingBases + ReadUtils.getLastInsertionOffset(read) + absEventLength;

                if (DEBUG)
                    System.out.format("orig Start:%d orig stop: %d\n", startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes);

                int readLength = read.getReadLength()-numStartSoftClippedBases-numEndSoftClippedBases;

                if (startLocationInRefForHaplotypes < ref.getWindow().getStart()) {
                    startLocationInRefForHaplotypes = ref.getWindow().getStart();                                       // read starts before haplotype: read will have to be cut numStartSoftClippedBases += ref.getWindow().getStart() - startLocationInRefForHaplotypes;
                }
                else if (startLocationInRefForHaplotypes > ref.getWindow().getStop()) {
                    startLocationInRefForHaplotypes = ref.getWindow().getStop();                                        // read starts after haplotype: read will have to be clipped completely;
                }

                // candidate haplotype cannot go beyond reference context
                if (stopLocationInRefForHaplotypes > ref.getWindow().getStop()) {
                    stopLocationInRefForHaplotypes = ref.getWindow().getStop();                                         // check also if end of read will go beyond reference context
                }

                if (stopLocationInRefForHaplotypes <= startLocationInRefForHaplotypes + readLength) {
                    stopLocationInRefForHaplotypes = startLocationInRefForHaplotypes + readLength-1;                    // if there's an insertion in the read, the read stop position will be less than start + read legnth, but we want to compute likelihoods in the whole region that a read might overlap
                }

                // ok, we now figured out the total number of clipped bases on both ends.
                // Figure out where we want to place the haplotype to score read against

                if (DEBUG)
                    System.out.format("numStartSoftClippedBases: %d numEndSoftClippedBases: %d WinStart:%d WinStop:%d start: %d stop: %d readLength: %d\n",
                            numStartSoftClippedBases, numEndSoftClippedBases, ref.getWindow().getStart(), ref.getWindow().getStop(), startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes, read.getReadLength());

               // LinkedHashMap<Allele,Double> readEl = new LinkedHashMap<Allele,Double>();

                /**
                 * Check if we'll end up with an empty read once all clipping is done
                 */
                if (numStartSoftClippedBases + numEndSoftClippedBases >= unclippedReadBases.length) {
                    int j=0;
                    for (Allele a: haplotypeMap.keySet()) {
                        perReadAlleleLikelihoodMap.add(p,a,0.0);
                        readLikelihoods[readIdx][j++] = 0.0;
                    }
                }
                else {
                    final int endOfCopy = unclippedReadBases.length - numEndSoftClippedBases;
                    final byte[] readBases = Arrays.copyOfRange(unclippedReadBases, numStartSoftClippedBases, endOfCopy);
                    final byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals, numStartSoftClippedBases, endOfCopy);

                    int j=0;

                    final byte[] contextLogGapOpenProbabilities = new byte[readBases.length];
                    final byte[] contextLogGapContinuationProbabilities  = new byte[readBases.length];

                    // get homopolymer length profile for current haplotype
                    final int[] hrunProfile = new int[readBases.length];
                    getContextHomopolymerLength(readBases,hrunProfile);
                    fillGapProbabilities(hrunProfile, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities);

                    // get the base insertion and deletion qualities to use
                    final byte[] baseInsertionQualities, baseDeletionQualities;
                    if ( read.hasBaseIndelQualities() ) {
                        baseInsertionQualities = Arrays.copyOfRange(read.getBaseInsertionQualities(), numStartSoftClippedBases, endOfCopy);
                        baseDeletionQualities = Arrays.copyOfRange(read.getBaseDeletionQualities(), numStartSoftClippedBases, endOfCopy);
                    } else {
                        baseInsertionQualities = contextLogGapOpenProbabilities;
                        baseDeletionQualities = contextLogGapOpenProbabilities;
                    }

                    // Create a new read based on the current one, but with trimmed bases/quals, for use in the HMM
                    final GATKSAMRecord processedRead = GATKSAMRecord.createQualityModifiedRead(read, readBases, readQuals, baseInsertionQualities, baseDeletionQualities);

                    // Pack the shortened read and its associated gap-continuation-penalty array into a map, as required by PairHMM
                    final Map<GATKSAMRecord,byte[]> readGCPArrayMap = Collections.singletonMap(processedRead,contextLogGapContinuationProbabilities);

                    // Create a map of alleles to a new set of haplotypes, whose bases have been trimmed to the appropriate genomic locations
                    final Map<Allele, Haplotype> trimmedHaplotypeMap = trimHaplotypes(haplotypeMap, startLocationInRefForHaplotypes, stopLocationInRefForHaplotypes, ref);

                    // Apparently more than one allele can map to the same haplotype after trimming
                    final Set<Haplotype> distinctHaplotypesSet = new LinkedHashSet<>(trimmedHaplotypeMap.values());
                    final AlleleList<Haplotype> distinctHaplotypesList = new IndexedAlleleList<>(distinctHaplotypesSet.toArray(new Haplotype[distinctHaplotypesSet.size()]));
                    // Get the likelihoods for our clipped read against each of our trimmed haplotypes.
                    final ReadLikelihoods<Haplotype> rl = new ReadLikelihoods<>(
                            new IndexedSampleList(Collections.singletonList("DUMMY_SAMPLE")),distinctHaplotypesList,Collections.singletonMap("DUMMY_SAMPLE",Collections.singletonList(processedRead)));

                    final ReadLikelihoods.Matrix<Haplotype> dummySampleLikelihoods = rl.sampleMatrix(0);
                    pairHMM.computeLikelihoods(rl.sampleMatrix(0), Collections.singletonList(processedRead), readGCPArrayMap);

                    // Pack the original pilup element, each allele, and each associated log10 likelihood into a final map, and add each likelihood to the array
                    for (final Allele a: trimmedHaplotypeMap.keySet()){
                        final Haplotype h = trimmedHaplotypeMap.get(a);
                        final int hIndex = rl.alleleIndex(h);
                        final double readLikelihood = dummySampleLikelihoods.get(hIndex,0);
                        readLikelihoods[readIdx][j++] = readLikelihood;
                        perReadAlleleLikelihoodMap.add(p,a,readLikelihood);
                    }
                }
            }
            readIdx++;
        }

        if (DEBUG) {
            System.out.println("\nLikelihood summary");
            for (readIdx=0; readIdx < pileup.getNumberOfElements(); readIdx++) {
                System.out.format("Read Index: %d ",readIdx);
                for (int i=0; i < readLikelihoods[readIdx].length; i++)
                    System.out.format("L%d: %f ",i,readLikelihoods[readIdx][i]);
                System.out.println();
            }

        }

        return readLikelihoods;
    }

    private boolean useSoftClippedBases(GATKSAMRecord read, long eventStartPos, int eventLength) {
        return !((read.getAlignmentStart() >= eventStartPos-eventLength && read.getAlignmentStart() <= eventStartPos+1) || (read.getAlignmentEnd() >= eventStartPos && read.getAlignmentEnd() <= eventStartPos + eventLength));
    }

//    private int computeFirstDifferingPosition(byte[] b1, byte[] b2) {
//        if (b1.length != b2.length)
//            return 0; // sanity check
//
//        for (int i=0; i < b1.length; i++ ){
//            if ( b1[i]!= b2[i] )
//                return i;
//        }
//        return b1.length;
//    }

    private static double[] getDiploidHaplotypeLikelihoods(final int numHaplotypes, final double readLikelihoods[][]) {
        final double[][] haplotypeLikehoodMatrix = new double[numHaplotypes][numHaplotypes];

        // todo: MAD 09/26/11 -- I'm almost certain this calculation can be simplified to just a single loop without the intermediate NxN matrix
        for (int i=0; i < numHaplotypes; i++) {
            for (int j=i; j < numHaplotypes; j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                for (int readIdx = 0; readIdx < readLikelihoods.length; readIdx++) {
                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    if (Double.isInfinite(readLikelihoods[readIdx][i]) && Double.isInfinite(readLikelihoods[readIdx][j]))
                        continue;
                    final double li = readLikelihoods[readIdx][i];
                    final double lj = readLikelihoods[readIdx][j];
                    haplotypeLikehoodMatrix[i][j] += MathUtils.approximateLog10SumLog10(li, lj) + MathUtils.LOG_ONE_HALF;
                }
            }
        }

        final double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int k=0;
        for (int j=0; j < numHaplotypes; j++) {
            for (int i=0; i <= j; i++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize so that max element is zero.
        return MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
    }
}
