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

import java.util.*;

/**
 * generic utility function that error corrects kmers based on counts
 *
 * This class provides a generic facility for remapping kmers (byte[] of constant size)
 * that occur infrequently to those that occur frequently, based on their simple edit distance
 * as measured by mismatches.
 *
 * The overall workflow of using this class is simple.  First, you create the class with
 * parameters determining how the error correction should proceed.  Next, you provide all
 * of the kmers you see in your data.  Once all kmers have been added, you call computeErrorCorrectionMap
 * to tell this class that all kmers have been added and its time to determine error correcting
 * mapping from observed kmers to corrected kmers.  This correction looks for low-count (as determined
 * by maxCountToCorrect) kmers and chooses the best kmer (minimizing mismatches) among those
 * with at least minCountOfKmerToBeCorrection occurrences to error correct the kmer to.  If
 * there is no kmer with less than maxMismatchesToCorrect then the kmer will be mapped to
 * null, indicating the kmer should not be used.
 *
 * TODO -- for ease of implementation this class uses strings instead of byte[] as those cannot
 * TODO -- be added to hashmaps (more specifically, those don't implement .equals).  A more efficient
 * TODO -- version would use the byte[] directly
 *
 * User: depristo
 * Date: 3/8/13
 * Time: 1:16 PM
 */
public class KMerErrorCorrector {
    /**
     * A map of for each kmer to its num occurrences in addKmers
     */
    Map<String, Integer> countsByKMer = new HashMap<String, Integer>();

    /**
     * A map from raw kmer -> error corrected kmer
     */
    Map<String, String> rawToErrorCorrectedMap = null;

    final int kmerLength;
    final int maxCountToCorrect;
    final int maxMismatchesToCorrect;
    final int minCountOfKmerToBeCorrection;

    /**
     * Create a new kmer corrector
     *
     * @param kmerLength the length of kmers we'll be counting to error correct, must be >= 1
     * @param maxCountToCorrect kmers with < maxCountToCorrect will try to be error corrected to another kmer, must be >= 0
     * @param maxMismatchesToCorrect the maximum number of mismatches between a to-be-corrected kmer and its
     *                               best match that we attempt to error correct.  If no sufficiently similar
     *                               kmer exists, it will be remapped to null.  Must be >= 1
     * @param minCountOfKmerToBeCorrection the minimum count of a kmer to be considered a target for correction.
     *                                     That is, kmers that need correction will only be matched with kmers
     *                                     with at least minCountOfKmerToBeCorrection occurrences.  Must be >= 1
     */
    public KMerErrorCorrector(final int kmerLength,
                              final int maxCountToCorrect,
                              final int maxMismatchesToCorrect,
                              final int minCountOfKmerToBeCorrection) {
        if ( kmerLength < 1 ) throw new IllegalArgumentException("kmerLength must be > 0 but got " + kmerLength);
        if ( maxCountToCorrect < 0 ) throw new IllegalArgumentException("maxCountToCorrect must be >= 0 but got " + maxCountToCorrect);
        if ( maxMismatchesToCorrect < 1 ) throw new IllegalArgumentException("maxMismatchesToCorrect must be >= 1 but got " + maxMismatchesToCorrect);
        if ( minCountOfKmerToBeCorrection < 1 ) throw new IllegalArgumentException("minCountOfKmerToBeCorrection must be >= 1 but got " + minCountOfKmerToBeCorrection);

        this.kmerLength = kmerLength;
        this.maxCountToCorrect = maxCountToCorrect;
        this.maxMismatchesToCorrect = maxMismatchesToCorrect;
        this.minCountOfKmerToBeCorrection = minCountOfKmerToBeCorrection;
    }

    /**
     * For testing purposes
     *
     * @param kmers
     */
    protected void addKmers(final String ... kmers) {
        for ( final String kmer : kmers )
            addKmer(kmer, 1);
        computeErrorCorrectionMap();
    }

    /**
     * Add a kmer that occurred kmerCount times
     *
     * @param rawKmer a kmer
     * @param kmerCount the number of occurrences
     */
    public void addKmer(final byte[] rawKmer, final int kmerCount) {
        addKmer(new String(rawKmer), kmerCount);
    }


    /**
     * Get the error corrected kmer for rawKmer
     *
     * @param rawKmer a kmer that was already added that we want to get an error corrected version for
     * @return an error corrected kmer to use instead of rawKmer.  May be == rawKmer if no error correction
     *         is not necessary.  May be null, indicating the rawKmer shouldn't be used at all
     */
    public byte[] getErrorCorrectedKmer(final byte[] rawKmer) {
        final String result = getErrorCorrectedKmer(new String(rawKmer));
        return result == null ? null : result.getBytes();
    }

    /**
     * Indicate that no more kmers will be added to the kmer error corrector, so that the
     * error correction data structure should be computed from the added kmers.  Enabled calls
     * to getErrorCorrectedKmer, and disable calls to addKmer.
     */
    public void computeErrorCorrectionMap() {
        if ( countsByKMer == null )
            throw new IllegalStateException("computeErrorCorrectionMap can only be called once");

        final LinkedList<String> needsCorrection = new LinkedList<String>();
        final LinkedList<String> goodKmers = new LinkedList<String>();

        rawToErrorCorrectedMap = new HashMap<String, String>();
        for ( Map.Entry<String, Integer> kmerCounts: countsByKMer.entrySet() ) {
            if ( kmerCounts.getValue() <= maxCountToCorrect )
                needsCorrection.add(kmerCounts.getKey());
            else {
                // todo -- optimization could make not in map mean ==
                rawToErrorCorrectedMap.put(kmerCounts.getKey(), kmerCounts.getKey());

                // only allow corrections to kmers with at least this count
                if ( kmerCounts.getValue() >= minCountOfKmerToBeCorrection )
                    goodKmers.add(kmerCounts.getKey());
            }
        }

        for ( final String toCorrect : needsCorrection ) {
            final String corrected = findClosestKMer(toCorrect, goodKmers);
            rawToErrorCorrectedMap.put(toCorrect, corrected);
        }

        // cleanup memory -- we don't need the counts for each kmer any longer
        countsByKMer = null;
    }

    protected void addKmer(final String rawKmer, final int kmerCount) {
        if ( rawKmer.length() != kmerLength ) throw new IllegalArgumentException("bad kmer length " + rawKmer + " expected size " + kmerLength);
        if ( kmerCount < 0 ) throw new IllegalArgumentException("bad kmerCount " + kmerCount);
        if ( countsByKMer == null ) throw new IllegalStateException("Cannot add kmers to an already finalized error corrector");

        final Integer countFromMap = countsByKMer.get(rawKmer);
        final int count = countFromMap == null ? 0 : countFromMap;
        countsByKMer.put(rawKmer, count + kmerCount);
    }

    protected String findClosestKMer(final String kmer, final Collection<String> goodKmers) {
        String bestMatch = null;
        int minMismatches = Integer.MAX_VALUE;

        for ( final String goodKmer : goodKmers ) {
            final int mismatches = countMismatches(kmer, goodKmer);
            if ( mismatches < minMismatches ) {
                minMismatches = mismatches;
                bestMatch = goodKmer;
            }
        }

        return minMismatches > maxMismatchesToCorrect ? null : bestMatch;
    }

    protected int countMismatches(final String one, final String two) {
        int mismatches = 0;
        for ( int i = 0; i < one.length(); i++ )
            mismatches += one.charAt(i) == two.charAt(i) ? 0 : 1;
        return mismatches;
    }

    protected String getErrorCorrectedKmer(final String rawKmer) {
        if ( rawToErrorCorrectedMap == null ) throw new IllegalStateException("Cannot get error corrected kmers until after computeErrorCorrectionMap has been called");
        if ( rawKmer.length() != kmerLength ) throw new IllegalArgumentException("bad kmer length " + rawKmer + " expected size " + kmerLength);
        return rawToErrorCorrectedMap.get(rawKmer);
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("KMerErrorCorrector{");
        for ( Map.Entry<String, String> toCorrect : rawToErrorCorrectedMap.entrySet() ) {
            final boolean correcting = ! toCorrect.getKey().equals(toCorrect.getValue());
            if ( correcting )
                b.append(String.format("%n\t%s / %d -> %s / %d [correcting? %b]",
                        toCorrect.getKey(), getCounts(toCorrect.getKey()),
                        toCorrect.getValue(), getCounts(toCorrect.getValue()),
                        correcting));
        }
        b.append("\n}");
        return b.toString();
    }

    /**
     * Get a simple count estimate for printing for kmer
     * @param kmer the kmer
     * @return an integer count for kmer
     */
    private int getCounts(final String kmer) {
        if ( kmer == null ) return 0;
        final Integer count = countsByKMer == null ? -1 : countsByKMer.get(kmer);
        if ( count == null )
            throw new IllegalArgumentException("kmer not found in counts -- bug " + kmer);
        return count;
    }
}
