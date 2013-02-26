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

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * The statistics calculator for a specific sample given the interval
 */
class SampleStatistics {
    private final GenomeLoc interval;
    private final ArrayList<LocusStatistics> loci;

    private int[] preSortedDepths = null;
    private int preComputedTotalCoverage = -1;         // avoids re-calculating the total sum (-1 means we haven't pre-computed it yet)

    private int nReads = -1;
    private int nBadMates = -1;

    private SampleStatistics(GenomeLoc interval, ArrayList<LocusStatistics> loci) {
        this.interval = interval;
        this.loci = loci;
        nReads = 0;
        nBadMates = 0;
    }

    public SampleStatistics(GenomeLoc interval) {
        this(interval, new ArrayList<LocusStatistics>(interval.size()));

        // Initialize every loci (this way we don't have to worry about non-existent loci in the object
        for (int i = 0; i < interval.size(); i++)
            this.loci.add(new LocusStatistics());

    }

    public long totalCoverage() {
        if (preComputedTotalCoverage < 0)
            calculateTotalCoverage();
        return preComputedTotalCoverage;
    }

    public double averageCoverage() {
        if (preComputedTotalCoverage < 0)
            calculateTotalCoverage();
        return (double) preComputedTotalCoverage / loci.size();
    }

    /**
     * Calculates the callable statuses of the entire sample
     *
     * @param thresholds the class contains the statistical threshold for making calls
     * @return the callable statuses of the entire sample
     */
    public Set<CallableStatus> getCallableStatuses(ThresHolder thresholds) {
        // We check if reads are present ot prevent div / 0 exceptions
        if (nReads == 0) {
            return Collections.singleton(CallableStatus.NO_READS);
        }

        Set<CallableStatus> output = new HashSet<CallableStatus>();
        Map<CallableStatus, Double> totals = new HashMap<CallableStatus, Double>(CallableStatus.values().length);

        // initialize map
        for (CallableStatus status : CallableStatus.values())
            totals.put(status, 0.0);

        // sum up all the callable statuses for each locus
        for (int i = 0; i < interval.size(); i++) {
            for (CallableStatus status : callableStatus(i, thresholds)) {
                double count = totals.get(status);

                totals.put(status, count + 1);
            }
        }

        double intervalSize = interval.size();

        if (((double) nBadMates / nReads) >= thresholds.getBadMateStatusThreshold())
            output.add(CallableStatus.BAD_MATE);

        if ((totals.get(CallableStatus.COVERAGE_GAPS) / intervalSize) >= thresholds.getCoverageStatusThreshold())
            output.add(CallableStatus.COVERAGE_GAPS);

        if ((totals.get(CallableStatus.LOW_COVERAGE) / intervalSize) >= thresholds.getCoverageStatusThreshold())
            output.add(CallableStatus.LOW_COVERAGE);

        if ((totals.get(CallableStatus.EXCESSIVE_COVERAGE) / intervalSize) >= thresholds.getExcessiveCoverageThreshold())
            output.add(CallableStatus.EXCESSIVE_COVERAGE);

        if ((totals.get(CallableStatus.POOR_QUALITY) / intervalSize) >= thresholds.getQualityStatusThreshold())
            output.add(CallableStatus.POOR_QUALITY);

        if (totals.get(CallableStatus.REF_N) > 0)
            output.add(CallableStatus.REF_N);


        if (output.isEmpty()) {
            output.add(CallableStatus.PASS);
        }

        return output;
    }

    /**
     * Adds a locus to the interval wide stats
     *
     * @param locus      The locus given as a GenomeLoc
     * @param pileup     The pileup of that locus, this exclusively contains the sample
     * @param thresholds the class contains the statistical threshold for making calls
     */
    public void addLocus(GenomeLoc locus, ReadBackedPileup pileup, ThresHolder thresholds) {
        if (!interval.containsP(locus))
            throw new ReviewedStingException(String.format("Locus %s is not part of the Interval %s", locus, interval));

        // a null pileup means there nothing ot add
        if (pileup != null) {

            int locusIndex = locus.getStart() - interval.getStart();

            int rawCoverage = pileup.depthOfCoverage();
            int coverage = thresholds.getFilteredCoverage(pileup);

            LocusStatistics locusData = new LocusStatistics(coverage, rawCoverage);

            loci.set(locusIndex, locusData);

            for (GATKSAMRecord read : pileup.getReads())
                processRead(read, thresholds);
        }
    }

    private void processRead(GATKSAMRecord read, ThresHolder thresholds) {
        // Was this read already processed?
        if (read.getTemporaryAttribute("checkedBadMate") == null) {
            nReads++;
            if (!hasValidMate(read, thresholds))
                nBadMates++;
            read.setTemporaryAttribute("checkedBadMate", true);
        }
    }

    /**
     * returns the callable status of a given locus without taking the reference base into account.
     *
     * @param locusIndex location in the genome to inquire (only one locus)
     * @param thresholds the class contains the statistical threshold for making calls
     * @return the callable status of a locus
     */
    private Set<CallableStatus> callableStatus(int locusIndex, ThresHolder thresholds) {
        LocusStatistics locus = loci.get(locusIndex);

        return locus.callableStatuses(thresholds);
    }

    private void calculateTotalCoverage() {
        preComputedTotalCoverage = 0;
        for (LocusStatistics locus : loci)
            preComputedTotalCoverage += locus.getCoverage();
    }

    public double getQuantileDepth(double percentage) {
        if (preSortedDepths == null)
            getDepthsAsSortedArray();

        return getQuartile(preSortedDepths, percentage);
    }

    static double getQuartile(int[] data, double percentage) {
        int size = data.length;
        if (size == 1)
            return (double) data[0];

        if (percentage == 0.5) {
            return getMedian(data);
        }

        double position = (size - 1.0) / 2;
        if (percentage == 0.25) {
            // if the position is a whole number
            return getMedian(Arrays.copyOfRange(data, 0, (int) position + 1));

        }
        if (percentage == 0.75) {
            if (position % 1 == 0) {
                return getMedian(Arrays.copyOfRange(data, (int) position, size));
            } else {
                return getMedian(Arrays.copyOfRange(data, (int) position + 1, size));
            }
        }
        return -1;
    }

    // Assumes data is sorted
    private static double getMedian(int[] data) {
        double size = (double) data.length;
        if (size == 1)
            return (double) data[0];

        double position = (size - 1.0) / 2;

        if (position % 1 == 0)
            return (double) data[(int) position];

        else {
            double high = (double) data[(int) Math.ceil(position)];
            double low = (double) data[(int) Math.floor(position)];

            return (high + low) / 2;

        }

    }

    private void getDepthsAsSortedArray() {
        preSortedDepths = new int[loci.size()];

        for (int i = 0; i < loci.size(); i++)
            preSortedDepths[i] = loci.get(i).getCoverage();

        Arrays.sort(preSortedDepths);
    }

    boolean hasValidMate(GATKSAMRecord read, ThresHolder thresholds) {
        /** Check the following
         * Does it have a pair?
         * reasonable insert size?
         * inverted?
         * same orientation?
         * same contig?
         * is pair mapped?
         * todo - is forced mate?
         *
         */

        // has NO pair
        if (!read.getReadPairedFlag())
            return false;

        // different contigs
        if (!read.getMateReferenceIndex().equals(read.getReferenceIndex()))
            return false;

        // unmapped
        if (read.getMateUnmappedFlag() || read.getReadUnmappedFlag())
            return false;

        // same orientation
        if (read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())
            return false;

        // inverted
        if (read.getReadNegativeStrandFlag() ==
                read.getAlignmentStart() < read.getMateAlignmentStart())
            return false;

        // TODO note: IGV uses a different algorithm for insert size, there should be a common util class that does this for you
        // mates are too far apart
        if (Math.abs(read.getAlignmentStart() - read.getMateAlignmentStart()) > thresholds.getMaximumInsertSize())
            return false;

        return true;
    }

    public int getnReads() {
        return nReads;
    }

    public int getnBadMates() {
        return nBadMates;
    }
}
