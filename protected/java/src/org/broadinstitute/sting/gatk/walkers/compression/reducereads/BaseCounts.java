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

package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;


/**
 * An object to keep track of the number of occurrences of each base and it's quality.
 *
 * User: depristo
 * Date: 4/8/11
 * Time: 2:55 PM
 */

 public class BaseCounts {
    public final static BaseIndex MAX_BASE_INDEX_WITH_NO_COUNTS = BaseIndex.N;
    public final static byte MAX_BASE_WITH_NO_COUNTS = MAX_BASE_INDEX_WITH_NO_COUNTS.getByte();

    private final int[] counts;       // keeps track of the base counts
    private final long[] sumQuals;    // keeps track of the quals of each base
    private int totalCount = 0;       // keeps track of total count since this is requested so often

    public BaseCounts() {
        counts = new int[BaseIndex.values().length];
        sumQuals = new long[BaseIndex.values().length];
        // Java primitive arrays comes zero-filled, so no need to do it explicitly.
    }

    public static BaseCounts createWithCounts(int[] countsACGT) {
        BaseCounts baseCounts = new BaseCounts();
        baseCounts.counts[BaseIndex.A.index] = countsACGT[0];
        baseCounts.counts[BaseIndex.C.index] = countsACGT[1];
        baseCounts.counts[BaseIndex.G.index] = countsACGT[2];
        baseCounts.counts[BaseIndex.T.index] = countsACGT[3];
        baseCounts.totalCount = countsACGT[0] + countsACGT[1] + countsACGT[2] + countsACGT[3];
        return baseCounts;
    }

    @Requires("other != null")
    public void add(final BaseCounts other) {
        for (final BaseIndex i : BaseIndex.values()) {
            final int otherCount = other.counts[i.index];
            counts[i.index] += otherCount;
            totalCount += otherCount;
        }
    }

    @Requires("other != null")
    public void sub(final BaseCounts other) {
        for (final BaseIndex i : BaseIndex.values()) {
            final int otherCount = other.counts[i.index];
            counts[i.index] -= otherCount;
            totalCount -= otherCount;
        }
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(final byte base) {
        final BaseIndex i = BaseIndex.byteToBase(base);
        counts[i.index]++;
        totalCount++;
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(final BaseIndex base, final byte qual) {
        counts[base.index]++;
        totalCount++;
        sumQuals[base.index] += qual;
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) - 1")
    public void decr(final byte base) {
        final BaseIndex i = BaseIndex.byteToBase(base);
        counts[i.index]--;
        totalCount--;
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) - 1")
    public void decr(final BaseIndex base, final byte qual) {
        counts[base.index]--;
        totalCount--;
        sumQuals[base.index] -= qual;
    }

    @Ensures("result >= 0")
    public long getSumQuals(final byte base) {
        return getSumQuals(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public long getSumQuals(final BaseIndex base) {
        return sumQuals[base.index];
    }

    @Ensures("result >= 0")
    public byte averageQuals(final byte base) {
        return averageQuals(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public byte averageQuals(final BaseIndex base) {
        return (byte) (getSumQuals(base) / countOfBase(base));
    }

    @Ensures("result >= 0")
    public int countOfBase(final byte base) {
        return countOfBase(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public int countOfBase(final BaseIndex base) {
        return counts[base.index];
    }

    @Ensures("result >= 0")
    public long sumQualsOfBase(final BaseIndex base) {
        return sumQuals[base.index];
    }

    @Ensures("result >= 0")
    public byte averageQualsOfBase(final BaseIndex base) {
        return (byte) (sumQualsOfBase(base) / countOfBase(base));
    }


    @Ensures("result >= 0")
    public int totalCount() {
        return totalCount;
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     *
     * @param base     base
     * @return the proportion of this base over all other bases
     */
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportion(final byte base) {
        return baseCountProportion(BaseIndex.byteToBase(base));
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     *
     * @param baseIndex    base
     * @return the proportion of this base over all other bases
     */
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportion(final BaseIndex baseIndex) {
        return (totalCount == 0) ? 0.0 : (double)counts[baseIndex.index] / (double)totalCount;
    }

    @Ensures("result != null")
    public String toString() {
        StringBuilder b = new StringBuilder();
        for (final BaseIndex i : BaseIndex.values()) {
            b.append(i.toString()).append("=").append(counts[i.index]).append(",");
        }
        return b.toString();
    }

    public byte baseWithMostCounts() {
        return baseIndexWithMostCounts().getByte();
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostCounts() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (counts[i.index] > counts[maxI.index])
                maxI = i;
        }
        return maxI;
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostCountsWithoutIndels() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (i.isNucleotide() && counts[i.index] > counts[maxI.index])
                maxI = i;
        }
        return maxI;
    }

    public byte baseWithMostProbability() {
        return baseIndexWithMostProbability().getByte();
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostProbability() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (sumQuals[i.index] > sumQuals[maxI.index])
                maxI = i;
        }
        return (sumQuals[maxI.index] > 0L ? maxI : baseIndexWithMostCounts());
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostProbabilityWithoutIndels() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (i.isNucleotide() && sumQuals[i.index] > sumQuals[maxI.index])
                maxI = i;
        }
        return (sumQuals[maxI.index] > 0L ? maxI : baseIndexWithMostCountsWithoutIndels());
    }

    @Ensures("result >=0")
    public int totalCountWithoutIndels() {
        return totalCount - counts[BaseIndex.D.index] - counts[BaseIndex.I.index];
    }

    /**
     * Calculates the proportional count of a base compared to all other bases except indels (I and D)
     *
     * @param base      base
     * @return the proportion of this base over all other bases except indels
     */
    @Requires("base.isNucleotide()")
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportionWithoutIndels(final BaseIndex base) {
        final int total = totalCountWithoutIndels();
        return (total == 0) ? 0.0 : (double)counts[base.index] / (double)total;
    }

    public int[] countsArray() {
        return counts.clone();
    }
}
