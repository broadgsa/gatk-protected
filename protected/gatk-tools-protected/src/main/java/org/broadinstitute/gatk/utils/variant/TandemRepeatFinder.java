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

package org.broadinstitute.gatk.utils.variant;

import org.broadinstitute.gatk.utils.Utils;

import java.util.Arrays;

/**
 * Utility to find and quantify tandem repeat units in a byte array.
 *
 * <p>
 *     This class provide a more efficient implementation of deprecated
 *     {@link GATKVariantContextUtils#findNumberOfRepetitions(byte[], byte[], boolean)}
 *     and RepeatCovariate which are proven to be inefficient and buggy.
 * </p>
 *
 * <p>
 *     For now it does not change the logic of those methods in order to preserve current behaviour, but this
 *     needs to be revisited at some point with the proper re-evaluation.
 *
 *     Example.
 *
 *     ttcttcttCtgca
 *
 *     Where the current offset is in the capital C, will result in the STR unit returned to be TGCA with only one repeat.
 *     whereas the logical choice is TTC with 3 repeats.
 *
 *     And for further proof, a small modification and its effect:
 *
 *     ttcttcttCttca
 *
 *     Unit T, repeated 2.
 *
 *     I would say it should be 4 TTC instead.
 *
 *     I think we might well be failing to model the actual PCR artifact appropriately:
 *
 *     <a>http://nar.oxfordjournals.org/content/24/14/2807.full</a>
 *     <a>http://www.ncbi.nlm.nih.gov/pubmed/12560493</a>
 *
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TandemRepeatFinder {


    private final byte[] bases;
    private final int maxRepeatCount;
    private final int maxUnitLength;

    public TandemRepeatFinder(final byte[] bases, final int maxUnitLength, final int maxRepeatLength) {
        if (bases == null) throw new IllegalArgumentException();
        if (maxRepeatLength < 0) throw new IllegalArgumentException();
        if (maxUnitLength < 0) throw new IllegalArgumentException();
        this.maxRepeatCount = maxRepeatLength;
        this.maxUnitLength = maxUnitLength;
        this.bases = bases;
    }

    /**
     * Calculates the number of repeated units of certain length starting at a position.
     *
     * <p>
     * The repeat unit is determined by the original byte array passed to this tandem repeat finder and the input
     * offset and length passed to this method based on the following pseudo-code:
     *
     * <pre>
     *     if (length > 0) {
     *         unit = bytes[offset .. (offset + length - 1)]
     *     } else if (length < 0) {
     *         unit = bytes[offset + length + 1  .. offset]
     *     } else { // length == 0
     *         throw IllegalArgumentException() // not allowed.
     *     }
     * </pre>
     * </p>
     *
     * <p>
     *     0 will be returned if given the offset and length, part of the unit falls outside the byte array.
     * </p>
     * <p>
     *     Otherwise,
     *     this method will return the number of repeats (minimum 1 indicating that there is no duplicates) only looking into
     *     a single direction: if <code>length > 0</code> forward in the byte array <code>byte[offset .. END]</code>,
     *     if <code>length < 0</code> then backward in the array <code>byte[0 .. offset]</code>.
     * </p>
     *
     * @param offset the offset in the bases byte for which to start
     * @param length the unit length, a negative indicates a backward unit.
     * @return the number of repeats.
     * @throws IllegalArgumentException if {@code length} is 0 or {@code offset} is outside boundaries: (0 .. bases.length - 1)
     *  where bases is the array passed to this finder at construction.
     */
    protected int numberOfRepeats(final int offset, final int length) {
        if (length == 0) throw new IllegalArgumentException();
        if (offset < 0 || offset >= bases.length) throw new IllegalArgumentException();
        int from = offset;
        int to = offset + length;
        if (to > bases.length || to < -1) return 0;
        final int increment = length < 0 ? -1 : 1;
        final int stop = length < 0? -1 : bases.length;
        int totalLength = 0;
        while (to != stop) {
           if (bases[to] != bases[from]) break;
           to += increment;
           from += increment;
           totalLength++;
        }
        return 1 + totalLength / Math.abs(length);
    }

    public final class Result {

        private final int unitLength;
        private final int unitOffset;
        private final int repeatCount;

        private Result(final int unitOffset, final int unitLength, final int repeatCount) {
            this.unitOffset = unitOffset;
            this.unitLength = unitLength;
            this.repeatCount = repeatCount;
        }

        /**
         * Returns the repeated unit byte sequence.
         * @return never {@code null}.
         */
        public byte[] getUnit() {
            return Arrays.copyOfRange(bases,unitOffset, unitOffset + unitLength);
        }

        /**
         * Returns the original search bases.
         *
         * @return never {@code null}.
         */
        public byte[] getBases() {
            return bases;
        }

        /**
         * Returns the unit offset.
         *
         * @return 0 to {@link #getBases().length - 1}
         */
        public int getUnitOffset() {
            return unitOffset;
        }

        /**
         * Returns the unit length.
         *
         * @return 0 to {@link #getBases().length - 1}
         */
        public int getUnitLength() {
            return unitLength;
        }

        /**
         * Returns the number of repeats of the unit in the input sequence.
         * @return 0 or greater.
         */
        public int getRepeatCount() {
            return repeatCount;
        }

    }

    /**
     * Re-implements {@link RepeatCovariate#findTandemRepeatUnits(byte[], int)}.
     *
     * @param offset search offset.
     * @return never {@code null}.
     */
    public Result findMostRelevantTandemRepeatUnitAt(final int offset) {

        // Notice that this code is not very nice and is rather long but is just a copy of the existing one implemented
        // in RepeatCovariate, eventually this should be improved.

        // first we look forward for a repeat.

        // first we find the best backward
        int bestBWRepeatCount = 0;
        int bestBWOffset = offset;
        int bestBWLength = 1;
        for (int str = 1; str <= maxUnitLength; str++) {
            final int repeatCount = numberOfRepeats(offset, -str);
            if (repeatCount == 0) {
                break;
            } else if ((bestBWRepeatCount = repeatCount) > 1) {
                bestBWOffset = offset - str + 1;
                bestBWLength = str;
                break;
            }
        }

        // The best forward:
        final int bestFWOffset = offset + 1;
        int bestFWLength = 1;
        int bestFWRepeatCount = 0;
        for (int str = 1; str <= maxUnitLength; str++) {
            final int repeatCount = numberOfRepeats(bestFWOffset, str);
            if (repeatCount == 0) {
                break;
            } else if ((bestFWRepeatCount = repeatCount) > 1) {
                bestFWLength = str;
                break;
            }
        }

        // And we combine forward and backwards results; if different forward repeat has priority:
        if (bestFWLength == bestBWLength && Utils.equalRange(bases, bestFWOffset, bases, bestBWOffset, bestFWLength)) {
            return new Result(bestBWOffset, bestBWLength, Math.min(maxRepeatCount, bestBWRepeatCount + bestFWRepeatCount));
        }
          else {
            final int bestFWBackwardRepeatCount = numberOfRepeats(bestFWOffset + bestFWLength - 1, - bestFWLength) - 1;
            return new Result(bestFWOffset, bestFWLength, Math.min(maxRepeatCount, bestFWRepeatCount + bestFWBackwardRepeatCount));
        }
    }
}
