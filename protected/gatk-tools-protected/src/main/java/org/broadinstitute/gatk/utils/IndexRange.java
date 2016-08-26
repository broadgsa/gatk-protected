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

package org.broadinstitute.gatk.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.*;

/**
 * Represents 0-based integer index range.
 *
 * <p>
 * It represents an integer index range as the pair values:
 * <dl>
 *     <dt>{@link #from}</dt>
 *     <dd>- index of the first element in range (i.e. inclusive).</dd>
 *     <dt>{@link #to}</dt>
 *     <dd>- index of the element following the last element in range (i.e. exclusive).</dd>
 * </dl>
 * </p>
 *
 * <p>
 *     This class is intended to specify a valid index range in arrays or ordered collections.
 * </p>
 *
 * <p>
 *     All instances are constraint so that neither <code>from</code> nor <code>to</code> can
 *     be negative nor <code>from</code> can be larger than <code>to</code>.
 * </p>
 *
 * <p>
 *     You can use {@link #isValidLength(int) isValidFor(length)} to verify that a range instance represents a valid
 *     range for an 0-based indexed object with {@code length} elements.
 * </p>
 */
public final class IndexRange {

    /**
     * First index in the range.
     * <p>
     *     It won't ever be negative nor greater than {@link #to}.
     * </p>
     */
    public final int from;

    /**
     * Index following the last index included in the range.
     *
     * <p>
     *     It won't ever be negative nor less than {@link #from}.
     * </p>
     */
    public final int to;

    /**
     * Creates a new range given its {@code from} and {@code to} indices.
     *
     * @param fromIndex the {@code from} index value.
     * @param toIndex   the {@code to} index value.
     * @throws IllegalArgumentException if {@code fromIndex} is larger than {@code toIndex} or either is
     *                                  negative.
     */
    public IndexRange(final int fromIndex, final int toIndex) {
        Utils.validateArg(fromIndex <= toIndex, "the range size cannot be negative");
        Utils.validateArg(fromIndex >= 0, "the range cannot contain negative indices");
        from = fromIndex;
        to = toIndex;
    }

    /**
     * Checks whether this range is valid for a collection or array of a given size.
     *
     * <p>
     *     It assume that 0 is the first valid index for target indexed object which is true
     *     for Java Arrays and mainstream collections.
     * </p>
     *
     * <p>
     *     If the input length is less than 0, thus incorrect, this method is guaranteed to return
     *     {@code false}. No exception is thrown.
     * </p>
     *
     *
     * @param length the targeted collection or array length.
     * @return {@code true} if this range is valid for that {@code length}, {@code false} otherwise.
     */
    public boolean isValidLength(final int length) {
        return to <= length;
    }

    /**
     * Returns number indexes expanded by this range.
     *
     * @return 0 or greater.
     */
    public int size() {
        return to - from;
    }

    /**
     * Iterate through all indexes in the range in ascending order to be processed by the
     * provided {@link IntConsumer integer consumer} lambda function.
     *
     * <p>
     *     Exceptions thrown by the execution of the index consumer {@code lambda}
     *     will be propagated to the caller immediately thus stopping early and preventing
     *     further indexes to be processed.
     * </p>
     * @param lambda the index consumer lambda.
     * @throws IllegalArgumentException if {@code lambda} is {@code null}.
     * @throws RuntimeException if thrown by {@code lambda} for some index.
     * @throws Error if thrown by {@code lambda} for some index.
     */
    public void forEach(final IntConsumer lambda) {
        Utils.nonNull(lambda, "the lambda function cannot be null");
        for (int i = from; i < to; i++) {
            lambda.accept(i);
        }
    }

    /**
     * Apply an int -> double function to this range, producing a double[]
     *
     * @param lambda the int -> double function
     */
    public double[] mapToDouble(final IntToDoubleFunction lambda) {
        Utils.nonNull(lambda, "the lambda function cannot be null");
        final double[] result = new double[size()];
        for (int i = from; i < to; i++) {
            result[i - from] = lambda.applyAsDouble(i);
        }
        return result;
    }

    /**
     * Sums the values of an int -> double function applied to this range
     *
     * @param lambda the int -> double function
     */
    public double sum(final IntToDoubleFunction lambda) {
        Utils.nonNull(lambda, "the lambda function cannot be null");
        double result = 0;
        for (int i = from; i < to; i++) {
            result += lambda.applyAsDouble(i);
        }
        return result;
    }

    /**
     * Apply an int -> int function to this range, producing an int[]
     *
     * @param lambda the int -> int function
     */
    public int[] mapToInteger(final IntUnaryOperator lambda) {
        Utils.nonNull(lambda, "the lambda function cannot be null");
        final int[] result = new int[size()];
        for (int i = from; i < to; i++) {
            result[i - from] = lambda.applyAsInt(i);
        }
        return result;
    }

    /**
     * Find the elements of this range for which an int -> boolean predicate is true
     *
     * @param predicate the int -> boolean predicate
     * @return
     */
    public List<Integer> filter(final IntPredicate predicate) {
        Utils.nonNull(predicate, "predicate may not be null");
        final List<Integer> result = new ArrayList<>();
        forEach(i -> {if (predicate.test(i)) result.add(i); } );
        return result;
    }

    @Override
    public boolean equals(final Object other) {
        if (other == this) {
            return true;
        } else if (!(other instanceof IndexRange)) {
            return false;
        } else {
            final IndexRange otherCasted = (IndexRange) other;
            return otherCasted.from == this.from && otherCasted.to == this.to;
        }
    }

    @Override
    public int hashCode() {
        // Inspired on {@link Arrays#hashCode(Object[])}.
        return (( 31 + Integer.hashCode(from) ) * 31 ) + Integer.hashCode(to);
    }

    @Override
    public String toString() {
        return String.format("%d-%d",from,to);
    }
}
