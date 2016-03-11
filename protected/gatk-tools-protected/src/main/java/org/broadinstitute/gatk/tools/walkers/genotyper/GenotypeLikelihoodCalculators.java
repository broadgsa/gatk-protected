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

package org.broadinstitute.gatk.tools.walkers.genotyper;

import java.util.Arrays;

/**
 * Genotype likelihood calculator utility.
 *
 * <p>
 *     This class provide genotype likelihood calculators with any number of alleles able given an arbitrary ploidy and allele
 *     count (number of distinct alleles).
 * </p>
 *
 * <p>
 *     This class is thread-safe.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeLikelihoodCalculators {

    /**
     * Maximum possible number of genotypes that this calculator can handle.
     */
    public static final int MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY = 1000;

    /**
     * Mark to indicate genotype-count overflow due to a large number of allele and ploidy;
     */
    protected static final int GENOTYPE_COUNT_OVERFLOW = -1;

    /**
     * The current maximum allele index supported by the tables.
     * <p>
     *     Its initial value indicates the initial capacity of the shared {@link #alleleFirstGenotypeOffsetByPloidy} table.
     *     Feel free to change it to anything reasonable that is non-negative.
     * </p>
     */
    private static int maximumAllele = 1; // its initial value is the initial capacity of the shared tables.

    /**
     * The current maximum ploidy supported by the tables.
     * <p>
     *     Its initial value indicates the initial capacity of the shared {@link #genotypeTableByPloidy}. Feel free
     *     to change it to anything reasonable that is non-negative.
     * </p>
     */
    private static int maximumPloidy = 2; // its initial value is the initial capacity of the shared tables.

    /**
     * Shared copy of the offset table as described in {@link #buildGenotypeAlleleCountsTable(int, int, int[][])}.
     *
     * This reference holds the largest requested so far in terms of maximum-allele and maximum-ploidy.
     */
    private volatile static int[][] alleleFirstGenotypeOffsetByPloidy =
            buildAlleleFirstGenotypeOffsetTable(maximumPloidy, maximumAllele);


    /**
     * Shared table of genotypes give the ploidy sorted by their index in the likelihood array.
     *
     * <p>
     *  Its format is described in {@link #buildGenotypeAlleleCountsTable(int, int, int[][])}.
     * </p>
     */
    private volatile static GenotypeAlleleCounts[][] genotypeTableByPloidy =
            buildGenotypeAlleleCountsTable(maximumPloidy,maximumAllele,alleleFirstGenotypeOffsetByPloidy);



    /**
     * Build the table with the genotype offsets based on ploidy and the maximum allele index with representation
     * in the genotype.
     * <p>
     * The result is a matrix containing the offset of the first genotype that contain a particular allele
     * stratified by ploidy.
     * <p>
     *     Row (first dimension) represent the ploidy, whereas
     *     the second dimension represents the allele.
     * </p>
     *
     * <p>
     *     Thus the value a position <i>[p][a]</i> indicates how many genotypes of ploidy <i>p</i> there are before the first
     *     one that contains allele <i>a</i>. <br/>
     *
     *     For example, considering ploidy 3 and alleles A, B, C, D, etc ... (indexed 0, 1, 2, ... respectively):
     *     <br/>
     *     [3][A] == [3][0] == 0 as the first genotype AAA contains A.
     *     <br/>
     *     [3][C] == [3][2] == 4 as the first genotype that contains C, AAC follows: AAA AAB ABB BBB
     *     <br/>
     *     [4][D] == [4][3] == 14  as the first genotype that contains D, AAAD follows: AAAA AAAB AABB ABBB BBBB AAAC
     *     AABC ABBC BBBC AACC ABCC BBCC ACCC BCCC CCCC.
     *
     * </p>
     *
     * <p>
     *     This value are calculated recursively as follows:
     * </p>
     * <pre>
     *
     *     Offset[p][a] := Offset[p-1][a] + Offset[p][a-1] when a > 0, p > 0
     *                     0                               when a == 0
     *                     1                               otherwise
     *
     *
     *         0 1 1  1  1  1   1 ...
     *         0 1 2  3  4  5   6 ...
     *         0 1 3  6 10 15  21 ...
     *         0 1 4 10 20 35  56 ...
     *         0 1 5 15 35 70 126 ...
     *         0 ..................
     * </pre>
     *
     * <p>
     *    Note: if someone can come with a close form computable 0(1) (respect to ploidy and allele count)
     *     please let the author know.
     * </p>
     *
     * <p>
     *     The matrix is guaranteed to have as many rows as indicated by {@code maximumPloidy} + 1; the first
     *     row refers to the special case of ploidy == 0, the second row to ploidy 1 and so forth. Thus the ploidy
     *     matches the index.
     * </p>
     * <p>
     *     The matrix is guaranteed to have as many columns as indicate by {@code maximumAllele} + 1. In this case however
     *     the first allele index 0 is a sense allele (typically the reference allele). The reason to have at least the total
     *     genotype count up to allele count {@link @alleleCapacity} that is equal to the offset of the first genotype
     *     of the following allele; thus we need an extra one.
     * </p>
     *
     * <p>
     *     Although it might seem non-sense to have genotypes of ploidy 0. The values in the first row are used when
     *     filling up values in row 1 and so forth so it is present for programmatic convenience.
     *     Offsets in this row are 0 for the first column and 1 for any others.
     * </p>
     *
     * @param maximumPloidy maximum supported ploidy.
     * @param maximumAllele maximum supported allele index.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative.
     *
     * @return never {@code null}, the matrix described with enough information to address
     *       problems concerning up to the requested maximum allele index and ploidy.
     */
    private static int[][] buildAlleleFirstGenotypeOffsetTable(final int maximumPloidy, final int maximumAllele) {
        checkPloidyAndMaximumAllele(maximumPloidy, maximumAllele);
        final int rowCount = maximumPloidy + 1;
        final int colCount = maximumAllele + 1;
        final int[][] result = new int[rowCount][colCount];

        // Ploidy 0 array must be { 0, 1, 1, ...., 1}
        Arrays.fill(result[0],1,colCount,1);
        // Now we take care of the rest of ploidies.
        // We leave the first allele offset to it correct value 0 by starting with allele := 1.
        for (int ploidy = 1; ploidy < rowCount; ploidy++)
            for (int allele = 1; allele < colCount; allele++) {
                result[ploidy][allele] = result[ploidy][allele - 1] + result[ploidy - 1][allele];
                if (result[ploidy][allele] < result[ploidy][allele - 1])
                    result[ploidy][allele] = GENOTYPE_COUNT_OVERFLOW;
            }
        return result;
    }

    /**
     * Composes a table with the lists of all possible genotype allele counts given the the ploidy and maximum allele index.
     * <p>
     *     The resulting matrix has at least as many rows as {@code maximumPloidy } + 1 as the first row with index 0 correspond
     *     to ploidy == 0. Each row array has as many positions as necessary to contain all possible genotype-allele-counts in increasing order.
     *     This quantity varies with the ploidy.
     * </p>
     *
     * <p>
     *     Therefore <code>result[3][4]</code> would contain the 5th genotype with ploidy 3, and <code>result[4].length</code>
     *     would be equal to the count of possible genotypes for ploidy 4.
     * </p>
     *
     * @param maximumPloidy maximum ploidy to use in queries to the resulting table.
     * @param maximumAllele maximum allele index to use in queries to the resulting table.
     * @param offsetTable an allele first genotype offset table as constructed using {@link #buildAlleleFirstGenotypeOffsetTable(int, int)}
     *                    that supports at least up to {@code maximumAllele} and {@code maximumPloidy}.
     *
     * @throws IllegalArgumentException if {@code maximumPloidy} or {@code maximumAllele} is negative, or {@code offsetTable} is {@code null},
     *   or it does not have the capacity to handle the requested maximum ploidy or allele index.
     *
     * @return never {@code null}.
     */
    private static GenotypeAlleleCounts[][] buildGenotypeAlleleCountsTable(final int maximumPloidy, final int maximumAllele, final int[][] offsetTable) {
        checkPloidyAndMaximumAllele(maximumPloidy, maximumAllele);
        checkOffsetTableCapacity(offsetTable,maximumPloidy,maximumAllele);
        final int rowCount = maximumPloidy + 1;
        final GenotypeAlleleCounts[][] result = new GenotypeAlleleCounts[rowCount][]; // each row has a different number of columns.

        for (int ploidy = 0; ploidy <= maximumPloidy; ploidy++)
            result[ploidy] = buildGenotypeAlleleCountsArray(ploidy, maximumAllele, offsetTable);

        return result;
    }

    /**
     * Builds a genotype-allele-counts array given the genotype ploidy and how many genotype you need.
     * <p>
     *     The result is guarantee to have exactly {@code length} positions and the elements are sorted
     *     in agreement with the standard way to display genotypes following the VCF standard.
     * </p>
     *
     * <p> Notice that is possible to request ploidy ==0. In that case the resulting array will have repetitions
     * of the empty genotype allele count.
     * </p>
     *
     * <p>
     *     For example,
     *
     *     <pre>
     *         ploidy = 1, length = 5 : [ {A}, {B}, {C}, {D}, {E} ]
     *         ploidy = 2, length = 7 : [ {AA}, {AB}, {BB}, {AC}, {BC}, {CC}, {AD}
     *         ploidy = 3, length = 10 : [ {AAA}, {AAB}, {ABB}, {BBB}, {AAC}, {ABC}, {BBC}, {BCC}, {CCC}, {AAD} ]
     *     </pre>
     * </p>
     *
     * @param ploidy requested ploidy.
     * @param alleleCount number of different alleles that the genotype table must support.
     * @param genotypeOffsetTable table with the offset of the first genotype that contain an allele given
     *                            the ploidy and its index.
     *
     * @throws IllegalArgumentException if {@code ploidy} or {@code length} is negative.
     *
     * @return never {@code null}, follows the specification above.
     */
    private static GenotypeAlleleCounts[] buildGenotypeAlleleCountsArray(final int ploidy, final int alleleCount, final int[][] genotypeOffsetTable) {
        if (ploidy < 0)
            throw new IllegalArgumentException("the requested ploidy cannot be negative: " + ploidy);
        if (alleleCount < 0)
            throw new IllegalArgumentException("the requested maximum allele cannot be negative: " + alleleCount);
        final int length = genotypeOffsetTable[ploidy][alleleCount];
        final int strongRefLength = length == GENOTYPE_COUNT_OVERFLOW ? MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY : Math.min(length, MAXIMUM_STRONG_REF_GENOTYPE_PER_PLOIDY);
        final GenotypeAlleleCounts[] result = new GenotypeAlleleCounts[strongRefLength];
        result[0] = GenotypeAlleleCounts.first(ploidy);
        for (int genotypeIndex = 1; genotypeIndex < strongRefLength; genotypeIndex++)
            result[genotypeIndex] = result[genotypeIndex-1].next();
        return result;
    }

    /**
     * Cached log10 values for the first integer up to the maximum ploidy requested thus far.
     */
    private volatile static double[] ploidyLog10;

    // Initialize {@link #ploidyLog10}.
    static {
        ploidyLog10 = new double[maximumPloidy + 1];
        for (int i = 0; i <= maximumPloidy; i++)
            ploidyLog10[i] = Math.log10(i);
    }

    /**
     * Returns an instance given its ploidy and the number of alleles.
     *
     * @param alleleCount the required allele-count.
     * @param ploidy the required ploidy-count.
     *
     * @throws IllegalArgumentException if either {@code ploidy} or {@code alleleCount} is {@code null}, or
     *      the resulting number of genotypes is too large.
     *
     * @return never {@code null}.
     */
    public static GenotypeLikelihoodCalculator getInstance(final int ploidy,
                                                   final int alleleCount) {
        checkPloidyAndMaximumAllele(ploidy, alleleCount);

        // Non-thread safe (fast) check on tables capacities,
        // if not enough capacity we expand the tables in a thread-safe manner:
        if (alleleCount > maximumAllele || ploidy > maximumPloidy)
            ensureCapacity(alleleCount, ploidy);

        // At this point the tables must have at least the requested capacity, likely to be much more.
        return new GenotypeLikelihoodCalculator(ploidy,alleleCount,alleleFirstGenotypeOffsetByPloidy,genotypeTableByPloidy,ploidyLog10);
    }

    /**
     * Thread safe update of shared tables
     *
     * @param requestedMaximumAllele the new requested maximum allele maximum.
     * @param requestedMaximumPloidy the new requested ploidy maximum.
     */
    private synchronized static void ensureCapacity(final int requestedMaximumAllele, final int requestedMaximumPloidy) {

        final boolean needsToExpandAlleleCapacity = requestedMaximumAllele > maximumAllele;
        final boolean needsToExpandPloidyCapacity = requestedMaximumPloidy > maximumPloidy;

        // Double check with the lock on to avoid double work.
        if (!needsToExpandAlleleCapacity && !needsToExpandPloidyCapacity)
            return;

        final int newMaximumPloidy = Math.max(maximumPloidy,requestedMaximumPloidy);
        final int newMaximumAllele = Math.max(maximumAllele,requestedMaximumAllele);

        // Update tables first.
        alleleFirstGenotypeOffsetByPloidy = buildAlleleFirstGenotypeOffsetTable(newMaximumPloidy,newMaximumAllele);
        genotypeTableByPloidy = buildGenotypeAlleleCountsTable(newMaximumPloidy,newMaximumAllele,alleleFirstGenotypeOffsetByPloidy);

        if (needsToExpandPloidyCapacity)
            ploidyLog10 = ploidyLog10Extension(newMaximumPloidy);

        // Since tables are volatile fields, it is guaranteed that tables changes will be seen before
        // than any change on ploidyCapacity and alleleCapacity ensuring that the non-thread safe
        // capacity verification test in {@link #getInstance} wont ever allow a thread
        // to proceed to use a table without the required capacity.
        // Just after updating tables update the capacity fields:

        if (needsToExpandAlleleCapacity)
            maximumAllele = requestedMaximumAllele;
        if (needsToExpandPloidyCapacity)
            maximumPloidy = requestedMaximumPloidy;
    }

    /**
     * Extends the existing {@link #ploidyLog10} with more log10 as needed by maximum-ploidy expansion.
     * @param newMaximumPloidy the new maximum ploidy.
     *
     * @return never code {@code null}.
     */
    private static double[] ploidyLog10Extension(final int newMaximumPloidy) {
        final int start = ploidyLog10.length;
        final double[] result = Arrays.copyOf(ploidyLog10,newMaximumPloidy + 1);
        for (int i = start; i < result.length; i++)
            result[i] = Math.log10(i);
        return result;
    }

    /**
     * Perform value checks on maximumPloidy and allele passed to diverse methods in this class.
     * <p>
     *     Throws an exception if there is any issues.
     * </p>
     *
     * @param ploidy the maximum ploidy value.
     * @param maximumAllele the maximum allele value.
     *
     * @throws IllegalArgumentException if either value is negative.
     */
    private static void checkPloidyAndMaximumAllele(final int ploidy, final int maximumAllele) {
        if (ploidy < 0)
            throw new IllegalArgumentException("the ploidy provided cannot be negative: " + ploidy);
        if (maximumAllele < 0)
            throw new IllegalArgumentException("the maximum allele index provided cannot be negative: " + maximumAllele);
    }

    private static void checkOffsetTableCapacity(final int[][] offsetTable, final int maximumPloidy, final int maximumAllele) {
        if (offsetTable == null)
            throw new IllegalArgumentException("the allele first genotype offset table provided cannot be null");
        if (offsetTable.length <= maximumPloidy )
            throw new IllegalArgumentException("the allele first genotype offset table provided does not have enough " +
                    "capacity for requested maximum ploidy: " + maximumPloidy);
        if (offsetTable[0].length < maximumAllele)
            throw new IllegalArgumentException("the allele first genotype offset table provided does not have enough " +
                    "capacity for requested maximum allele index: " + maximumAllele);
    }


    /**
     * Returns the number of possible genotypes given the ploidy and number of different alleles.
     * @param ploidy the requested ploidy.
     * @param alleleCount the requested number of alleles.
     *
     * @throws IllegalArgumentException if {@code ploidy} or {@code alleleCount} is negative.
     *
     * @return 0 or greater.
     */
    public final static int genotypeCount(final int ploidy, final int alleleCount) {
        checkPloidyAndMaximumAllele(ploidy, alleleCount);
        if (ploidy > maximumPloidy || alleleCount > maximumAllele)
            ensureCapacity(alleleCount,ploidy);
        return alleleFirstGenotypeOffsetByPloidy[ploidy][alleleCount];
    }
}