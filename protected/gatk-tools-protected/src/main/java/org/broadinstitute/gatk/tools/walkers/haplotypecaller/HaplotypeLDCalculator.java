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
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.genotyper.AlleleList;
import org.broadinstitute.gatk.utils.genotyper.AlleleListUtils;
import org.broadinstitute.gatk.utils.genotyper.SampleListUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;

import java.util.*;

/**
 * Computes the likelihood based probability that haplotypes for first and second variant contexts
 * only appear in their fully linked form (x11 and x22) given a set of haplotypes where they might occur
 * and read likelihoods per sample
 *
 * User: depristo
 * Date: 3/29/13
 * Time: 9:23 AM
 */
public class HaplotypeLDCalculator {
    private final List<Haplotype> haplotypes;
    private final ReadLikelihoods<Haplotype> readLikelihoods;
    private List<Map<Haplotype, Double>> haplotypeLikelihoodsPerSample = null;

    // linear contigency table with table[0] == [0][0], table[1] = [0][1], table[2] = [1][0], table[3] = [1][1]
    private final double[] table = new double[4];

    /**
     * For testing
     */
    @SuppressWarnings("unchecked")
    protected HaplotypeLDCalculator() {
        haplotypes = Collections.emptyList();
        final AlleleList<Haplotype> alleleList = AlleleListUtils.emptyList();
        readLikelihoods = new ReadLikelihoods<>(SampleListUtils.emptyList(),
                alleleList, Collections.EMPTY_MAP);
    }

    public HaplotypeLDCalculator(final List<Haplotype> haplotypes, final ReadLikelihoods<Haplotype> haplotypeReadMap) {
        this.haplotypes = haplotypes;
        this.readLikelihoods = haplotypeReadMap;
    }

    /**
     * Construct the cached list of summed haplotype likelihoods per sample if it
     * hasn't already been computed.  This data structure is lazy created but only
     * needs to be made once when we make 1 merge decision as the data doesn't change
     * no matter how many calls to computeProbOfBeingPhased
     */
    private void buildHaplotypeLikelihoodsPerSampleIfNecessary() {
        if ( haplotypeLikelihoodsPerSample == null ) {
            // do the lazy computation
            final Set<String> samples = new LinkedHashSet<>(readLikelihoods.samples());
            haplotypeLikelihoodsPerSample = new LinkedList<>();
            for( final String sample : samples ) {
                final Map<Haplotype, Double> map = new HashMap<>(haplotypes.size());
                for( final Haplotype h : haplotypes ) {
                    // count up the co-occurrences of the events for the R^2 calculation
                    final double haplotypeLikelihood = PairHMMLikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(sample, readLikelihoods, Collections.singletonList(h), false)[0][0];
                    map.put(h, haplotypeLikelihood);
                }
                haplotypeLikelihoodsPerSample.add(map);
            }
        }
    }

    /**
     * Compute the likelihood based probability that that haplotypes for first and second are only x11 and x22
     *
     * As opposed to the hypothesis that all four haplotypes (x11, x12, x21, and x22) exist in the population
     *
     * @param first a non-null VariantContext
     * @param second a non-null VariantContext
     * @return the probability that only x11 and x22 exist among the samples
     */
    protected double computeProbOfBeingPhased(final VariantContext first, final VariantContext second) {
        buildHaplotypeLikelihoodsPerSampleIfNecessary();

        Arrays.fill(table, Double.NEGATIVE_INFINITY);

        for ( final Map<Haplotype, Double> entry : haplotypeLikelihoodsPerSample ) {
            for ( final Map.Entry<Haplotype, Double> haplotypeLikelihood : entry.entrySet() ) {
                final Haplotype h = haplotypeLikelihood.getKey();
                // count up the co-occurrences of the events for the R^2 calculation
                final VariantContext thisHapVC = h.getEventMap().get(first.getStart());
                final VariantContext nextHapVC = h.getEventMap().get(second.getStart()); // TODO -- add function to take a VC
                final int i = thisHapVC == null ? 0 : 1;
                final int j = nextHapVC == null ? 0 : 1;
                final int index = 2 * i + j;
                table[index] = MathUtils.approximateLog10SumLog10(table[index], haplotypeLikelihood.getValue());
            }
        }

        return pPhased(table);
    }

    /**
     * Compute probability that two variants are in phase with each other and that no
     * compound hets exist in the population.
     *
     * Implemented as a likelihood ratio test of the hypothesis:
     *
     * x11 and x22 are the only haplotypes in the populations
     *
     * vs.
     *
     * all four haplotype combinations (x11, x12, x21, and x22) all exist in the population.
     *
     * Now, since we have to have both variants in the population, we exclude the x11 & x11 state.  So the
     * p of having just x11 and x22 is P(x11 & x22) + p(x22 & x22).
     *
     * Alternatively, we might have any configuration that gives us both 1 and 2 alts, which are:
     *
     * - P(x11 & x12 & x21) -- we have hom-ref and both hets
     * - P(x22 & x12 & x21) -- we have hom-alt and both hets
     * - P(x22 & x12) -- one haplotype is 22 and the other is het 12
     * - P(x22 & x21) -- one haplotype is 22 and the other is het 21
     *
     * The probability is just p11_22 / (p11_22 + p hets)
     *
     * @param table linear contigency table with table[0] == [0][0], table[1] = [0][1], table[2] = [1][0], table[3] = [1][1]
     *      doesn't have to be normalized as this function does the normalization internally
     * @return the real space probability that the data is phased
     */
    @Requires("table.length == 4")
    protected double pPhased( double[] table ) {
        final double[] normTable = MathUtils.normalizeFromLog10(table, true);

        final double x11 = normTable[0], x12 = normTable[1], x21 = normTable[2], x22 = normTable[3];

        // probability that we are only x11 && x22
        final double p11_22 = MathUtils.approximateLog10SumLog10(x11 + x22, x22 + x22);

        // probability of having any of the other pairs
        final double p11_12_21 = MathUtils.approximateLog10SumLog10(x11 + x12, x11 + x21, x12 + x21);
        final double p22_12_21 = MathUtils.approximateLog10SumLog10(x22 + x12, x22 + x21, x12 + x21);
        final double p22_12 = x22 + x12;
        final double p22_21 = x22 + x21;
        final double pOthers = MathUtils.approximateLog10SumLog10(new double[]{p11_12_21, p22_12_21, p22_12, p22_21});

        // probability of being phases is the ratio of p11_22 / pOthers which in log space is just a substraction
        final double log10phased = p11_22 - (MathUtils.approximateLog10SumLog10(p11_22, pOthers));

        return Math.pow(10.0, log10phased);
    }

    protected double pPhasedTest( final double x11, final double x12, final double x21, final double x22 ) {
        return pPhased(new double[]{x11, x12, x21, x22});
    }
}
