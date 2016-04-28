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

package org.broadinstitute.gatk.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Uses the Exact calculation of Heng Li
 */
abstract class ExactAFCalculator extends AFCalculator {

    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first

    // useful so that we don't keep printing out the same warning messages
    protected static boolean printedMaxAltAllelesWarning = false;

    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_SUM_COMPARATOR = new Comparator<LikelihoodSum>() {

        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            return - Double.compare(o1.sum,o2.sum);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first but make sure that
     * NON_REF alleles are place are last.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR = new Comparator<LikelihoodSum>() {
        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            if (o1.allele == GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)
                return 1;
            else if (o2.allele == GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)
                return -1;
            else
                return o1.compareTo(o2);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with lower alternative allele index are first regardless of
     * the likelihood sum.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_INDEX_COMPARATOR = new Comparator<LikelihoodSum>() {
        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            return Integer.compare(o1.index, o2.index);
        }
    };

    protected ExactAFCalculator() {

    }

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public final Allele allele;
        public final int index;

        public LikelihoodSum(final Allele allele, final int index) { this.allele = allele; this.index = index; }

        public int compareTo(LikelihoodSum other) {
            final double diff = sum - other.sum;
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        return getGLs(GLs, includeDummy, false);
    }

    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @param keepUninformative Don't filter out uninformative genotype likelihoods (i.e. all log likelihoods near 0)
     *                          This is useful for VariantContexts with a NON_REF allele
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy, final boolean keepUninformative) {
        final ArrayList<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);

        if ( includeDummy ) genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();

                if ( MathUtils.sum(gls) < GATKVariantContextUtils.SUM_GL_THRESH_NOCALL || keepUninformative )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }

    @Override
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        // don't try to genotype too many alternate alleles
        final List<Allele> inputAltAlleles = vc.getAlternateAlleles();
        final List<Allele> outputAltAlleles = reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles);

        // only if output allele has reduced from the input alt allele set size we should care.
        final int altAlleleReduction = inputAltAlleles.size() - outputAltAlleles.size();

        if (altAlleleReduction == 0)
            return vc;

        final String message = String.format("This tool is currently set to genotype at most %d " +
                        "alternate alleles in a given context, but the context at %s: %d has %d " +
                        "alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument",
                maximumAlternativeAlleles, vc.getContig(), vc.getStart(), vc.getAlternateAlleles().size());

        if ( !printedMaxAltAllelesWarning ) {
            printedMaxAltAllelesWarning = true;
            logger.warn(message + ". Unless the DEBUG logging level is used, this warning message is output just once per run and further warnings are suppressed.");
        } else {
            logger.debug(message);
        }

        final List<Allele> alleles = new ArrayList<>(maximumAlternativeAlleles + 1);
        alleles.add(vc.getReference());
        alleles.addAll(reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles));
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.alleles(alleles);
        builder.genotypes(reduceScopeGenotypes(vc, defaultPloidy, alleles));
        if (altAlleleReduction < 0)
            throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!: " + - altAlleleReduction + " " + vc + " " + builder.make());
        return builder.make();
    }

    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative allele to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {

        // Look  for the <NON_REF> allele to exclude it from the pruning if present.
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();

        final int nonRefAltAlleleIndex = GATKVariantContextUtils.indexOfAltAllele(vc,
                GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final boolean nonRefAltAllelePresent = nonRefAltAlleleIndex >= 0;

        // <NON_REF> should not be considered in the downsizing, so we need to count it out when
        // considering if alt. allele downsizing is required.
        final int numProperOriginalAltAlleles = numOriginalAltAlleles - (nonRefAltAllelePresent ? 1 : 0);

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= numProperOriginalAltAlleles)
            return vc.getAlternateAlleles();

        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ ) {
            final Allele allele = vc.getAlternateAllele(i);
            likelihoodSums[i] = new LikelihoodSum(allele,i);
        }

        // Calculate the allele likelihood sums.
        reduceScopeCalculateLikelihoodSums(vc, defaultPloidy, likelihoodSums);

        // sort them by probability mass and choose the best ones
        // Make sure that the <NON_REF> allele is last if present.
        Collections.sort(Arrays.asList(likelihoodSums), nonRefAltAllelePresent ? LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR : LIKELIHOOD_SUM_COMPARATOR);

        // We need to return the best likelihood alleles in the original alternative allele index order.
        // This heap will keep track of that index order.
        final PriorityQueue<LikelihoodSum> mostLikelyAllelesHeapByIndex = new PriorityQueue<>(numOriginalAltAlleles, LIKELIHOOD_INDEX_COMPARATOR);

        for ( int i = 0; i < numAllelesToChoose; i++ )
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[i]);

        // guaranteed no to have been added at this point thanks for checking on whether reduction was
        // needed in the first place.
        if (nonRefAltAllelePresent)
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[nonRefAltAlleleIndex]);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<>(numAllelesToChoose);

        while (!mostLikelyAllelesHeapByIndex.isEmpty())
            orderedBestAlleles.add(mostLikelyAllelesHeapByIndex.remove().allele);

        return orderedBestAlleles;
    }

    protected static final int PL_INDEX_OF_HOM_REF = 0;

    /**
     * Update the likelihood sums with using the variant context genotype likelihoods.
     * @param vc source variant context.
     * @param likelihoodSums where to update the likelihood sums.
     */
    protected abstract void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums);

    /**
     * Transforms the genotypes of the variant context according to the new subset of possible alleles.
     *
     * @param vc original variant-context.
     * @param allelesToUse possible alleles.
     * @return never {@code null}, the new set of genotype calls for the reduced scope.
     */
    protected abstract GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse);
}