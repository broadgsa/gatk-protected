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

import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

public abstract class DiploidExactAFCalculator extends ExactAFCalculator {

    private static final double LOG10_OF_2 = MathUtils.Log10Cache.get(2);

    public DiploidExactAFCalculator() {
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                               final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        final int numAlternateAlleles = vc.getNAlleles() - 1;

        final ArrayList<double[]> genotypeLikelihoods = getGLs(vc.getGenotypes(), true, vc.hasAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE));
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        // queue of AC conformations to process
        final LinkedList<ExactACset> ACqueue = new LinkedList<>();

        // mapping of ExactACset indexes to the objects
        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<>(numChr+1);

        // add AC=0 to the queue
        final int[] zeroCounts = new int[numAlternateAlleles];
        ExactACset zeroSet = new ExactACset(numSamples+1, new ExactACcounts(zeroCounts));
        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.getACcounts(), zeroSet);

        while ( !ACqueue.isEmpty() ) {
            stateTracker.incNEvaluations(); // keep track of the number of evaluations

            // compute log10Likelihoods
            final ExactACset set = ACqueue.remove();

            calculateAlleleCountConformation(set, genotypeLikelihoods, numChr, ACqueue,
                    indexesToACset, log10AlleleFrequencyPriors,stateTracker);

            // clean up memory
            indexesToACset.remove(set.getACcounts());
            //if ( DEBUG )
            //    System.out.printf(" *** removing used set=%s%n", set.ACcounts);
        }

        return getResultFromFinalState(vc, log10AlleleFrequencyPriors, stateTracker);
    }


    @Override
    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        return GATKVariantContextUtils.subsetDiploidAlleles(vc, allelesToUse, GATKVariantContextUtils.GenotypeAssignmentMethod.SET_TO_NO_CALL);
    }

    @Override
    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes(), true);
        for ( final double[] likelihoods : GLs ) {
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PL_INDEX_OF_HOM_REF ) {
                final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindexOfBestGL);
                final int alleleLikelihoodIndex1 = alleles.alleleIndex1 - 1;
                final int alleleLikelihoodIndex2 = alleles.alleleIndex2 - 1;
                if ( alleles.alleleIndex1 != 0 )
                    likelihoodSums[alleleLikelihoodIndex1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
                // don't double-count it
                if ( alleles.alleleIndex2 != 0 && alleles.alleleIndex2 != alleles.alleleIndex1 )
                    likelihoodSums[alleleLikelihoodIndex2].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
            }
        }
    }

    private static final class DependentSet {
        public final int[] ACcounts;
        public final int PLindex;

        public DependentSet(final int[] ACcounts, final int PLindex) {
            this.ACcounts = ACcounts;
            this.PLindex = PLindex;
        }
    }

    private double calculateAlleleCountConformation(final ExactACset set,
                                                    final ArrayList<double[]> genotypeLikelihoods,
                                                    final int numChr,
                                                    final LinkedList<ExactACset> ACqueue,
                                                    final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                                    final double[] log10AlleleFrequencyPriors,
                                                    final StateTracker stateTracker) {

        //if ( DEBUG )
        //    System.out.printf(" *** computing LofK for set=%s%n", set.ACcounts);

        // compute the log10Likelihoods
        computeLofK(set, genotypeLikelihoods, log10AlleleFrequencyPriors, stateTracker);

        final double log10LofK = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];

        // can we abort early because the log10Likelihoods are so small?
        if ( stateTracker.abort(log10LofK, set.getACcounts(), true, false) ) {
            //if ( DEBUG )
            //    System.out.printf(" *** breaking early set=%s log10L=%.2f maxLog10L=%.2f%n", set.ACcounts, log10LofK, maxLog10L);
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        final int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;

        final int numAltAlleles = set.getACcounts().getCounts().length;

        // add conformations for the k+1 case
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // to get to this conformation, a sample would need to be AB (remember that ref=0)
            final int PLindex = GenotypeLikelihoods.calculatePLindex(0, allele+1);
            updateACset(ACcountsClone, numChr, set, PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            final ArrayList<DependentSet> differentAlleles = new ArrayList<>(numAltAlleles * numAltAlleles);
            final ArrayList<DependentSet> sameAlleles = new ArrayList<>(numAltAlleles);

            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    final int[] ACcountsClone = set.getACcounts().getCounts().clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;

                    // to get to this conformation, a sample would need to be BB or BC (remember that ref=0, so add one to the index)
                    final int PLindex = GenotypeLikelihoods.calculatePLindex(allele_i+1, allele_j+1);
                    if ( allele_i == allele_j )
                        sameAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    else
                        differentAlleles.add(new DependentSet(ACcountsClone, PLindex));
                }
            }

            // IMPORTANT: we must first add the cases where the 2 new alleles are different so that the queue maintains its ordering
            for ( DependentSet dependent : differentAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            for ( DependentSet dependent : sameAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        return log10LofK;
    }

    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also pushes its value to the given callingSetIndex.
    private void updateACset(final int[] newSetCounts,
                             final int numChr,
                             final ExactACset dependentSet,
                             final int PLsetIndex,
                             final Queue<ExactACset> ACqueue,
                             final HashMap<ExactACcounts, ExactACset> indexesToACset,
                             final ArrayList<double[]> genotypeLikelihoods) {
        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            ExactACset set = new ExactACset(numChr/2 +1, index);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // push data from the dependency to the new set
        //if ( DEBUG )
        //    System.out.println(" *** pushing data from " + index + " to " + dependencySet.ACcounts);
        pushData(indexesToACset.get(index), dependentSet, PLsetIndex, genotypeLikelihoods);
    }

    private void computeLofK(final ExactACset set,
                             final ArrayList<double[]> genotypeLikelihoods,
                             final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {

        set.getLog10Likelihoods()[0] = 0.0; // the zero case
        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            for ( int j = 1; j < set.getLog10Likelihoods().length; j++ )
                set.getLog10Likelihoods()[j] = set.getLog10Likelihoods()[j-1] + genotypeLikelihoods.get(j)[HOM_REF_INDEX];

            final double log10Lof0 = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];
            stateTracker.setLog10LikelihoodOfAFzero(log10Lof0);
            stateTracker.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return;
        }

        // if we got here, then k > 0 for at least one k.
        // the non-AA possible conformations were already dealt with by pushes from dependent sets;
        // now deal with the AA case (which depends on previous cells in this column) and then update the L(j,k) value
        for ( int j = 1; j < set.getLog10Likelihoods().length; j++ ) {

            if ( totalK < 2*j-1 ) {
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue = MathUtils.Log10Cache.get(2*j-totalK) + MathUtils.Log10Cache.get(2*j-totalK-1) + set.getLog10Likelihoods()[j-1] + gl[HOM_REF_INDEX];
                set.getLog10Likelihoods()[j] = MathUtils.approximateLog10SumLog10(set.getLog10Likelihoods()[j], conformationValue);
            }

            final double logDenominator = MathUtils.Log10Cache.get(2*j) + MathUtils.Log10Cache.get(2*j-1);
            set.getLog10Likelihoods()[j] = set.getLog10Likelihoods()[j] - logDenominator;
        }

        double log10LofK = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];

        // update the MLE if necessary
        stateTracker.updateMLEifNeeded(log10LofK, set.getACcounts().getCounts());

        // apply the priors over each alternate allele
        for ( final int ACcount : set.getACcounts().getCounts() ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }

        stateTracker.updateMAPifNeeded(log10LofK, set.getACcounts().getCounts());
    }

    private void pushData(final ExactACset targetSet,
                          final ExactACset dependentSet,
                          final int PLsetIndex,
                          final ArrayList<double[]> genotypeLikelihoods) {
        final int totalK = targetSet.getACsum();
        final double[] targetLog10Likelihoods = targetSet.getLog10Likelihoods();
        final double[] dependentLog10Likelihoods = dependentSet.getLog10Likelihoods();
        final int[] targetACcounts = targetSet.getACcounts().getCounts();

        // skip impossible conformations, namely those for which there aren't enough possible chromosomes (2*j in the if loop below)
        // to fill up the totalK alternate alleles needed; we do want to ensure that there's at least one sample included (hence the Math.max)
        final int firstIndex = Math.max(1, (totalK + 1) / 2);

        // find the 2 alleles that are represented by this PL index
        final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLsetIndex);
        // if neither allele is reference then the coefficient is constant throughout this invocation
        final Double constantCoefficient = (alleles.alleleIndex1 == 0) ? null : determineCoefficient(alleles, 0, targetACcounts, totalK);

        for ( int j = firstIndex; j < targetLog10Likelihoods.length; j++ ) {
            final double[] gl = genotypeLikelihoods.get(j);
            final double coefficient = (constantCoefficient != null) ? constantCoefficient : determineCoefficient(alleles, j, targetACcounts, totalK);

            final double conformationValue = coefficient + dependentLog10Likelihoods[j-1] + gl[PLsetIndex];
            targetLog10Likelihoods[j] = MathUtils.approximateLog10SumLog10(targetLog10Likelihoods[j], conformationValue);
        }
    }

    private double determineCoefficient(final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles, final int j, final int[] ACcounts, final int totalK) {
        // the closed form representation generalized for multiple alleles is as follows:
        // AA: (2j - totalK) * (2j - totalK - 1)
        // AB: 2k_b * (2j - totalK)
        // AC: 2k_c * (2j - totalK)
        // BB: k_b * (k_b - 1)
        // BC: 2 * k_b * k_c
        // CC: k_c * (k_c - 1)

        // *** note that throughout this method we subtract one from the alleleIndex because ACcounts ***
        // *** doesn't consider the reference allele whereas the GenotypeLikelihoods PL cache does.   ***

        // the AX het case
        if ( alleles.alleleIndex1 == 0 )
            return MathUtils.Log10Cache.get(2*ACcounts[alleles.alleleIndex2-1]) + MathUtils.Log10Cache.get(2*j-totalK);

        final int k_i = ACcounts[alleles.alleleIndex1-1];

        // the hom var case (e.g. BB, CC, DD)
        final double coeff;
        if ( alleles.alleleIndex1 == alleles.alleleIndex2 ) {
            coeff = MathUtils.Log10Cache.get(k_i) + MathUtils.Log10Cache.get(k_i - 1);
        }
        // the het non-ref case (e.g. BC, BD, CD)
        else {
            final int k_j = ACcounts[alleles.alleleIndex2-1];
            coeff = LOG10_OF_2 + MathUtils.Log10Cache.get(k_i) + MathUtils.Log10Cache.get(k_j);
        }

        return coeff;
    }

    @Override
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                          final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {
        if (defaultPloidy != 2)
            throw new IllegalArgumentException("cannot support ploidy different than 2 and the default ploidy is " + defaultPloidy);
        return allelesToUse.size() == 1
                ? GATKVariantContextUtils.subsetToRefOnly(vc, 2)
                : GATKVariantContextUtils.subsetDiploidAlleles(vc, allelesToUse,
                     assignGenotypes ? GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN : GATKVariantContextUtils.GenotypeAssignmentMethod.SET_TO_NO_CALL);
    }
}
