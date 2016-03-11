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

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.fragments.FragmentUtils;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;

import java.util.List;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Stable, error checking version of the Bayesian genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
public class DiploidSNPGenotypeLikelihoods implements Cloneable {

    public final static double DEFAULT_PCR_ERROR_RATE = FragmentUtils.DEFAULT_PCR_ERROR_RATE;

    protected final static int FIXED_PLOIDY = 2;
    protected final static int MAX_PLOIDY = FIXED_PLOIDY + 1;
    protected final static double ploidyAdjustment = log10(FIXED_PLOIDY);
    protected final static double log10_3 = log10(3.0);

    protected boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelihoods object
    //
    protected double[] log10Likelihoods = null;

    // TODO: don't calculate this each time through
    protected double log10_PCR_error_3;
    protected double log10_1_minus_PCR_error;

    /**
     * Create a new GenotypeLikelhoods object with given PCR error rate for each diploid genotype
     *
     * @param PCR_error_rate  the PCR error rate
     */
    public DiploidSNPGenotypeLikelihoods(double PCR_error_rate) {
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - PCR_error_rate);
        setToZero();
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        DiploidSNPGenotypeLikelihoods c = (DiploidSNPGenotypeLikelihoods)super.clone();
        c.log10Likelihoods = log10Likelihoods.clone();
        return c;
    }

    protected void setToZero() {
        log10Likelihoods = genotypeZeros.clone();                 // likelihoods are all zeros
    }

    /**
     * Returns an array of log10 likelihoods for each genotype, indexed by DiploidGenotype.ordinal values()
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }

    // -------------------------------------------------------------------------------------
    //
    // add() routines.  These are the workhorse routines for calculating the overall genotype
    // likelihoods given observed bases and reads.  Includes high-level operators all the
    // way down to single base and qual functions.
    //
    // -------------------------------------------------------------------------------------

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @param minBaseQual               the minimum base quality at which to consider a base valid
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        int n = 0;

        // for each fragment, add to the likelihoods
        FragmentCollection<PileupElement> fpile = pileup.toFragments();

        for ( PileupElement p : fpile.getSingletonReads() )
            n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        for ( List<PileupElement> overlappingPair : fpile.getOverlappingPairs() )
            n += add(overlappingPair, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        return n;
    }

    public int add(PileupElement elt, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        byte obsBase = elt.getBase();
        byte qual = qualToUse(elt, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        if ( qual == 0 )
            return 0;

        return add(obsBase, qual, (byte)0, (byte)0, 1);
    }

    public int add(List<PileupElement> overlappingPair, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        final PileupElement p1 = overlappingPair.get(0);
        final PileupElement p2 = overlappingPair.get(1);

        final byte observedBase1 = p1.getBase();
        final byte qualityScore1 = qualToUse(p1, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        final byte observedBase2 = p2.getBase();
        final byte qualityScore2 = qualToUse(p2, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        if ( qualityScore1 == 0 ) {
            if ( qualityScore2 == 0 ) // abort early if we didn't see any good bases
                return 0;
            else {
                return add(observedBase2, qualityScore2, (byte)0, (byte)0);
            }
        } else {
            return add(observedBase1, qualityScore1, observedBase2, qualityScore2);
        }
    }

    /**
     *
     * @param obsBase1    first observed base
     * @param qual1       base qual of first observed base
     * @param obsBase2    second observed base
     * @param qual2       base qual of second observed base; can be 0, indicating no second base was observed for this fragment
     * @param nObs        the number of times this quad of values was seen.  Generally 1, but reduced reads can have nObs > 1 for synthetic reads
     * @return 0 if the base is bad, 1 otherwise
     */
    private int add(byte obsBase1, byte qual1, byte obsBase2, byte qual2, int nObs) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.

        // Just look up the cached result if it's available, or compute and store it
        DiploidSNPGenotypeLikelihoods gl;
        if ( ! inCache(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY) ) {
            gl = calculateCachedGenotypeLikelihoods(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY);
        } else {
            gl = getCachedGenotypeLikelihoods(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY);
        }

        // for bad bases, there are no likelihoods
        if ( gl == null )
            return 0;

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];
            log10Likelihoods[g.ordinal()] += likelihood * nObs;
        }

        return 1;
    }

    private int add(byte obsBase1, byte qual1, byte obsBase2, byte qual2) {
        return add(obsBase1, qual1, obsBase2, qual2, 1);
    }

    // -------------------------------------------------------------------------------------
    //
    // Dealing with the cache routines
    //
    // -------------------------------------------------------------------------------------

    static DiploidSNPGenotypeLikelihoods[][][][][] CACHE = new DiploidSNPGenotypeLikelihoods[BaseUtils.BASES.length][QualityUtils.MAX_SAM_QUAL_SCORE +1][BaseUtils.BASES.length+1][QualityUtils.MAX_SAM_QUAL_SCORE +1][MAX_PLOIDY];

    protected boolean inCache(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        return getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy) != null;
    }

    protected DiploidSNPGenotypeLikelihoods getCachedGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        DiploidSNPGenotypeLikelihoods gl = getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy);
        if ( gl == null )
            throw new RuntimeException(String.format("BUG: trying to fetch an unset cached genotype likelihood at base1=%c, qual1=%d, base2=%c, qual2=%d, ploidy=%d",
                    observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy));
        return gl;
    }

    protected DiploidSNPGenotypeLikelihoods calculateCachedGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        DiploidSNPGenotypeLikelihoods gl = calculateGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);
        setCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy, gl);
        return gl;
    }

    protected void setCache( DiploidSNPGenotypeLikelihoods[][][][][] cache,
                             byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy,
                             DiploidSNPGenotypeLikelihoods val ) {
        int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
        int j = qualityScore1;
        int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
        int l = qualityScore2;
        int m = ploidy;

        cache[i][j][k][l][m] = val;
    }

    protected DiploidSNPGenotypeLikelihoods getCache(DiploidSNPGenotypeLikelihoods[][][][][] cache,
                                            byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
        int j = qualityScore1;
        int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
        int l = qualityScore2;
        int m = ploidy;
        return cache[i][j][k][l][m];
    }

    protected DiploidSNPGenotypeLikelihoods calculateGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
        double[] log10FourBaseLikelihoods = computeLog10Likelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);

        try {

            DiploidSNPGenotypeLikelihoods gl = (DiploidSNPGenotypeLikelihoods)this.clone();
            gl.setToZero();

            // we need to adjust for ploidy.  We take the raw p(obs | chrom) / ploidy, which is -log10(ploidy) in log space
            for ( DiploidGenotype g : DiploidGenotype.values() ) {

                // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
                double p_base = 0.0;
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - ploidyAdjustment);
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - ploidyAdjustment);

                final double likelihood = log10(p_base);
                gl.log10Likelihoods[g.ordinal()] += likelihood;
            }

            if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
            for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
            }

            return gl;

         } catch ( CloneNotSupportedException e ) {
             throw new RuntimeException(e);
         }
    }

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase1  the base observed on the 1st read of the fragment
     * @param qualityScore1  the qual of the base on the 1st read of the fragment, or zero if NA
     * @param observedBase2  the base observed on the 2nd read of the fragment
     * @param qualityScore2  the qual of the base on the 2nd read of the fragment, or zero if NA
     * @return likelihoods for this observation or null if the base was not considered good enough to add to the likelihoods (Q0 or 'N', for example)
     */
    protected double[] computeLog10Likelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
        double[] log10FourBaseLikelihoods = baseZeros.clone();

        for ( byte trueBase : BaseUtils.BASES ) {
            double likelihood = 0.0;

            for ( byte fragmentBase : BaseUtils.BASES ) {
                double log10FragmentLikelihood = (trueBase == fragmentBase ? log10_1_minus_PCR_error : log10_PCR_error_3);
                if ( qualityScore1 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase1, fragmentBase, qualityScore1);
                }
                if ( qualityScore2 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase2, fragmentBase, qualityScore2);
                }

                //if ( VERBOSE ) {
                //    System.out.printf("  L(%c | b=%s, Q=%d) = %f / %f%n",
                //            observedBase, trueBase, qualityScore, pow(10,likelihood) * 100, likelihood);
                //}

                likelihood += pow(10, log10FragmentLikelihood);
            }

            log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(trueBase)] = log10(likelihood);
        }

        return log10FourBaseLikelihoods;
    }

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param qual         base quality
     * @return log10 likelihood
     */
    protected double log10PofObservingBaseGivenChromosome(byte observedBase, byte chromBase, byte qual) {

        double logP;

        if ( observedBase == chromBase ) {
            // the base is consistent with the chromosome -- it's 1 - e
            //logP = oneMinusData[qual];
            double e = pow(10, (qual / -10.0));
            logP = log10(1.0 - e);
        } else {
            // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
            logP = qual / -10.0 + (-log10_3);
        }

        //System.out.printf("%c %c %d => %f%n", observedBase, chromBase, qual, logP);
        return logP;
    }

    /**
     * Helper function that returns the phred-scaled base quality score we should use for calculating
     * likelihoods for a pileup element.  May return 0 to indicate that the observation is bad, and may
     * cap the quality score by the mapping quality of the read itself.
     *
     * @param p                           Pileup element
     * @param ignoreBadBases              Should we ignore bad bases?
     * @param capBaseQualsAtMappingQual   Should we cap the base qualities at the mapping quality of the read?
     * @param minBaseQual                 Minimum allowed base quality
     * @return the actual base quality to use
     */
    private static byte qualToUse(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        if ( ignoreBadBases && !BaseUtils.isRegularBase( p.getBase() ) )
            return 0;

        byte qual = p.getQual();

        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MisencodedBAM(p.getRead(), "we encountered an extremely high quality score (" + (int)qual + ")");
        if ( capBaseQualsAtMappingQual )
            qual = (byte) Math.min( 0xff & qual, p.getMappingQual());
        if ( (int)qual < minBaseQual )
            qual = (byte)0;

        return qual;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        double sum = 0;
        StringBuilder s = new StringBuilder();
        for (DiploidGenotype g : DiploidGenotype.values()) {
            s.append(String.format("%s %.10f ", g, log10Likelihoods[g.ordinal()]));
			sum += Math.pow(10,log10Likelihoods[g.ordinal()]);
        }
		s.append(String.format(" %f", sum));
        return s.toString();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Validation routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        try {
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                String bad = null;

                int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(log10Likelihoods[i]) || ! MathUtils.isNegativeOrZero(log10Likelihoods[i]) ) {
                    bad = String.format("Likelihood %f is badly formed", log10Likelihoods[i]);
                }

                if ( bad != null ) {
                    throw new IllegalStateException(String.format("At %s: %s", g.toString(), bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

    //
    // Constant static data
    //
    private final static double[] genotypeZeros = new double[DiploidGenotype.values().length];
    private final static double[] baseZeros = new double[BaseUtils.BASES.length];

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            genotypeZeros[g.ordinal()] = 0.0;
        }
        for ( byte base : BaseUtils.BASES ) {
            baseZeros[BaseUtils.simpleBaseToBaseIndex(base)] = 0.0;
        }
    }
}
