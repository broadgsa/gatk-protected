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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFConstants;

import java.util.*;

public class PosteriorLikelihoodsUtils {

    public static VariantContext calculatePosteriorGLs(final VariantContext vc1,
                                                       final Collection<VariantContext> resources,
                                                       final int numRefSamplesFromMissingResources,
                                                       final double globalFrequencyPriorDirichlet,
                                                       final boolean useInputSamples,
                                                       final boolean useEM,
                                                       final boolean useAC) {
        if ( useEM )
            throw new IllegalArgumentException("EM loop for posterior GLs not yet implemented");

        final Map<Allele,Integer> totalAlleleCounts = new HashMap<>();
        for ( final VariantContext resource : resources ) {
            addAlleleCounts(totalAlleleCounts,resource,useAC);
        }

        if ( useInputSamples ) {
            addAlleleCounts(totalAlleleCounts,vc1,useAC);
        }

        totalAlleleCounts.put(vc1.getReference(),totalAlleleCounts.get(vc1.getReference())+numRefSamplesFromMissingResources);

        // now extract the counts of the alleles present within vc1, and in order
        final double[] alleleCounts = new double[vc1.getNAlleles()];
        int alleleIndex = 0;
        for ( final Allele allele : vc1.getAlleles() ) {

            alleleCounts[alleleIndex++] = globalFrequencyPriorDirichlet + ( totalAlleleCounts.containsKey(allele) ?
                    totalAlleleCounts.get(allele) : 0 );
        }

        final List<double[]> likelihoods = new ArrayList<>(vc1.getNSamples());
        for ( final Genotype genotype : vc1.getGenotypes() ) {
            likelihoods.add(genotype.hasLikelihoods() ? genotype.getLikelihoods().getAsVector() : null );
        }

        final List<double[]> posteriors = calculatePosteriorGLs(likelihoods,alleleCounts,vc1.getMaxPloidy(2));

        final GenotypesContext newContext = GenotypesContext.create();
        for ( int genoIdx = 0; genoIdx < vc1.getNSamples(); genoIdx ++ ) {
            final GenotypeBuilder builder = new GenotypeBuilder(vc1.getGenotype(genoIdx));
            if ( posteriors.get(genoIdx) != null ) {
                GATKVariantContextUtils.updateGenotypeAfterSubsetting(vc1.getAlleles(), builder,
                        GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, posteriors.get(genoIdx), vc1.getAlleles());
                builder.attribute(VCFConstants.GENOTYPE_POSTERIORS_KEY,
                        Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(posteriors.get(genoIdx)).getAsPLs()));

            }
            newContext.add(builder.make());
        }

        final List<Integer> priors = Utils.listFromPrimitives(
                GenotypeLikelihoods.fromLog10Likelihoods(getDirichletPrior(alleleCounts, vc1.getMaxPloidy(2))).getAsPLs());

        return new VariantContextBuilder(vc1).genotypes(newContext).attribute("PG",priors).make();
    }

    /**
     * Given genotype likelihoods and known allele counts, calculate the posterior likelihoods
     * over the genotype states
     * @param genotypeLikelihoods - the genotype likelihoods for the individual
     * @param knownAlleleCountsByAllele - the known allele counts in the population. For AC=2 AN=12 site, this is {10,2}
     * @param ploidy - the ploidy to assume
     * @return - the posterior genotype likelihoods
     */
    protected static List<double[]> calculatePosteriorGLs(final List<double[]> genotypeLikelihoods,
                                                       final double[] knownAlleleCountsByAllele,
                                                       final int ploidy) {
        if ( ploidy != 2 ) {
            throw new IllegalStateException("Genotype posteriors not yet implemented for ploidy != 2");
        }

        final double[] genotypePriorByAllele = getDirichletPrior(knownAlleleCountsByAllele,ploidy);
        final List<double[]> posteriors = new ArrayList<>(genotypeLikelihoods.size());
        for ( final double[] likelihoods : genotypeLikelihoods ) {
            double[] posteriorLikelihoods = null;

            if ( likelihoods != null ) {
                if ( likelihoods.length != genotypePriorByAllele.length ) {
                    throw new IllegalStateException(String.format("Likelihoods not of correct size: expected %d, observed %d",
                            knownAlleleCountsByAllele.length*(knownAlleleCountsByAllele.length+1)/2,likelihoods.length));
                }

                posteriorLikelihoods = new double[genotypePriorByAllele.length];
                for ( int genoIdx = 0; genoIdx < likelihoods.length; genoIdx ++ ) {
                    posteriorLikelihoods[genoIdx] = likelihoods[genoIdx] + genotypePriorByAllele[genoIdx];
                }

                posteriorLikelihoods = MathUtils.toLog10(MathUtils.normalizeFromLog10(posteriorLikelihoods));

            }

            posteriors.add(posteriorLikelihoods);
        }

        return posteriors;
    }

    // convenience function for a single genotypelikelihoods array. Just wraps.
    protected static double[] calculatePosteriorGLs(final double[] genotypeLikelihoods,
                                                 final double[] knownAlleleCountsByAllele,
                                                 final int ploidy) {
        return calculatePosteriorGLs(Arrays.asList(genotypeLikelihoods),knownAlleleCountsByAllele,ploidy).get(0);
    }


    /**
     * Given known allele counts (whether external, from the sample, or both), calculate the prior distribution
     * over genotype states. This assumes
     *   1) Random sampling of alleles (known counts are unbiased, and frequency estimate is Dirichlet)
     *   2) Genotype states are independent (Hardy-Weinberg)
     * These assumptions give rise to a Dirichlet-Multinomial distribution of genotype states as a prior
     * (the "number of trials" for the multinomial is simply the ploidy)
     * @param knownCountsByAllele - the known counts per allele. For an AC=2, AN=12 site this is {10,2}
     * @param ploidy - the number of chromosomes in the sample. For now restricted to 2.
     * @return - the Dirichlet-Multinomial distribution over genotype states
     */
    protected static double[] getDirichletPrior(final double[] knownCountsByAllele, final int ploidy) {
        if ( ploidy != 2 ) {
            throw new IllegalStateException("Genotype priors not yet implemented for ploidy != 2");
        }

        // multi-allelic format is
        // AA AB BB AC BC CC AD BD CD DD ...
        final double sumOfKnownCounts = MathUtils.sum(knownCountsByAllele);
        final double[] priors = new double[knownCountsByAllele.length*(knownCountsByAllele.length+1)/2];
        int priorIndex = 0;
        for ( int allele2 = 0; allele2 < knownCountsByAllele.length; allele2++ ) {
            for ( int allele1 = 0; allele1 <= allele2; allele1++) {
                final int[] counts = new int[knownCountsByAllele.length];
                counts[allele1] += 1;
                counts[allele2] += 1;
                priors[priorIndex++] = MathUtils.dirichletMultinomial(knownCountsByAllele,sumOfKnownCounts,counts,ploidy);
            }
        }

        return priors;
    }

    private static void addAlleleCounts(final Map<Allele,Integer> counts, final VariantContext context, final boolean useAC) {
        final int[] ac;
        if ( context.hasAttribute(VCFConstants.MLE_ALLELE_COUNT_KEY) && ! useAC ) {
            ac = extractInts(context.getAttribute(VCFConstants.MLE_ALLELE_COUNT_KEY));
        } else if ( context.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            ac = extractInts(context.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
        } else {
            ac = new int[context.getAlternateAlleles().size()];
            int idx = 0;
            for ( final Allele allele : context.getAlternateAlleles() ) {
                ac[idx++] = context.getCalledChrCount(allele);
            }
        }

        for ( final Allele allele : context.getAlleles() ) {
            final int count;
            if ( allele.isReference() ) {
                if ( context.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
                    count = context.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY,-1) - (int) MathUtils.sum(ac);
                } else {
                    count = context.getCalledChrCount() - (int) MathUtils.sum(ac);
                }
            } else {
                count = ac[context.getAlternateAlleles().indexOf(allele)];
            }
            if ( ! counts.containsKey(allele) ) {
                counts.put(allele,0);
            }
            counts.put(allele,count + counts.get(allele));
        }
    }

    public static int[] extractInts(final Object integerListContainingVCField) {
        List<Integer> mleList = null;
        if ( integerListContainingVCField instanceof List ) {
            if ( ((List) integerListContainingVCField).get(0) instanceof String ) {
                mleList = new ArrayList<>(((List) integerListContainingVCField).size());
                for ( Object s : ((List)integerListContainingVCField)) {
                    mleList.add(Integer.parseInt((String) s));
                }
            } else {
                mleList = (List<Integer>) integerListContainingVCField;
            }
        } else if ( integerListContainingVCField instanceof Integer ) {
            mleList = Arrays.asList((Integer) integerListContainingVCField);
        } else if ( integerListContainingVCField instanceof String ) {
            mleList = Arrays.asList(Integer.parseInt((String)integerListContainingVCField));
        }
        if ( mleList == null )
            throw new IllegalArgumentException(String.format("VCF does not have properly formatted "+
                    VCFConstants.MLE_ALLELE_COUNT_KEY+" or "+VCFConstants.ALLELE_COUNT_KEY));

        final int[] mle = new int[mleList.size()];

        if ( ! ( mleList.get(0) instanceof Integer ) ) {
            throw new IllegalStateException("BUG: The AC values should be an Integer, but was "+mleList.get(0).getClass().getCanonicalName());
        }

        for ( int idx = 0; idx < mle.length; idx++) {
            mle[idx] = mleList.get(idx);
        }

        return mle;
    }
}
