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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;

import java.util.*;

public class PosteriorLikelihoodsUtils {

    public static VariantContext calculatePosteriorGLs(final VariantContext vc1,
                                                       final Collection<VariantContext> resources,
                                                       final int numRefSamplesFromMissingResources,
                                                       final double globalFrequencyPriorDirichlet,
                                                       final boolean useInputSamples,
                                                       final boolean useAC,
                                                       final boolean useACoff) {

        final Map<Allele,Integer> totalAlleleCounts = new HashMap<>();
        boolean nonSNPprior = false;
        if (vc1 == null) throw new IllegalArgumentException("VariantContext vc1 is null");
        final boolean nonSNPeval = !vc1.isSNP();
        final double[] alleleCounts = new double[vc1.getNAlleles()];
        //only use discovered allele count if there are at least 10 samples
        final boolean useDiscoveredAC = !useACoff && vc1.getNSamples() >= 10;

        if(vc1.isSNP())
        {
            //store the allele counts for each allele in the variant priors
            for ( final VariantContext resource : resources ) {
                if( !resource.isSNP()) nonSNPprior = true;
                addAlleleCounts(totalAlleleCounts,resource,useAC);
            }

            //add the allele counts from the input samples (if applicable)
            if ( useInputSamples ) {
                addAlleleCounts(totalAlleleCounts,vc1,useAC);
            }

            //add zero allele counts for any reference alleles not seen in priors (if applicable)
            int existingRefCounts = 0;
            if (totalAlleleCounts.containsKey(vc1.getReference()))
                existingRefCounts += totalAlleleCounts.get(vc1.getReference());
            totalAlleleCounts.put(vc1.getReference(),existingRefCounts+numRefSamplesFromMissingResources);
        }

        // now extract the counts of the alleles present within vc1, and in order
        int alleleIndex = 0;
        for ( final Allele allele : vc1.getAlleles() ) {

            alleleCounts[alleleIndex++] = globalFrequencyPriorDirichlet + ( totalAlleleCounts.containsKey(allele) ?
                    totalAlleleCounts.get(allele) : 0 );
        }


        //parse the likelihoods for each sample's genotype
        final List<double[]> likelihoods = new ArrayList<>(vc1.getNSamples());
        for ( final Genotype genotype : vc1.getGenotypes() ) {
            if (!genotype.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY)){
                likelihoods.add(genotype.hasLikelihoods() ? genotype.getLikelihoods().getAsVector() : null );

            }
            else {
                Object PPfromVCF = genotype.getExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
                //parse the PPs into a vector of probabilities
                if (PPfromVCF instanceof String) {
                    final String PPstring = (String)PPfromVCF;
                    if (PPstring.charAt(0)=='.')  //samples not in trios will have PP tag like ".,.,." if family priors are applied
                        likelihoods.add(genotype.hasLikelihoods() ? genotype.getLikelihoods().getAsVector() : null );
                    else {
                        final String[] likelihoodsAsStringVector = PPstring.split(",");
                        double[] likelihoodsAsVector = new double[likelihoodsAsStringVector.length];
                        for ( int i = 0; i < likelihoodsAsStringVector.length; i++ ) {
                            likelihoodsAsVector[i] = Double.parseDouble(likelihoodsAsStringVector[i])/-10.0;
                        }
                        likelihoods.add(likelihoodsAsVector);
                    }
                }
                else {
                    int[] likelihoodsAsInts = extractInts(PPfromVCF);
                    double[] likelihoodsAsVector = new double[likelihoodsAsInts.length];
                    for ( int i = 0; i < likelihoodsAsInts.length; i++ ) {
                        likelihoodsAsVector[i] = likelihoodsAsInts[i]/-10.0;
                    }
                    likelihoods.add(likelihoodsAsVector);
                }
            }

        }

        //TODO: for now just use priors that are SNPs because indel priors will bias SNP calls
        final boolean useFlatPriors = nonSNPeval || nonSNPprior || (resources.isEmpty() && !useDiscoveredAC);

        final List<double[]> posteriors = calculatePosteriorGLs(likelihoods,alleleCounts,vc1.getMaxPloidy(2), useFlatPriors);

        final GenotypesContext newContext = GenotypesContext.create();
        for ( int genoIdx = 0; genoIdx < vc1.getNSamples(); genoIdx ++ ) {
            final GenotypeBuilder builder = new GenotypeBuilder(vc1.getGenotype(genoIdx));
            builder.phased(vc1.getGenotype(genoIdx).isPhased());
            if ( posteriors.get(genoIdx) != null ) {
                GATKVariantContextUtils.updateGenotypeAfterSubsetting(vc1.getAlleles(), builder,
                        GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, posteriors.get(genoIdx), vc1.getAlleles());
                builder.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY,
                        Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(posteriors.get(genoIdx)).getAsPLs()));
            }
            newContext.add(builder.make());
        }

        final List<Integer> priors = Utils.listFromPrimitives(
                GenotypeLikelihoods.fromLog10Likelihoods(getDirichletPrior(alleleCounts, vc1.getMaxPloidy(2),useFlatPriors)).getAsPLs());

        final VariantContextBuilder builder = new VariantContextBuilder(vc1).genotypes(newContext).attribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY, priors);
        // add in the AC, AF, and AN attributes
        VariantContextUtils.calculateChromosomeCounts(builder, true);
        return builder.make();
    }

    /**
     * Given genotype likelihoods and known allele counts, calculate the posterior likelihoods
     * over the genotype states
     * @param genotypeLikelihoods - the genotype likelihoods for the individual
     * @param knownAlleleCountsByAllele - the known allele counts in the population. For AC=2 AN=12 site, this is {10,2}
     * @param ploidy - the ploidy to assume
     * @param useFlatPriors - if true, apply flat priors to likelihoods in order to calculate posterior probabilities
     * @return - the posterior genotype likelihoods
     */
    protected static List<double[]> calculatePosteriorGLs(final List<double[]> genotypeLikelihoods,
                                                          final double[] knownAlleleCountsByAllele,
                                                          final int ploidy,
                                                          final boolean useFlatPriors) {
        if ( ploidy != 2 ) {
            throw new IllegalStateException("Genotype posteriors not yet implemented for ploidy != 2");
        }

        final double[] genotypePriorByAllele = getDirichletPrior(knownAlleleCountsByAllele,ploidy, useFlatPriors);
        final List<double[]> posteriors = new ArrayList<>(genotypeLikelihoods.size());
        for ( final double[] likelihoods : genotypeLikelihoods ) {
            double[] posteriorProbabilities = null;

            if ( likelihoods != null ) {
                if ( likelihoods.length != genotypePriorByAllele.length ) {
                    throw new IllegalStateException(String.format("Likelihoods not of correct size: expected %d, observed %d",
                            knownAlleleCountsByAllele.length*(knownAlleleCountsByAllele.length+1)/2,likelihoods.length));
                }

                posteriorProbabilities = new double[genotypePriorByAllele.length];
                for ( int genoIdx = 0; genoIdx < likelihoods.length; genoIdx ++ ) {
                    posteriorProbabilities[genoIdx] = likelihoods[genoIdx] + genotypePriorByAllele[genoIdx];
                }

                posteriorProbabilities = MathUtils.normalizeFromLog10(posteriorProbabilities, true);

            }

            posteriors.add(posteriorProbabilities);
        }

        return posteriors;
    }

    // convenience function for a single genotypelikelihoods array. Just wraps.
    protected static double[] calculatePosteriorGLs(final double[] genotypeLikelihoods,
                                                 final double[] knownAlleleCountsByAllele,
                                                 final int ploidy,
                                                 final boolean useFlatPriors) {
        return calculatePosteriorGLs(Arrays.asList(genotypeLikelihoods),knownAlleleCountsByAllele,ploidy, useFlatPriors).get(0);
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
    protected static double[] getDirichletPrior(final double[] knownCountsByAllele, final int ploidy, final boolean useFlatPrior) {
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
                if (useFlatPrior)
                    priors[priorIndex++] = 1.0;
                else {
                    final int[] counts = new int[knownCountsByAllele.length];
                    counts[allele1] += 1;
                    counts[allele2] += 1;
                    priors[priorIndex++] = MathUtils.dirichletMultinomial(knownCountsByAllele,sumOfKnownCounts,counts,ploidy);
                }
            }
        }

        return priors;
    }

    /**
     * Parse counts for each allele
     * @param counts - Map to store and return data
     * @param context - line to be parsed from the input VCF file
     * @param useAC - use allele count annotation value from VariantContext (vs. MLEAC)
     */
    private static void addAlleleCounts(final Map<Allele,Integer> counts, final VariantContext context, final boolean useAC) {
        final int[] ac;
        //use MLEAC value...
        if ( context.hasAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) && ! useAC ) {
            ac = getAlleleCounts(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, context);
        }
        //...unless specified by the user in useAC or unless MLEAC is absent
        else if ( context.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            ac = getAlleleCounts(VCFConstants.ALLELE_COUNT_KEY, context);
        }
        //if VariantContext annotation doesn't contain AC or MLEAC then get the data from direct evaluation
        else {
            ac = new int[context.getAlternateAlleles().size()];
            int idx = 0;
            for ( final Allele allele : context.getAlternateAlleles() ) {
                ac[idx++] = context.getCalledChrCount(allele);
            }
        }

        //since the allele count for the reference allele is not given in the VCF format,
        //calculate it from the allele number minus the total counts for alternate alleles
        for ( final Allele allele : context.getAlleles() ) {
            final int count;
            if ( allele.isReference() ) {
                if ( context.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
                    count = Math.max(context.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY,-1) - (int) MathUtils.sum(ac),0); //occasionally an MLEAC value will sneak in that's greater than the AN
                } else {
                    count = Math.max(context.getCalledChrCount() - (int) MathUtils.sum(ac),0);
                }
            } else {
                count = ac[context.getAlternateAlleles().indexOf(allele)];
            }
            //if this allele isn't in the map yet, add it
            if ( ! counts.containsKey(allele) ) {
                counts.put(allele,0);
            }
            //add the count for the current allele to the existing value in the map
            counts.put(allele,count + counts.get(allele));
        }
    }

    /**
     * Retrieve allele count data from VariantContext using VCFkey, checks for correct number of values in VCF
     * @param VCFkey VariantContext annotation tag of interest (should be AC or MLEAC)
     * @param context VariantContext from which to extract the data
     * @return int[] with allele count data
     */
    private static int[] getAlleleCounts(final String VCFkey, final VariantContext context) {
        final Object alleleCountsFromVCF = context.getAttribute(VCFkey);
        if ( alleleCountsFromVCF instanceof List ) {
            if ( ((List) alleleCountsFromVCF).size() != context.getAlternateAlleles().size() )
                throw new UserException(String.format("Variant does not contain the same number of MLE allele counts as alternate alleles for record at %s:%d", context.getChr(), context.getStart()));
        }
        else if ( alleleCountsFromVCF instanceof String || alleleCountsFromVCF instanceof Integer) {//here length is 1
            if (context.getAlternateAlleles().size() != 1)
                throw new UserException(String.format("Variant does not contain the same number of MLE allele counts as alternate alleles for record at %s:%d", context.getChr(), context.getStart()));
        }
        return extractInts(alleleCountsFromVCF);
    }

    /**
     * Check the formatting on the Object returned by a call to VariantContext::getAttribute() and parse appropriately
     * @param integerListContainingVCField - Object returned by a call to VariantContext::getAttribute()
     * @return - array of ints
     */
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
                    GATKVCFConstants.MLE_ALLELE_COUNT_KEY+" or "+VCFConstants.ALLELE_COUNT_KEY));

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
