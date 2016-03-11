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

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.MathUtils;
import java.util.HashMap;
import java.util.Map;

/**
 * A class containing utility methods used in the calculation of annotations related to cohort heterozygosity, e.g. InbreedingCoefficient and ExcessHet
 * Stores sample count to make sure we never have to iterate the genotypes more than once
 * Should be reinitialized for each VariantContext
 */
public class HeterozygosityUtils {

    final public static int REF_INDEX = 0;
    final public static int HET_INDEX = 1;
    final public static int VAR_INDEX = 2;

    protected int sampleCount = -1;
    private Map<Allele, Double> hetCounts;
    private Map<Allele, Double> alleleCounts;
    boolean returnRounded = false;

    /**
     * Create a new HeterozygosityUtils -- a new class should be instantiated for each VariantContext to store data for that VC
     * @param returnRounded round the likelihoods to return integer numbers of counts (as doubles)
     */
    protected HeterozygosityUtils(final boolean returnRounded) {
        this.returnRounded = returnRounded;
    }

    /**
     * Get the genotype counts for A/A, A/B, and B/B where A is the reference and B is any alternate allele
     * @param vc
     * @param genotypes may be subset to just founders if a pedigree file is provided
     * @return may be null, otherwise length-3 double[] representing homRef, het, and homVar counts
     */
    protected double[] getGenotypeCountsForRefVsAllAlts(final VariantContext vc, final GenotypesContext genotypes) {
        if (genotypes == null || !vc.isVariant())
            return null;

        final boolean doMultiallelicMapping = !vc.isBiallelic();

        int idxAA = 0, idxAB = 1, idxBB = 2;

        double refCount = 0;
        double hetCount = 0;
        double homCount = 0;

        sampleCount = 0;
        for (final Genotype g : genotypes) {
            if (g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2)  // only work for diploid samples
                sampleCount++;
            else
                continue;

            //Need to round the likelihoods to deal with small numerical deviations due to normalizing
            final double[] normalizedLikelihoodsUnrounded = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
            double[] normalizedLikelihoods = new double[normalizedLikelihoodsUnrounded.length];
            if (returnRounded) {
                for (int i = 0; i < normalizedLikelihoodsUnrounded.length; i++) {
                    normalizedLikelihoods[i] = Math.round(normalizedLikelihoodsUnrounded[i]);
                }
            } else {
                normalizedLikelihoods = normalizedLikelihoodsUnrounded;
            }

            if (doMultiallelicMapping) {
                if (g.isHetNonRef()) {
                    //all likelihoods go to homCount
                    homCount++;
                    continue;
                }

                if (!g.isHomRef()) {
                    //get alternate allele for each sample
                    final Allele a1 = g.getAllele(0);
                    final Allele a2 = g.getAllele(1);
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2.isNonReference() ? a2 : a1);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
            }

            refCount += normalizedLikelihoods[idxAA];
            hetCount += normalizedLikelihoods[idxAB];
            homCount += normalizedLikelihoods[idxBB];
        }
        return new double[]{refCount, hetCount, homCount};
    }

    /**
     * Get the count of heterozygotes in vc for a specific altAllele (both reference and non-reference hets, e.g. 1/2)
     * @param vc
     */
    protected void doGenotypeCalculations(final VariantContext vc) {
        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || !vc.isVariant())
            return;

        final int numAlleles = vc.getNAlleles();

        sampleCount = 0;
        if (hetCounts == null && alleleCounts == null) {
            hetCounts = new HashMap<>();
            alleleCounts = new HashMap<>();
            for (final Allele a : vc.getAlleles()) {
                if (a.isNonReference())
                    hetCounts.put(a, 0.0);
                alleleCounts.put(a, 0.0);
            }

            int idxAB;

            //for each sample
            for (final Genotype g : genotypes) {
                if (g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2)  // only work for diploid samples
                    sampleCount++;
                else
                    continue;

                int altIndex = 0;
                for(final Allele a : vc.getAlternateAlleles()) {
                    //for each alt allele index from 1 to N
                    altIndex++;

                        final double[] normalizedLikelihoodsUnrounded = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
                        double[] normalizedLikelihoods = new double[normalizedLikelihoodsUnrounded.length];
                        if (returnRounded) {
                            for (int i = 0; i < normalizedLikelihoodsUnrounded.length; i++) {
                                normalizedLikelihoods[i] = Math.round(normalizedLikelihoodsUnrounded[i]);
                            }
                        } else {
                            normalizedLikelihoods = normalizedLikelihoodsUnrounded;
                        }
                        //iterate over the other alleles
                        for (int i = 0; i < numAlleles; i++) {
                            //only add homozygotes to alleleCounts, not hetCounts
                            if (i == altIndex) {
                                final double currentAlleleCounts = alleleCounts.get(a);
                                alleleCounts.put(a, currentAlleleCounts + 2*normalizedLikelihoods[GenotypeLikelihoods.calculatePLindex(altIndex,altIndex)]);
                                continue;
                            }
                            //pull out the heterozygote PL index, ensuring that the first allele index < second allele index
                            idxAB = GenotypeLikelihoods.calculatePLindex(Math.min(i,altIndex),Math.max(i,altIndex));
                            final double aHetCounts = hetCounts.get(a);
                            hetCounts.put(a, aHetCounts + normalizedLikelihoods[idxAB]);
                            final double currentAlleleCounts = alleleCounts.get(a);
                            //these are guaranteed to be hets
                            alleleCounts.put(a, currentAlleleCounts + normalizedLikelihoods[idxAB]);
                            final double refAlleleCounts = alleleCounts.get(vc.getReference());
                            alleleCounts.put(vc.getReference(), refAlleleCounts + normalizedLikelihoods[idxAB]);
                        }
                    //add in ref/ref likelihood
                    final double refAlleleCounts = alleleCounts.get(vc.getReference());
                    alleleCounts.put(vc.getReference(), refAlleleCounts + 2*normalizedLikelihoods[0]);
                }

            }
        }
    }

    /**
     * Get the count of heterozygotes in vc for a specific altAllele (both reference and non-reference hets, e.g. 1/2)
     * @param vc
     * @param altAllele the alternate allele of interest
     * @return number of hets
     */
    protected double getHetCount(final VariantContext vc, final Allele altAllele) {
        if (hetCounts == null)
            doGenotypeCalculations(vc);
        return hetCounts.containsKey(altAllele)? hetCounts.get(altAllele) : 0;
    }

    protected double getAlleleCount(final VariantContext vc, final Allele allele) {
        if (alleleCounts == null)
            doGenotypeCalculations(vc);
        return alleleCounts.containsKey(allele)? alleleCounts.get(allele) : 0;
    }

    protected int getSampleCount() {
        return sampleCount;
    }
}
