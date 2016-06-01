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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.MathUtils;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.*;


/**
 * Phred-scaled p-value for exact test of excess heterozygosity
 *
 * <p>This annotation estimates excess heterozygosity in a population of samples. It is related to but distinct from InbreedingCoeff, which estimates evidence for inbreeding in a population. ExcessHet scales more reliably to large cohort sizes.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>This annotation is a one-sided phred-scaled p-value using an exact test of the Hardy-Weinberg Equilibrium. The null hypothesis is that the number of heterozygotes follows the Hardy-Weinberg Equilibrium. The p-value is the probability of getting the same or more heterozygotes as was observed, given the null hypothesis. </p>
 * <p>The implementation used is adapted from Wigginton JE, Cutler DJ, Abecasis GR. A Note on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 2005;76(5):887-893.</p>
 * <p>The p-value is calculated exactly by using the Levene-Haldane distribution. This implementation also uses a mid-p correction as described by Graffelman, J. & Moreno, V. (2013). The mid p-value in exact tests for Hardy-Weinberg equilibrium. Statistical Applications in Genetics and Molecular Biology, 12(4), pp. 433-448. </p>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The annotation is not accurate for very small p-values. Beyond 1.0E-16 there is no guarantee that the p-value is accurate, just that it is in fact smaller than 1.0E-16. </li>
 * <li>For multiallelic sites, all non-reference alleles are treated as a single alternate allele.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_InbreedingCoeff.php">InbreedingCoeff</a></b> estimates whether there is evidence of inbreeding in a population</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_AS_InbreedingCoeff.php">AS_InbreedingCoeff</a></b> outputs an allele-specific version of the InbreedingCoeff annotation.</li>
 * </ul>
 *
 */
public class ExcessHet extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
    private final static Logger logger = Logger.getLogger(ExcessHet.class);
    private final double minNeededValue = 1.0E-16;
    private Set<String> founderIds;
    private final boolean RETURN_ROUNDED = true;
    private int sampleCount = -1;

    @Override
    public void initialize ( AnnotatorCompatible walker, GenomeAnalysisEngine toolkit, Set<VCFHeaderLine> headerLines ) {
        //If available, get the founder IDs and cache them. The ExcessHet value will only be computed on founders then.
        //excessHet respects pedigree files, but doesn't require a minimum number of samples
        if(founderIds == null && walker != null) {
            founderIds = ((Walker) walker).getSampleDB().getFounderIds();
        }

    }

    @Override
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        return makeEHAnnotation(vc);
    }

    protected double calculateEH(final VariantContext vc, final GenotypesContext genotypes) {
        HeterozygosityUtils heterozygosityUtils = new HeterozygosityUtils(RETURN_ROUNDED);
        final double[] genotypeCountsDoubles = heterozygosityUtils.getGenotypeCountsForRefVsAllAlts(vc, genotypes);
        sampleCount = heterozygosityUtils.getSampleCount();
        final int[] genotypeCounts = new int[genotypeCountsDoubles.length];
        for(int i = 0; i < genotypeCountsDoubles.length; i++) {
            genotypeCounts[i] = (int)genotypeCountsDoubles[i];
        }

        double pval = exactTest(genotypeCounts);

        //If the actual phredPval would be infinity we will probably still filter out just a very large number
        if (pval == 0) {
            return Integer.MAX_VALUE;
        }
        double phredPval = -10.0 * Math.log10(pval);

        return phredPval;
    }

    /**
     * Note that this method is not accurate for very small p-values. Beyond 1.0E-16 there is no guarantee that the
     * p-value is accurate, just that it is in fact smaller than 1.0E-16 (and therefore we should filter it). It would
     * be more computationally expensive to calculate accuracy beyond a given threshold. Here we have enough accuracy
     * to filter anything below a p-value of 10E-6.
     *
     * @param genotypeCounts Number of observed genotypes (n_aa, n_ab, n_bb)
     * @return Right sided p-value or the probability of getting the observed or higher number of hets given the sample
     * size (N) and the observed number of allele a (rareCopies)
     */
    protected double exactTest(final int[] genotypeCounts) {
        if (genotypeCounts.length != 3) {
            throw new IllegalStateException("Input genotype counts must be length 3 for the number of genotypes with {2, 1, 0} ref alleles.");
        }
        final int REF_INDEX = 0;
        final int HET_INDEX = 1;
        final int VAR_INDEX = 2;

        final int refCount = genotypeCounts[REF_INDEX];
        final int hetCount = genotypeCounts[HET_INDEX];
        final int homCount = genotypeCounts[VAR_INDEX];

        if (hetCount < 0 || refCount < 0 || homCount < 0) {
            throw new IllegalArgumentException("Genotype counts cannot be less than 0");
        }

        //Split into observed common allele and rare allele
        final int obsHomR;
        final int obsHomC;
        if (refCount < homCount) {
            obsHomR = refCount;
            obsHomC = homCount;
        } else {
            obsHomR = homCount;
            obsHomC = refCount;
        }

        final int rareCopies = 2 * obsHomR + hetCount;
        final int N = hetCount + obsHomC + obsHomR;

        //If the probability distribution has only 1 point, then the mid p-value is .5
        if (rareCopies <= 1) {
            return .5;
        }

        double[] probs = new double[rareCopies + 1];

        //Find (something close to the) mode for the midpoint
        int mid = (int) Math.floor(((double) rareCopies * (2.0 * (double) N - (double) rareCopies)) / (2.0 * (double) N - 1.0));
        if ((mid % 2) != (rareCopies % 2)) {
            mid++;
        }

        probs[mid] = 1.0;
        double mysum = 1.0;

        //Calculate probabilities from midpoint down
        int currHets = mid;
        int currHomR = (rareCopies - mid) / 2;
        int currHomC = N - currHets - currHomR;

        while (currHets >= 2) {
            double potentialProb = probs[currHets] * (double) currHets * ((double) currHets - 1.0) / (4.0 * ((double) currHomR + 1.0) * ((double) currHomC + 1.0));
            if (potentialProb < minNeededValue) {
                break;
            }

            probs[currHets - 2] = potentialProb;
            mysum = mysum + probs[currHets - 2];

            //2 fewer hets means one additional homR and homC each
            currHets = currHets - 2;
            currHomR = currHomR + 1;
            currHomC = currHomC + 1;
        }

        //Calculate probabilities from midpoint up
        currHets = mid;
        currHomR = (rareCopies - mid) / 2;
        currHomC = N - currHets - currHomR;

        while (currHets <= rareCopies - 2) {
            double potentialProb = probs[currHets] * 4.0 * (double) currHomR * (double) currHomC / (((double) currHets + 2.0) * ((double) currHets + 1.0));
            if (potentialProb < minNeededValue) {
                break;
            }

            probs[currHets + 2] = potentialProb;
            mysum = mysum + probs[currHets + 2];

            //2 more hets means 1 fewer homR and homC each
            currHets = currHets + 2;
            currHomR = currHomR - 1;
            currHomC = currHomC - 1;
        }

        double rightPval = probs[hetCount] / (2.0 * mysum);
        //Check if we observed the highest possible number of hets
        if (hetCount == rareCopies) {
            return rightPval;
        }
        rightPval = rightPval + StatUtils.sum(Arrays.copyOfRange(probs, hetCount + 1, probs.length)) / mysum;

        return (rightPval);
    }

    protected Map<String, Object> makeEHAnnotation(final VariantContext vc) {
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
        if (genotypes == null || !vc.isVariant())
            return null;
        double EH = calculateEH(vc, genotypes);
        if (sampleCount < 1)
            return null;
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", EH));
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.EXCESS_HET_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

}
