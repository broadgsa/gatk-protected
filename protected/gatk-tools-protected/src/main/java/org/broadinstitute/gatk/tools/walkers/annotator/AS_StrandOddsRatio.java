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
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.*;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Allele-specific strand bias estimated by the Symmetric Odds Ratio test
 *
 * <p>Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other. </p>
 *
 * <p>The AS_StrandOddsRatio annotation is one of several methods that aims to evaluate whether there is strand bias in the data. It is an updated form of the Fisher Strand Test that is better at taking into account large amounts of data in high coverage situations. It is used to determine if there is strand bias between forward and reverse strands for the reference or alternate allele. It does so separately for each allele. The reported value is ln-scaled.</p>
 *
 * <h3>Statistical notes</h3>
 * <p> Odds Ratios in the 2x2 contingency table below are</p>
 *
 * $$ R = \frac{X[0][0] * X[1][1]}{X[0][1] * X[1][0]} $$
 *
 * <p>and its inverse:</p>
 *
 * <table>
 *      <tr><td>&nbsp;</td><td>+ strand </td><td>- strand</td></tr>
 *      <tr><td>REF;</td><td>X[0][0]</td><td>X[0][1]</td></tr>
 *      <tr><td>ALT;</td><td>X[1][0]</td><td>X[1][1]</td></tr>
 * </table>
 *
 * <p>The sum R + 1/R is used to detect a difference in strand bias for REF and for ALT (the sum makes it symmetric). A high value is indicative of large difference where one entry is very small compared to the others. A scale factor of refRatio/altRatio where</p>
 *
 * $$ refRatio = \frac{max(X[0][0], X[0][1])}{min(X[0][0], X[0][1} $$
 *
 * <p>and </p>
 *
 * $$ altRatio = \frac{max(X[1][0], X[1][1])}{min(X[1][0], X[1][1]} $$
 *
 * <p>ensures that the annotation value is large only. </p>
 *
 * <p>See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this statistical test.</p>
 *
 * <h3>Caveat</h3>
 * <p>
 * The name AS_StrandOddsRatio is not entirely appropriate because the implementation was changed somewhere between the start of development and release of this annotation. Now SOR isn't really an odds ratio anymore. The goal was to separate certain cases of data without penalizing variants that occur at the ends of exons because they tend to only be covered by reads in one direction (depending on which end of the exon they're on), so if a variant has 10 ref reads in the + direction, 1 ref read in the - direction, 9 alt reads in the + direction and 2 alt reads in the - direction, it's actually not strand biased, but the FS score is pretty bad. The implementation that resulted derived in part from empirically testing some read count tables of various sizes with various ratios and deciding from there.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php">StrandOddsRatio</a></b> outputs a version of this annotation that includes all alternate alleles in a single calculation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandBiasBySample.php">StrandBiasBySample</a></b> outputs counts of read depth per allele for each strand orientation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_FisherStrand.php">FisherStrand</a></b> uses Fisher's Exact Test to evaluate strand bias.</li>
 * </ul>
 *
 */
public class AS_StrandOddsRatio extends AS_StrandBiasTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoodMap(Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                                       final VariantContext vc){
        // either SNP with no alignment context, or indels: per-read likelihood map needed
        final int[][] table = getContingencyTable(stratifiedPerReadAlleleLikelihoodMap, vc, MIN_COUNT);
        final double ratio = calculateSOR(table);
        return Collections.singletonMap(getKeyNames().get(0), (Object)String.format("%.3f",ratio));
    }

    @Override
    protected Map<Allele,Double> calculateReducedData(AlleleSpecificAnnotationData<List<Integer>> combinedData) {
        final Map<Allele,Double> annotationMap = new HashMap<>();
        final Map<Allele, List<Integer>> perAlleleData = combinedData.getAttributeMap();
        final List<Integer> refStrandCounts = perAlleleData.get(combinedData.getRefAllele());
        for (final Allele a : perAlleleData.keySet()) {
            List<Integer> altStrandCounts = perAlleleData.get(a);
            int[][] refAltTable = new int[][] {new int[]{refStrandCounts.get(0),refStrandCounts.get(1)},new int[]{altStrandCounts.get(0),altStrandCounts.get(1)}};
            annotationMap.put(a,calculateSOR(refAltTable));
        }
        return annotationMap;
    }

    /**
     * Computes the SOR value of a table after augmentation (adding pseudocounts). Based on the symmetric odds ratio but modified to take on
     * low values when the reference +/- read count ratio is skewed but the alt count ratio is not.  Natural log is taken
     * to keep values within roughly the same range as other annotations.
     *
     * Adding pseudocounts prevent divide-by-zero.
     *
     * @param originalTable The table before augmentation
     * @return the SOR annotation value
     */
    final protected double calculateSOR(final int[][] originalTable) {
        final double[][] augmentedTable = StrandBiasTableUtils.augmentContingencyTable(originalTable);

        double ratio = 0;

        ratio += (augmentedTable[0][0] / augmentedTable[0][1]) * (augmentedTable[1][1] / augmentedTable[1][0]);
        ratio += (augmentedTable[0][1] / augmentedTable[0][0]) * (augmentedTable[1][0] / augmentedTable[1][1]);

        final double refRatio = (Math.min(augmentedTable[0][0], augmentedTable[0][1])/Math.max(augmentedTable[0][0], augmentedTable[0][1]));
        final double altRatio = (Math.min(augmentedTable[1][0], augmentedTable[1][1])/Math.max(augmentedTable[1][0], augmentedTable[1][1]));

        ratio = ratio*refRatio/altRatio;

        return Math.log(ratio);
    }
}
