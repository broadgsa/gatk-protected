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

package org.broadinstitute.gatk.engine.arguments;

import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;

import java.util.Collections;
import java.util.List;

public class GenotypeCalculationArgumentCollection implements Cloneable{


    public static final String MAX_ALTERNATE_ALLELES_SHORT_NAME = "maxAltAlleles";

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles
     * being sent on for genotyping. Using this argument instructs the genotyper to annotate (in the INFO field) the
     * number of alternate alleles that were originally discovered (but not necessarily genotyped) at the site.
     */
    @Argument(fullName = "annotateNDA", shortName = "nda", doc = "Annotate number of alleles observed", required = false)
    public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

    /**
     * This activates a model for calculating QUAL that was introduced in version 3.7 (November 2016). We expect this
     * model will become the default in future versions.
     */
    @Argument(fullName = "useNewAFCalculator", shortName = "newQual", doc = "Use new AF model instead of the so-called exact model", required = false)
    public boolean USE_NEW_AF_CALCULATOR = false;

    /**
     * The expected heterozygosity value used to compute prior probability that a locus is non-reference. See
     * https://software.broadinstitute.org/gatk/documentation/article?id=8603 for more details.
     */
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double snpHeterozygosity = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * This argument informs the prior probability of having an indel at a site.
     */
    @Argument(fullName = "indel_heterozygosity", shortName = "indelHeterozygosity", doc = "Heterozygosity for indel calling", required = false)
    public double indelHeterozygosity = HomoSapiensConstants.INDEL_HETEROZYGOSITY;

    /**
     * The standard deviation of the distribution of alt allele fractions. The above heterozygosity parameters give
     * the *mean* of this distribution; this parameter gives its spread.
     */
    @Argument(fullName = "heterozygosity_stdev", shortName = "heterozygosityStandardDeviation", doc = "Standard deviation of eterozygosity for SNP and indel calling.", required = false)
    public double heterozygosityStandardDeviation = 0.01;

    /**
     * The minimum phred-scaled confidence threshold at which variants should be called. Only variant sites with QUAL equal
     * or greater than this threshold will be called. Note that since version 3.7, we no longer differentiate high confidence
     * from low confidence calls at the calling step. The default call confidence threshold is set low intentionally to achieve
     * high sensitivity, which will allow false positive calls as a side effect. Be sure to perform some kind of filtering after
     * calling to reduce the amount of false positives in your final callset. Note that when HaplotypeCaller is used in GVCF mode
     * (using either -ERC GVCF or -ERC BP_RESOLUTION) the call threshold is automatically set to zero. Call confidence thresholding
     * will then be performed in the subsequent GenotypeGVCFs command.
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_calling", shortName = "stand_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be called", required = false)
    public double STANDARD_CONFIDENCE_FOR_CALLING = 10.0;

    /**
     * This argument allows you to emit low quality calls as filtered records.
     */
    @Hidden
    @Deprecated
    @Argument(fullName = "standard_min_confidence_threshold_for_emitting", shortName = "stand_emit_conf",
            doc = "This argument is no longer used in GATK versions 3.7 and newer. Please see the online documentation for the latest usage recommendations.", required = false)
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or
     * GENOTYPE_GIVEN_ALLELES), then only this many alleles will be used.  Note that genotyping sites with many
     * alternate alleles is both CPU and memory intensive and it scales exponentially based on the number of alternate
     * alleles.  Unless there is a good reason to change the default value, we highly recommend that you not play around
     * with this parameter.
     *
     * See also {@link #MAX_GENOTYPE_COUNT}.
     */
    @Advanced
    @Argument(fullName = "max_alternate_alleles", shortName = MAX_ALTERNATE_ALLELES_SHORT_NAME, doc = "Maximum number of alternate alleles to genotype", required = false)
    public int MAX_ALTERNATE_ALLELES = 6;

    /**
     * If there are more than this number of genotypes at a locus presented to the genotyper, then only this many
     * genotypes will be used. This is intended to deal with sites where the combination of high ploidy and high alt
     * allele count can lead to an explosion in the number of possible genotypes, with extreme adverse effects on
     * runtime performance.
     *
     * How does it work? The possible genotypes are simply different ways of partitioning alleles given a specific
     * ploidy assumption. Therefore, we remove genotypes from consideration by removing alternate alleles that are the
     * least well supported. The estimate of allele support is based on the ranking of the candidate haplotypes coming
     * out of the graph building step. Note however that the reference allele is always kept.
     *
     * The maximum number of alternative alleles used in the genotyping step will be the lesser of the two:
     * 1. the largest number of alt alleles, given ploidy, that yields a genotype count no higher than {@link #MAX_GENOTYPE_COUNT}
     * 2. the value of {@link #MAX_ALTERNATE_ALLELES}
     *
     * As noted above, genotyping sites with large genotype counts is both CPU and memory intensive. Unless you have
     * a good reason to change the default value, we highly recommend that you not play around with this parameter.
     *
     * See also {@link #MAX_ALTERNATE_ALLELES}.
     */
    @Advanced
    @Argument(fullName = "max_genotype_count", shortName = "maxGT", doc = "Maximum number of genotypes to consider at any site", required = false)
    public int MAX_GENOTYPE_COUNT = 1024;

    /**
     * Determines the maximum number of PL values that will be logged in the output.  If the number of genotypes
     * (which is determined by the ploidy and the number of alleles) exceeds the value provided by this argument,
     * then output of all of the PL values will be suppressed.
     */
    @Advanced
    @Argument(fullName = "max_num_PL_values", shortName = "maxNumPLValues", doc = "Maximum number of PL values to output", required = false)
    public int MAX_NUM_PL_VALUES = AFCalculator.MAX_NUM_PL_VALUES_DEFAULT;

    /**
     * By default, the prior specified with the argument --heterozygosity/-hets is used for variant discovery at a
     * particular locus, using an infinite sites model (see e.g. Waterson, 1975 or Tajima, 1996). This model asserts that
     * the probability of having a population of k variant sites in N chromosomes is proportional to theta/k, for 1=1:N.
     * However, there are instances where using this prior might not be desirable, e.g. for population studies where prior
     * might not be appropriate, as for example when the ancestral status of the reference allele is not known.
     *
     * This argument allows you to manually specify a list of probabilities for each AC>0 to be used as
     * priors for genotyping, with the following restrictions: only diploid calls are supported; you must specify 2 *
     * N values where N is the number of samples; probability values must be positive and specified in Double format,
     * in linear space (not log10 space nor Phred-scale); and all values must sume to 1.
     *
     * For completely flat priors,  specify the same value (=1/(2*N+1)) 2*N times, e.g.
     *      -inputPrior 0.33 -inputPrior 0.33
     * for the single-sample diploid case.
     */
    @Advanced
    @Argument(fullName = "input_prior", shortName = "inputPrior", doc = "Input prior for calls", required = false)
    public List<Double> inputPrior = Collections.emptyList();

    /**
     *   Sample ploidy - equivalent to number of chromosome copies per pool. For pooled experiments this should be set to
     *   the number of samples in pool multiplied by individual sample ploidy.
     */
    @Argument(shortName="ploidy", fullName="sample_ploidy", doc="Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).", required=false)
    public int samplePloidy = HomoSapiensConstants.DEFAULT_PLOIDY;

    /**
     * Creates a copy of this configuration.
     * @return never {@code null}.
     */
    @Override
    public GenotypeCalculationArgumentCollection clone() {
        try {
            return (GenotypeCalculationArgumentCollection) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException("unreachable code");
        }
    }
}
