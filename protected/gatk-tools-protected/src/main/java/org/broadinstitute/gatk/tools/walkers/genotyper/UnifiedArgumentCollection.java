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

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.gatk.utils.pairhmm.PairHMM;
import htsjdk.variant.variantcontext.VariantContext;

public class UnifiedArgumentCollection extends StandardCallerArgumentCollection {

    @Argument(fullName = "genotype_likelihoods_model", shortName = "glm", doc = "Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together", required = false)
    public GenotypeLikelihoodsCalculationModel.Model GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;

    /**
     * The PCR error rate is independent of the sequencing error rate, which is necessary because we cannot necessarily
     * distinguish between PCR errors vs. sequencing errors.  The practical implication for this value is that it
     * effectively acts as a cap on the base qualities.
     */
    @Argument(fullName = "pcr_error_rate", shortName = "pcr_error", doc = "The PCR error rate to be used for computing fragment-based likelihoods", required = false)
    public Double PCR_error = DiploidSNPGenotypeLikelihoods.DEFAULT_PCR_ERROR_RATE;

    /**
     * Note that calculating the SLOD increases the runtime by an appreciable amount.
     */
    @Argument(fullName = "computeSLOD", shortName = "slod", doc = "If provided, we will calculate the SLOD (SB annotation)", required = false)
    public boolean COMPUTE_SLOD = false;

    /**
     * The PairHMM implementation to use for -glm INDEL genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for -glm INDEL genotype likelihood calculations", required = false)
    public PairHMM.HMM_IMPLEMENTATION pairHMM = PairHMM.HMM_IMPLEMENTATION.LOGLESS_CACHING;

    /**
     * The minimum confidence needed in a given base for it to be used in variant calling.  Note that the base quality of a base
     * is capped by the mapping quality so that bases on reads with low mapping quality may get filtered out depending on this value.
     * Note too that this argument is ignored in indel calling. In indel calling, low-quality ends of reads are clipped off (with fixed threshold of Q20).
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 17;

    /**
     * If the fraction of reads with deletions spanning a locus is greater than this value, the site will not be considered callable and will be skipped.
     * To disable the use of this parameter, set its value to >1.
     */
    @Argument(fullName = "max_deletion_fraction", shortName = "deletions", doc = "Maximum fraction of reads with deletions spanning this locus for it to be callable", required = false)
    public Double MAX_DELETION_FRACTION = 0.05;

    // indel-related arguments
    /**
     * A candidate indel is genotyped (and potentially called) if there are this number of reads with a consensus indel at a site.
     * Decreasing this value will increase sensitivity but at the cost of larger calling time and a larger number of false positives.
     */
    @Argument(fullName = "min_indel_count_for_genotyping", shortName = "minIndelCnt", doc = "Minimum number of consensus indels required to trigger genotyping run", required = false)
    public int MIN_INDEL_COUNT_FOR_GENOTYPING = 5;

    /**
     * Complementary argument to minIndelCnt.  Only samples with at least this fraction of indel-containing reads will contribute
     * to counting and overcoming the threshold minIndelCnt.  This parameter ensures that in deep data you don't end
     * up summing lots of super rare errors up to overcome the 5 read default threshold.  Should work equally well for
     * low-coverage and high-coverage samples, as low coverage samples with any indel containing reads should easily over
     * come this threshold.
     */
    @Argument(fullName = "min_indel_fraction_per_sample", shortName = "minIndelFrac", doc = "Minimum fraction of all reads at a locus that must contain an indel (of any allele) for that sample to contribute to the indel count for alleles", required = false)
    public double MIN_INDEL_FRACTION_PER_SAMPLE = 0.25;

    @Advanced
    @Argument(fullName = "indelGapContinuationPenalty", shortName = "indelGCP", doc = "Indel gap continuation penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10", required = false)
    public byte INDEL_GAP_CONTINUATION_PENALTY = 10;

    @Advanced
    @Argument(fullName = "indelGapOpenPenalty", shortName = "indelGOP", doc = "Indel gap open penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10", required = false)
    public byte INDEL_GAP_OPEN_PENALTY = 45;

    @Hidden
    @Argument(fullName = "indelHaplotypeSize", shortName = "indelHSize", doc = "Indel haplotype size", required = false)
    public int INDEL_HAPLOTYPE_SIZE = 80;

    @Hidden
    @Argument(fullName = "indelDebug", shortName = "indelDebug", doc = "Output indel debug info", required = false)
    public boolean OUTPUT_DEBUG_INDEL_INFO = false;

    @Hidden
    @Argument(fullName = "ignoreSNPAlleles", shortName = "ignoreSNPAlleles", doc = "expt", required = false)
    public boolean IGNORE_SNP_ALLELES = false;

    /*
        Generalized ploidy argument (debug only): squash all reads into a single pileup without considering sample info
     */
    @Hidden
    @Argument(fullName = "allReadsSP", shortName = "dl", doc = "expt", required = false)
    public boolean TREAT_ALL_READS_AS_SINGLE_POOL = false;

    /*
       Generalized ploidy argument (debug only): When building site error models, ignore lane information and build only
       sample-level error model
     */
    @Hidden
    @Argument(fullName = "ignoreLaneInfo", shortName = "ignoreLane", doc = "Ignore lane when building error model, error model is then per-site", required = false)
    public boolean IGNORE_LANE_INFO = false;

    /*
        Generalized ploidy argument: VCF file that contains truth calls for reference sample. If a reference sample is included through argument -refsample,
        then this argument is required.
     */
    @Hidden
    @Input(fullName="reference_sample_calls", shortName = "referenceCalls", doc="VCF file with the truth callset for the reference sample", required=false)
    RodBinding<VariantContext> referenceSampleRod;

    /*
        Reference sample name: if included, a site-specific error model will be built in order to improve calling quality. This requires ideally
        that a bar-coded reference sample be included with the polyploid/pooled data in a sequencing experimental design.
        If argument is absent, no per-site error model is included and calling is done with a generalization of traditional statistical calling.
     */
    @Hidden
    @Argument(shortName="refsample", fullName="reference_sample_name", doc="Reference sample name.", required=false)
    String referenceSampleName;

    /**
     * The following argument are for debug-only tweaks when running generalized ploidy with a reference sample
     */
    @Hidden
    @Argument(shortName="minqs", fullName="min_quality_score", doc="Min quality score to consider. Smaller numbers process faster. Default: Q1.", required=false)
    byte minQualityScore= 1;

    @Hidden
    @Argument(shortName="maxqs", fullName="max_quality_score", doc="Max quality score to consider. Smaller numbers process faster. Default: Q40.", required=false)
    byte maxQualityScore= 40;

    @Hidden
    @Argument(shortName="site_prior", fullName="site_quality_prior", doc="Phred-Scaled prior quality of the site. Default: Q20.", required=false)
    byte phredScaledPrior = 20;

    @Hidden
    @Argument(shortName = "min_call_power", fullName = "min_power_threshold_for_calling", doc="The minimum confidence in the error model to make a call. Number should be between 0 (no power requirement) and 1 (maximum power required).", required = false)
    double minPower = 0.95;

    /**
     * Create a new UAC with defaults for all UAC arguments
     */
    public UnifiedArgumentCollection() {
        super();
    }

    @Override
    public UnifiedArgumentCollection clone() {
        return (UnifiedArgumentCollection) super.clone();
    }
}
