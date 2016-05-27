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


import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.samples.Trio;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;

import java.util.*;

/**
 * Calculate genotype posterior likelihoods given panel data
 *
 * <p>
 * Given a VCF with genotype likelihoods from the HaplotypeCaller, UnifiedGenotyper, or another source which provides
 * <b>unbiased</b> genotype likelihoods, calculate the posterior genotype state and likelihood given allele frequency
 * information from both the samples themselves and input VCFs describing allele frequencies in related populations.</p>
 *
 * <p>The AF field will not be used in this calculation as it does not provide a way to estimate the confidence interval
 * or uncertainty around the allele frequency, while AN provides this necessary information. This uncertainty is
 * modeled by a Dirichlet distribution: that is, the frequency is known up to a Dirichlet distribution with
 * parameters AC1+q,AC2+q,...,(AN-AC1-AC2-...)+q, where "q" is the global frequency prior (typically q << 1). The
 * genotype priors applied then follow a Dirichlet-Multinomial distribution, where 2 alleles per sample are drawn
 * independently. This assumption of independent draws is the assumption Hardy-Weinberg Equilibrium. Thus, HWE is
 * imposed on the likelihoods as a result of CalculateGenotypePosteriors.</p>
 *
 * <h3>Input</h3>
 * <p>
 *     <ul>
 *         <li>A VCF with genotype likelihoods, and optionally genotypes, AC/AN fields, or MLEAC/AN fields</li>
 *         <li>(Optional) A PED pedigree file containing the description of the individuals relationships.</li>
 *     </ul>
 * </p>
 *
 * <p>
 * A collection of VCFs to use for informing allele frequency priors. Each VCF must have one of
 * </p>
 * <ul>
 *     <li>AC field and AN field</li>
 *     <li>MLEAC field and AN field</li>
 *     <li>genotypes</li>
 * </ul>
 * </p>
 *
 * <h3>Output</h3>
 * <p>A new VCF with:</p>
 * <ul>
 *     <li>Genotype posteriors added to the genotype fields ("PP")</li>
 *     <li>Genotypes and GQ assigned according to these posteriors</li>
 *     <li>Per-site genotype priors added to the INFO field ("PG")</li>
 *     <li>(Optional) Per-site, per-trio joint likelihoods (JL) and joint posteriors (JL) given as Phred-scaled probability
 *  of all genotypes in the trio being correct based on the PLs for JL and the PPs for JP. These annotations are added to
 *  the genotype fields.</li>
 * </ul>
 *
 * <h3>Notes</h3>
 * <p>
 * Using the default behavior, priors will only be applied for each variants (provided each variant has at least 10
 * called samples.) SNP sites in the input callset that have a SNP at the matching site in the supporting VCF will have
 * priors applied based on the AC from the supporting samples and the input callset (unless the --ignoreInputSamples
 * flag is used). If the site is not called in the supporting VCF, priors will be applied using the discovered AC from
 * the input samples (unless the --discoveredACpriorsOff flag is used). Flat priors are applied for any non-SNP sites in
 * the input callset.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <h4>Inform the genotype assignment of NA12878 using the 1000G Euro panel</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CalculateGenotypePosteriors \
 *   -R reference.fasta \
 *   -V NA12878.wgs.HC.vcf \
 *   -supporting 1000G_EUR.genotypes.combined.vcf \
 *   -o NA12878.wgs.HC.posteriors.vcf 
 * </pre>
 *
 * <h4>Refine the genotypes of a large panel based on the discovered allele frequency</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CalculateGenotypePosteriors \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.withPosteriors.vcf
 * </pre>
 *
 * <h4>Apply frequency and HWE-based priors to the genotypes of a family without including the family allele counts
 * in the allele frequency estimates the genotypes of a large panel based on the discovered allele frequency</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CalculateGenotypePosteriors \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.withPosteriors.vcf \
 *   --ignoreInputSamples
 * </pre>
 *
 * <h4>Calculate the posterior genotypes of a callset, and impose that a variant *not seen* in the external panel
 * is tantamount to being AC=0, AN=100 within that panel</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CalculateGenotypePosteriors \
 *   -R reference.fasta \
 *   -supporting external.panel.vcf \
 *   -V input.vcf \
 *   -o output.withPosteriors.vcf \
 *   --numRefSamplesIfNoCall 100
 * </pre>
 *
 * <h4>Apply only family priors to a callset</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CalculateGenotypePosteriors \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   --skipPopulationPriors \
 *   -ped family.ped \
 *   -o output.withPosteriors.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
public class CalculateGenotypePosteriors extends RodWalker<Integer,Integer> {

    /**
     * The input VCF (posteriors will be calculated for these samples, and written to the output)
     */
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * Supporting external panels. Allele counts from these panels (taken from AC,AN or MLEAC,AN or raw genotypes) will
     * be used to inform the frequency distribution underying the genotype priors. These files must be VCF 4.2 spec or later.
     */
    @Input(fullName="supporting", shortName = "supporting", doc="Other callsets to use in generating genotype posteriors", required=false)
    public List<RodBinding<VariantContext>> supportVariants = new ArrayList<>();

    /**
     * The global prior of a variant site -- i.e. the expected allele frequency distribution knowing only that N alleles
     * exist, and having observed none of them. This is the "typical" 1/x trend, modeled here as not varying
     * across alleles. The calculation for this parameter is (Effective population size) * (steady state mutation rate)
     *
     */
     @Argument(fullName="globalPrior",shortName="G",doc="The global Dirichlet prior parameters for the allele frequency",required=false)
     public double globalPrior = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * The mutation prior -- i.e. the probability that a new mutation occurs. Sensitivity analysis on known de novo 
     * mutations suggests a default value of 10^-6.
     *
     */
    @Argument(fullName="deNovoPrior",shortName="DNP",doc="The de novo mutation prior",required=false)
    public double deNovoPrior = 1e-6;

    /**
     * When a variant is not seen in a panel, whether to infer (and with what effective strength) that only reference
     * alleles were ascertained at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000". This is
     * applied across all external panels, so if numRefIsMissing = 10, and the variant is absent in two panels, this
     * confers evidence of AC=0,AN=20
     */
    @Argument(fullName="numRefSamplesIfNoCall",shortName="nrs",doc="The number of homozygous reference to infer were " +
            "seen at a position where an \"other callset\" contains no site or genotype information",required=false)
    public int numRefIfMissing = 0;

    /**
     * Rather than looking for the MLEAC field first, and then falling back to AC; first look for the AC field and then
     * fall back to MLEAC or raw genotypes
     */
    @Argument(fullName="defaultToAC",shortName="useAC",doc="Use the AC field as opposed to MLEAC. Does nothing if VCF lacks MLEAC field",required=false)
    public boolean defaultToAC = false;

    /**
     * Do not use the [MLE] allele count from the input samples (the ones for which you're calculating posteriors)
     * in the site frequency distribution; only use the AC and AN calculated from external sources.
     */
    @Argument(fullName="ignoreInputSamples",shortName="ext",doc="Use external information only; do not inform genotype priors by "+
          "the discovered allele frequency in the callset whose posteriors are being calculated. Useful for callsets containing "+
          "related individuals.",required=false)
    public boolean ignoreInputSamples = false;

    /**
     * Calculate priors for missing external variants from sample data -- default behavior is to apply flat priors
     */
    @Argument(fullName="discoveredACpriorsOff",shortName="useACoff",doc="Do not use discovered allele count in the input callset " +
            "for variants that do not appear in the external callset. ", required=false)
    public boolean useACoff = false;

    /**
     * Skip application of population-based priors
     */
    @Argument(fullName="skipPopulationPriors",shortName="skipPop",doc="Skip application of population-based priors", required=false)
    public boolean skipPopulationPriors = false;

    /**
     * Skip application of family-based priors.  Note: if pedigree file is absent, family-based priors will be skipped.
     */
    @Argument(fullName="skipFamilyPriors",shortName="skipFam",doc="Skip application of family-based priors", required=false)
    public boolean skipFamilyPriors = false;

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    private FamilyLikelihoodsUtils famUtils = new FamilyLikelihoodsUtils();

    public void initialize() {
        // Get list of samples to include in the output
        final List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        final Map<String,VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);

        final Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //Get the trios from the families passed as ped
        if (!skipFamilyPriors){
            final Set<Trio> trios = getSampleDB().getTrios();
            if(trios.size()<1) {
                logger.info("No PED file passed or no *non-skipped* trios found in PED file. Skipping family priors.");
                skipFamilyPriors = true;
            }
        }

        if ( vcfRods.size() > 1 )
            throw new IllegalStateException("Somehow more than one variant was bound?");

        final VCFHeader header = new ArrayList<>(vcfRods.values()).get(0); // pure laziness

        if ( ! header.hasGenotypingData() ) {
            throw new UserException("VCF has no genotypes");
        }

        if ( header.hasInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) ) {
            final VCFInfoHeaderLine mleLine = header.getInfoHeaderLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
            if ( mleLine.getCountType() != VCFHeaderLineCount.A ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the count type should be \"A\"");
            }

            if ( mleLine.getType() != VCFHeaderLineType.Integer ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the field type should be \"Integer\"");
            }
        }

        // Initialize VCF header
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
        if (!skipFamilyPriors) {
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.JOINT_LIKELIHOOD_TAG_NAME));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.JOINT_POSTERIOR_TAG_NAME));
        }
        headerLines.add(new VCFHeaderLine("source", "CalculateGenotypePosteriors"));

        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        Map<String,Set<Sample>> families = this.getSampleDB().getFamilies(vcfSamples);
        famUtils.initialize(deNovoPrior, vcfSamples, families);
    }

    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || context == null || ref == null ) {
            return 0;
        }

        final Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, ref.getLocus());

        final Collection<VariantContext> otherVCs = tracker.getValues(supportVariants, context.getLocation());

        final int missing = supportVariants.size() - otherVCs.size();

            for ( final VariantContext vc : vcs ) {
                VariantContext vc_familyPriors,vc_bothPriors;

                //do family priors first (if applicable)
                final VariantContextBuilder builder = new VariantContextBuilder(vc);
                //only compute family priors for biallelelic sites
                if (!skipFamilyPriors && vc.isBiallelic()){
                    GenotypesContext gc = famUtils.calculatePosteriorGLs(vc);
                    builder.genotypes(gc);
                }
                VariantContextUtils.calculateChromosomeCounts(builder, false);
                vc_familyPriors = builder.make();

                if (!skipPopulationPriors)
                    vc_bothPriors = PosteriorLikelihoodsUtils.calculatePosteriorGLs(vc_familyPriors, otherVCs, missing * numRefIfMissing, globalPrior, !ignoreInputSamples, defaultToAC, useACoff);
                else {
                    final VariantContextBuilder builder2 = new VariantContextBuilder(vc_familyPriors);
                    VariantContextUtils.calculateChromosomeCounts(builder, false);
                    vc_bothPriors = builder2.make();
                }
                vcfWriter.add(vc_bothPriors);
            }

        return 1;
    }

    public Integer reduce(Integer l, Integer r) { return r + l; }
}

