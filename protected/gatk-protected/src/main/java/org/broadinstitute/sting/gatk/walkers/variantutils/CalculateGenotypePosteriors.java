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

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.variant.vcf.*;
import java.util.*;

/**
 * Calculates genotype posterior likelihoods given panel data
 *
 * <p>
 * Given a VCF with genotype likelihoods from the HaplotypeCaller, UnifiedGenotyper, or another source which provides
 * -unbiased- GLs, calculate the posterior genotype state and likelihood given allele frequency information from
 * both the samples themselves and input VCFs describing allele frequencies in related populations.
 *
 * VCFs to use for informing the genotype likelihoods (e.g. a population-specific VCF from 1000 genomes) should have
 * at least one of:
 * - AC field and AN field
 * - MLEAC field and AN field
 * - genotypes
 *
 * The AF field will not be used in this calculation as it does not provide a way to estimate the confidence interval
 * or uncertainty around the allele frequency, while AN provides this necessary information. This uncertainty is
 * modeled by a Dirichlet distribution: that is, the frequency is known up to a Dirichlet distribution with
 * parameters AC1+q,AC2+q,...,(AN-AC1-AC2-...)+q, where "q" is the global frequency prior (typically q << 1). The
 * genotype priors applied then follow a Dirichlet-Multinomial distribution, where 2 alleles per sample are drawn
 * independently. This assumption of independent draws is the assumption Hardy-Weinberg Equilibrium. Thus, HWE is
 * imposed on the likelihoods as a result of CalculateGenotypePosteriors.
 *
 * <h3>Input</h3>
 * <p>
 * A VCF with genotype likelihoods, and optionally genotypes, AC/AN fields, or MLEAC/AN fields
 * </p>
 *
 * <p>
 * A collection of VCFs to use for informing allele frequency priors. Each VCF must have one of
 * - AC field and AN field
 * - MLEAC field and AN field
 * - genotypes
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A new VCF with:
 *  1) Genotype posteriors added to the genotype fields ("GP")
 *  2) Genotypes and GQ assigned according to these posteriors
 *  3) Per-site genotype priors added to the INFO field ("PG")
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * Inform the genotype assignment of NA12878 using the 1000G Euro panel
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CalculateGenotypePosteriors \
 *   -V NA12878.wgs.HC.vcf \
 *   -VV 1000G_EUR.genotypes.combined.vcf \
 *   -o NA12878.wgs.HC.posteriors.vcf \
 *
 * Refine the genotypes of a large panel based on the discovered allele frequency
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CalculateGenotypePosteriors \
 *   -V input.vcf \
 *   -o output.withPosteriors.vcf
 *
 * Apply frequency and HWE-based priors to the genotypes of a family without including the family allele counts
 * in the allele frequency estimates
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CalculateGenotypePosteriors \
 *   -V input.vcf \
 *   -o output.withPosteriors.vcf \
 *   --ignoreInputSamples
 *
 * Calculate the posterior genotypes of a callset, and impose that a variant *not seen* in the external panel
 * is tantamount to being AC=0, AN=100 within that panel
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CalculateGenotypePosteriors \
 *   -VV external.panel.vcf \
 *   -V input.vcf \
 *   -o output.withPosteriors.vcf
 *   --numRefSamplesIfNoCall 100
 *
 * </pre>
 *
 */
public class CalculateGenotypePosteriors extends RodWalker<Integer,Integer> {

    /**
     * The input VCF (posteriors will be calculated for these samples, and written to the output)
     */
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * Supporting external panels. Allele counts from these panels (taken from AC,AN or MLEAC,AN or raw genotypes) will
     * be used to inform the frequency distribution underying the genotype priors.
     */
    @Input(fullName="supporting", shortName = "VV", doc="Other callsets to use in generating genotype posteriors", required=false)
    public List<RodBinding<VariantContext>> supportVariants = new ArrayList<RodBinding<VariantContext>>();

    /**
     * The global prior of a variant site -- i.e. the expected allele frequency distribution knowing only that N alleles
     * exist, and having observed none of them. This is the "typical" 1/x trend, modeled here as not varying
     * across alleles. The calculation for this parameter is (Effective population size) * (steady state mutation rate)
     *
     */
     @Argument(fullName="globalPrior",shortName="G",doc="The global Dirichlet prior parameters for the allele frequency",required=false)
     public double globalPrior = UnifiedGenotyperEngine.HUMAN_SNP_HETEROZYGOSITY;

    /**
     * When a variant is not seen in a panel, whether to infer (and with what effective strength) that only reference
     * alleles were ascertained at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000". This is
     * applied across all external panels, so if numRefIsMissing = 10, and the variant is absent in two panels, this
     * confers evidence of AC=0,AN=20
     */
    @Argument(fullName="numRefSamplesIfNoCall",shortName="nrs",doc="The number of homozygous reference to infer were " +
            "seen at a position where an \"other callset\" contains no site or genotype information",required=false)
    public int numRefIfMissing = 1;

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

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    private final boolean NO_EM = false;

    public void initialize() {
        // Get list of samples to include in the output
        final List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        final Map<String,VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);

        if ( vcfRods.size() > 1 )
            throw new IllegalStateException("Somehow more than one variant was bound?");

        final VCFHeader header = new ArrayList<>(vcfRods.values()).get(0); // pure laziness

        if ( ! header.hasGenotypingData() ) {
            throw new UserException("VCF has no genotypes");
        }

        if ( header.hasInfoLine(VCFConstants.MLE_ALLELE_COUNT_KEY) ) {
            final VCFInfoHeaderLine mleLine = header.getInfoHeaderLine(VCFConstants.MLE_ALLELE_COUNT_KEY);
            if ( mleLine.getCountType() != VCFHeaderLineCount.A ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the count type should be \"A\"");
            }

            if ( mleLine.getType() != VCFHeaderLineType.Integer ) {
                throw new UserException("VCF does not have a properly formatted MLEAC field: the field type should be \"Integer\"");
            }
        }

        final TreeSet<String> vcfSamples = new TreeSet<>(SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

        // Initialize VCF header
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_POSTERIORS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Posterior Genotype Likelihoods"));
        headerLines.add(new VCFInfoHeaderLine("PG", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Genotype Likelihood Prior"));
        headerLines.add(new VCFHeaderLine("source", "CalculateGenotypePosteriors"));

        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
    }

    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || context == null || ref == null ) {
            return 0;
        }

        final Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, ref.getLocus());

        final Collection<VariantContext> otherVCs = tracker.getValues(supportVariants, context.getLocation());

        final int missing = supportVariants.size() - otherVCs.size();

        for ( VariantContext vc : vcs ) {
            vcfWriter.add(PosteriorLikelihoodsUtils.calculatePosteriorGLs(vc, otherVCs, missing * numRefIfMissing, globalPrior, !ignoreInputSamples, NO_EM, defaultToAC));
        }

        return 1;
    }

    public Integer reduce(Integer l, Integer r) { return r + l; }
}
