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

package org.broadinstitute.gatk.tools.walkers.cancer;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;

import java.util.*;

/**
 * Assigns somatic status to a set of calls
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class}, enable = false )
public class AssignSomaticStatus extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName="n", fullName="normalSample", required=true, doc="The normal sample")
    public String normalSample;

    @Argument(shortName="t", fullName="tumorSample", required=true, doc="The tumor sample")
    public String tumorSample;

    @Argument(shortName="somaticPriorQ", fullName="somaticPriorQ", required=false, doc="Phred-scaled probability that a site is a somatic mutation")
    public byte somaticPriorQ = 60;

    @Argument(shortName="somaticMinLOD", fullName="somaticMinLOD", required=false, doc="Phred-scaled min probability that a site should be called somatic mutation")
    public byte somaticMinLOD = 1;

    @Argument(shortName="minimalVCF", fullName="minimalVCF", required=false, doc="If provided, the attributes of the output VCF will only contain the somatic status fields")
    public boolean minimalVCF = false;

    @Output
    protected VariantContextWriter vcfWriter = null;

    private final String SOMATIC_LOD_TAG_NAME = "SOMATIC_LOD";
    private final String SOMATIC_AC_TAG_NAME = "SOMATIC_AC";
    private final String SOMATIC_NONREF_TAG_NAME = "SOMATIC_NNR";

    private final Set<String> samples = new HashSet<String>(2);

    /**
     * Parse the familial relationship specification, and initialize VCF writer
     */
    public void initialize() {
        List<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        // set up tumor and normal samples
        if ( !vcfSamples.contains(normalSample) )
            throw new UserException.BadArgumentValue("--normalSample", "the normal sample " + normalSample + " doesn't match any sample from the input VCF");
        if ( !vcfSamples.contains(tumorSample) )
            throw new UserException.BadArgumentValue("--tumorSample", "the tumor sample " + tumorSample + " doesn't match any sample from the input VCF");

        logger.info("Normal sample: " + normalSample);
        logger.info("Tumor  sample: " + tumorSample);

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(GATKVCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(new VCFInfoHeaderLine(VCFConstants.SOMATIC_KEY, 0, VCFHeaderLineType.Flag, "Is this a confidently called somatic mutation"));
        headerLines.add(new VCFInfoHeaderLine(SOMATIC_LOD_TAG_NAME, 1, VCFHeaderLineType.Float, "log10 probability that the site is a somatic mutation"));
        headerLines.add(new VCFInfoHeaderLine(SOMATIC_AC_TAG_NAME, 1, VCFHeaderLineType.Integer, "Allele count of samples with somatic event"));
        headerLines.add(new VCFInfoHeaderLine(SOMATIC_NONREF_TAG_NAME, 1, VCFHeaderLineType.Integer, "Number of samples with somatic event"));

        samples.add(normalSample);
        samples.add(tumorSample);
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    private double log10pNonRefInSamples(final VariantContext vc, final String sample) {
        return log10PLFromSamples(vc, sample, false);
     }

    private double log10pRefInSamples(final VariantContext vc, final String sample) {
        return log10PLFromSamples(vc, sample, true);
    }

    private double log10PLFromSamples(final VariantContext vc, final String sample, boolean calcRefP) {

        Genotype g = vc.getGenotype(sample);
        double log10pSample = -1000;
        if ( ! g.isNoCall() ) {
            final double[] gLikelihoods = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
            log10pSample = Math.log10(calcRefP ? gLikelihoods[0] : 1 - gLikelihoods[0]);
            log10pSample = Double.isInfinite(log10pSample) ? -10000 : log10pSample;
        }
        return log10pSample;
    }

    private int calculateTumorAC(final VariantContext vc) {
        int ac = 0;
        switch ( vc.getGenotype(tumorSample).getType() ) {
            case HET:       ac += 1; break;
            case HOM_VAR:   ac += 2; break;
            case NO_CALL: case UNAVAILABLE: case HOM_REF: break;
        }
        return ac;
    }

    private int calculateTumorNNR(final VariantContext vc) {
        int nnr = 0;
        switch ( vc.getGenotype(tumorSample).getType() ) {
            case HET: case HOM_VAR: nnr += 1; break;
            case NO_CALL: case UNAVAILABLE: case HOM_REF: break;
        }
        return nnr;
    }

    /**
     * P(somatic | D)
     *   = P(somatic) * P(D | somatic)
     *   = P(somatic) * P(D | normals are ref) * P(D | tumors are non-ref)
     *
     * P(! somatic | D)
     *   = P(! somatic) * P(D | ! somatic)
     *   = P(! somatic) *
     *      * (  P(D | normals are non-ref) * P(D | tumors are non-ref) [germline]
     *         + P(D | normals are ref) * P(D | tumors are ref)) [no-variant at all]
     *
     * @param vc
     * @return
     */
    private double calcLog10pSomatic(final VariantContext vc) {
        // walk over tumors
        double log10pNonRefInTumors = log10pNonRefInSamples(vc, tumorSample);
        double log10pRefInTumors = log10pRefInSamples(vc, tumorSample);

        // walk over normals
        double log10pNonRefInNormals = log10pNonRefInSamples(vc, normalSample);
        double log10pRefInNormals = log10pRefInSamples(vc, normalSample);

        // priors
        double log10pSomaticPrior = QualityUtils.qualToErrorProbLog10(somaticPriorQ);
        double log10pNotSomaticPrior = Math.log10(1 - QualityUtils.qualToErrorProb(somaticPriorQ));

        double log10pNotSomaticGermline = log10pNonRefInNormals + log10pNonRefInTumors;
        double log10pNotSomaticNoVariant = log10pRefInNormals + log10pRefInTumors;

        double log10pNotSomatic = log10pNotSomaticPrior + MathUtils.log10sumLog10(new double[]{log10pNotSomaticGermline, log10pNotSomaticNoVariant});
        double log10pSomatic = log10pSomaticPrior + log10pNonRefInTumors + log10pRefInNormals;
        double lod = log10pSomatic - log10pNotSomatic;

        return Double.isInfinite(lod) ? -10000 : lod;
    }

    /**
     * For each variant in the file, determine the phasing for the child and replace the child's genotype with the trio's genotype
     *
     * @param tracker  the reference meta-data tracker
     * @param ref      the reference context
     * @param context  the alignment context
     * @return null
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            for ( VariantContext vc : tracker.getValues(variantCollection.variants, context.getLocation()) ) {
                vc = vc.subContextFromSamples(samples);
                if ( !vc.isPolymorphicInSamples() )
                    continue;

                double log10pSomatic = calcLog10pSomatic(vc);

                // write in the somatic status probability
                Map<String, Object> attrs = new HashMap<String, Object>(); // vc.getAttributes());
                if ( ! minimalVCF ) attrs.putAll(vc.getAttributes());
                attrs.put(SOMATIC_LOD_TAG_NAME, log10pSomatic);
                if ( log10pSomatic > somaticMinLOD ) {
                    attrs.put(VCFConstants.SOMATIC_KEY, true);
                    attrs.put(SOMATIC_NONREF_TAG_NAME, calculateTumorNNR(vc));
                    attrs.put(SOMATIC_AC_TAG_NAME, calculateTumorAC(vc));

                }
                final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(attrs);
                VariantContextUtils.calculateChromosomeCounts(builder, false);
                VariantContext newvc = builder.make();

                vcfWriter.add(newvc);
            }

            return null;
        }

        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    @Override
    public Integer treeReduce(Integer sum1, Integer sum2) {
        return reduce(sum1, sum2);
    }
}
