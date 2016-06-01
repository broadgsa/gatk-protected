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

package org.broadinstitute.gatk.tools.walkers.validation.validationsiteselector;

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.util.*;


/**
 * Randomly select variant records according to specified options
 *
 * <p>
 * This tool is intended for use in experiments where we sample data randomly from a set of variants, for example
 * in order to choose sites for a follow-up validation study.</p>
 *
 * <p>Sites are selected randomly but within certain restrictions. There are two main sources of restrictions:</p>
 * <ul>
 *     <li><b>Sample restrictions:</b> A user can specify a set of samples, and we will only consider sites which are
 *     polymorphic within the given sample subset. These sample restrictions can be given as a set of individual
 *     samples, a text file (each line containing a sample name), or a regular expression. A user can additionally
 *     specify whether samples will be considered based on their genotypes (a non-reference genotype means that the
 *     sample is polymorphic in that variant, and hence that variant will be considered for inclusion in set), or
 *     based on their PLs.</li>
 *     <li><b>Sampling methods:</b>
 *          <ol>
 *              <li>Uniform sampling will just sample uniformly from variants that are polymorphic in selected samples</li>
 *              <li>Sampling based on Allele Frequency spectrum will ensure that output sites have the same AF distribution as the input set</li>
 *          </ol>
 *     </li>
 *     <li>Variant type (SNP, Indel, etc.)</li>
 * </ul>
 *
 * <h3>Input</h3>
 * <p>
 * One or more variant sets to choose from.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A sites-only VCF with the desired number of randomly selected sites.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ValidationSiteSelectorWalker \
 *   -R reference.fasta \
 *   -V input1.vcf \
 *   -V input2.vcf \
 *   -sn NA12878 \
 *   -o output.vcf \
 *   --numValidationSites 200   \
 *   -sampleMode POLY_BASED_ON_GT \
 *   -freqMode KEEP_AF_SPECTRUM
 * </pre>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ValidationSiteSelectorWalker \
 *   -R reference.fasta \
 *   -V:foo input1.vcf \
 *   -V:bar input2.vcf \
 *   --numValidationSites 200 \
 *   -sf samples.txt \
 *   -o output.vcf \
 *   -sampleMode  POLY_BASED_ON_GT \
  *   -freqMode UNIFORM \
 *   -selectType INDEL
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class ValidationSiteSelector extends RodWalker<Integer, Integer> {

    public enum AF_COMPUTATION_MODE {
        KEEP_AF_SPECTRUM,
        UNIFORM
    }

    public enum SAMPLE_SELECTION_MODE {
        NONE,
        POLY_BASED_ON_GT,
        POLY_BASED_ON_GL
    }

    /**
     * The input VCF file
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file, can be specified multiple times", required=true)
    public List<RodBinding<VariantContext>> variants;

    /**
     * The output VCF file
     */
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    /**
     * Sample name(s) to subset the input VCF to, prior to selecting variants. -sn A -sn B subsets to samples A and B.
     */
    @Argument(fullName="sample_name", shortName="sn", doc="Include genotypes from this sample. Can be specified multiple times", required=false)
    public Set<String> sampleNames = new HashSet<String>(0);

    /**
     * Sample regexps to subset the input VCF to, prior to selecting variants. -sn NA12* subsets to all samples with prefix NA12
     */
    @Argument(fullName="sample_expressions", shortName="se", doc="Regular expression to select many samples from the ROD tracks provided. Can be specified multiple times", required=false)
    public Set<String> sampleExpressions ;

    /**
     * File containing a list of sample names to subset the input vcf to. Equivalent to specifying the contents of the file separately with -sn
     */
    @Input(fullName="sample_file", shortName="sf", doc="File containing a list of samples (one per line) to include. Can be specified multiple times", required=false)
    public Set<File> sampleFiles;

    /**
     * A mode for selecting sites based on sample-level data. See the wiki documentation for more information.
     */
    @Argument(fullName="sampleMode", shortName="sampleMode", doc="Sample selection mode", required=false)
    private SAMPLE_SELECTION_MODE sampleMode = SAMPLE_SELECTION_MODE.NONE;

    /**
     * An P[nonref] threshold for SAMPLE_SELECTION_MODE=POLY_BASED_ON_GL. See the wiki documentation for more information.
     */
    @Argument(shortName="samplePNonref",fullName="samplePNonref", doc="GL-based selection mode only: the probability" +
            " that a site is non-reference in the samples for which to include the site",required=false)
    private double samplePNonref = 0.99;

    /**
     * The number of sites in your validation set
     */
    @Argument(fullName="numValidationSites", shortName="numSites", doc="Number of output validation sites", required=true)
    private int numValidationSites;

    /**
     * Do not exclude filtered sites (e.g. not PASS or .) from consideration for validation
     */
    @Argument(fullName="includeFilteredSites", shortName="ifs", doc="If true, will include filtered sites in set to choose variants from", required=false)
    private boolean INCLUDE_FILTERED_SITES = false;

    /**
     * Argument for the frequency selection mode. (AC/AF/AN) are taken from VCF info field, not recalculated. Typically specified for sites-only VCFs that still have AC/AF/AN information.
     */
    @Argument(fullName="ignoreGenotypes", shortName="ignoreGenotypes", doc="If true, will ignore genotypes in VCF, will take AC,AF from annotations and will make no sample selection", required=false)
    private boolean IGNORE_GENOTYPES = false;

    /**
     * Argument for the frequency selection mode. Allows reference (non-polymorphic) sites to be included in the validation set.
     */
    @Argument(fullName="ignorePolymorphicStatus", shortName="ignorePolymorphicStatus", doc="If true, will ignore polymorphic status in VCF, and will take VCF record directly without pre-selection", required=false)
    private boolean IGNORE_POLYMORPHIC = false;

    @Hidden
    @Argument(fullName="numFrequencyBins", shortName="numBins", doc="Number of frequency bins if we're to match AF distribution", required=false)
    private int numFrequencyBins = 20;

    /**
      * This argument selects allele frequency selection mode. See the wiki for more information.
      */
    @Argument(fullName="frequencySelectionMode", shortName="freqMode", doc="Allele Frequency selection mode", required=false)
    private AF_COMPUTATION_MODE freqMode = AF_COMPUTATION_MODE.KEEP_AF_SPECTRUM;

    /**
      * This argument selects particular kinds of variants (i.e. SNP, INDEL) out of a list. If left unspecified, all types are considered.
      */
     @Argument(fullName="selectTypeToInclude", shortName="selectType", doc="Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times", required=false)
     private List<VariantContext.Type> TYPES_TO_INCLUDE = new ArrayList<VariantContext.Type>();


    private TreeSet<String> samples = new TreeSet<String>();
    SampleSelector sampleSelector = null;
    FrequencyModeSelector frequencyModeSelector = null;
    private ArrayList<VariantContext.Type> selectedTypes = new ArrayList<VariantContext.Type>();

    public void initialize() {
         // Get list of samples to include in the output
         Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit());
         TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

         Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
         Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

         // first, add any requested samples
         samples.addAll(samplesFromFile);
         samples.addAll(samplesFromExpressions);
         samples.addAll(sampleNames);

         // if none were requested, we want all of them
         if ( samples.isEmpty() ) {
             samples.addAll(vcfSamples);

         }

         sampleSelector = getSampleSelectorObject(sampleMode, samples);

        // initialize frequency mode selector
        frequencyModeSelector = getFrequencyModeSelectorObject(freqMode, getToolkit().getGenomeLocParser());

        // if user specified types to include, add these, otherwise, add all possible variant context types to list of vc types to include
        if (TYPES_TO_INCLUDE.isEmpty()) {

            for (VariantContext.Type t : VariantContext.Type.values())
                selectedTypes.add(t);

        }
        else {
            for (VariantContext.Type t : TYPES_TO_INCLUDE)
                selectedTypes.add(t);

        }

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(new VCFHeaderLine("source", "ValidationSiteSelector"));
        vcfWriter.writeHeader(new VCFHeader(headerLines));

    }


    @Override
     public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
         if ( tracker == null )
             return 0;

        Collection<VariantContext> vcs = tracker.getValues(variants, context.getLocation());

         if ( vcs == null || vcs.size() == 0) {
             return 0;
         }


        for (VariantContext vc : vcs) {
            if (!selectedTypes.contains(vc.getType()))
                continue;

            // skip if site isn't polymorphic and if user didn't request to ignore polymorphic status
            if (!vc.isPolymorphicInSamples() && !IGNORE_POLYMORPHIC)
                continue;

            if (!INCLUDE_FILTERED_SITES && vc.filtersWereApplied() && vc.isFiltered())
                continue;


            // does this site pass the criteria for the samples we are interested in?
            boolean passesSampleSelectionCriteria;
            if (samples.isEmpty())
                passesSampleSelectionCriteria = true;
            else
                passesSampleSelectionCriteria = sampleSelector.selectSiteInSamples(vc);

            frequencyModeSelector.logCurrentSiteData(vc,passesSampleSelectionCriteria,IGNORE_GENOTYPES,IGNORE_POLYMORPHIC);
        }
        return 1;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    public void onTraversalDone(Integer result) {
        logger.info("Outputting validation sites...");
        ArrayList<VariantContext> selectedSites = frequencyModeSelector.selectValidationSites(numValidationSites);

        for (VariantContext vc : selectedSites) {
            vcfWriter.add(vc);
        }
        logger.info(result + " records processed.");

    }

    private SampleSelector getSampleSelectorObject(SAMPLE_SELECTION_MODE sampleMode, TreeSet<String> samples) {
        SampleSelector sm;
         switch ( sampleMode ) {
             case POLY_BASED_ON_GL:
                 sm = new GLBasedSampleSelector(samples, Math.log10(1.0-samplePNonref));
                 break;
             case POLY_BASED_ON_GT:
                 sm = new GTBasedSampleSelector(samples);
                 break;
             case NONE:
                 sm = new NullSampleSelector(samples);
                 break;
             default:
                 throw new IllegalArgumentException("Unsupported Sample Selection Mode: " + sampleMode);
         }

         return sm;
    }

    private FrequencyModeSelector getFrequencyModeSelectorObject (AF_COMPUTATION_MODE freqMode, GenomeLocParser parser) {
        FrequencyModeSelector fm;

        switch (freqMode) {
            case KEEP_AF_SPECTRUM:
                fm = new KeepAFSpectrumFrequencySelector(numFrequencyBins, parser);
                break;
            case UNIFORM:
                fm = new UniformSamplingFrequencySelector(parser);
                break;
            default: throw new IllegalArgumentException("Unexpected Frequency Selection Mode: "+ freqMode);

        }
        return fm;
    }
}
