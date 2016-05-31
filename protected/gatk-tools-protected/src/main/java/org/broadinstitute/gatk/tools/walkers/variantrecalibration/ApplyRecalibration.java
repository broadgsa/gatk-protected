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

package org.broadinstitute.gatk.tools.walkers.variantrecalibration;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.PartitionBy;
import org.broadinstitute.gatk.engine.walkers.PartitionType;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Apply a score cutoff to filter variants based on a recalibration table
 *
 * <p>
 * This tool performs the second pass in a two-stage process called VQSR; the first pass is performed by the
 * <a href='https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php'>VariantRecalibrator</a> tool.
 * In brief, the first pass consists of creating a Gaussian mixture model by looking at the distribution of annotation
 * values over a high quality subset of the input call set, and then scoring all input variants according to the model.
 * The second pass consists of filtering variants based on score cutoffs identified in the first pass.
 *</p>
 *
 * <p>
 * Using the tranche file and recalibration table generated by the previous step, the ApplyRecalibration tool looks at each variant's VQSLOD value
 * and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level
 * have their FILTER field annotated with the corresponding tranche level. This will result in a call set that is filtered
 * to the desired level but retains the information necessary to increase sensitivity if needed.</p>
 *
 * <p>To be clear, please note that by "filtered", we mean that variants failing the requested tranche cutoff are <b>marked
 * as filtered</b> in the output VCF; they are <b>not discarded</b>.</p>
 *
 * <p>VQSR is probably the hardest part of the Best Practices to get right, so be sure to read the
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=39'>method documentation</a>,
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=1259'>parameter recommendations</a> and
 * <a href='https://www.broadinstitute.org/gatk/guide/article?id=2805'>tutorial</a> to really understand what these
 * tools and how to use them for best results on your own data.</p>
 *
 * <h3>Input</h3>
 * <ul>
 * <li>The raw input variants to be filtered.</li>
 * <li>The recalibration table file that was generated by the VariantRecalibrator tool.</li>
 * <li>The tranches file that was generated by the VariantRecalibrator tool.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 * <li>A recalibrated VCF file in which each variant of the requested type is annotated with its VQSLOD and marked as filtered if the score is below the desired quality level.</li>
 * </ul>
 *
 * <h3>Usage example for filtering SNPs</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ApplyRecalibration \
 *   -R reference.fasta \
 *   -input raw_variants.vcf \
 *   --ts_filter_level 99.0 \
 *   -tranchesFile output.tranches \
 *   -recalFile output.recal \
 *   -mode SNP \
 *   -o path/to/output.recalibrated.filtered.vcf
 * </pre>
 *
 * <h3>Allele-specific usage</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ApplyRecalibration \
 *   -R reference.fasta \
 *   -input raw_variants.withASannotations.vcf \
 *   -AS \
 *   --ts_filter_level 99.0 \
 *   -tranchesFile output.AS.tranches \
 *   -recalFile output.AS.recal \
 *   -mode SNP \
 *   -o path/to/output.recalibrated.ASfiltered.vcf
 * </pre>
 * Each allele will be annotated by its corresponding entry in the AS_FilterStatus INFO field annotation.  Allele-specific VQSLOD and culprit are also carried through from VariantRecalibrator and stored in the AS_VQSLOD and AS_culprit INFO fields, respectively.
 * The site-level filter is set to the most lenient of any of the allele filters.  That is, if one allele passes, the whole site will be PASS.  If no alleles pass, the site-level filter will be set to the lowest sensitivity tranche among all the alleles.
 *
 * Note that the .tranches and .recal files should be derived from an allele-specific run of VariantRecalibrator
 * Also note that the AS_culprit, AS_FilterStatus, and AS_VQSLOD fields will have placeholder values (NA or NaN) for alleles of a type that have not yet been processed by ApplyRecalibration
 * The spanning deletion allele (*) will not be recalibrated because it represents missing data. Its VQSLOD will remain NaN and it's culprit and FilterStatus will be NA.
 *
 * <h3>Caveats</h3>
 *
 * <ul>
 * <li>The tranche values used in the example above are only meant to be a general example. You should determine the level of sensitivity
 * that is appropriate for your specific project. Remember that higher sensitivity (more power to detect variants, yay!) comes
 * at the cost of specificity (more false negatives, boo!). You have to choose at what point you want to set the tradeoff.</li>
 * <li>In order to create the tranche reporting plots (which are only generated for SNPs, not indels!) Rscript needs to be
 * in your environment PATH (this is the scripting version of R, not the interactive version).</li>
 * </ul>
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
public class ApplyRecalibration extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

    public static final String LOW_VQSLOD_FILTER_NAME = "LOW_VQSLOD";
    private final double DEFAULT_VQSLOD_CUTOFF = 0.0;

    boolean foundSNPTranches = false;
    boolean foundINDELTranches = false;

    /////////////////////////////
    // Inputs
    /////////////////////////////
    /**
     * These calls should be unfiltered and annotated with the error covariates that are intended to use for modeling.
     */
    @Input(fullName="input", shortName = "input", doc="The raw input variants to be recalibrated", required=true)
    public List<RodBinding<VariantContext>> input;
    @Input(fullName="recal_file", shortName="recalFile", doc="The input recal file used by ApplyRecalibration", required=true)
    protected RodBinding<VariantContext> recal;
    @Input(fullName="tranches_file", shortName="tranchesFile", doc="The input tranches file describing where to cut the data", required=false)
    protected File TRANCHES_FILE;

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Output( doc="The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value")
    private VariantContextWriter vcfWriter = null;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="ts_filter_level", shortName="ts_filter_level", doc="The truth sensitivity level at which to start filtering", required=false)
    protected Double TS_FILTER_LEVEL = null;

    /**
     *  Filter the input file based on allele-specific recalibration data.  See tool docs for site-level and allele-level filtering details.
     *  Requires a .recal file produced using an allele-specific run of VariantRecalibrator
     */
    @Argument(fullName="useAlleleSpecificAnnotations", shortName="AS", doc="If specified, the tool will attempt to apply a filter to each allele based on the input tranches and allele-specific .recal file.", required=false)
    public boolean useASannotations = false;

    @Advanced
    @Argument(fullName="lodCutoff", shortName="lodCutoff", doc="The VQSLOD score below which to start filtering", required=false)
    protected Double VQSLOD_CUTOFF = null;

    /**
     * For this to work properly, the -ignoreFilter argument should also be applied to the VariantRecalibration command.
     */
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified, the recalibration will be applied to variants marked as filtered by the specified filter name in the input VCF file", required=false)
    private String[] IGNORE_INPUT_FILTERS = null;
    @Argument(fullName="ignore_all_filters", shortName="ignoreAllFilters", doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.", required=false)
    private boolean IGNORE_ALL_FILTERS = false;
    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't output filtered loci after applying the recalibration", required=false)
    protected boolean EXCLUDE_FILTERED = false;
    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ: 1.) SNP for recalibrating only SNPs (emitting indels untouched in the output VCF); 2.) INDEL for indels; and 3.) BOTH for recalibrating both SNPs and indels simultaneously.", required = false)
    public VariantRecalibratorArgumentCollection.Mode MODE = VariantRecalibratorArgumentCollection.Mode.SNP;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    final private List<Tranche> tranches = new ArrayList<>();
    final private Set<String> inputNames = new HashSet<>();
    final private Set<String> ignoreInputFilterSet = new TreeSet<>();
    final static private String listPrintSeparator = ",";
    final static private String trancheFilterString = "VQSRTranche";
    final static private String arrayParseRegex = "[\\[\\]\\s]";
    final static private String emptyStringValue = "NA";
    final static private String emptyFloatValue = "NaN";


    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if( TS_FILTER_LEVEL != null ) {
            for ( final Tranche t : Tranche.readTranches(TRANCHES_FILE) ) {
                if ( t.ts >= TS_FILTER_LEVEL ) {
                    tranches.add(t);
                }
                logger.info(String.format("Read tranche " + t));
            }
            Collections.reverse(tranches); // this algorithm wants the tranches ordered from best (lowest truth sensitivity) to worst (highest truth sensitivity)
        }

        for( final RodBinding rod : input ) {
            inputNames.add( rod.getName() );
        }

        if( IGNORE_INPUT_FILTERS != null ) {
            ignoreInputFilterSet.addAll( Arrays.asList(IGNORE_INPUT_FILTERS) );
        }

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<>();
        final Set<VCFHeaderLine> inputHeaders = GATKVCFUtils.getHeaderFields(getToolkit(), inputNames);
        hInfo.addAll(inputHeaders);
        addVQSRStandardHeaderLines(hInfo);
        if (useASannotations)
            addAlleleSpecificVQSRHeaderLines(hInfo);

        checkForPreviousApplyRecalRun(Collections.unmodifiableSet(inputHeaders));

        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames));

        //generate headers from tranches file
        //TODO: throw away old tranche headers if we're ignoring filters
        if( TS_FILTER_LEVEL != null ) {
            // if the user specifies both ts_filter_level and lodCutoff then throw a user error
            if( VQSLOD_CUTOFF != null ) {
                throw new UserException("Arguments --ts_filter_level and --lodCutoff are mutually exclusive. Please only specify one option.");
            }

            if( tranches.size() >= 2 ) {
                for( int iii = 0; iii < tranches.size() - 1; iii++ ) {
                    final Tranche t = tranches.get(iii);
                    hInfo.add(new VCFFilterHeaderLine(t.name, String.format("Truth sensitivity tranche level for " + t.model.toString() + " model at VQS Lod: " + t.minVQSLod + " <= x < " + tranches.get(iii+1).minVQSLod)));
                }
            }
            if( tranches.size() >= 1 ) {
                hInfo.add(new VCFFilterHeaderLine(tranches.get(0).name + "+", String.format("Truth sensitivity tranche level for " + tranches.get(0).model.toString() + " model at VQS Lod < " + tranches.get(0).minVQSLod)));
            } else {
                throw new UserException("No tranches were found in the file or were above the truth sensitivity filter level " + TS_FILTER_LEVEL);
            }

            logger.info("Keeping all variants in tranche " + tranches.get(tranches.size()-1));
        } else {
            if( VQSLOD_CUTOFF == null ) {
                VQSLOD_CUTOFF = DEFAULT_VQSLOD_CUTOFF;
            }
            hInfo.add(new VCFFilterHeaderLine(LOW_VQSLOD_FILTER_NAME, "VQSLOD < " + VQSLOD_CUTOFF));
            logger.info("Keeping all variants with VQSLOD >= " + VQSLOD_CUTOFF);
        }

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    private boolean trancheIntervalIsValid(final String sensitivityLimits) {
        final String[] vals = sensitivityLimits.split("to");
        if(vals.length != 2)
            return false;
        try {
            double lowerLimit = Double.parseDouble(vals[0]);
            double upperLimit = Double.parseDouble(vals[1].replace("+",""));    //why does our last tranche end with 100+? Is there anything greater than 100 percent?  Really???
        }
        catch(NumberFormatException e) {
            throw new UserException("Poorly formatted tranche filter name does not contain two sensitivity interval end points.");
        }
        return true;
    }

    public static void addVQSRStandardHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VQS_LOD_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CULPRIT_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NEGATIVE_LABEL_KEY));
    }

    public static void addAlleleSpecificVQSRHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_FILTER_STATUS_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_CULPRIT_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VQS_LOD_KEY));
    }

    /**
     * Check the filter declarations in the input VCF header to see if any ApplyRecalibration mode has been run
     * Here we assume that the tranches are named with a specific format: VQSRTranche[SNP|INDEL][lowerLimit]to[upperLimit]
     * @param inputHeaders
     */
    private void checkForPreviousApplyRecalRun(final Set<VCFHeaderLine> inputHeaders) {
        for(final VCFHeaderLine header : inputHeaders) {
            if(header instanceof VCFFilterHeaderLine) {
                final String filterName = ((VCFFilterHeaderLine)header).getID();
                if(filterName.length() < 12 || !filterName.substring(0, 11).equalsIgnoreCase(trancheFilterString)) {
                    continue;
                }
                if(filterName.charAt(11) == 'S') {
                    //for SNP tranches, get sensitivity limit
                    final String sensitivityLimits = filterName.substring(14);
                    if(trancheIntervalIsValid(sensitivityLimits))
                        foundSNPTranches = true;
                }
                else if(filterName.charAt(11) == 'I') {
                    //for INDEL tranches, get sensitivity limit
                    final String sensitivityLimits = filterName.substring(16);
                    if(trancheIntervalIsValid(sensitivityLimits))
                        foundINDELTranches = true;
                }
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return 1;
        }

        final List<VariantContext> VCs =  tracker.getValues(input, context.getLocation());
        final List<VariantContext> recals =  tracker.getValues(recal, context.getLocation());

        for( final VariantContext vc : VCs ) {

            final boolean evaluateThisVariant = useASannotations || VariantDataManager.checkVariationClass( vc, MODE );
            final boolean variantIsNotFiltered = IGNORE_ALL_FILTERS || vc.isNotFiltered() || (!ignoreInputFilterSet.isEmpty() && ignoreInputFilterSet.containsAll(vc.getFilters()));     //vc.isNotFiltered is true for PASS; vc.filtersHaveBeenApplied covers PASS and filters
            if( evaluateThisVariant && variantIsNotFiltered) {


                String filterString;
                final VariantContextBuilder builder = new VariantContextBuilder(vc);
                if (!useASannotations) {
                    filterString = doSiteSpecificFiltering(vc, recals, builder);
                }
                else {  //allele-specific mode
                    filterString = doAlleleSpecificFiltering(vc, recals, builder);
                }

                //for both non-AS and AS modes:

                if( filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                    builder.passFilters();
                } else if(filterString.equals(VCFConstants.UNFILTERED)) {
                    builder.unfiltered();
                } else {
                    builder.filters(filterString);
                }

                final VariantContext outputVC = builder.make();
                if( !EXCLUDE_FILTERED || outputVC.isNotFiltered() ) {
                    vcfWriter.add( outputVC );
                }
            } else { // valid VC but not compatible with this mode, so just emit the variant untouched
                vcfWriter.add( vc );
            }
        }

        return 1; // This value isn't used for anything
    }

    public double parseFilterLowerLimit(final String trancheFilter) {
        final Pattern pattern = Pattern.compile("VQSRTranche\\S+(\\d+\\.\\d+)to(\\d+\\.\\d+)");
        final Matcher m = pattern.matcher(trancheFilter);
        return m.find() ? Double.parseDouble(m.group(1)) : -1;
    }

    /**
     * Generate the VCF filter string for this record based on the ApplyRecalibration modes run so far
     * @param vc the input VariantContext (with at least one ApplyRecalibration mode already run)
     * @param bestLod best LOD from the alleles we've seen in this recalibration mode
     * @return the String to use as the VCF filter field
     */
    protected String generateFilterStringFromAlleles(final VariantContext vc, final double bestLod) {
        String filterString = ".";

        final boolean bothModesWereRun = (MODE == VariantRecalibratorArgumentCollection.Mode.SNP && foundINDELTranches) || (MODE == VariantRecalibratorArgumentCollection.Mode.INDEL && foundSNPTranches);
        final boolean onlyOneModeNeeded = !vc.isMixed() && VariantDataManager.checkVariationClass( vc, MODE );

        //if both SNP and INDEL modes have not yet been run (and need to be), leave this variant as unfiltered and add the filters for the alleles in this mode to the INFO field
        if (!bothModesWereRun && !onlyOneModeNeeded) {
            return VCFConstants.UNFILTERED;
        }

        //if both SNP and INDEL modes have been run or the site is not mixed, generate a filter string for this site based on both models

        //pull out the allele filter status from the info field (there may be more than one entry in the list if there were multiple snp/indel alleles assessed in the other mode)
        final String prevFilterStatus = vc.getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, null);

        //if this site hasn't had a filter applied yet
        if (prevFilterStatus != null && !prevFilterStatus.equals(VCFConstants.UNFILTERED)) {
            final String prevAllelesFilterStatusString = vc.getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, null);
            final String[] prevAllelesFilterStatusList = prevAllelesFilterStatusString.split(listPrintSeparator);
            //start with the current best allele filter as the most lenient filter across all modes and all alleles
            String mostLenientFilterName = generateFilterString(bestLod);
            //if the current mode's best allele passes the tranche filter, then let the whole site pass
            if (mostLenientFilterName.equals(VCFConstants.PASSES_FILTERS_v4)) {
                filterString = mostLenientFilterName;
            }
            //if the current mode's best allele does not pass the tranche filter, compare the most lenient filter of this mode with those from the previous mode
            else {
                double mostLenientSensitivityLowerLimit = parseFilterLowerLimit(mostLenientFilterName);
                for (int i = 0; i < prevAllelesFilterStatusList.length; i++) {
                    final String alleleFilterString = prevAllelesFilterStatusList[i].replaceAll(arrayParseRegex, "").trim();
                    //if any allele from the previous mode passed the tranche filter, then let the whole site pass
                    if (alleleFilterString.equals(VCFConstants.PASSES_FILTERS_v4)) { //this allele is PASS
                        mostLenientFilterName = alleleFilterString;
                        break;
                    }
                    //if there's no PASS, then we need to parse the filters to find out how lenient they are
                    else {
                        final double alleleLowerLimit = parseFilterLowerLimit(alleleFilterString);
                        if (alleleLowerLimit == -1)
                            continue;
                        if (alleleLowerLimit < mostLenientSensitivityLowerLimit) {
                            mostLenientSensitivityLowerLimit = alleleLowerLimit;
                            mostLenientFilterName = alleleFilterString;
                        }
                    }

                }
                filterString = mostLenientFilterName;
            }
        }

        //if both modes have been run, but the previous mode didn't apply a filter, use the current mode's best allele VQSLOD filter (shouldn't get run, but just in case)
        else {
            filterString = generateFilterString(bestLod);
        }

        return filterString;
    }

    /**
     * Generate the VCF filter string for this record based on the provided lod score
     * @param lod non-null double
     * @return the String to use as the VCF filter field
     */
    protected String generateFilterString( final double lod ) {
        String filterString = null;
        if( TS_FILTER_LEVEL != null ) {
            for( int i = tranches.size() - 1; i >= 0; i-- ) {
                final Tranche tranche = tranches.get(i);
                if( lod >= tranche.minVQSLod ) {
                    if( i == tranches.size() - 1 ) {
                        filterString = VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        filterString = tranche.name;
                    }
                    break;
                }
            }

            if( filterString == null ) {
                filterString = tranches.get(0).name+"+";
            }
        } else {
            filterString = ( lod < VQSLOD_CUTOFF ? LOW_VQSLOD_FILTER_NAME : VCFConstants.PASSES_FILTERS_v4 );
        }

        return filterString;
    }

    private VariantContext getMatchingRecalVC(final VariantContext target, final List<VariantContext> recalVCs, final Allele allele) {
        for( final VariantContext recalVC : recalVCs ) {
            if ( target.getEnd() == recalVC.getEnd() ) {
                if (!useASannotations)
                    return recalVC;
                else if (allele.equals(recalVC.getAlternateAllele(0)))
                    return recalVC;
            }
        }
        return null;
    }

    /**
     *
     * @param altIndex current alt allele
     * @param prevCulpritList culprits from previous ApplyRecalibration run
     * @param prevLodList lods from previous ApplyRecalibration run
     * @param prevASfiltersList AS_filters from previous ApplyRecalibration run
     * @param culpritString
     * @param lodString
     * @param AS_filterString
     */
    private void updateAnnotationsWithoutRecalibrating(final int altIndex, final String[] prevCulpritList, final String[] prevLodList, final String[] prevASfiltersList,
                                                       final List<String> culpritString, final List<String> lodString, final List<String> AS_filterString) {
        if (foundINDELTranches || foundSNPTranches) {
            if (altIndex < prevCulpritList.length) {
                culpritString.add(prevCulpritList[altIndex].replaceAll(arrayParseRegex, "").trim());
                lodString.add(prevLodList[altIndex].replaceAll(arrayParseRegex, "").trim());
                AS_filterString.add(prevASfiltersList[altIndex].replaceAll(arrayParseRegex, "").trim());
            }
        } else { //if the other allele type hasn't been processed yet, make sure there are enough entries
            culpritString.add(emptyStringValue);
            lodString.add(emptyFloatValue);
            AS_filterString.add(emptyStringValue);
        }
    }

    /**
     * Calculate the allele-specific filter status of vc
     * @param vc
     * @param recals
     * @param builder   is modified by adding attributes
     * @return a String with the filter status for this site
     */
    private String doAlleleSpecificFiltering(final VariantContext vc, final List<VariantContext> recals, final VariantContextBuilder builder) {
        double bestLod = VariantRecalibratorEngine.MIN_ACCEPTABLE_LOD_SCORE;
        final List<String> culpritStrings = new ArrayList<>();
        final List<String> lodStrings = new ArrayList<>();
        final List<String> AS_filterStrings = new ArrayList<>();

        String[] prevCulpritList = null;
        String[] prevLodList = null;
        String[] prevASfiltersList = null;

        //get VQSR annotations from previous run of ApplyRecalibration, if applicable
        if(foundINDELTranches || foundSNPTranches) {
            final String prevCulprits = vc.getAttributeAsString(GATKVCFConstants.AS_CULPRIT_KEY,"");
            prevCulpritList = prevCulprits.isEmpty()? new String[0] : prevCulprits.split(listPrintSeparator);
            final String prevLodString = vc.getAttributeAsString(GATKVCFConstants.AS_VQS_LOD_KEY,"");
            prevLodList = prevLodString.isEmpty()? new String[0] : prevLodString.split(listPrintSeparator);
            final String prevASfilters = vc.getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY,"");
            prevASfiltersList = prevASfilters.isEmpty()? new String[0] : prevASfilters.split(listPrintSeparator);
        }

        //for each allele in the current VariantContext
        for (int altIndex = 0; altIndex < vc.getNAlleles()-1; altIndex++) {
            final Allele allele = vc.getAlternateAllele(altIndex);

            //if the current allele is not part of this recalibration mode, add its annotations to the list and go to the next allele
            if (!VariantDataManager.checkVariationClass(vc, allele, MODE)) {
                updateAnnotationsWithoutRecalibrating(altIndex, prevCulpritList, prevLodList, prevASfiltersList, culpritStrings, lodStrings, AS_filterStrings);
                continue;
            }

            //if the current allele does need to have recalibration applied...

            //initialize allele-specific VQSR annotation data with values for spanning deletion
            String alleleLodString = emptyFloatValue;
            String alleleFilterString = emptyStringValue;
            String alleleCulpritString = emptyStringValue;

            //if it's not a spanning deletion, replace those allele strings with the real values
            if (!allele.equals(Allele.SPAN_DEL)) {
                VariantContext recalDatum = getMatchingRecalVC(vc, recals, allele);
                if (recalDatum == null) {
                    throw new UserException("Encountered input allele which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyRecalibration were run on the same set of input variants with flag -AS. First seen at: " + vc);
                }

                //compare VQSLODs for all alleles in the current mode for filtering later
                final double lod = recalDatum.getAttributeAsDouble(GATKVCFConstants.VQS_LOD_KEY, VariantRecalibratorEngine.MIN_ACCEPTABLE_LOD_SCORE);
                if (lod > bestLod)
                    bestLod = lod;

                alleleLodString = String.format("%.4f", lod);
                alleleFilterString = generateFilterString(lod);
                alleleCulpritString = recalDatum.getAttributeAsString(GATKVCFConstants.CULPRIT_KEY, ".");

                if(recalDatum != null) {
                    if (recalDatum.hasAttribute(GATKVCFConstants.POSITIVE_LABEL_KEY))
                        builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
                    if (recalDatum.hasAttribute(GATKVCFConstants.NEGATIVE_LABEL_KEY))
                        builder.attribute(GATKVCFConstants.NEGATIVE_LABEL_KEY, true);
                }
            }

            //append per-allele VQSR annotations
            lodStrings.add(alleleLodString);
            AS_filterStrings.add(alleleFilterString);
            culpritStrings.add(alleleCulpritString);
        }

        // Annotate the new record with its VQSLOD, AS_FilterStatus, and the worst performing annotation
        if(!AS_filterStrings.isEmpty() )
            builder.attribute(GATKVCFConstants.AS_FILTER_STATUS_KEY, AnnotationUtils.encodeStringList(AS_filterStrings));
        if(!lodStrings.isEmpty())
            builder.attribute(GATKVCFConstants.AS_VQS_LOD_KEY, AnnotationUtils.encodeStringList(lodStrings));
        if(!culpritStrings.isEmpty())
            builder.attribute(GATKVCFConstants.AS_CULPRIT_KEY, AnnotationUtils.encodeStringList(culpritStrings));

        return generateFilterStringFromAlleles(vc, bestLod);
    }

    /**
     * Calculate the filter status for a given VariantContext using the combined data from all alleles at a site
     * @param vc
     * @param recals
     * @param builder   is modified by adding attributes
     * @return a String with the filter status for this site
     */
    private String doSiteSpecificFiltering(final VariantContext vc, final List<VariantContext> recals, final VariantContextBuilder builder) {
        VariantContext recalDatum = getMatchingRecalVC(vc, recals, null);
        if( recalDatum == null ) {
            throw new UserException("Encountered input variant which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyRecalibration were run on the same set of input variants. First seen at: " + vc );
        }

        final String lodString = recalDatum.getAttributeAsString(GATKVCFConstants.VQS_LOD_KEY, null);
        if( lodString == null ) {
            throw new UserException("Encountered a malformed record in the input recal file. There is no lod for the record at: " + vc );
        }
        final double lod;
        try {
            lod = Double.valueOf(lodString);
        } catch (NumberFormatException e) {
            throw new UserException("Encountered a malformed record in the input recal file. The lod is unreadable for the record at: " + vc );
        }

        builder.attribute(GATKVCFConstants.VQS_LOD_KEY, lod);
        builder.attribute(GATKVCFConstants.CULPRIT_KEY, recalDatum.getAttribute(GATKVCFConstants.CULPRIT_KEY));
        if(recalDatum != null) {
            if (recalDatum.hasAttribute(GATKVCFConstants.POSITIVE_LABEL_KEY))
                builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
            if (recalDatum.hasAttribute(GATKVCFConstants.NEGATIVE_LABEL_KEY))
                builder.attribute(GATKVCFConstants.NEGATIVE_LABEL_KEY, true);
        }

        return generateFilterString(lod);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer reduceInit() {
        return 1; // This value isn't used for anything
    }

    public Integer reduce( final Integer mapValue, final Integer reduceSum ) {
        return 1; // This value isn't used for anything
    }

    public Integer treeReduce( final Integer lhs, final Integer rhs ) {
        return 1; // This value isn't used for anything
    }

    public void onTraversalDone( final Integer reduceSum ) {
    }
}

