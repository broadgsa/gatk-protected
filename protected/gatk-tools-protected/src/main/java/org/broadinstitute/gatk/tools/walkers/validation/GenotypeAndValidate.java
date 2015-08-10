/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.validation;

import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.genotyper.*;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.util.Map;
import java.util.Set;

import static org.broadinstitute.gatk.utils.IndelUtils.isInsideExtendedIndel;

/**
 * Genotype and validate a dataset and the calls of another dataset using the Unified Genotyper
 *
 *  <h4>Note that this is an old tool that makes use of the UnifiedGenotyper, which has since been
 *  deprecated in favor of the HaplotypeCaller.</h4>
 *  <p>
 *     Genotype and Validate is a tool to evaluate the quality of a dataset for calling SNPs
 *     and Indels given a secondary (validation) data source. The data sources are BAM or VCF
 *     files. You can use them interchangeably (i.e. a BAM to validate calls in a VCF or a VCF
 *     to validate calls on a BAM).
 *  </p>
 *
 *  <p>
 *     The simplest scenario is when you have a VCF of hand annotated SNPs and Indels, and you
 *     want to know how well a particular technology performs calling these snps. With a
 *     dataset (BAM file) generated by the technology in test, and the hand annotated VCF, you
 *     can run GenotypeAndValidate to asses the accuracy of the calls with the new technology's
 *     dataset.
 *  </p>
 *
 *  <p>
 *     Another option is to validate the calls on a VCF file, using a deep coverage BAM file
 *     that you trust the calls on. The GenotypeAndValidate walker will make calls using the
 *     reads in the BAM file and take them as truth, then compare to the calls in the VCF file
 *     and produce a truth table.
 *  </p>
 *
 *
 * <h3>Input</h3>
 *  <p>
 *      A BAM file to make calls on and a VCF file to use as truth validation dataset.
 *
 *      You also have the option to invert the roles of the files using the command line options listed below.
 *  </p>
 *
 * <h3>Output</h3>
 *  <p>
 *      GenotypeAndValidate has two outputs. The truth table and the optional VCF file. The truth table is a
 *      2x2 table correlating what was called in the dataset with the truth of the call (whether it's a true
 *      positive or a false positive). The table should look like this:
 *  </p>
 * <center>
 * <table id="description-table">
 *     <tr>
 *         <th></th>
 *         <th>ALT</th>
 *         <th>REF</th>
 *         <th>Predictive Value</th>
 *     </tr>
 *     <tr>
 *         <td><b>called alt</b></td>
 *         <td>True Positive (TP)</td>
 *         <td>False Positive (FP)</td>
 *         <td>Positive PV</td>
 *     </tr>
 *     <tr>
 *         <td><b>called ref</b></td>
 *         <td>False Negative (FN)</td>
 *         <td>True Negative (TN)</td>
 *         <td>Negative PV</td>
 *     </tr>
 * </table>
 * </center>
 *
 *  <p>
 *      The <b>positive predictive value (PPV)</b> is the proportion of subjects with positive test results
 *      who are correctly diagnosed.
 *  </p>
 *  <p>
 *      The <b>negative predictive value (NPV)</b> is the proportion of subjects with a negative test result
 *      who are correctly diagnosed.
 *  </p>
 *  <p>
 *      The VCF file will contain only the variants that were called or not called, excluding the ones that
 *      were uncovered or didn't pass the filters. This file is useful if you are trying to compare
 *      the PPV and NPV of two different technologies on the exact same sites (so you can compare apples to
 *      apples).
 *  </p>
 *
 *  <p>
 *      Here is an example of an annotated VCF file (info field clipped for clarity)
 *
 * <pre>
 * #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  NA12878
 * 1   20568807    .   C   T   0    HapMapHet        AC=1;AF=0.50;AN=2;DP=0;GV=T  GT  0/1
 * 1   22359922    .   T   C   282  WG-CG-HiSeq      AC=2;AF=0.50;GV=T;AN=4;DP=42 GT:AD:DP:GL:GQ  1/0 ./. 0/1:20,22:39:-72.79,-11.75,-67.94:99    ./.
 * 13  102391461   .   G   A   341  Indel;SnpCluster AC=1;GV=F;AF=0.50;AN=2;DP=45 GT:AD:DP:GL:GQ  ./. ./. 0/1:32,13:45:-50.99,-13.56,-112.17:99   ./.
 * 1   175516757   .   C   G   655  SnpCluster,WG    AC=1;AF=0.50;AN=2;GV=F;DP=74 GT:AD:DP:GL:GQ  ./. ./. 0/1:52,22:67:-89.02,-20.20,-191.27:99   ./.
 * </pre>
 *
 *  </p>
 *
 *  <h3>Additional Details</h3>
 *  <ul>
 *      <li>
 *          You should always use -L on your VCF track, so that the GATK only looks at the sites on the VCF file.
 *          This speeds up the process a lot.
 *      </li>
 *      <li>
 *          The total number of visited bases may be greater than the number of variants in the original
 *          VCF file because of extended indels, as they trigger one call per new insertion or deletion.
 *          (i.e. ACTG/- will count as 4 genotyper calls, but it's only one line in the VCF).
 *      </li>
 *  </ul>
 *
 * <h3>Usage examples</h3>
 * <h4>Genotypes BAM file from new technology using the VCF as a truth dataset</h4>
 * <pre>
 *  java
 *      -jar GenomeAnalysisTK.jar \
 *      -T  GenotypeAndValidate \
 *      -R reference.fasta \
 *      -I myNewTechReads.bam \
 *      -alleles handAnnotatedVCF.vcf \
 *      -L handAnnotatedVCF.vcf \
 *      -o output.vcf
 * </pre>
 *
 * <h4>Genotypes BAM file from new technology a BAM file as the truth dataset</h4>
 * <pre>
 *  java
 *      -jar GenomeAnalysisTK.jar \
 *      -T  GenotypeAndValidate \
 *      -R reference.fasta \
 *      -I myTruthDataset.bam \
 *      -alleles callsToValidate.vcf \
 *      -L callsToValidate.vcf \
 *      -bt \
 *      -o output.vcf
 * </pre>
 *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VALIDATION, extraDocs = {CommandLineGATK.class} )
@Requires(value={DataSource.READS, DataSource.REFERENCE})
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@By(DataSource.REFERENCE)
@Reference(window=@Window(start=-200,stop=200))
public class GenotypeAndValidate extends RodWalker<GenotypeAndValidate.CountedData, GenotypeAndValidate.CountedData> implements TreeReducible<GenotypeAndValidate.CountedData> {

    /**
     * The optional output file that will have all the variants used in the Genotype and Validation essay.
     * The new annotation `callStatus` will carry the value called in the validation VCF or BAM file."
     */
    @Output(doc="Output VCF file with annotated variants", required=false)
    protected VariantContextWriter vcfWriter = null;

    /**
     * The callset to be used as truth (default) or validated (if BAM file is set to truth).
     */
    @Input(fullName="alleles", shortName = "alleles", doc="The set of alleles at which to genotype", required=true)
    public RodBinding<VariantContext> alleles;

    /**
     * Makes the Unified Genotyper calls to the BAM file the truth dataset and validates the alleles ROD binding callset.
     */
    @Argument(fullName ="set_bam_truth", shortName ="bt", doc="Use the calls on the reads (bam file) as the truth dataset and validate the calls on the VCF", required=false)
    private boolean bamIsTruth = false;

    /**
     * The minimum base quality score necessary for a base to be considered when calling a genotype. This argument is passed to the Unified Genotyper.
     */
    @Argument(fullName="minimum_base_quality_score", shortName="mbq", doc="Minimum base quality score for calling a genotype", required=false)
    private int mbq = -1;

    /**
     * The maximum deletion fraction allowed in a site for calling a genotype. This argument is passed to the Unified Genotyper.
     */
    @Argument(fullName="maximum_deletion_fraction", shortName="deletions", doc="Maximum deletion fraction for calling a genotype", required=false)
    private double deletions = -1;

    /**
     * the minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls. This argument is passed to the Unified Genotyper.
     */
    @Argument(fullName="standard_min_confidence_threshold_for_calling", shortName="stand_call_conf", doc="the minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls", required=false)
    private double callConf = -1;

    /**
     * the minimum phred-scaled Qscore threshold to emit low confidence calls. This argument is passed to the Unified Genotyper.
     */
    @Argument(fullName="standard_min_confidence_threshold_for_emitting", shortName="stand_emit_conf", doc="the minimum phred-scaled Qscore threshold to emit low confidence calls", required=false)
    private double emitConf = -1;

    /**
     * Only validate sites that have at least a given depth
     */
    @Argument(fullName="condition_on_depth", shortName="depth", doc="Condition validation on a minimum depth of coverage by the reads", required=false)
    private int minDepth = -1;

   /**
     * Print out discordance sites to standard out.
     */
    @Hidden
    @Argument(fullName ="print_interesting_sites", shortName ="print_interesting", doc="Print out interesting sites to standard out", required=false)
    private boolean printInterestingSites = false;

    private UnifiedGenotypingEngine snpEngine;
    private UnifiedGenotypingEngine indelEngine;
    private Set<String> samples;

    private enum GVstatus {
        T, F, NONE
    }

    public static class CountedData {
        private long nAltCalledAlt = 0L;
        private long nAltCalledRef = 0L;
        private long nAltNotCalled = 0L;
        private long nRefCalledAlt = 0L;
        private long nRefCalledRef = 0L;
        private long nRefNotCalled = 0L;
        private long nNoStatusCalledAlt = 0L;
        private long nNoStatusCalledRef = 0L;
        private long nNoStatusNotCalled = 0L;
        private long nNotConfidentCalls = 0L;
        private long nUncovered = 0L;

        /**
         * Adds the values of other to this, returning this
         * @param other the other object
         */
        public void add(CountedData other) {
            nAltCalledAlt += other.nAltCalledAlt;
            nAltCalledRef += other.nAltCalledRef;
            nAltNotCalled += other.nAltNotCalled;
            nRefCalledAlt += other.nRefCalledAlt;
            nRefCalledRef += other.nRefCalledRef;
            nRefNotCalled += other.nRefNotCalled;
            nNoStatusCalledAlt += other.nNoStatusCalledAlt;
            nNoStatusCalledRef += other.nNoStatusCalledRef;
            nNoStatusNotCalled += other.nNoStatusNotCalled;
            nUncovered += other.nUncovered;
            nNotConfidentCalls += other.nNotConfidentCalls;
        }
    }



    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        // Initialize VCF header
        if (vcfWriter != null) {
            Map<String, VCFHeader> header = GATKVCFUtils.getVCFHeadersFromRodPrefix(getToolkit(), alleles.getName());
            samples = SampleUtils.getSampleList(header, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
            Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(header.values(), true);
            headerLines.add(new VCFHeaderLine("source", "GenotypeAndValidate"));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.GENOTYPE_AND_VALIDATE_STATUS_KEY));
            vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
        }

        // Filling in SNP calling arguments for UG
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.outputMode = OutputMode.EMIT_ALL_SITES;
        uac.alleles = alleles;

        // TODO -- if we change this tool to actually validate against the called allele, then this if statement is needed;
        // TODO -- for now, though, we need to be able to validate the right allele (because we only test isVariant below) [EB]
        //if (!bamIsTruth)
        uac.genotypingOutputMode = GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES;

        if (mbq >= 0) uac.MIN_BASE_QUALTY_SCORE = mbq;
        if (deletions >= 0)
            uac.MAX_DELETION_FRACTION = deletions;
        else
            uac.MAX_DELETION_FRACTION = 1.0;
        if (emitConf >= 0) uac.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING = emitConf;
        if (callConf >= 0) uac.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = callConf;

        final GenomeAnalysisEngine toolkit = getToolkit();
        uac.GLmodel = GenotypeLikelihoodsCalculationModel.Model.SNP;
        snpEngine = new UnifiedGenotypingEngine(uac,
                FixedAFCalculatorProvider.createThreadSafeProvider(toolkit, uac, logger),toolkit);


        // Adding the INDEL calling arguments for UG
        UnifiedArgumentCollection uac_indel = uac.clone();
        uac_indel.GLmodel = GenotypeLikelihoodsCalculationModel.Model.INDEL;
        indelEngine = new UnifiedGenotypingEngine(uac_indel,
                FixedAFCalculatorProvider.createThreadSafeProvider(toolkit, uac, logger),toolkit);

        // make sure we have callConf set to the threshold set by the UAC so we can use it later.
        callConf = uac.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public CountedData map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        final CountedData counter = new CountedData();

        // For some reason RodWalkers get map calls with null trackers
        if( tracker == null )
            return counter;

        VariantContext vcComp = tracker.getFirstValue(alleles);
        if( vcComp == null )
            return counter;

        //todo - not sure I want this, may be misleading to filter extended indel events.
        if (isInsideExtendedIndel(vcComp,  ref))
            return counter;

        // Do not operate on variants that are not covered to the optional minimum depth
        if (!context.hasReads() || (minDepth > 0 && context.getBasePileup().getBases().length < minDepth)) {
            counter.nUncovered = 1L;
            final GVstatus status = getGVstatus(vcComp);
            if ( status == GVstatus.T )
                counter.nAltNotCalled = 1L;
            else if ( status == GVstatus.F )
                counter.nRefNotCalled = 1L;
            else
                counter.nNoStatusNotCalled = 1L;

            return counter;
        }

        VariantCallContext call;
        if ( vcComp.isSNP() ) {
            call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context).get(0);
        } else if ( vcComp.isIndel() ) {
            call = indelEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context).get(0);
        } else if ( bamIsTruth ) {
            // assume it's a SNP if no variation is present; this is necessary so that we can test supposed monomorphic sites against the truth bam
            call = snpEngine.calculateLikelihoodsAndGenotypes(tracker, ref, context).get(0);
        } else {
            logger.info("Not SNP or INDEL " + vcComp.getChr() + ":" + vcComp.getStart() + " " + vcComp.getAlleles());
            return counter;
        }


        boolean writeVariant = true;

        if (bamIsTruth) {
            if (call.confidentlyCalled) {
                // If truth is a confident REF call
                if (call.isVariant()) {
                    if (vcComp.isVariant())
                        counter.nAltCalledAlt = 1L;
                    else {
                        counter.nAltCalledRef = 1L;
                        if ( printInterestingSites )
                            System.out.println("Truth=ALT Call=REF at " + call.getChr() + ":" + call.getStart());
                    }
                }
                // If truth is a confident ALT call
                else {
                    if (vcComp.isVariant()) {
                        counter.nRefCalledAlt = 1L;
                        if ( printInterestingSites )
                            System.out.println("Truth=REF Call=ALT at " + call.getChr() + ":" + call.getStart());
                    } else
                        counter.nRefCalledRef = 1L;
                }
            }
            else {
                counter.nNotConfidentCalls = 1L;
                if ( printInterestingSites )
                    System.out.println("Truth is not confident at " + call.getChr() + ":" + call.getStart());
                writeVariant = false;
            }
        }
        else {
//            if (!vcComp.hasExtendedAttribute("GV"))
//                throw new UserException.BadInput("Variant has no GV annotation in the INFO field. " + vcComp.getChr() + ":" + vcComp.getStart());

            final GVstatus status = getGVstatus(vcComp);
            if (call.isCalledAlt(callConf)) {
                if ( status == GVstatus.T )
                    counter.nAltCalledAlt = 1L;
                else if ( status == GVstatus.F ) {
                    counter.nRefCalledAlt = 1L;
                    if ( printInterestingSites )
                        System.out.println("Truth=REF Call=ALT at " + call.getChr() + ":" + call.getStart());
                }
                else
                    counter.nNoStatusCalledAlt = 1L;
            }
            else if (call.isCalledRef(callConf)) {
                if ( status == GVstatus.T ) {
                    counter.nAltCalledRef = 1L;
                    if ( printInterestingSites )
                        System.out.println("Truth=ALT Call=REF at " + call.getChr() + ":" + call.getStart());
                }
                else if ( status == GVstatus.F )
                    counter.nRefCalledRef = 1L;

                else
                    counter.nNoStatusCalledRef = 1L;
            }
            else {
                counter.nNotConfidentCalls = 1L;
                if ( status == GVstatus.T )
                    counter.nAltNotCalled = 1L;
                else if ( status == GVstatus.F )
                    counter.nRefNotCalled = 1L;
                else
                    counter.nNoStatusNotCalled = 1L;

                if ( printInterestingSites )
                    System.out.println("Truth is not confident at " + call.getChr() + ":" + call.getStart());
                writeVariant = false;
            }
        }

        if (vcfWriter != null && writeVariant) {
            if (!vcComp.hasAttribute(GATKVCFConstants.GENOTYPE_AND_VALIDATE_STATUS_KEY)) {
                vcfWriter.add(new VariantContextBuilder(vcComp).attribute(GATKVCFConstants.GENOTYPE_AND_VALIDATE_STATUS_KEY, call.isCalledAlt(callConf) ? "ALT" : "REF").make());
            }
            else
                vcfWriter.add(vcComp);
        }
        return counter;
    }

    private GVstatus getGVstatus(final VariantContext vc) {
        return ( !vc.hasAttribute("GV") ) ? GVstatus.NONE : (vc.getAttribute("GV").equals("T") ? GVstatus.T : GVstatus.F);
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public CountedData reduceInit() {
        return new CountedData();
    }

    public CountedData treeReduce( final CountedData sum1, final CountedData sum2) {
        sum2.add(sum1);
        return sum2;
    }

    public CountedData reduce( final CountedData mapValue, final CountedData reduceSum ) {
        reduceSum.add(mapValue);
        return reduceSum;
    }

    public void onTraversalDone( CountedData reduceSum ) {
        double ppv = 100 * ((double) reduceSum.nAltCalledAlt /( reduceSum.nAltCalledAlt + reduceSum.nRefCalledAlt));
        double npv = 100 * ((double) reduceSum.nRefCalledRef /( reduceSum.nRefCalledRef + reduceSum.nAltCalledRef));
        double sensitivity = 100 * ((double) reduceSum.nAltCalledAlt /( reduceSum.nAltCalledAlt + reduceSum.nAltCalledRef));
        double specificity = (reduceSum.nRefCalledRef + reduceSum.nRefCalledAlt > 0) ? 100 * ((double) reduceSum.nRefCalledRef /( reduceSum.nRefCalledRef + reduceSum.nRefCalledAlt)) : 100;
        logger.info(String.format("Resulting Truth Table Output\n\n" +
                                  "------------------------------------------------------------------\n" +
                                  "\t\t|\tALT\t|\tREF\t|\tNo Status\n"  +
                                  "------------------------------------------------------------------\n" +
                                  "called alt\t|\t%d\t|\t%d\t|\t%d\n" +
                                  "called ref\t|\t%d\t|\t%d\t|\t%d\n" +
                                  "not called\t|\t%d\t|\t%d\t|\t%d\n" +
                                  "------------------------------------------------------------------\n" +
                                  "positive predictive value: %f%%\n" +
                                  "negative predictive value: %f%%\n" +
                                  "------------------------------------------------------------------\n" +
                                  "sensitivity: %f%%\n" +
                                  "specificity: %f%%\n" +
                                  "------------------------------------------------------------------\n" +
                                  "not confident: %d\n" +
                                  "not covered: %d\n" +
                                  "------------------------------------------------------------------\n", reduceSum.nAltCalledAlt, reduceSum.nRefCalledAlt, reduceSum.nNoStatusCalledAlt, reduceSum.nAltCalledRef, reduceSum.nRefCalledRef, reduceSum.nNoStatusCalledRef, reduceSum.nAltNotCalled, reduceSum.nRefNotCalled, reduceSum.nNoStatusNotCalled, ppv, npv, sensitivity, specificity, reduceSum.nNotConfidentCalls, reduceSum.nUncovered));
    }
}
