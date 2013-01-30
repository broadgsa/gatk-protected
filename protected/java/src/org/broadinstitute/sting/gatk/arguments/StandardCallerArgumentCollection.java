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

package org.broadinstitute.sting.gatk.arguments;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.afcalc.AFCalcFactory;
import org.broadinstitute.sting.utils.collections.DefaultHashMap;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 8/20/12
 * A collection of arguments that are common to the various callers.
 * This is pulled out so that every caller isn't exposed to the arguments from every other caller.
 */

public class StandardCallerArgumentCollection {
    /**
     * The expected heterozygosity value used to compute prior likelihoods for any locus. The default priors are:
     * het = 1e-3, P(hom-ref genotype) = 1 - 3 * het / 2, P(het genotype) = het, P(hom-var genotype) = het / 2
     */
    @Argument(fullName = "heterozygosity", shortName = "hets", doc = "Heterozygosity value used to compute prior likelihoods for any locus", required = false)
    public Double heterozygosity = UnifiedGenotyperEngine.HUMAN_SNP_HETEROZYGOSITY;

    @Argument(fullName = "genotyping_mode", shortName = "gt_mode", doc = "Specifies how to determine the alternate alleles to use for genotyping", required = false)
    public GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY;

    @Argument(fullName = "output_mode", shortName = "out_mode", doc = "Specifies which type of calls we should output", required = false)
    public UnifiedGenotyperEngine.OUTPUT_MODE OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;

    /**
     * The minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls. Only genotypes with
     * confidence >= this threshold are emitted as called sites. A reasonable threshold is 30 for high-pass calling (this
     * is the default).
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_calling", shortName = "stand_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be called", required = false)
    public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;

    /**
     * This argument allows you to emit low quality calls as filtered records.
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_emitting", shortName = "stand_emit_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)", required = false)
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    /**
     * When the UnifiedGenotyper is put into GENOTYPE_GIVEN_ALLELES mode it will genotype the samples using only the alleles provide in this rod binding
     */
    @Input(fullName="alleles", shortName = "alleles", doc="The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES", required=false)
    public RodBinding<VariantContext> alleles;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     *
     * As of GATK 2.2 the genotyper can handle a very large number of events, so the default maximum has been increased to 6.
     */
    @Advanced
    @Argument(fullName = "max_alternate_alleles", shortName = "maxAltAlleles", doc = "Maximum number of alternate alleles to genotype", required = false)
    public int MAX_ALTERNATE_ALLELES = 6;

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    @Argument(fullName = "contamination_fraction_to_filter", shortName = "contamination", doc = "Fraction of contamination in sequencing data (for all samples) to aggressively remove", required = false)
    public double CONTAMINATION_FRACTION = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.05;

    /**
     *  This argument specifies a file with two columns "sample" and "contamination" specifying the contamination level for those samples.
     *  Samples that do not appear in this file will be processed with CONTAMINATION_FRACTION
     **/
    @Advanced
    @Argument(fullName = "contamination_fraction_per_sample_file", shortName = "contaminationFile", doc = "Tab-separated File containing fraction of contamination in sequencing data (per sample) to aggressively remove. Format should be \"<SampleID><TAB><Contamination>\" (Contamination is double) per line; No header.", required = false)
    public File CONTAMINATION_FRACTION_FILE = null;

    /**
     *
     * @return an _Immutable_ copy of the Sample-Contamination Map, defaulting to CONTAMINATION_FRACTION so that if the sample isn't in the map map(sample)==CONTAMINATION_FRACTION
     */
    public Map<String,Double> getSampleContamination(){
        //make sure that the default value is set up right
        sampleContamination.setDefaultValue(CONTAMINATION_FRACTION);
        return Collections.unmodifiableMap(sampleContamination);
    }

    public void setSampleContamination(DefaultHashMap<String, Double> sampleContamination) {
        this.sampleContamination.clear();
        this.sampleContamination.putAll(sampleContamination);
        this.sampleContamination.setDefaultValue(CONTAMINATION_FRACTION);
    }

    //Needs to be here because it uses CONTAMINATION_FRACTION
    private DefaultHashMap<String,Double> sampleContamination = new DefaultHashMap<String,Double>(CONTAMINATION_FRACTION);

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
    @Hidden
    @Argument(fullName = "p_nonref_model", shortName = "pnrm", doc = "Non-reference probability calculation model to employ", required = false)
    public AFCalcFactory.Calculation AFmodel = AFCalcFactory.Calculation.getDefaultModel();

    @Hidden
    @Argument(fullName = "logRemovedReadsFromContaminationFiltering", shortName="contaminationLog", required=false)
    public PrintStream contaminationLog = null;

    @Hidden
    @Argument(shortName = "logExactCalls", doc="x", required=false)
    public File exactCallsLog = null;

    public StandardCallerArgumentCollection() { }

    // Developers must remember to add any newly added arguments to the list here as well otherwise they won't get changed from their default value!
    public StandardCallerArgumentCollection(final StandardCallerArgumentCollection SCAC) {
        this.alleles = SCAC.alleles;
        this.GenotypingMode = SCAC.GenotypingMode;
        this.heterozygosity = SCAC.heterozygosity;
        this.MAX_ALTERNATE_ALLELES = SCAC.MAX_ALTERNATE_ALLELES;
        this.OutputMode = SCAC.OutputMode;
        this.STANDARD_CONFIDENCE_FOR_CALLING = SCAC.STANDARD_CONFIDENCE_FOR_CALLING;
        this.STANDARD_CONFIDENCE_FOR_EMITTING = SCAC.STANDARD_CONFIDENCE_FOR_EMITTING;
        this.CONTAMINATION_FRACTION = SCAC.CONTAMINATION_FRACTION;
        this.CONTAMINATION_FRACTION_FILE=SCAC.CONTAMINATION_FRACTION_FILE;
        this.contaminationLog = SCAC.contaminationLog;
        this.exactCallsLog = SCAC.exactCallsLog;
        this.sampleContamination=SCAC.sampleContamination;
        this.AFmodel = SCAC.AFmodel;
    }
}
