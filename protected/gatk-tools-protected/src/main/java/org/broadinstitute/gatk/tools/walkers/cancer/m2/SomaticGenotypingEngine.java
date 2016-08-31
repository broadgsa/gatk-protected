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

package org.broadinstitute.gatk.tools.walkers.cancer.m2;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

public class SomaticGenotypingEngine extends HaplotypeCallerGenotypingEngine {

    private final M2ArgumentCollection MTAC;
    private final TumorPowerCalculator strandArtifactPowerCalculator;

    private final String tumorSampleName;
    private final String matchedNormalSampleName;
    private final String DEBUG_READ_NAME;

    //Mutect2 does not run in GGA mode
    private static final List<VariantContext> NO_GIVEN_ALLELES = Collections.EMPTY_LIST;

    // {@link GenotypingEngine} requires a non-null {@link AFCalculatorProvider} but this class doesn't need it.  Thus we make a dummy
    private static AFCalculatorProvider DUMMY_AF_CALCULATOR_PROVIDER = new AFCalculatorProvider() {
        public AFCalculator getInstance(final int ploidy, final int maximumAltAlleles) { return null; }
    };

    private final static Logger logger = Logger.getLogger(SomaticGenotypingEngine.class);

    public SomaticGenotypingEngine(final M2ArgumentCollection configuration,
                                   final SampleList samples,
                                   final GenomeLocParser genomeLocParser,
                                   final boolean doPhysicalPhasing,
                                   final M2ArgumentCollection MTAC,
                                   final String tumorSampleName,
                                   final String matchedNormalSampleName,
                                   final String DEBUG_READ_NAME) {
        super(configuration, samples, genomeLocParser, DUMMY_AF_CALCULATOR_PROVIDER, doPhysicalPhasing);
        this.MTAC = MTAC;
        this.tumorSampleName = tumorSampleName;
        this.matchedNormalSampleName = matchedNormalSampleName;
        this.DEBUG_READ_NAME = DEBUG_READ_NAME;

        // coverage related initialization
        //TODO: in GATK4, use a QualityUtils method
        final double errorProbability = Math.pow(10, -MTAC.POWER_CONSTANT_QSCORE/10);
        strandArtifactPowerCalculator = new TumorPowerCalculator(errorProbability, MTAC.STRAND_ARTIFACT_LOD_THRESHOLD, 0.0f);
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param readLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     *
     */
    public CalledHaplotypes callMutations (
            final ReadLikelihoods<Haplotype> readLikelihoods,
            final Map<String, Integer> originalNormalReadQualities,
            final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList,
            final byte[] ref,
            final GenomeLoc refLoc,
            final GenomeLoc activeRegionWindow,
            final RefMetaDataTracker tracker) {
        //TODO: in GATK4 use Utils.nonNull
        if (readLikelihoods == null || readLikelihoods.sampleCount() == 0) throw new IllegalArgumentException("readLikelihoods input should be non-empty and non-null, got "+readLikelihoods);
        if (ref == null || ref.length == 0 ) throw new IllegalArgumentException("ref bytes input should be non-empty and non-null, got "+ref);
        if (refLoc == null || refLoc.size() != ref.length) throw new IllegalArgumentException(" refLoc must be non-null and length must match ref bytes, got "+refLoc);
        if (activeRegionWindow == null ) throw new IllegalArgumentException("activeRegionWindow must be non-null, got "+activeRegionWindow);

        final List<Haplotype> haplotypes = readLikelihoods.alleles();

        // Somatic Tumor/Normal Sample Handling
        if (!readLikelihoods.samples().contains(tumorSampleName)) {
            throw new IllegalArgumentException("readLikelihoods does not contain the tumor sample " + tumorSampleName);
        }
        final boolean hasNormal = matchedNormalSampleName != null;

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final TreeSet<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, readLikelihoods, ref, refLoc, NO_GIVEN_ALLELES);

        // Walk along each position in the key set and create each event to be outputted
        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        for( final int loc : startPosKeySet ) {
            if( loc < activeRegionWindow.getStart() || loc > activeRegionWindow.getStop() ) {
                continue;
            }

            final List<VariantContext> eventsAtThisLoc = getVCsAtThisLocation(haplotypes, loc, NO_GIVEN_ALLELES);

            if( eventsAtThisLoc.isEmpty() ) { continue; }

            // Create the event mapping object which maps the original haplotype events to the events present at just this locus
            final Map<Event, List<Haplotype>> eventMapper = createEventMapper(loc, eventsAtThisLoc, haplotypes);

            // TODO: priorityList is not sorted by priority, might as well just use eventsAtThisLoc.map(VariantContext::getSource)
            final List<String> priorityList = makePriorityList(eventsAtThisLoc);

            // merge variant contexts from multiple haplotypes into one variant context
            // TODO: we should use haplotypes if possible, but that may have to wait for GATK4
            VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(eventsAtThisLoc, priorityList,
                    GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                    GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

            if( mergedVC == null ) { continue; }

            // TODO: this varaible needs a descriptive name
            final Map<VariantContext, Allele> mergeMap = new LinkedHashMap<>();

            mergeMap.put(null, mergedVC.getReference()); // the reference event (null) --> the reference allele
            for(int i = 0; i < eventsAtThisLoc.size(); i++) {
                // TODO: as noted below, this operation seems dangerous. Understand how things can go wrong.
                mergeMap.put(eventsAtThisLoc.get(i), mergedVC.getAlternateAllele(i)); // BUGBUG: This is assuming that the order of alleles is the same as the priority list given to simpleMerge function
            }

            /** TODO: the code in the for loop up to here needs refactor. The goal, as far as I can tell, is to create two things: alleleMapper and mergedVC
             * alleleMapper maps alleles to haplotypes, and we need this to create readAlleleLikelihoods.
             * To make alleleMapper we make mergeMap (of type VC -> Allele) and eventMapper (of type Event -> List(Haplotypes), where Event is essentialy Variant Context)
             * If we just want a map of Alleles to Haplotypes, we should be able to do so directly; no need for intermediate maps, which just complicates the code.
             **/

            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergeMap, eventMapper);

            // converting ReadLikelihoods<Haplotype> to ReadLikeliHoods<Allele>
            ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper, genomeLocParser.createPaddedGenomeLoc(genomeLocParser.createGenomeLoc(mergedVC), ALLELE_EXTENSION));

            //LDG: do we want to do this before or after pulling out overlapping reads?
            if (MTAC.isSampleContaminationPresent()) {
                readAlleleLikelihoods.contaminationDownsampling(MTAC.getSampleContamination());
            }

            // TODO: this is a good break point for a new method
            // TODO: replace PRALM with ReadLikelihoods
            final PerReadAlleleLikelihoodMap tumorPRALM = readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(readAlleleLikelihoods.sampleIndex(tumorSampleName));
            filterPRALMForOverlappingReads(tumorPRALM, mergedVC.getReference(), loc, false);
            MuTect2.logReadInfo(DEBUG_READ_NAME, tumorPRALM.getLikelihoodReadMap().keySet(), "Present in Tumor PRALM after filtering for overlapping reads");
            // extend to multiple samples

            // compute tumor LOD for each alternate allele
            // TODO: somewhere we have to ensure that the all the alleles in the variant context is in alleleFractions passed to getHetGenotypeLogLikelihoods. getHetGenotypeLogLikelihoods will not check that for you
            final PerAlleleCollection<Double> altAlleleFractions = estimateAlleleFraction(mergedVC, tumorPRALM, false);
            final PerAlleleCollection<Double> tumorHetGenotypeLLs = getHetGenotypeLogLikelihoods(mergedVC, tumorPRALM, originalNormalReadQualities, altAlleleFractions);

            final PerAlleleCollection<Double> tumorLods = PerAlleleCollection.createPerAltAlleleCollection();
            for (final Allele altAllele : mergedVC.getAlternateAlleles()){
                tumorLods.set(altAllele, tumorHetGenotypeLLs.get(altAllele) - tumorHetGenotypeLLs.getRef());
            }

            // TODO: another good breakpoint e.g. compute normal LOD/set thresholds
            // TODO: anything related to normal should be encapsulated in Optional

            // A variant candidate whose normal LOD is below this threshold will be filtered as 'germline_risk'
            // This is a more stringent threshold than normalLodThresholdForVCF
            double normalLodFilterThreshold = -Double.MAX_VALUE;
            PerReadAlleleLikelihoodMap normalPRALM = null;
            final PerAlleleCollection<Double> normalLods = PerAlleleCollection.createPerAltAlleleCollection();

            // if normal bam is available, compute normal LOD
            // TODO: this if statement should be a standalone method for computing normal LOD
            // TODO: then we can do something like normalLodThreshold = hasNormal ? thisMethod() : Optional.empty()
            if (hasNormal) {
                normalPRALM = readAlleleLikelihoods.toPerReadAlleleLikelihoodMap(readAlleleLikelihoods.sampleIndex(matchedNormalSampleName));
                filterPRALMForOverlappingReads(normalPRALM, mergedVC.getReference(), loc, true);
                MuTect2.logReadInfo(DEBUG_READ_NAME, normalPRALM.getLikelihoodReadMap().keySet(), "Present after in Nomral PRALM filtering for overlapping reads");

                final GenomeLoc eventGenomeLoc = genomeLocParser.createGenomeLoc(activeRegionWindow.getContig(), loc);
                final Collection<VariantContext> cosmicVC = tracker.getValues(MTAC.cosmicRod, eventGenomeLoc);
                final Collection<VariantContext> dbsnpVC = tracker.getValues(MTAC.dbsnp.dbsnp, eventGenomeLoc);
                final boolean germlineAtRisk = !dbsnpVC.isEmpty() && cosmicVC.isEmpty();

                normalLodFilterThreshold = germlineAtRisk ? MTAC.NORMAL_DBSNP_LOD_THRESHOLD : MTAC.NORMAL_LOD_THRESHOLD;

                // compute normal LOD = LL(X|REF)/LL(X|ALT) where REF is the diploid HET with AF = 0.5
                // note normal LOD is REF over ALT, the reciprocal of the tumor LOD
                final PerAlleleCollection<Double> diploidHetAlleleFractions = PerAlleleCollection.createPerRefAndAltAlleleCollection();
                for (final Allele allele : mergedVC.getAlternateAlleles()){
                    diploidHetAlleleFractions.setAlt(allele, 0.5);
                }

                final PerAlleleCollection<Double> normalGenotypeLLs = getHetGenotypeLogLikelihoods(mergedVC, normalPRALM, originalNormalReadQualities, diploidHetAlleleFractions);

                for (final Allele altAllele : mergedVC.getAlternateAlleles()){
                    normalLods.setAlt(altAllele, normalGenotypeLLs.getRef() - normalGenotypeLLs.getAlt(altAllele));
                }
            }

            int numPassingAlts = 0;
            final Set<Allele> allelesThatPassThreshold = new HashSet<>();
            Allele alleleWithHighestTumorLOD = null;

            for (final Allele altAllele : mergedVC.getAlternateAlleles()) {
                final boolean passesTumorLodThreshold = tumorLods.getAlt(altAllele) >= MTAC.INITIAL_TUMOR_LOD_THRESHOLD;
                final boolean passesNormalLodThreshold = hasNormal ? normalLods.getAlt(altAllele) >= MTAC.INITIAL_NORMAL_LOD_THRESHOLD : true;
                if (passesTumorLodThreshold && passesNormalLodThreshold) {
                    numPassingAlts++;
                    allelesThatPassThreshold.add(altAllele);
                    if (alleleWithHighestTumorLOD == null || tumorLods.getAlt(altAllele) > tumorLods.getAlt(alleleWithHighestTumorLOD)){
                        alleleWithHighestTumorLOD = altAllele;
                    }
                }
            }

            if (numPassingAlts == 0) {
                continue;
            }

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC);
            final int haplotypeCount = alleleMapper.get(alleleWithHighestTumorLOD).size();
            callVcb.attribute(GATKVCFConstants.HAPLOTYPE_COUNT_KEY, haplotypeCount);
            callVcb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, tumorLods.getAlt(alleleWithHighestTumorLOD));

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY, normalLods.getAlt(alleleWithHighestTumorLOD));
                if (normalLods.getAlt(alleleWithHighestTumorLOD) < normalLodFilterThreshold) {
                    callVcb.filter(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
                }
            }

            // TODO: this should be a separate method
            // TODO: move code to MuTect2::calculateFilters()
            if (MTAC.ENABLE_STRAND_ARTIFACT_FILTER && numPassingAlts == 1) {
                final PerReadAlleleLikelihoodMap forwardPRALM = new PerReadAlleleLikelihoodMap();
                final PerReadAlleleLikelihoodMap reversePRALM = new PerReadAlleleLikelihoodMap();
                splitPRALMintoForwardAndReverseReads(tumorPRALM, forwardPRALM, reversePRALM);

                MuTect2.logReadInfo(DEBUG_READ_NAME, tumorPRALM.getLikelihoodReadMap().keySet(), "Present in tumor PRALM after PRALM is split");
                MuTect2.logReadInfo(DEBUG_READ_NAME, forwardPRALM.getLikelihoodReadMap().keySet(), "Present in forward PRALM after PRALM is split");
                MuTect2.logReadInfo(DEBUG_READ_NAME, reversePRALM.getLikelihoodReadMap().keySet(), "Present in reverse PRALM after PRALM is split");

                // TODO: build a new type for probability, likelihood, and log_likelihood. e.g. f_fwd :: probability[], tumorGLs_fwd :: likelihood[]
                // TODO: don't want to call getHetGenotypeLogLikelihoods on more than one alternate alelle. May need to overload it to take a scalar f_fwd.
                final PerAlleleCollection<Double> alleleFractionsForward = estimateAlleleFraction(mergedVC, forwardPRALM, true);
                final PerAlleleCollection<Double> tumorGenotypeLLForward = getHetGenotypeLogLikelihoods(mergedVC, forwardPRALM, originalNormalReadQualities, alleleFractionsForward);

                final PerAlleleCollection<Double> alleleFractionsReverse = estimateAlleleFraction(mergedVC, reversePRALM, true);
                final PerAlleleCollection<Double> tumorGenotypeLLReverse = getHetGenotypeLogLikelihoods(mergedVC, reversePRALM, originalNormalReadQualities, alleleFractionsReverse);

                final double tumorLod_fwd = tumorGenotypeLLForward.getAlt(alleleWithHighestTumorLOD) - tumorGenotypeLLForward.getRef();
                final double tumorLod_rev = tumorGenotypeLLReverse.getAlt(alleleWithHighestTumorLOD) - tumorGenotypeLLReverse.getRef();

                // Note that we use the observed combined (+ and -) allele fraction for power calculation in either direction
                final double tumorSBpower_fwd = strandArtifactPowerCalculator.cachedPowerCalculation(forwardPRALM.getNumberOfStoredElements(), altAlleleFractions.getAlt(alleleWithHighestTumorLOD));
                final double tumorSBpower_rev = strandArtifactPowerCalculator.cachedPowerCalculation(reversePRALM.getNumberOfStoredElements(), altAlleleFractions.getAlt(alleleWithHighestTumorLOD));

                callVcb.attribute(GATKVCFConstants.TLOD_FWD_KEY, tumorLod_fwd);
                callVcb.attribute(GATKVCFConstants.TLOD_REV_KEY, tumorLod_rev);
                callVcb.attribute(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY, tumorSBpower_fwd);
                callVcb.attribute(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY, tumorSBpower_rev);

                if ((tumorSBpower_fwd > MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && tumorLod_fwd < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                        (tumorSBpower_rev > MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && tumorLod_rev < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD))
                    callVcb.filter(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);
            }

            // TODO: this probably belongs in M2::calculateFilters()
            if (numPassingAlts > 1) {
                callVcb.filter(GATKVCFConstants.TRIALLELIC_SITE_FILTER_NAME);
            }

            // build genotypes TODO: this part needs review and refactor
            final List<Allele> tumorAlleles = Arrays.asList(mergedVC.getReference(), alleleWithHighestTumorLOD);
            // TODO: estimateAlleleFraction should not repeat counting allele depths
            final PerAlleleCollection<Integer> tumorAlleleDepths = getRefAltCount(mergedVC, tumorPRALM, false);
            final int tumorRefAlleleDepth = tumorAlleleDepths.getRef();
            final int tumorAltAlleleDepth = tumorAlleleDepths.getAlt(alleleWithHighestTumorLOD);
            final Genotype tumorGenotype = new GenotypeBuilder(tumorSampleName, tumorAlleles)
                    .AD(new int[] { tumorRefAlleleDepth, tumorAltAlleleDepth })
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, altAlleleFractions.getAlt(alleleWithHighestTumorLOD))
                    .make();

            final List<Genotype> genotypes = new ArrayList<>();
            genotypes.add(tumorGenotype);

            // We assume that the genotype in the normal is 0/0
            // TODO: is normal always homozygous reference?
            final List<Allele> homRefAllelesforNormalGenotype = Collections.nCopies(2, mergedVC.getReference());

            // if we are calling with a normal, build the genotype for the sample to appear in vcf
            if (hasNormal) {
                final PerAlleleCollection<Integer> normalAlleleDepths = getRefAltCount(mergedVC, normalPRALM, false);
                final int normalRefAlleleDepth = normalAlleleDepths.getRef();
                final int normalAltAlleleDepth = normalAlleleDepths.getAlt(alleleWithHighestTumorLOD);
                final double normalAlleleFraction = (double) normalAltAlleleDepth / ( normalRefAlleleDepth + normalAltAlleleDepth);

                final Genotype normalGenotype = new GenotypeBuilder(matchedNormalSampleName, homRefAllelesforNormalGenotype)
                        .AD(new int[] { normalRefAlleleDepth, normalAltAlleleDepth })
                        .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, normalAlleleFraction)
                        .make();
                genotypes.add(normalGenotype);
            }

            final VariantContext call = new VariantContextBuilder(callVcb).alleles(tumorAlleles).genotypes(genotypes).make();

            // how should we be making use of _perSampleFilteredReadList_?
            readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
                    genomeLocParser, false, alleleMapper, readAlleleLikelihoods, call);

            final ReferenceContext referenceContext = new ReferenceContext(genomeLocParser, genomeLocParser.createGenomeLoc(mergedVC.getChr(), mergedVC.getStart(), mergedVC.getEnd()), refLoc, ref);
            VariantContext annotatedCall = annotationEngine.annotateContextForActiveRegion(referenceContext, tracker, readAlleleLikelihoods, call, false);

            if( call.getAlleles().size() != mergedVC.getAlleles().size() )
                annotatedCall = GATKVariantContextUtils.reverseTrimAlleles(annotatedCall);

            // maintain the set of all called haplotypes
            call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            returnCalls.add( annotatedCall );
        }

        // TODO: understand effect of enabling this for somatic calling...
        final List<VariantContext> outputCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        return new CalledHaplotypes(outputCalls, calledHaplotypes);
    }

    /** Calculate the likelihoods of hom ref and each het genotype of the form ref/alt
     *
     * @param mergedVC                              input VC
     * @param tumorPRALM                            read likelihoods
     * @param originalNormalMQs                     original MQs, before boosting normals to avoid qual capping
     * @param alleleFractions                       allele fraction(s) for alternate allele(s)
     *
     * @return                                      genotype likelihoods for homRef and het for each alternate allele
     */
    private PerAlleleCollection<Double> getHetGenotypeLogLikelihoods(final VariantContext mergedVC,
                                                                     final PerReadAlleleLikelihoodMap tumorPRALM,
                                                                     final Map<String, Integer> originalNormalMQs,
                                                                     final PerAlleleCollection<Double> alleleFractions) {
        // make sure that alleles in alleleFraction are a subset of alleles in the variant context
        if (! mergedVC.getAlternateAlleles().containsAll(alleleFractions.getAltAlleles()) ){
            throw new IllegalArgumentException("alleleFractions has alleles that are not in the variant context");
        }

        final PerAlleleCollection<MutableDouble> genotypeLogLikelihoods = PerAlleleCollection.createPerRefAndAltAlleleCollection();
        mergedVC.getAlleles().forEach(a -> genotypeLogLikelihoods.set(a, new MutableDouble(0)));

        final Allele refAllele = mergedVC.getReference();
        for(Map.Entry<GATKSAMRecord,Map<Allele, Double>> readAlleleLikelihoodMap : tumorPRALM.getLikelihoodReadMap().entrySet()) {
            final Map<Allele, Double> alleleLikelihoodMap = readAlleleLikelihoodMap.getValue();
            if (originalNormalMQs.get(readAlleleLikelihoodMap.getKey().getReadName()) == 0) {
                continue;
            }

            final double readRefLogLikelihood = alleleLikelihoodMap.get(refAllele);
            genotypeLogLikelihoods.getRef().add(readRefLogLikelihood);

            for (final Allele altAllele : alleleFractions.getAltAlleles()) {
                final double readAltLogLikelihood = alleleLikelihoodMap.get(altAllele);
                final double adjustedReadAltLL = Math.log10(
                        Math.pow(10, readRefLogLikelihood) * (1 - alleleFractions.getAlt(altAllele)) +
                                Math.pow(10, readAltLogLikelihood) * alleleFractions.getAlt(altAllele)
                );
                genotypeLogLikelihoods.get(altAllele).add(adjustedReadAltLL);
            }

        }

        final PerAlleleCollection<Double> result = PerAlleleCollection.createPerRefAndAltAlleleCollection();
        mergedVC.getAlleles().stream().forEach(a -> result.set(a,genotypeLogLikelihoods.get(a).toDouble()));

        return result;
    }

    /**
     * Find the allele fractions for each alternate allele
     *
     * @param vc                        input VC, for alleles
     * @param pralm                     read likelihoods
     * @return                          estimated AF for each alt
     */
    // FIXME: calculate using the uncertainty rather than this cheap approach
    private PerAlleleCollection<Double> estimateAlleleFraction(final VariantContext vc,
                                                               final PerReadAlleleLikelihoodMap pralm,
                                                               final boolean oneStrandOnly) {
        final PerAlleleCollection<Integer> alleleCounts = getRefAltCount(vc, pralm, oneStrandOnly);
        final PerAlleleCollection<Double> alleleFractions = PerAlleleCollection.createPerAltAlleleCollection();

        final int refCount = alleleCounts.getRef();
        for ( final Allele altAllele : vc.getAlternateAlleles() ) {
            final int altCount = alleleCounts.getAlt(altAllele);
            double alleleFraction = (double) altCount / (refCount + altCount);
            // weird case, but I've seen it happen in one strand cases
            if (refCount == 0 && altCount == refCount ) {
                alleleFraction = 0;
            }
            alleleFractions.setAlt(altAllele, alleleFraction);
            // logger.info("Counted " + refCount + " ref and " + altCount + " alt " );
        }

        return alleleFractions;
    }

    /**
     *  Go through the PRALM and tally the most likely allele in each read. Only count informative reads.
     *
     * @param vc                      input VC, for alleles
     * @param pralm                         read likelihoods
     * @return                              an array giving the read counts for the ref and each alt allele
     */
    private PerAlleleCollection<Integer> getRefAltCount(final VariantContext vc,
                                                        final PerReadAlleleLikelihoodMap pralm,
                                                        final boolean oneStrandOnly) {
        // Check that the alleles in Variant Context are in PRALM
        // Skip the check for strand-conscious PRALM; + reads may not have alleles in - reads, for example.
        final Set<Allele> vcAlleles = new HashSet<>(vc.getAlleles());
        if ( ! oneStrandOnly && ! pralm.getAllelesSet().containsAll( vcAlleles ) ) {
            StringBuilder message = new StringBuilder();
            message.append("At Locus chr" + vc.getContig() + ":" + vc.getStart() + ", we detected that variant context had alleles that not in PRALM. ");
            message.append("VC alleles = " + vcAlleles + ", PRALM alleles = " + pralm.getAllelesSet());
            logger.warn(message);
        }


        final PerAlleleCollection<MutableInt> alleleCounts = PerAlleleCollection.createPerRefAndAltAlleleCollection();
        vcAlleles.stream().forEach(a -> alleleCounts.set(a, new MutableInt(0)));

        for (final Map.Entry<GATKSAMRecord, Map<Allele, Double>> readAlleleLikelihoodMap : pralm.getLikelihoodReadMap().entrySet()) {
            final GATKSAMRecord read = readAlleleLikelihoodMap.getKey();
            final Map<Allele, Double> alleleLikelihoodMap = readAlleleLikelihoodMap.getValue();
            final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(alleleLikelihoodMap, vcAlleles);

            if (read.getMappingQuality() > 0 && mostLikelyAllele.isInformative()) {
                alleleCounts.get(mostLikelyAllele.getMostLikelyAllele()).increment();
            }

        }

        final PerAlleleCollection<Integer> result = PerAlleleCollection.createPerRefAndAltAlleleCollection();
        vc.getAlleles().stream().forEach(a -> result.set(a, alleleCounts.get(a).toInteger()));

        return(result);
    }

    private void logM2Debug(String s) {
        if (MTAC.M2_DEBUG) {
            logger.info(s);
        }
    }

    private void filterPRALMForOverlappingReads(final PerReadAlleleLikelihoodMap pralm, final Allele ref, final int location, final boolean retainMismatches) {
        final Map<GATKSAMRecord, Map<Allele, Double>> m = pralm.getLikelihoodReadMap();

        // iterate through the reads, if the name has been seen before we have overlapping (potentially) fragments, so handle them
        final Map<String, GATKSAMRecord> nameToRead = new HashMap<>();
        final Set<GATKSAMRecord> readsToKeep = new HashSet<>();

        for(final GATKSAMRecord rec : m.keySet()) {
            // if we haven't seen it... just record the name and add it to the list of reads to keep
            final GATKSAMRecord existing = nameToRead.get(rec.getReadName());
            if (existing == null) {
                nameToRead.put(rec.getReadName(), rec);
                readsToKeep.add(rec);
            } else {
                logM2Debug("Found a paired read for " + rec.getReadName());

                // NOTE: Can we use FragmentUtils to do all of this processing (to find overlapping pairs?)
                // seems like maybe, but it has some requirements about the order of the reads supplied which may be painful to meet
                // TODO: CHECK IF THE READS BOTH OVERLAP THE POSITION!!!!
                if ( ReadUtils.isInsideRead(existing, location) && ReadUtils.isInsideRead(rec, location) ) {

                    final MostLikelyAllele existingMLA = PerReadAlleleLikelihoodMap.getMostLikelyAllele(pralm.getLikelihoodReadMap().get(existing));
                    final Allele existingAllele = existingMLA.getMostLikelyAllele();

                    final MostLikelyAllele recMLA = PerReadAlleleLikelihoodMap.getMostLikelyAllele(pralm.getLikelihoodReadMap().get(rec));
                    final Allele recAllele = recMLA.getMostLikelyAllele();

                    // if the reads disagree at this position...
                    if (!existingAllele.equals(recAllele)) {
                        //... and we're not retaining mismatches, throw them both out
                        if (!retainMismatches) {
                            logM2Debug("Discarding read-pair due to disagreement" + rec.getReadName() + " and allele " + existingAllele);
                            readsToKeep.remove(existing);

                            //... and we are retaining mismatches, keep the mismatching one
                        } else {
                            if (existingAllele.equals(ref)) {
                                logM2Debug("Discarding read to keep mismatching " + rec.getReadName() + " and allele " + existingAllele);
                                readsToKeep.remove(existing);
                                readsToKeep.add(rec);
                            }
                        }
                        // Otherwise, keep the element with the higher quality score
                    } else {
                        logM2Debug("Discarding lower quality read of overlapping pair " + rec.getReadName() + " and allele " + existingAllele);
                        if (existingMLA.getLog10LikelihoodOfMostLikely() < recMLA.getLog10LikelihoodOfMostLikely()) {
                            readsToKeep.remove(existing);
                            readsToKeep.add(rec);
                        }
                    }
                } else {
                    // although these are overlapping fragments, they don't overlap at the position in question
                    // so keep the read
                    readsToKeep.add(rec);
                }
            }

        }

        // perhaps moved into PRALM
        final Iterator<Map.Entry<GATKSAMRecord, Map<Allele, Double>>> it = m.entrySet().iterator();
        while ( it.hasNext() ) {
            final Map.Entry<GATKSAMRecord, Map<Allele, Double>> record = it.next();
            if(!readsToKeep.contains(record.getKey())) {
                it.remove();
                logM2Debug("Dropping read " + record.getKey() + " due to overlapping read fragment rules");
            }
        }
    }

    private void splitPRALMintoForwardAndReverseReads(final PerReadAlleleLikelihoodMap originalPRALM, final PerReadAlleleLikelihoodMap forwardPRALM, final PerReadAlleleLikelihoodMap reversePRALM) {
        final Map<GATKSAMRecord, Map<Allele, Double>> origReadAlleleLikelihoodMap = originalPRALM.getLikelihoodReadMap();
        for (final GATKSAMRecord read : origReadAlleleLikelihoodMap.keySet()) {
            if (read.isStrandless())
                continue;

            for (final Map.Entry<Allele, Double> alleleLikelihoodMap : origReadAlleleLikelihoodMap.get(read).entrySet()) {
                final Allele allele = alleleLikelihoodMap.getKey();
                final Double likelihood = alleleLikelihoodMap.getValue();
                if (read.getReadNegativeStrandFlag())
                    reversePRALM.add(read, allele, likelihood);
                else
                    forwardPRALM.add(read, allele, likelihood);
            }
        }
    }
}
