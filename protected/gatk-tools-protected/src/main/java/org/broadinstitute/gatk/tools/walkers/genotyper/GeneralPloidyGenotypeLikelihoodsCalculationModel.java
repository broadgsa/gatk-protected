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

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

public abstract class GeneralPloidyGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {

    final protected UnifiedArgumentCollection UAC;

    protected GeneralPloidyGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC,logger);
        this.UAC = UAC;

    }


    /*
       Get vc with alleles from reference sample. Can be null if there's no ref sample call or no ref sample coverage at this site.
    */
    protected VariantContext getTrueAlleles(final RefMetaDataTracker tracker,
                                            final ReferenceContext ref,
                                            Map<String,AlignmentContext> contexts) {
        // Get reference base from VCF or Reference
        if (UAC.referenceSampleName == null)
            return null;

        AlignmentContext context = contexts.get(UAC.referenceSampleName);
        ArrayList<Allele> trueReferenceAlleles = new ArrayList<Allele>();

        VariantContext referenceSampleVC;

        if (tracker != null && context != null)
            referenceSampleVC = tracker.getFirstValue(UAC.referenceSampleRod, context.getLocation());
        else
            return null;

        if (referenceSampleVC == null) {
            trueReferenceAlleles.add(Allele.create(ref.getBase(),true));
            return new VariantContextBuilder("pc",ref.getLocus().getContig(), ref.getLocus().getStart(), ref.getLocus().getStop(),trueReferenceAlleles).make();

        }
        else {
            Genotype referenceGenotype = referenceSampleVC.getGenotype(UAC.referenceSampleName);
            List<Allele> referenceAlleles = referenceGenotype.getAlleles();

            return new VariantContextBuilder("pc",referenceSampleVC.getChr(), referenceSampleVC.getStart(), referenceSampleVC.getEnd(),
                    referenceSampleVC.getAlleles())
                    .genotypes(new GenotypeBuilder(UAC.referenceSampleName, referenceAlleles).GQ(referenceGenotype.getGQ()).make())
                    .make();
        }
    }


    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     * This function returns the list of lane identifiers.
     *
     * @param readGroups readGroups A collection of read group strings (obtained from the alignment context pileup)
     * @return a collection of lane ids.
     */
    public static Set<String> parseLaneIDs(Collection<String> readGroups) {
        HashSet<String> result = new HashSet<String>();
        for (String readGroup : readGroups) {
            result.add(getLaneIDFromReadGroupString(readGroup));
        }
        return result;
    }

    /**
     * GATK Engine creates readgroups of the form XXX.Y.Z
     * XXX.Y is the unique lane identifier.
     *     Z is the id of the sample to make the read group id unique
     *
     * @param readGroupID the read group id string
     * @return just the lane id (the XXX.Y string)
     */
    public static String getLaneIDFromReadGroupString(String readGroupID) {
        // System.out.println(readGroupID);
        String [] parsedID = readGroupID.split("\\.");
        if (parsedID.length > 1)
            return parsedID[0] + "." + parsedID[1];
        else
            return parsedID[0] + ".0";
    }


    /** Wrapper class that encapsulates likelihood object and sample name
     *
     */
    protected static class PoolGenotypeData {

        public final String name;
        public final GeneralPloidyGenotypeLikelihoods GL;
        public final int depth;
        public final List<Allele> alleles;

        public PoolGenotypeData(final String name, final GeneralPloidyGenotypeLikelihoods GL, final int depth, final List<Allele> alleles) {
            this.name = name;
            this.GL = GL;
            this.depth = depth;
            this.alleles = alleles;
        }
    }

    // determines the alleles to use
    protected List<Allele> determineAlternateAlleles(final List<PoolGenotypeData> sampleDataList) {

        if (sampleDataList.isEmpty())
            return Collections.emptyList();

        final int REFERENCE_IDX = 0;
        final List<Allele> allAlleles = sampleDataList.get(0).GL.getAlleles();
        double[] likelihoodSums = new double[allAlleles.size()];

        // based on the GLs, find the alternate alleles with enough probability
        for ( PoolGenotypeData sampleData : sampleDataList ) {
            final Pair<int[],Double> mlACPair = sampleData.GL.getMostLikelyACCount();
            final double topLogGL = mlACPair.second;

            if (sampleData.GL.getAlleles().size() != allAlleles.size())
                throw new ReviewedGATKException("BUG: inconsistent size of alleles!");

            // ref allele is always first in array list
            if (sampleData.GL.alleles.get(0).isNonReference())
                throw new ReviewedGATKException("BUG: first allele in list is not reference!");

            double refGL = sampleData.GL.getLikelihoods()[REFERENCE_IDX];

            // check if maximum likelihood AC is all-ref for current pool. If so, skip
            if (mlACPair.first[REFERENCE_IDX] == sampleData.GL.numChromosomes)
                continue;

            // most likely AC is not all-ref: for all non-ref alleles, add difference of max likelihood and all-ref likelihood
            for (int i=0; i < mlACPair.first.length; i++) {
                if (i==REFERENCE_IDX) continue;

                if (mlACPair.first[i] > 0)
                    likelihoodSums[i] += topLogGL - refGL;

            }
        }

        final List<Allele> allelesToUse = new ArrayList<Allele>();
        for ( int i = 0; i < likelihoodSums.length; i++ ) {
            if ( likelihoodSums[i] > 0.0 )
                allelesToUse.add(allAlleles.get(i));
        }

        return allelesToUse;
    }


    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> allAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser,
                                         final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        HashMap<String, ErrorModel> perLaneErrorModels = getPerLaneErrorModels(tracker, ref, contexts);
        if (perLaneErrorModels == null && UAC.referenceSampleName != null)
            return null;

        if (UAC.TREAT_ALL_READS_AS_SINGLE_POOL) {
            AlignmentContext mergedContext = AlignmentContextUtils.joinContexts(contexts.values());
            Map<String,AlignmentContext> newContext = new HashMap<String,AlignmentContext>();
            newContext.put(DUMMY_SAMPLE_NAME,mergedContext);
            contexts = newContext;
        }
        if (contextType == AlignmentContextUtils.ReadOrientation.COMPLETE) {
            // starting a new site: clear allele list
            perReadAlleleLikelihoodMap.clear(); // clean mapping sample-> per read, per allele likelihoods
        }
            // get initial alleles to genotype
        final List<Allele> allAlleles = new ArrayList<Allele>();
        if (allAllelesToUse == null || allAllelesToUse.isEmpty())
            allAlleles.addAll(getInitialAllelesToUse(tracker, ref,contexts,contextType,locParser, allAllelesToUse));
        else
            allAlleles.addAll(allAllelesToUse);

        if (allAlleles.isEmpty())
            return null;

        final ArrayList<PoolGenotypeData> GLs = new ArrayList<PoolGenotypeData>(contexts.size());

        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            // skip reference sample
            if (UAC.referenceSampleName != null && sample.getKey().equals(UAC.referenceSampleName))
                continue;

            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup();
            if (!perReadAlleleLikelihoodMap.containsKey(sample.getKey())){
                // no likelihoods have been computed for this sample at this site
                perReadAlleleLikelihoodMap.put(sample.getKey(), new PerReadAlleleLikelihoodMap());
            }

            // create the GenotypeLikelihoods object
            final GeneralPloidyGenotypeLikelihoods GL = getPoolGenotypeLikelihoodObject(allAlleles, null, UAC.genotypeArgs.samplePloidy, perLaneErrorModels, useBAQedPileup, ref, UAC.IGNORE_LANE_INFO, perReadAlleleLikelihoodMap.get(sample.getKey()));
            // actually compute likelihoods
            final int nGoodBases = GL.add(pileup, UAC);
            if ( nGoodBases > 0 )
                // create wrapper object for likelihoods and add to list
                GLs.add(new PoolGenotypeData(sample.getKey(), GL, getFilteredDepth(pileup), allAlleles));
        }

        // find the alternate allele(s) that we should be using
        final List<Allele> alleles = getFinalAllelesToUse(tracker, ref, allAllelesToUse, GLs);
        if (alleles == null || alleles.isEmpty() || (alleles.size() == 1 && alleles.get(0).isReference()))
            return null;
        // start making the VariantContext
        final GenomeLoc loc = ref.getLocus();
        final int endLoc = getEndLocation(tracker, ref, alleles);

        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), endLoc, alleles);
        builder.alleles(alleles);

        final HashMap<String, Object> attributes = new HashMap<String, Object>();

        if (UAC.referenceSampleName != null && perLaneErrorModels != null)
            attributes.put(GATKVCFConstants.REFSAMPLE_DEPTH_KEY, ErrorModel.getTotalReferenceDepth(perLaneErrorModels));

        builder.attributes(attributes);
        // create the genotypes; no-call everyone for now
        final GenotypesContext genotypes = GenotypesContext.create();
        final int ploidy = UAC.genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy);

        for ( PoolGenotypeData sampleData : GLs ) {
            // extract from multidimensional array
            final double[] myLikelihoods = GeneralPloidyGenotypeLikelihoods.subsetToAlleles(sampleData.GL.getLikelihoods(), sampleData.GL.numChromosomes,
                    allAlleles, alleles);

            // normalize in log space so that max element is zero.
            final GenotypeBuilder gb = new GenotypeBuilder(sampleData.name, noCall);
            gb.DP(sampleData.depth);
            gb.PL(MathUtils.normalizeFromLog10(myLikelihoods, false, true));
            genotypes.add(gb.make());
        }

        return builder.genotypes(genotypes).make();

    }


    protected HashMap<String, ErrorModel> getPerLaneErrorModels(final RefMetaDataTracker tracker,
                                                                final ReferenceContext ref,
                                                                Map<String, AlignmentContext> contexts) {
        VariantContext refVC =  getTrueAlleles(tracker, ref, contexts);


        // Build error model for site based on reference sample, and keep stratified for each lane.
        AlignmentContext refContext = null;
        if (UAC.referenceSampleName != null)
            refContext = contexts.get(UAC.referenceSampleName);

        ReadBackedPileup refPileup = null;
        if (refContext != null) {
            HashMap<String, ErrorModel> perLaneErrorModels = new HashMap<String, ErrorModel>();
            refPileup = refContext.getBasePileup();

            Set<String> laneIDs = new TreeSet<String>();
            if (UAC.TREAT_ALL_READS_AS_SINGLE_POOL || UAC.IGNORE_LANE_INFO)
                laneIDs.add(DUMMY_LANE);
            else
                laneIDs = parseLaneIDs(refPileup.getReadGroups());
            // build per-lane error model for all lanes present in ref sample
            for (String laneID : laneIDs) {
                // get reference pileup for this lane
                ReadBackedPileup refLanePileup = refPileup;
                // subset for this lane
                if (refPileup != null && !(UAC.TREAT_ALL_READS_AS_SINGLE_POOL || UAC.IGNORE_LANE_INFO))
                    refLanePileup = refPileup.getPileupForLane(laneID);

                //ReferenceSample referenceSample = new ReferenceSample(UAC.referenceSampleName, refLanePileup, trueReferenceAlleles);
                perLaneErrorModels.put(laneID, new ErrorModel(UAC,  refLanePileup, refVC, ref));
            }
            return perLaneErrorModels;

        }
        else
            return null;

    }

    /*
       Abstract methods - must be implemented in derived classes
    */

    protected abstract GeneralPloidyGenotypeLikelihoods getPoolGenotypeLikelihoodObject(final List<Allele> alleles,
                                                                               final double[] logLikelihoods,
                                                                               final int ploidy,
                                                                               final HashMap<String, ErrorModel> perLaneErrorModels,
                                                                               final boolean useBQAedPileup,
                                                                               final ReferenceContext ref,
                                                                               final boolean ignoreLaneInformation,
                                                                               final org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap);

    protected abstract List<Allele> getInitialAllelesToUse(final RefMetaDataTracker tracker,
                                                           final ReferenceContext ref,
                                                           Map<String, AlignmentContext> contexts,
                                                           final AlignmentContextUtils.ReadOrientation contextType,
                                                           final GenomeLocParser locParser,
                                                           final List<Allele> allAllelesToUse);

    protected abstract List<Allele> getFinalAllelesToUse(final RefMetaDataTracker tracker,
                                                         final ReferenceContext ref,
                                                         final List<Allele> allAllelesToUse,
                                                         final ArrayList<PoolGenotypeData> GLs);

    protected abstract int getEndLocation(final RefMetaDataTracker tracker,
                                          final ReferenceContext ref,
                                          final List<Allele> alternateAllelesToUse);
}
