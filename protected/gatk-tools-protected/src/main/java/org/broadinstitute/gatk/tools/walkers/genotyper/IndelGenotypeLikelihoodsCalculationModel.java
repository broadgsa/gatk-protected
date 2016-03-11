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
import org.broadinstitute.gatk.tools.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

public class IndelGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {
    private static final int HAPLOTYPE_SIZE = 80;

    private boolean DEBUG = false;
    private boolean ignoreSNPAllelesWhenGenotypingIndels = false;
    private PairHMMIndelErrorModel pairModel;


    private LinkedHashMap<Allele, Haplotype> haplotypeMap;

    private List<Allele> alleleList = new ArrayList<Allele>();


    protected IndelGenotypeLikelihoodsCalculationModel(final UnifiedArgumentCollection UAC,
                                                       final Logger logger) {
        super(UAC, logger);
        pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.pairHMM);
        DEBUG = UAC.OUTPUT_DEBUG_INDEL_INFO;
        haplotypeMap = new LinkedHashMap<Allele, Haplotype>();
        ignoreSNPAllelesWhenGenotypingIndels = UAC.IGNORE_SNP_ALLELES;
    }

    protected static List<Allele> computeConsensusAlleles(final ReferenceContext ref,
                                                 final Map<String, AlignmentContext> contexts,
                                                 final AlignmentContextUtils.ReadOrientation contextType,
                                                 final UnifiedArgumentCollection UAC) {
        ConsensusAlleleCounter counter = new ConsensusAlleleCounter(true, UAC.MIN_INDEL_COUNT_FOR_GENOTYPING, UAC.MIN_INDEL_FRACTION_PER_SAMPLE);
        return counter.computeConsensusAlleles(ref, contexts, contextType);
    }

    private final static EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.INDEL, VariantContext.Type.MIXED);


    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> allAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser,
                                         final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        GenomeLoc loc = ref.getLocus();
//        if (!ref.getLocus().equals(lastSiteVisited)) {
        if (contextType == AlignmentContextUtils.ReadOrientation.COMPLETE) {
            // starting a new site: clear allele list
            haplotypeMap.clear();
            perReadAlleleLikelihoodMap.clear(); // clean mapping sample-> per read, per allele likelihoods
            alleleList = getInitialAlleleList(tracker, ref, contexts, contextType, UAC, ignoreSNPAllelesWhenGenotypingIndels);
            if (alleleList.isEmpty())
                return null;
        }

        getHaplotypeMapFromAlleles(alleleList, ref, loc, haplotypeMap); // will update haplotypeMap adding elements
        if (haplotypeMap == null || haplotypeMap.isEmpty())
            return null;

        // start making the VariantContext
        // For all non-snp VC types, VC end location is just startLocation + length of ref allele including padding base.
        final int endLoc = loc.getStart() + alleleList.get(0).length() - 1;
        final int eventLength = getEventLength(alleleList);

        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), endLoc, alleleList);

        // create the genotypes; no-call everyone for now
        GenotypesContext genotypes = GenotypesContext.create();
        final int ploidy = UAC.genotypeArgs.samplePloidy;
        final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy);

        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            if (!perReadAlleleLikelihoodMap.containsKey(sample.getKey())){
                // no likelihoods have been computed for this sample at this site
                perReadAlleleLikelihoodMap.put(sample.getKey(), new PerReadAlleleLikelihoodMap());
            }
            final ReadBackedPileup pileup = context.getBasePileup();
            if (pileup != null) {
                final GenotypeBuilder b = new GenotypeBuilder(sample.getKey());
                final double[] genotypeLikelihoods = pairModel.computeDiploidReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, perReadAlleleLikelihoodMap.get(sample.getKey()), UAC.getSampleContamination().get(sample.getKey()));
                b.PL(genotypeLikelihoods);
                b.alleles(noCall);
                b.DP(getFilteredDepth(pileup));
                genotypes.add(b.make());

                if (DEBUG) {
                    System.out.format("Sample:%s Alleles:%s GL:", sample.getKey(), alleleList.toString());
                    for (int k = 0; k < genotypeLikelihoods.length; k++)
                        System.out.format("%1.4f ", genotypeLikelihoods[k]);
                    System.out.println();
                }
            }
        }

        return builder.genotypes(genotypes).make();
    }

    public static void getHaplotypeMapFromAlleles(final List<Allele> alleleList,
                                                 final ReferenceContext ref,
                                                 final GenomeLoc loc,
                                                 final LinkedHashMap<Allele, Haplotype> haplotypeMap) {
        // protect against having an indel too close to the edge of a contig
        if (loc.getStart() <= HAPLOTYPE_SIZE)
            haplotypeMap.clear();
        // check if there is enough reference window to create haplotypes (can be an issue at end of contigs)
        else if (ref.getWindow().getStop() < loc.getStop() + HAPLOTYPE_SIZE)
            haplotypeMap.clear();
        else if (alleleList.isEmpty())
            haplotypeMap.clear();
        else {
            final int eventLength = getEventLength(alleleList);
            final int hsize = ref.getWindow().size() - Math.abs(eventLength) - 1;
            final int numPrefBases = ref.getLocus().getStart() - ref.getWindow().getStart() + 1;

            if (hsize <= 0)  // protect against event lengths larger than ref window sizes
                haplotypeMap.clear();
            else
                haplotypeMap.putAll(Haplotype.makeHaplotypeListFromAlleles(alleleList, loc.getStart(),
                    ref, hsize, numPrefBases));
        }
    }

    public static int getEventLength(List<Allele> alleleList) {
        Allele refAllele = alleleList.get(0);
        Allele altAllele = alleleList.get(1);
        // look for alt allele that has biggest length distance to ref allele
        int maxLenDiff = 0;
        for (Allele a : alleleList) {
            if (a.isNonReference()) {
                int lenDiff = Math.abs(a.getBaseString().length() - refAllele.getBaseString().length());
                if (lenDiff > maxLenDiff) {
                    maxLenDiff = lenDiff;
                    altAllele = a;
                }
            }
        }

        return altAllele.getBaseString().length() - refAllele.getBaseString().length();

    }
    
    public static List<Allele> getInitialAlleleList(final RefMetaDataTracker tracker,
                                                    final ReferenceContext ref,
                                                    final Map<String, AlignmentContext> contexts,
                                                    final AlignmentContextUtils.ReadOrientation contextType,
                                                    final UnifiedArgumentCollection UAC,
                                                    final boolean ignoreSNPAllelesWhenGenotypingIndels) {

        List<Allele> alleles = new ArrayList<Allele>();
        if (UAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) {
            VariantContext vc = null;
            for (final VariantContext vc_input : tracker.getValues(UAC.alleles, ref.getLocus())) {
                if (vc_input != null &&
                        allowableTypes.contains(vc_input.getType()) &&
                        ref.getLocus().getStart() == vc_input.getStart()) {
                    vc = vc_input;
                    break;
                }
            }
           // ignore places where we don't have a variant
            if (vc == null)
                return alleles;

            if (ignoreSNPAllelesWhenGenotypingIndels) {
                // if there's an allele that has same length as the reference (i.e. a SNP or MNP), ignore it and don't genotype it
                for (Allele a : vc.getAlleles())
                    if (a.isNonReference() && a.getBases().length == vc.getReference().getBases().length)
                        continue;
                    else
                        alleles.add(a);

            } else {
                alleles.addAll(vc.getAlleles());
            }

        } else {
            alleles = computeConsensusAlleles(ref, contexts, contextType, UAC);
        }
        return alleles;
    }

    // Overload function in GenotypeLikelihoodsCalculationModel so that, for an indel case, we consider a deletion as part of the pileup,
    // so that per-sample DP will include deletions covering the event.
    protected int getFilteredDepth(ReadBackedPileup pileup) {
        int count = 0;
        for (PileupElement p : pileup) {
            if (p.isDeletion() || BaseUtils.isRegularBase(p.getBase()))
                count++;
        }

        return count;
    }

}
