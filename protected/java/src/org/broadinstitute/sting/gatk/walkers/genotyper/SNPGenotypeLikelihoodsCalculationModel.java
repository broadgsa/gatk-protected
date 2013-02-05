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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

public class SNPGenotypeLikelihoodsCalculationModel extends GenotypeLikelihoodsCalculationModel {

    private final boolean useAlleleFromVCF;

    private final double[] likelihoodSums = new double[4];

    private final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap;

    protected SNPGenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        super(UAC, logger);
        useAlleleFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
        perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
    }

    public VariantContext getLikelihoods(final RefMetaDataTracker tracker,
                                         final ReferenceContext ref,
                                         final Map<String, AlignmentContext> contexts,
                                         final AlignmentContextUtils.ReadOrientation contextType,
                                         final List<Allele> allAllelesToUse,
                                         final boolean useBAQedPileup,
                                         final GenomeLocParser locParser,
                                         final Map<String, PerReadAlleleLikelihoodMap> sampleLikelihoodMap) {

        sampleLikelihoodMap.clear(); // not used in SNP model, sanity check to delete any older data

        final byte refBase = ref.getBase();
        final int indexOfRefBase = BaseUtils.simpleBaseToBaseIndex(refBase);
        // handle non-standard reference bases
        if ( indexOfRefBase == -1 )
            return null;
        final Allele refAllele = Allele.create(refBase, true);

        // calculate the GLs
        ArrayList<SampleGenotypeData> GLs = new ArrayList<SampleGenotypeData>(contexts.size());
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
        // Down-sample with bias according to the contamination level (global or per file)
            ReadBackedPileup pileup = AlignmentContextUtils.stratify(sample.getValue(), contextType).getBasePileup();
            final Double contamination =  UAC.getSampleContamination().get(sample.getKey());
            if( contamination > 0.0 ) //no need to enter if no contamination reduction
                pileup = perReadAlleleLikelihoodMap.createPerAlleleDownsampledBasePileup(pileup,contamination, UAC.contaminationLog);
            if ( useBAQedPileup )
                pileup = createBAQedPileup(pileup);

            // create the GenotypeLikelihoods object
            final DiploidSNPGenotypeLikelihoods GL = new DiploidSNPGenotypeLikelihoods(UAC.PCR_error);
            final int nGoodBases = GL.add(pileup, true, true, UAC.MIN_BASE_QUALTY_SCORE);
            if ( nGoodBases > 0 )
                GLs.add(new SampleGenotypeData(sample.getKey(), GL, getFilteredDepth(pileup)));
        }

        // start making the VariantContext
        final GenomeLoc loc = ref.getLocus();
        final List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(refAllele);


        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", loc.getContig(), loc.getStart(), loc.getStop(), alleles);
        // find the alternate allele(s) that we should be using
        if ( allAllelesToUse != null ) {
            alleles.addAll(allAllelesToUse.subList(1,allAllelesToUse.size()));   // this includes ref allele
        } else if ( useAlleleFromVCF ) {
            final VariantContext vc = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, ref.getLocus(), true, logger, UAC.alleles);

            // ignore places where we don't have a SNP
            if ( vc == null || !vc.isSNP() )
                return null;

            // make sure a user isn't passing the REF base in as an ALT
            if ( vc.hasAlternateAllele(refAllele, true) )
                throw new UserException.BadInput("Alternate allele '" + (char)refBase + "' passed in is the same as the reference at location " + vc.getChr() + ":" + vc.getStart());

            alleles.addAll(vc.getAlternateAlleles());
        } else {

            alleles.addAll(determineAlternateAlleles(refBase, GLs));

            // if there are no non-ref alleles...
            if ( alleles.size() == 1 ) {
                // if we only want variants, then we don't need to calculate genotype likelihoods
                if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    return builder.make();

                // otherwise, choose any alternate allele (it doesn't really matter)
                alleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(indexOfRefBase == 0 ? 1 : 0)));
             }
        }

        // create the alternate alleles and the allele ordering (the ordering is crucial for the GLs)
        final int numAlleles = alleles.size();
        final int numAltAlleles = numAlleles - 1;

        final int[] alleleOrdering = new int[numAlleles];
        int alleleOrderingIndex = 0;
        int numLikelihoods = 0;
        for ( Allele allele : alleles ) {
            alleleOrdering[alleleOrderingIndex++] = BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0]);
            numLikelihoods += alleleOrderingIndex;
        }
        builder.alleles(alleles);

        // create the PL ordering to use based on the allele ordering.
        final int[] PLordering = new int[numLikelihoods];
        for ( int i = 0; i <= numAltAlleles; i++ ) {
            for ( int j = i; j <= numAltAlleles; j++ ) {
                // As per the VCF spec: "the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
                // In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
                PLordering[(j * (j+1) / 2) + i] = DiploidGenotype.createDiploidGenotype(alleleOrdering[i], alleleOrdering[j]).ordinal();
            }
        }

        // create the genotypes; no-call everyone for now
        final GenotypesContext genotypes = GenotypesContext.create();

        for ( SampleGenotypeData sampleData : GLs ) {
            final double[] allLikelihoods = sampleData.GL.getLikelihoods();
            final double[] myLikelihoods = new double[numLikelihoods];

            for ( int i = 0; i < numLikelihoods; i++ )
                myLikelihoods[i] = allLikelihoods[PLordering[i]];

            // normalize in log space so that max element is zero.
            final GenotypeBuilder gb = new GenotypeBuilder(sampleData.name);
            final double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(myLikelihoods, false, true);
            gb.PL(genotypeLikelihoods);
            gb.DP(sampleData.depth);
            genotypes.add(gb.make());
        }

        return builder.genotypes(genotypes).make();
    }
    
    // determines the alleles to use
    protected List<Allele> determineAlternateAlleles(final byte ref, final List<SampleGenotypeData> sampleDataList) {

        final int baseIndexOfRef = BaseUtils.simpleBaseToBaseIndex(ref);
        final int PLindexOfRef = DiploidGenotype.createDiploidGenotype(ref, ref).ordinal();
        for ( int i = 0; i < 4; i++ )
            likelihoodSums[i] = 0.0;
        
        // based on the GLs, find the alternate alleles with enough probability
        for ( SampleGenotypeData sampleData : sampleDataList ) {
            final double[] likelihoods = sampleData.GL.getLikelihoods();
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PLindexOfRef ) {
                GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindexOfBestGL);
                if ( alleles.alleleIndex1 != baseIndexOfRef )
                    likelihoodSums[alleles.alleleIndex1] += likelihoods[PLindexOfBestGL] - likelihoods[PLindexOfRef];
                // don't double-count it
                if ( alleles.alleleIndex2 != baseIndexOfRef && alleles.alleleIndex2 != alleles.alleleIndex1 )
                    likelihoodSums[alleles.alleleIndex2] += likelihoods[PLindexOfBestGL] - likelihoods[PLindexOfRef];
            }
        }

        final List<Allele> allelesToUse = new ArrayList<Allele>(3);
        for ( int i = 0; i < 4; i++ ) {
            if ( likelihoodSums[i] > 0.0 )
                allelesToUse.add(Allele.create(BaseUtils.baseIndexToSimpleBase(i), false));
        }

        return allelesToUse;
    }

    public ReadBackedPileup createBAQedPileup( final ReadBackedPileup pileup ) {
        final List<PileupElement> BAQedElements = new ArrayList<PileupElement>();
        for( final PileupElement PE : pileup ) {
            final PileupElement newPE = new BAQedPileupElement( PE );
            BAQedElements.add( newPE );
        }
        return new ReadBackedPileupImpl( pileup.getLocation(), BAQedElements );
    }

    public static class BAQedPileupElement extends PileupElement {
        public BAQedPileupElement( final PileupElement PE ) {
            super(PE);
        }

        @Override
        public byte getQual() {
            if ( isDeletion() )
                return super.getQual();
            else
                return BAQ.calcBAQFromTag(getRead(), offset, true);
        }
    }

    private static class SampleGenotypeData {

        public final String name;
        public final DiploidSNPGenotypeLikelihoods GL;
        public final int depth;

        public SampleGenotypeData(final String name, final DiploidSNPGenotypeLikelihoods GL, final int depth) {
            this.name = name;
            this.GL = GL;
            this.depth = depth;
        }
    }
}
