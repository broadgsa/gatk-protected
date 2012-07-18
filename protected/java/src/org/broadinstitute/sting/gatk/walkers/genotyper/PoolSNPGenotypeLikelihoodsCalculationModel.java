package org.broadinstitute.sting.gatk.walkers.genotyper;
/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class PoolSNPGenotypeLikelihoodsCalculationModel extends PoolGenotypeLikelihoodsCalculationModel {


    protected PoolSNPGenotypeLikelihoodsCalculationModel( UnifiedArgumentCollection UAC,  Logger logger) {
        super(UAC, logger);

    }

    protected PoolGenotypeLikelihoods getPoolGenotypeLikelihoodObject(final List<Allele> alleles,
                                                                               final double[] logLikelihoods,
                                                                               final int ploidy,
                                                                               final HashMap<String, ErrorModel> perLaneErrorModels,
                                                                               final boolean useBQAedPileup,
                                                                               final ReferenceContext ref,
                                                                               final boolean ignoreLaneInformation) {
        return new PoolSNPGenotypeLikelihoods(alleles, null, UAC.samplePloidy, perLaneErrorModels, useBQAedPileup, UAC.IGNORE_LANE_INFO);
    }

    protected List<Allele> getInitialAllelesToUse(final RefMetaDataTracker tracker,
                                                  final ReferenceContext ref,
                                                  Map<String, AlignmentContext> contexts,
                                                  final AlignmentContextUtils.ReadOrientation contextType,
                                                  final GenomeLocParser locParser,
                                                  final List<Allele> allAllelesToUse) {

        if (allAllelesToUse != null)
            return allAllelesToUse;


        final byte refBase = ref.getBase();
        final List<Allele> allAlleles = new ArrayList<Allele>();
        // first add ref allele
        allAlleles.add(Allele.create(refBase, true));
        // add all possible alt alleles
        for (byte b: BaseUtils.BASES) {
            if (refBase != b)
                allAlleles.add(Allele.create(b));
        }

        return allAlleles;
    } 
    
    protected List<Allele> getFinalAllelesToUse(final RefMetaDataTracker tracker,
                                                final ReferenceContext ref,
                                                final List<Allele> allAllelesToUse,
                                                final ArrayList<PoolGenotypeData> GLs) {
        // find the alternate allele(s) that we should be using
        final List<Allele> alleles = new ArrayList<Allele>();
        if ( allAllelesToUse != null ) {
            alleles.addAll(allAllelesToUse);
        } else if ( UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext vc = UnifiedGenotyperEngine.getVCFromAllelesRod(tracker, ref, ref.getLocus(), true, logger, UAC.alleles);

            // ignore places where we don't have a SNP
            if ( vc == null || !vc.isSNP() )
                return null;

            alleles.addAll(vc.getAlleles());
        } else {

            alleles.add(Allele.create(ref.getBase(),true));
            alleles.addAll(determineAlternateAlleles( GLs));

            // if there are no non-ref alleles...
            if ( alleles.size() == 1 ) {
                final int indexOfRefBase = BaseUtils.simpleBaseToBaseIndex(ref.getBase());
                // if we only want variants, then we don't need to calculate genotype likelihoods
                if ( UAC.OutputMode != UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                    // otherwise, choose any alternate allele (it doesn't really matter)
                    alleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(indexOfRefBase == 0 ? 1 : 0)));
            }
        }
        return alleles;
    }

    /**
     * @param tracker           dummy parameter here
     * @param ref               Reference context
     * @param alternateAllelesToUse alt allele list
     * @return end location for vc to be created
      */
    protected int getEndLocation(final RefMetaDataTracker tracker,
                                 final ReferenceContext ref,
                                 final List<Allele> alternateAllelesToUse) {
        // for SNPs, end loc is is the same as start loc
        return ref.getLocus().getStart();

    }


}
