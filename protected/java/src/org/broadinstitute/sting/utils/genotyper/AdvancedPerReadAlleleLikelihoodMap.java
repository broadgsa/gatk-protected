/*
 * Copyright (c) 2011 The Broad Institute
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
package org.broadinstitute.sting.utils.genotyper;


import org.broadinstitute.sting.gatk.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.sting.utils.classloader.ProtectedPackageSource;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.io.PrintStream;
import java.util.*;

public class AdvancedPerReadAlleleLikelihoodMap extends StandardPerReadAlleleLikelihoodMap implements ProtectedPackageSource {

    public ReadBackedPileup createPerAlleleDownsampledBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction, final PrintStream log) {
        return AlleleBiasedDownsamplingUtils.createAlleleBiasedBasePileup(pileup, downsamplingFraction, log);
    }

    public void performPerAlleleDownsampling(final double downsamplingFraction, final PrintStream log) {
        // special case removal of all or no reads
        if ( downsamplingFraction <= 0.0 )
            return;
        if ( downsamplingFraction >= 1.0 ) {
            likelihoodReadMap.clear();
            return;
        }

        // start by stratifying the reads by the alleles they represent at this position
        final Map<Allele, List<GATKSAMRecord>> alleleReadMap = new HashMap<Allele, List<GATKSAMRecord>>(alleles.size());
        for ( Allele allele : alleles )
            alleleReadMap.put(allele, new ArrayList<GATKSAMRecord>());

        for ( Map.Entry<GATKSAMRecord, Map<Allele, Double>> entry : likelihoodReadMap.entrySet() ) {
            // do not remove reduced reads!
            if ( !entry.getKey().isReducedRead() ) {
                final Allele bestAllele = getMostLikelyAllele(entry.getValue());
                if ( bestAllele != Allele.NO_CALL )
                    alleleReadMap.get(bestAllele).add(entry.getKey());
            }
        }

        // compute the reads to remove and actually remove them
        final List<GATKSAMRecord> readsToRemove = AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(alleleReadMap, downsamplingFraction, log);
        for ( final GATKSAMRecord read : readsToRemove )
            likelihoodReadMap.remove(read);
    }
}
