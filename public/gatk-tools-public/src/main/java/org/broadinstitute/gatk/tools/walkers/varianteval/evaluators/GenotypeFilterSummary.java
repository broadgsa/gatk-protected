package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.Iterator;

/**
 * Created by bimber on 5/17/2017.
 */
@Analysis(description = "Counts called and filtered genotypes across samples")
public class GenotypeFilterSummary extends VariantEvaluator {
    @DataPoint(description = "Number of Called Genotypes", format = "%d")
    public long nCalledNotFiltered = 0;

    @DataPoint(description = "Number of No-Call Genotypes", format = "%d")
    public long nNoCallOrFiltered = 0;

    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Iterator<Genotype> it = eval.getGenotypes().iterator();
        while (it.hasNext()){
            Genotype g = it.next();
            if (g.isCalled() && !g.isFiltered()){
                nCalledNotFiltered++;
            }
            else if (g.isNoCall() || g.isFiltered()){
                nNoCallOrFiltered++;
            }
        }
    }
}
