package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.HashMap;
import java.util.Map;

/**
 * NOTE: unfortunately VariantEvaluators dont easily support flexible column numbers, so this will not work in the current form
 * Created by bimber on 5/17/2017.
 */
public class SiteFilterSummary extends VariantEvaluator {

    public Map<String, Long> filterCountMap = new HashMap<>();

    //@DataPoint(description = "Number of SNPs", format = "%d")
    //public long nSites = 0;

    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (eval.isFiltered()){
            eval.getFilters().forEach(x -> append(x));
        }
        else {
            append("PASS");
        }
    }

    private void append(String filterName){
        if (!filterCountMap.containsKey(filterName)){
            filterCountMap.put(filterName, 0L);
        }

        filterCountMap.put(filterName, filterCountMap.get(filterName) + 1);
    }
}
