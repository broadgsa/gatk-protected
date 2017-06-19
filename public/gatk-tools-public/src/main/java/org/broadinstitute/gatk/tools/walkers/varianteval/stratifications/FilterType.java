/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.*;

/**
 * Stratifies by the FILTER type(s) for each line, with PASS used for passing
 */
public class FilterType extends VariantStratifier {
    @Override
    public void initialize() {
        Set<String> filterNames = new HashSet<>();
        Map<String, VCFHeader> headers = GATKVCFUtils.getVCFHeadersFromRods(getVariantEvalWalker().getToolkit(), getVariantEvalWalker ().getEvals());
        for (VCFHeader header : headers.values()){
            for (VCFFilterHeaderLine line : header.getFilterLines()) {
                filterNames.add(line.getID());
            }
        }

        states.addAll(filterNames);
        states.add("PASS");
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval == null){
            return Collections.emptyList();
        }

        ArrayList<Object> relevantStates = new ArrayList<Object>();
        if (eval.isFiltered()){
            relevantStates.addAll(eval.getFilters());
        }
        else {
            relevantStates.add("PASS");
        }

        return relevantStates;
    }
}
