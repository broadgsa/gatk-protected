package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This is designed to add the VRT annotation required by dbSNP, which is effectively the same thing as VariantType
 *
 * Related to this thread:
 * http://gatkforums.broadinstitute.org/gatk/discussion/6405/format-vcf-for-dbsnp-submission
 *
 * Created by bimber on 5/4/2017.
 */
public class DbSnpVariationType extends InfoFieldAnnotation {
    private static final String KEY = "VRT";

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(KEY);
    }

    @Override
    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatible walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        Map<String, Object> ret = new HashMap<>();
        switch (vc.getType())
        {
            case SNP:
                ret.put(KEY, 1);
                break;
            case INDEL:
                ret.put(KEY, 2);
                break;
            case MIXED:
                ret.put(KEY, 7);
                break;
            case MNP:
                ret.put(KEY, 8);
                break;
        }

        return ret;
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(
            new VCFInfoHeaderLine(KEY, 1, VCFHeaderLineType.Integer, "Variation type,1 - SNV: single nucleotide variation,2 - DIV: deletion/insertion variation,3 - HETEROZYGOUS: variable, but undefined at nucleotide level,4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element,6 - NO VARIATON: sequence scanned for variation, but none observed,7 - MIXED: cluster contains submissions from 2 or more allelic classes (not used),8 - MNV: multiple nucleotide variation with alleles of common length greater than 1,9 - Exception"));
    }
}

