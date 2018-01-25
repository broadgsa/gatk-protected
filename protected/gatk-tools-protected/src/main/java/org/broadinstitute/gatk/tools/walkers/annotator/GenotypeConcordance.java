package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.Arrays;
import java.util.List;

/**
 * Created by bimber on 4/20/2017.
 */
public class GenotypeConcordance extends GenotypeAnnotation {
    private final static Logger logger = Logger.getLogger(GenotypeConcordance.class);
    public final static String BINDING_NAME = "GT_SOURCE";
    private boolean foundRodBinding = false;
    private RodBinding<VariantContext> rodBinding = null;
    private int warningsLogged = 0;

    public static final String KEY = "GTD";
    public static final String D_KEY = "REF_GT";

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(KEY, D_KEY);
    }

    @Override
    public void annotate(RefMetaDataTracker tracker, AnnotatorCompatible walker, ReferenceContext ref, AlignmentContext stratifiedContext, VariantContext vc, Genotype g, GenotypeBuilder gb, PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        if (rodBinding == null) {
            for (RodBinding<VariantContext> binding : walker.getResourceRodBindings()) {
                if (BINDING_NAME.equalsIgnoreCase(binding.getName())) {
                    foundRodBinding = true;
                    rodBinding = binding;

                    break;
                }
            }
        }

        if (vc.getContig().equalsIgnoreCase("chrUn")){
            return;
        }


        if (!foundRodBinding){
            throw new IllegalArgumentException("Must provide a VCF resource with the name: " + BINDING_NAME);
        }

        if (!g.isFiltered() && !g.isNoCall()) {
            List<VariantContext> list = tracker.getValues(rodBinding);
            if (list == null){
                if (warningsLogged < 10) {
                    logger.warn("position not found in reference VCF: " + vc.getContig() + ":" + vc.getStart());
                    warningsLogged++;
                }

                if (warningsLogged == 10){
                    logger.warn("future warnings will not be logged");
                }
                return;
            }

            if (!list.isEmpty()){
                for (VariantContext c : list){
                    Genotype refGenotype = c.getGenotype(g.getSampleName());
                    if (refGenotype != null && !refGenotype.isFiltered() && !refGenotype.isNoCall()) {
                        if (!refGenotype.sameGenotype(g)) {
                            gb.attribute(KEY, "1");
                            gb.attribute(D_KEY, refGenotype.getGenotypeString());
                        }
                        else {
                            gb.attribute(KEY, "0");
                        }

                    }
                }
            }
        }
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFFormatHeaderLine(KEY, 1, VCFHeaderLineType.Integer, "Flags genotypes (as 1) discordant with those from the same sample/position in the provided VCF file.  Concordant genotypes are flagged a 0.  Genotypes not called in either VCF are ignored."),
                new VCFFormatHeaderLine(D_KEY, 1, VCFHeaderLineType.String, "When comparing genotypes against an alternate VCF, this will store the genotype of this sample in that alternate VCF, if discordant.  Genotypes not called in either VCF are ignored.")
        );
    }
}
