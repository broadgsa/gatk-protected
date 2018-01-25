package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.RodRequiringAnnotation;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.*;

/**
 * Created by bimber on 4/20/2017.
 */
public class GenotypeConcordanceBySite extends InfoFieldAnnotation implements RodRequiringAnnotation {
    private final static Logger logger = Logger.getLogger(GenotypeConcordanceBySite.class);

    private boolean foundRodBinding = false;
    private RodBinding<VariantContext> rodBinding = null;
    private int warningsLogged = 0;
    private boolean resourceWarningsLogged = false;

    public static final String DISCORD_KEY = "GTD";
    public static final String CONCORD_KEY = "GTC";

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(DISCORD_KEY, CONCORD_KEY);
    }

    @Override
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if (rodBinding == null) {
            for (RodBinding<VariantContext> binding : walker.getResourceRodBindings()) {
                if (GenotypeConcordance.BINDING_NAME.equalsIgnoreCase(binding.getName())) {
                    foundRodBinding = true;
                    rodBinding = binding;

                    break;
                }
            }
        }

        if (!foundRodBinding){
            throw new IllegalArgumentException("Must provide a VCF resource with the name: " + GenotypeConcordance.BINDING_NAME);
        }

        List<VariantContext> list = tracker.getValues(rodBinding);
        if (list == null) {
            if (warningsLogged < 10) {
                logger.warn("position not found in reference VCF: " + vc.getContig() + ":" + vc.getStart());
                warningsLogged++;
            }

            if (warningsLogged == 10) {
                logger.warn("future warnings will not be logged");
            }

            return null;
        }
        else if (list.isEmpty()) {
            return null;
        }
        else if (list.size() > 1 && !resourceWarningsLogged) {
            logger.warn("more than one resource found with the name: " + GenotypeConcordance.BINDING_NAME);
            resourceWarningsLogged = true;
        }

        VariantContext refVC  = list.get(0);

        int discord = 0;
        int concord = 0;
        Iterator<Genotype> it = vc.getGenotypes().iterator();
        while (it.hasNext()) {
            Genotype g = it.next();
            if (!g.isFiltered() && !g.isNoCall()) {
                Genotype refGenotype = refVC.getGenotype(g.getSampleName());
                if (refGenotype != null && !refGenotype.isFiltered() && !refGenotype.isNoCall()) {
                    if (!refGenotype.sameGenotype(g)) {
                        discord++;
                    }
                    else {
                        concord++;
                    }
                }
            }
        }

        Map<String,Object> attributeMap = new HashMap<>(2);
        attributeMap.put(DISCORD_KEY, discord);
        attributeMap.put(CONCORD_KEY, concord);

        return attributeMap;
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFInfoHeaderLine(DISCORD_KEY, 1, VCFHeaderLineType.Integer, "The total number of genotypes discordant with those from the same sample/position in the provided VCF file.  Genotypes not called in either VCF are ignored."),
                new VCFInfoHeaderLine(CONCORD_KEY, 1, VCFHeaderLineType.Integer, "The total number of genotypes concordant with those from the same sample/position in the provided VCF file.  Genotypes not called in either VCF are ignored.")
        );
    }
}
