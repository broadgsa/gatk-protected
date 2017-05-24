package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.util.*;

/**
 * Created by bimber on 3/13/2017.
 */
public class MendelianViolationBySample extends GenotypeAnnotation {
    private final static Logger logger = Logger.getLogger(MendelianViolationBySample.class);
    private double minGenotypeQuality = 10.0;
    private SampleDB sampleDB = null;
    private Boolean isVariantAnnotator = null;
    private int nonAutosomePositions = 0;

    public static final String MV_KEY = "MV";

    @Override
    public void initialize(AnnotatorCompatible walker, GenomeAnalysisEngine toolkit, Set<VCFHeaderLine> headerLines) {
        super.initialize(walker, toolkit, headerLines);

        if (((VariantAnnotator)walker).minGenotypeQualityP > 0.0)
        {
            minGenotypeQuality = ((VariantAnnotator)walker).minGenotypeQualityP;
        }
        logger.info("MendelianViolationBySample minimum genotype quality: " + minGenotypeQuality);

        if ( sampleDB == null ) {
            sampleDB = ((Walker) walker).getSampleDB();
        }

    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(MV_KEY);
    }

    @Override
    public void annotate(RefMetaDataTracker tracker, AnnotatorCompatible walker, ReferenceContext ref, AlignmentContext stratifiedContext, VariantContext vc, Genotype g, GenotypeBuilder gb, PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        // Can only be called from VariantAnnotator
        isVariantAnnotator = MendelianViolationCount.logWalkerIdentityCheck(walker, isVariantAnnotator);
        if (!isVariantAnnotator){
            return;
        }

        for (String key : getKeyNames())
        {
            gb.attribute(key, null);
        }

        //not autosome
        if (vc.getContig().equalsIgnoreCase("chrUn")){
            nonAutosomePositions++;
            return;
        }

        if (sampleDB != null) {
            int totalViolations = MendelianViolationCount.countViolations(sampleDB.getSample(g.getSampleName()), vc, minGenotypeQuality);
            gb.attribute(MV_KEY, totalViolations);
        }
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(new VCFFormatHeaderLine(MV_KEY, 1, VCFHeaderLineType.Integer, "Number of mendelian violations observed for this sample."));
    }
}
