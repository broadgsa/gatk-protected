package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Created by bimber on 5/4/2017.
 */
@By(DataSource.REFERENCE)
public class RemoveAnnotations extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName="annotationToKeep", shortName="A", doc="List the specific INFO field annotations to retain (by their keys, such as AC, AF, etc).  If specified, all other INFO field annotations will be removed.", required=false)
    protected List<String> annotationToKeep = new ArrayList<>();

    @Argument(fullName="annotationToRemove", shortName="XA", doc="One or more specific INFO field annotations to remove (by their keys, such as AC, AF, etc).", required=false)
    protected List<String> annotationsToExclude = new ArrayList<>();

    @Argument(fullName="genotypeAnnotationToKeep", shortName="GA", doc="List the specific genotype (format field) annotations to retain", required=false)
    protected List<String> genotypeAnnotationToKeep = new ArrayList<>();

    @Argument(fullName="genotypeAnnotationToRemove", shortName="XGA", doc="One or more specific genotype (format field) annotations to remove.", required=false)
    protected List<String> genotypeAnnotationsToExclude = new ArrayList<>();

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered sites", required=false)
    protected boolean excludeFiltered = true;

    @Argument(fullName="retainExtraHeaderLines", shortName="rh", doc="If provided, additional header lines (metadata, etc) will be retained.", required=false)
    protected boolean retainExtraHeaderLines = false;

    @Argument(fullName="clearGenotypeFilter", shortName="cgf", doc="Clear the filter field on all genotypes.  ", required=false)
    protected boolean clearGTfilter = true;

    @Argument(fullName="sitesOnly", shortName="sitesOnly", doc="Omit samples and genotypes from the output VCF.  ", required=false)
    protected boolean sitesOnly = false;

    private VCFHeader header;
    private Set<String> allowableInfoKeys;
    private Set<String> allowableFormatKeys;

    @Override
    public void initialize() {
        super.initialize();

        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Collections.singletonList(variantCollection.variants.getName()));
        final Set<String> samples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        final VCFHeader initialHeader = new VCFHeader(VCFUtils.smartMergeHeaders(vcfRods.values(), true), samples);

        //retain only specified header lines
        Set<VCFHeaderLine> headerLines = new HashSet<>();
        int skippedHeaderLines = 0;

        //strip info annotations
        for (VCFInfoHeaderLine line : initialHeader.getInfoHeaderLines()) {
            skippedHeaderLines += inspectAnnotation(headerLines, line, annotationToKeep, annotationsToExclude);
        }

        //strip format fields
        for (VCFFormatHeaderLine line : initialHeader.getFormatHeaderLines()) {
            skippedHeaderLines += inspectAnnotation(headerLines, line, genotypeAnnotationToKeep, genotypeAnnotationsToExclude);
        }

        //strip filters, if selected
        if (!excludeFiltered){
            headerLines.addAll(initialHeader.getFilterLines());
        }

        //metadata
        if (retainExtraHeaderLines){
            headerLines.addAll(initialHeader.getOtherHeaderLines());
        }
        else {
            skippedHeaderLines += initialHeader.getOtherHeaderLines().size();
        }

        headerLines.addAll(initialHeader.getContigLines());

        logger.info("total header lines skipped: " + skippedHeaderLines);

        header = new VCFHeader(headerLines, (sitesOnly ? Collections.emptySet() : samples));

        allowableInfoKeys = new HashSet<>();
        //NOTE: if lambdas are used here, the walker will not be picked up by PluginManager
        // http://gatkforums.broadinstitute.org/gatk/discussion/comment/38892#Comment_38892
        for (VCFHeaderLine line : header.getInfoHeaderLines()){
            allowableInfoKeys.add(line.getKey());
        }

        allowableFormatKeys = new HashSet<>();
        for (VCFHeaderLine line : header.getFormatHeaderLines()){
            allowableFormatKeys.add(line.getKey());
        }

        vcfWriter.writeHeader(header);
    }

    private int inspectAnnotation(Set<VCFHeaderLine> headerLines, VCFHeaderLine line, List<String> annotationToKeep, List<String> annotationsToExclude){
        if (annotationToKeep != null && !annotationToKeep.contains(line.getKey())){
            return 1;
        }

        if (annotationsToExclude != null && annotationsToExclude.contains(line.getKey())){
            return 1;
        }

        headerLines.add(line);
        return 0;
    }

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return null;
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return 0;

        Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
        if (vcs == null || vcs.isEmpty())
            return 0;

        for (VariantContext vc : vcs){
            if (excludeFiltered && vc.isFiltered()){
                continue;
            }

            VariantContextBuilder vcb = new VariantContextBuilder(vc);

            //retain only allowable info fields
            Set<String> keys = new HashSet<>(vc.getAttributes().keySet());
            keys.removeAll(allowableInfoKeys);
            vcb.rmAttributes(new ArrayList<>(keys));

            //now genotypes:
            if (sitesOnly){
                vcb.noGenotypes();
            }
            else {
                GenotypesContext ctx = vc.getGenotypes();
                Set<String> sampleNames = ctx.getSampleNames();
                for (String sn : sampleNames){
                    Genotype g = ctx.get(sn);
                    GenotypeBuilder gb = new GenotypeBuilder(g);
                    if (clearGTfilter){
                        gb.unfiltered();
                    }

                    //retain only allowable info fields
                    Map<String, Object> gtAttributes = new HashMap<>(g.getExtendedAttributes());
                    gtAttributes.keySet().retainAll(allowableFormatKeys);
                    gb.attributes(gtAttributes);

                    ctx.replace(gb.make());
                }
                vcb.genotypes(ctx);
            }

            vcfWriter.add(vcb.make());
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }
}
