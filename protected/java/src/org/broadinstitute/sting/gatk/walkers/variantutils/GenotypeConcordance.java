package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * A simple walker for performing genotype concordance calculations between two callsets
 */
public class GenotypeConcordance extends RodWalker<Pair<VariantContext,VariantContext>,ConcordanceMetrics> {

    @Input(fullName="eval",shortName="eval",doc="The variants and genotypes to evaluate",required=true)
    RodBinding<VariantContext> evalBinding;

    @Input(fullName="comp",shortName="comp",doc="The variants and genotypes to compare against",required=true)
    RodBinding<VariantContext> compBinding;

    @Argument(fullName="ignoreFilters",doc="Filters will be ignored",required=false)
    boolean ignoreFilters = false;

    @Output
    PrintStream out;

    List<String> evalSamples;
    List<String> compSamples;

    // todo -- integration test coverage
    // todo -- deal with occurrences like:
    //     Eval: 20   4000     A     C
    //     Eval: 20   4000     A    AC
    //     Comp: 20   4000     A     C
    //  currently this results in a warning and skipping
    // todo -- extend to multiple eval, multiple comp
    // todo -- table with "proportion of overlapping sites" (not just eval/comp margins)


    public ConcordanceMetrics reduceInit() {
        Map<String,VCFHeader> headerMap = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(evalBinding,compBinding));
        VCFHeader evalHeader = headerMap.get(evalBinding.getName());
        evalSamples = evalHeader.getGenotypeSamples();
        VCFHeader compHeader = headerMap.get(compBinding.getName());
        compSamples = compHeader.getGenotypeSamples();
        return new ConcordanceMetrics(evalHeader,compHeader);
    }


    public Pair<VariantContext,VariantContext> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Pair<VariantContext,VariantContext> evalCompPair = null;
        if ( tracker != null && (
                tracker.getValues(evalBinding,ref.getLocus()).size() > 0 ||
                tracker.getValues(compBinding,ref.getLocus()).size() > 0 ) ) {

            List<VariantContext> eval = tracker.getValues(evalBinding,ref.getLocus());
            List<VariantContext> comp = tracker.getValues(compBinding,ref.getLocus());
            if ( eval.size() > 1 || comp.size() > 1 ) {
                logger.warn("Eval or Comp Rod at position "+ref.getLocus().toString()+" has multiple records. Site will be skipped.");
                return evalCompPair;
            }
            // if a rod is missing, explicitly create a variant context with 'missing' genotypes. Slow, but correct.
            // note that if there is no eval rod there must be a comp rod, and also the reverse
            VariantContext evalContext = eval.size() == 1 ? eval.get(0) : createEmptyContext(ref,comp.get(0),evalSamples);
            VariantContext compContext = comp.size() == 1 ? comp.get(0) : createEmptyContext(ref,eval.get(0),compSamples);
            evalContext = filterGenotypes(evalContext,ignoreFilters);
            compContext = filterGenotypes(compContext,ignoreFilters);
            evalCompPair = new Pair<VariantContext, VariantContext>(evalContext,compContext);
        }

        return evalCompPair;
    }

    public ConcordanceMetrics reduce(Pair<VariantContext,VariantContext> evalComp, ConcordanceMetrics metrics) {
        if ( evalComp != null )
            metrics.update(evalComp.getFirst(),evalComp.getSecond());
        return metrics;
    }

    public void onTraversalDone(ConcordanceMetrics metrics) {
        GATKReport report = new GATKReport();
        GATKReportTable concordanceCounts = new GATKReportTable("GenotypeConcordance_Counts","Per-sample concordance tables: comparison counts",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceEvalProportions = new GATKReportTable("GenotypeConcordance_EvalProportions", "Per-sample concordance tables: proportions of genotypes called in eval",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceCompProportions = new GATKReportTable("GenotypeConcordance_CompProportions", "Per-sample concordance tables: proportions of genotypes called in comp",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceSummary = new GATKReportTable("GenotypeConcordance_Summary","Per-sample summary statistics: NRS and NRD",2);
        GATKReportTable siteConcordance = new GATKReportTable("SiteConcordance_Summary","Site-level summary statistics",ConcordanceMetrics.SiteConcordanceType.values().length);
        concordanceCompProportions.addColumn("Sample","%s");
        concordanceCounts.addColumn("Sample","%s");
        concordanceEvalProportions.addColumn("Sample","%s");
        concordanceSummary.addColumn("Sample","%s");
        for ( GenotypeType evalType : GenotypeType.values() ) {
            for ( GenotypeType compType : GenotypeType.values() ) {
                String colKey = String.format("%s_%s", evalType.toString(), compType.toString());
                concordanceCounts.addColumn(colKey,"%d");
                if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                    concordanceEvalProportions.addColumn(colKey,"%.3f");
                if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                    concordanceCompProportions.addColumn(colKey,"%.3f");
            }
        }
        concordanceEvalProportions.addColumn("Mismatching_Alleles","%.3f");
        concordanceCompProportions.addColumn("Mismatching_Alleles","%.3f");
        concordanceCounts.addColumn("Mismatching_Alleles","%d");
        concordanceSummary.addColumn("Non-Reference Sensitivity","%.3f");
        concordanceSummary.addColumn("Non-Reference Discrepancy","%.3f");
        for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
            siteConcordance.addColumn(type.toString(),"%d");
        }

        for ( Map.Entry<String,ConcordanceMetrics.GenotypeConcordanceTable> entry : metrics.getPerSampleGenotypeConcordance().entrySet() ) {
            ConcordanceMetrics.GenotypeConcordanceTable table = entry.getValue();
            concordanceEvalProportions.set(entry.getKey(),"Sample",entry.getKey());
            concordanceCompProportions.set(entry.getKey(),"Sample",entry.getKey());
            concordanceCounts.set(entry.getKey(),"Sample",entry.getKey());
            for ( GenotypeType evalType : GenotypeType.values() ) {
                for ( GenotypeType compType : GenotypeType.values() ) {
                    String colKey = String.format("%s_%s",evalType.toString(),compType.toString());
                    int count = table.get(evalType, compType);
                    concordanceCounts.set(entry.getKey(),colKey,count);
                    if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                        concordanceEvalProportions.set(entry.getKey(),colKey,( (double) count)/table.getnEvalGenotypes(evalType));
                    if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                        concordanceCompProportions.set(entry.getKey(),colKey,( (double) count)/table.getnCompGenotypes(compType));
                }
            }
            concordanceEvalProportions.set(entry.getKey(),"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes());
            concordanceCompProportions.set(entry.getKey(),"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes());
            concordanceCounts.set(entry.getKey(),"Mismatching_Alleles",table.getnMismatchingAlt());
        }

        String rowKey = "ALL";
        concordanceCompProportions.set(rowKey,"Sample",rowKey);
        concordanceEvalProportions.set(rowKey,"Sample",rowKey);
        concordanceCounts.set(rowKey,"Sample",rowKey);
        ConcordanceMetrics.GenotypeConcordanceTable table = metrics.getOverallGenotypeConcordance();
        for ( GenotypeType evalType : GenotypeType.values() ) {
            for ( GenotypeType compType : GenotypeType.values() ) {
                String colKey = String.format("%s_%s",evalType.toString(),compType.toString());
                int count = table.get(evalType,compType);
                concordanceCounts.set(rowKey,colKey,count);
                if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                    concordanceEvalProportions.set(rowKey,colKey,( (double) count)/table.getnEvalGenotypes(evalType));
                if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                    concordanceCompProportions.set(rowKey,colKey,( (double) count)/table.getnCompGenotypes(compType));
            }
        }
        concordanceEvalProportions.set(rowKey,"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes());
        concordanceCompProportions.set(rowKey,"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes());
        concordanceCounts.set(rowKey,"Mismatching_Alleles",table.getnMismatchingAlt());

        for ( Map.Entry<String,Double> nrsEntry : metrics.getPerSampleNRS().entrySet() ) {
            concordanceSummary.set(nrsEntry.getKey(),"Sample",nrsEntry.getKey());
            concordanceSummary.set(nrsEntry.getKey(),"Non-Reference Sensitivity",nrsEntry.getValue());
        }
        for ( Map.Entry<String,Double> nrdEntry : metrics.getPerSampleNRD().entrySet() ) {
            concordanceSummary.set(nrdEntry.getKey(),"Non-Reference Discrepancy",nrdEntry.getValue());
        }
        concordanceSummary.set("ALL","Sample","ALL");
        concordanceSummary.set("ALL","Non-Reference Sensitivity",metrics.getOverallNRS());
        concordanceSummary.set("ALL","Non-Reference Discrepancy",metrics.getOverallNRD());

        for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
            siteConcordance.set("Comparison",type.toString(),metrics.getOverallSiteConcordance().get(type));
        }

        report.addTable(concordanceCompProportions);
        report.addTable(concordanceEvalProportions);
        report.addTable(concordanceCounts);
        report.addTable(concordanceSummary);
        report.addTable(siteConcordance);

        report.print(out);
    }

    public VariantContext createEmptyContext(ReferenceContext ref, VariantContext other, List<String> samples) {
        VariantContextBuilder builder = new VariantContextBuilder();
        // set the alleles to be the same
        builder.alleles(other.getAlleles());
        builder.loc(other.getChr(),other.getStart(),other.getEnd());
        // set all genotypes to empty
        List<Genotype> genotypes = new ArrayList<Genotype>(samples.size());
        for ( String sample : samples )
            genotypes.add(GenotypeBuilder.create(sample, new ArrayList<Allele>(0)));
        builder.genotypes(genotypes);
        return builder.make();
    }

    public VariantContext filterGenotypes(VariantContext context, boolean ignoreSiteFilter) {
        // placeholder method for genotype-level filtering. However if the site itself is filtered,
        // and such filters are not ignored, the genotype-level data should be altered to reflect this
        if ( ! context.isFiltered() || ignoreSiteFilter ) {
            // todo -- add genotype-level jexl filtering here
            return context;
        }
        VariantContextBuilder builder = new VariantContextBuilder();
        builder.alleles(Arrays.asList(context.getReference()));
        builder.loc(context.getChr(),context.getStart(),context.getEnd());
        List<Genotype> newGeno = new ArrayList<Genotype>(context.getNSamples());
        for ( Genotype g : context.getGenotypes().iterateInSampleNameOrder() ) {
            newGeno.add(GenotypeBuilder.create(g.getSampleName(),new ArrayList<Allele>()));
        }
        builder.genotypes(newGeno);
        return builder.make();
    }
}