package org.broadinstitute.gatk.tools.walkers.variantqc;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.utils.classloader.JVMUtils;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.report.GATKReportTable;
import org.broadinstitute.gatk.utils.report.GATKReportVersion;

import java.io.*;
import java.lang.reflect.Field;
import java.util.*;

/**
 * Created by bimber on 5/4/2017.
 */
@Reference(window=@Window(start=-50, stop=50))
@PartitionBy(PartitionType.NONE)
public class VariantQC extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output
    protected PrintStream out;

    protected VariantEvalWrapper[] wrappers = new VariantEvalWrapper[]{
            //TODO: add MendelianViolationEvaluator, SiteFilterSummary?
            new VariantEvalWrapper("Entire VCF", new String[]{"EvalRod"}, new String[]{"CountVariants", "IndelSummary", "TiTvVariantEvaluator", "GenotypeFilterSummary"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(),
                    BarPlotReportDescriptor.getVariantTypeBarPlot(),
            }),
            new VariantEvalWrapper("By Contig", new String[]{"Contig"}, new String[]{"CountVariants", "IndelSummary", "TiTvVariantEvaluator", "GenotypeFilterSummary"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(),
                    BarPlotReportDescriptor.getVariantTypeBarPlot()
            }),
            new VariantEvalWrapper("By Sample", new String[]{"Sample"}, new String[]{"CountVariants", "IndelSummary", "TiTvVariantEvaluator", "GenotypeFilterSummary"}, new ReportDescriptor[]{
                    TableReportDescriptor.getCountVariantsTable(),
                    BarPlotReportDescriptor.getVariantTypeBarPlot()
            })
            //TODO: FilterType?
    };

    @Override
    public void initialize() {
        super.initialize();

        //configure the child walkers
        for (VariantEvalWrapper wrapper : this.wrappers){
            configureWalker(wrapper);
            
            wrapper.walker.initialize();
        }
    }
    
    private void configureWalker(VariantEvalWrapper wrapper){
        wrapper.walker.setToolkit(getToolkit());

        //manually set arguments
        try
        {
            Field evalField = VariantEval.class.getDeclaredField("evals");
            JVMUtils.setFieldValue(evalField, wrapper.walker, Arrays.asList(variantCollection.variants));

            Field outField = VariantEval.class.getDeclaredField("out");
            JVMUtils.setFieldValue(outField, wrapper.walker, new PrintStream(wrapper.out));

            Field stratificationsToUseField = VariantEval.class.getDeclaredField("MODULES_TO_USE");
            JVMUtils.setFieldValue(stratificationsToUseField, wrapper.walker, wrapper.evaluationModules);

            Field noevField = VariantEval.class.getDeclaredField("NO_STANDARD_MODULES");
            JVMUtils.setFieldValue(noevField, wrapper.walker, true);

            Field modulesField = VariantEval.class.getDeclaredField("STRATIFICATIONS_TO_USE");
            JVMUtils.setFieldValue(modulesField, wrapper.walker, wrapper.stratifications);

            Field noSTField = VariantEval.class.getDeclaredField("NO_STANDARD_STRATIFICATIONS");
            JVMUtils.setFieldValue(noSTField, wrapper.walker, true);

            //TODO: set unbound?
            Field dbSnpField= VariantEval.class.getDeclaredField("dbsnp");
            dbSnpField.setAccessible(true);

            DbsnpArgumentCollection dbsnpArgumentCollection1 = (DbsnpArgumentCollection)dbSnpField.get(wrapper.walker);
            dbsnpArgumentCollection1.dbsnp = new RodBinding<>(VariantContext.class, "", "UNBOUND", "", new Tags());
            JVMUtils.setFieldValue(dbSnpField, wrapper.walker, dbsnpArgumentCollection1);
        }
        catch (Exception e)
        {
            throw new GATKException(e.getMessage(), e);
        }
    }

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return null;
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for (VariantEvalWrapper wrapper : this.wrappers) {
            wrapper.walker.map(tracker, ref, context);
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

    @Override
    public void onTraversalDone(Integer result) {
        super.onTraversalDone(result);

        for (VariantEvalWrapper wrapper : this.wrappers) {
            wrapper.walker.onTraversalDone(result);
        }

        //make classes like this to translate from the GATKReportTable into the config object we need in our HTML
        Map<String, SectionJsonDescriptor> translatorMap = new TreeMap<>();

        for (VariantEvalWrapper wrapper : this.wrappers) {
            try (BufferedReader sampleReader = new BufferedReader(new StringReader(new String(wrapper.out.toByteArray())))) {
                sampleReader.readLine(); //read first GATKReport line

                for (String evalModule : wrapper.evaluationModules){
                    List<ReportDescriptor> rds = wrapper.getReportsForModule(evalModule);
                    GATKReportTable table = new GATKReportTable(sampleReader, GATKReportVersion.V1_1);
                    if (!translatorMap.containsKey(wrapper.sectionLabel)){
                        translatorMap.put(wrapper.sectionLabel, new SectionJsonDescriptor(wrapper.sectionLabel, wrapper.stratifications));
                    }

                    for (ReportDescriptor rd : rds){
                        translatorMap.get(wrapper.sectionLabel).addReportDescriptor(rd, table);
                    }


                }
            } catch (IOException e) {
                throw new GATKException(e.getMessage(), e);
            }
        }

        try {
            List<SectionJsonDescriptor> sections = new ArrayList<>();

            //NOTE: if lambdas are used here, the walker will not be picked up by PluginManager
            // http://gatkforums.broadinstitute.org/gatk/discussion/comment/38892#Comment_38892
            for (String key : translatorMap.keySet()) {
                sections.add(translatorMap.get(key));
            }

            HtmlGenerator generator = new HtmlGenerator();
            generator.generateHtml(sections, out);
        }
        catch (IOException e){
            throw new GATKException(e.getMessage(), e);
        }
    }
    
    private static class VariantEvalWrapper {
        private VariantEval walker = new VariantEval();
        private ByteArrayOutputStream out = new ByteArrayOutputStream();
        String[] stratifications;
        String[] evaluationModules;
        String sectionLabel;
        ReportDescriptor[] reportDescriptors;

        public VariantEvalWrapper(String sectionLabel, String[] stratifications, String[] evaluationModules, ReportDescriptor[] reportDescriptors)
        {
            this.stratifications = stratifications;
            this.evaluationModules = evaluationModules;
            this.sectionLabel = sectionLabel;

            this.reportDescriptors = reportDescriptors;
        }

        public List<ReportDescriptor> getReportsForModule(String evalModule){
            List<ReportDescriptor> ret = new ArrayList<>();
            for (ReportDescriptor rd : reportDescriptors){
                if (rd.evaluatorModuleName.equals(evalModule)){
                    ret.add(rd);
                }
            }

            return ret;
        }
    }
}
