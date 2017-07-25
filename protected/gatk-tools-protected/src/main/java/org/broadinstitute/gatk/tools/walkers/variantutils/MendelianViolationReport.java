package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.tools.walkers.annotator.MendelianViolationCount;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by bimber on 6/21/2017.
 */
@By(DataSource.REFERENCE)
public class MendelianViolationReport extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output
    PrintStream out;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered sites", required=false)
    protected boolean excludeFiltered = true;

    @Argument(fullName="violationReportThreshold", shortName="rt", doc="Any sample with more than this many MVs will be reported", required=false)
    protected long violationReportThreshold = 500L;

    private Map<String, MVSummary> sampleMap;
    private SampleDB sampleDB = null;

    @Override
    public void initialize () {
        super.initialize();

        sampleMap = new HashMap<>();
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Collections.singletonList(variantCollection.variants.getName()));
        Set<String> samples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        for (String sample : samples){
            sampleMap.put(sample, new MVSummary());
        }

        if ( sampleDB == null ) {
            sampleDB = getSampleDB();
        }
    }

    private static class MVSummary {
        long violationsDad = 0L;
        long violationsMom = 0L;
        long violationsTogether = 0L;
        long totalViolations = 0L;
        long totalCalled = 0L;

        public void addMV(MendelianViolationCount.MV mv, Genotype g){
            if (g.isCalled()) {
                totalCalled++;
            }

            if (mv == null){
                return;
            }

            if (mv.motherIsViolation){
                violationsMom++;
            }

            if (mv.fatherIsViolation){
                violationsDad++;
            }

            if (mv.violationCombined){
                violationsTogether++;
            }

            if (mv.isViolation()){
                totalViolations++;
            }
        }
    }

    @Override
    public Integer treeReduce (Integer lhs, Integer rhs){
        return null;
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return 0;

        Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
        if (vcs == null || vcs.isEmpty())
            return 0;

        for (VariantContext vc : vcs) {
            if (excludeFiltered && vc.isFiltered()) {
                continue;
            }

            //not autosome
            if (vc.getContig().equalsIgnoreCase("chrUn")){
                return null;
            }

            for (String sample : vc.getSampleNames()){
                Genotype g = vc.getGenotype(sample);
                if (g != null){
                    Sample s = sampleDB.getSample(g);
                    if (s != null){
                        MendelianViolationCount.MV mv = MendelianViolationCount.getMendelianViolation(s, vc, -1.0);
                        sampleMap.get(sample).addMV(mv, g);
                    }
                }
            }
        }

        return 0;
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
        out.println(StringUtils.join(Arrays.asList("SampleName", "TotalCalled", "TotalViolations", "MotherInconsistent", "FatherInconsistent", "InconsistentCombined", "Mother", "MotherHasData", "MotherMVs", "Father", "FatherHasData", "FatherMVs"), "\t"));

        Set<String> samplesReported = new HashSet<>();
        Set<String> additionalSamplesToReport = new HashSet<>();
        for (String sn : sampleMap.keySet()){
            MVSummary summary = sampleMap.get(sn);
            if (summary.totalViolations > violationReportThreshold){
                samplesReported.add(sn);
                reportSample(sn, summary, additionalSamplesToReport);
            }
        }

        additionalSamplesToReport.removeAll(samplesReported);

        for (String sn : additionalSamplesToReport){
            MVSummary summary = sampleMap.get(sn);
            if (summary != null){
                reportSample(sn, summary, new HashSet<>());
            }
        }
    }

    private void reportSample(String sn, MVSummary summary, Set<String> additionalSamplesToReport){
        Sample sample = sampleDB.getSample(sn);
        List<String> line = new ArrayList<>(Arrays.asList(
                sn,
                String.valueOf(summary.totalCalled),
                String.valueOf(summary.totalViolations),
                String.valueOf(summary.violationsMom),
                String.valueOf(summary.violationsDad),
                String.valueOf(summary.violationsTogether)
        ));

        appendParentToLine(sample.getMother(), line, additionalSamplesToReport);
        appendParentToLine(sample.getFather(), line, additionalSamplesToReport);

        out.println(StringUtils.join(line, "\t"));
    }

    private void appendParentToLine(Sample parent, List<String> line, Set<String> additionalSamplesToReport){
        if (parent == null){
            line.add("Unknown");
            line.add("false");
            line.add("");
        }
        else {
            line.add(parent.getID());
            MVSummary summary = sampleMap.get(parent.getID());
            line.add(String.valueOf(summary != null));
            line.add(summary == null ? "" : String.valueOf(summary.totalViolations));
            additionalSamplesToReport.add(parent.getID());
        }
    }
}