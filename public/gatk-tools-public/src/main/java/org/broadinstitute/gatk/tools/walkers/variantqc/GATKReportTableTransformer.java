package org.broadinstitute.gatk.tools.walkers.variantqc;

import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import javax.annotation.Nullable;

/**
 * Created by bimber on 5/31/2017.
 */
public interface GATKReportTableTransformer {
    public String getEvalModuleName();

    public GATKReportTable transform(GATKReportTable table, @Nullable SampleDB sampleDB);
}
