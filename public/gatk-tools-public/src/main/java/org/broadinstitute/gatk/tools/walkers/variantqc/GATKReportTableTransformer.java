package org.broadinstitute.gatk.tools.walkers.variantqc;

import org.broadinstitute.gatk.utils.report.GATKReportTable;

/**
 * Created by bimber on 5/31/2017.
 */
public interface GATKReportTableTransformer {
    public String getEvalModuleName();

    public GATKReportTable transform(GATKReportTable table);
}
