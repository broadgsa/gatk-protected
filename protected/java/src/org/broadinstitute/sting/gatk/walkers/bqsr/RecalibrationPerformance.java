/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.recalibration.*;

import java.io.*;

/**
 */

@ReadFilters({MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class, UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class})
@PartitionBy(PartitionType.READ)
public class RecalibrationPerformance extends RodWalker<Integer, Integer> implements NanoSchedulable {

    @Output(doc="Write output to this file", required = true)
    public PrintStream out;

    @Input(fullName="recal", shortName="recal", required=false, doc="The input covariates table file")
    public File RECAL_FILE = null;

    public void initialize() {
        out.println("Cycle\tQrep\tQemp\tIsJoint\tObservations\tErrors");

        final GATKReport report = new GATKReport(RECAL_FILE);
        final GATKReportTable table = report.getTable(RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE);
        for ( int row = 0; row < table.getNumRows(); row++ ) {

            final int nObservations = (int)asDouble(table.get(row, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME));
            final int nErrors = (int)Math.round(asDouble(table.get(row, RecalUtils.NUMBER_ERRORS_COLUMN_NAME)));
            final double empiricalQuality = asDouble(table.get(row, RecalUtils.EMPIRICAL_QUALITY_COLUMN_NAME));

            final byte QReported = Byte.parseByte((String) table.get(row, RecalUtils.QUALITY_SCORE_COLUMN_NAME));

            final double jointEstimateQemp = RecalDatum.bayesianEstimateOfEmpiricalQuality(nObservations, nErrors, QReported);

            //if ( Math.abs((int)(jointEstimateQemp - empiricalQuality)) > 1 )
            //    System.out.println(String.format("Qreported = %f, nObservations = %f, nErrors = %f, point Qemp = %f, joint Qemp = %f", estimatedQReported, nObservations, nErrors, empiricalQuality, jointEstimateQemp));

            if ( table.get(row, RecalUtils.COVARIATE_NAME_COLUMN_NAME).equals("Cycle") &&
                    table.get(row, RecalUtils.EVENT_TYPE_COLUMN_NAME).equals("M") &&
                    table.get(row, RecalUtils.READGROUP_COLUMN_NAME).equals("20FUKAAXX100202.6") &&
                    (QReported == 6 || QReported == 10 || QReported == 20 || QReported == 30 || QReported == 45) ) {
                out.println(String.format("%s\t%d\t%d\t%s\t%d\t%d", table.get(row, RecalUtils.COVARIATE_VALUE_COLUMN_NAME), QReported, Math.round(empiricalQuality), "False", (int)nObservations, (int)nErrors));
                out.println(String.format("%s\t%d\t%d\t%s\t%d\t%d", table.get(row, RecalUtils.COVARIATE_VALUE_COLUMN_NAME), QReported, (int)jointEstimateQemp, "True", (int)nObservations, (int)nErrors));
            }
        }

    }

    @Override
    public boolean isDone() {
        return true;
    }

    private double asDouble(final Object o) {
        if ( o instanceof Double )
            return (Double)o;
        else if ( o instanceof Integer )
            return (Integer)o;
        else if ( o instanceof Long )
            return (Long)o;
        else
            throw new ReviewedStingException("Object " + o + " is expected to be either a double, long or integer but its not either: " + o.getClass());
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer counter, Integer sum) { return 0; }

    @Override
    public void onTraversalDone(Integer sum) {}
}