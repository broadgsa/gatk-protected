package org.broadinstitute.sting.gatk.walkers.bqsr;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.classloader.ProtectedPackageSource;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.LinkedList;
import java.util.List;

public class AdvancedRecalibrationEngine extends StandardRecalibrationEngine implements ProtectedPackageSource {
    private final static Logger logger = Logger.getLogger(AdvancedRecalibrationEngine.class);

    final List<NestedIntegerArray<RecalDatum>> allThreadLocalQualityScoreTables = new LinkedList<NestedIntegerArray<RecalDatum>>();
    private ThreadLocal<NestedIntegerArray<RecalDatum>> threadLocalQualityScoreTables = new ThreadLocal<NestedIntegerArray<RecalDatum>>() {
        @Override
        protected synchronized NestedIntegerArray<RecalDatum> initialValue() {
            final NestedIntegerArray<RecalDatum> table = recalibrationTables.makeQualityScoreTable();
            allThreadLocalQualityScoreTables.add(table);
            return table;
        }
    };

    @Override
    public void updateDataForRead( final ReadRecalibrationInfo recalInfo ) {
        final GATKSAMRecord read = recalInfo.getRead();
        final ReadCovariates readCovariates = recalInfo.getCovariatesValues();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = getThreadLocalQualityScoreTable();

        for( int offset = 0; offset < read.getReadBases().length; offset++ ) {
            if( ! recalInfo.skip(offset) ) {

                for (final EventType eventType : EventType.values()) {
                    final int[] keys = readCovariates.getKeySet(offset, eventType);
                    final int eventIndex = eventType.index;
                    final byte qual = recalInfo.getQual(eventType, offset);
                    final double isError = recalInfo.getErrorFraction(eventType, offset);

                    incrementDatumOrPutIfNecessary(qualityScoreTable, qual, isError, keys[0], keys[1], eventIndex);

                    for (int i = 2; i < covariates.length; i++) {
                        if (keys[i] < 0)
                            continue;

                        incrementDatumOrPutIfNecessary(recalibrationTables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
                    }
                }
            }
        }
    }

    /**
     * Get a NestedIntegerArray for a QualityScore table specific to this thread
     * @return a non-null NestedIntegerArray ready to be used to collect calibration info for the quality score covariate
     */
    private NestedIntegerArray<RecalDatum> getThreadLocalQualityScoreTable() {
        return threadLocalQualityScoreTables.get();
    }

    @Override
    public void finalizeData() {
        // merge in all of the thread local tables
        logger.info("Merging " + allThreadLocalQualityScoreTables.size() + " thread-local quality score tables");
        for ( final NestedIntegerArray<RecalDatum> localTable : allThreadLocalQualityScoreTables ) {
            recalibrationTables.combineQualityScoreTable(localTable);
        }
        allThreadLocalQualityScoreTables.clear(); // cleanup after ourselves

        super.finalizeData();
    }
}
