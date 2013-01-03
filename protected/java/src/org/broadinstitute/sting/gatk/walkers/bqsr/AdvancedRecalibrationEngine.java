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

import org.broadinstitute.sting.utils.classloader.ProtectedPackageSource;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.recalibration.RecalDatum;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class AdvancedRecalibrationEngine extends StandardRecalibrationEngine implements ProtectedPackageSource {
    @Override
    public void updateDataForRead( final ReadRecalibrationInfo recalInfo ) {
        final GATKSAMRecord read = recalInfo.getRead();
        final ReadCovariates readCovariates = recalInfo.getCovariatesValues();
        final RecalibrationTables tables = getRecalibrationTables();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = tables.getQualityScoreTable();

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

                        incrementDatumOrPutIfNecessary(tables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
                    }
                }
            }
        }
    }
}
