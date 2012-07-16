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

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.classloader.ProtectedPackageSource;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;

public class AdvancedRecalibrationEngine extends RecalibrationEngine implements ProtectedPackageSource {

    // optimizations: don't reallocate an array each time
    private byte[] tempQualArray;
    private boolean[] tempErrorArray;

    public void initialize(final Covariate[] covariates, final RecalibrationTables recalibrationTables) {
        super.initialize(covariates, recalibrationTables);
        tempQualArray = new byte[EventType.values().length];
        tempErrorArray = new boolean[EventType.values().length];
    }

    /**
     * Loop through the list of requested covariates and pick out the value from the read, offset, and reference
     * Using the list of covariate values as a key, pick out the RecalDatum and increment,
     * adding one to the number of observations and potentially one to the number of mismatches for all three
     * categories (mismatches, insertions and deletions).
     *
     * @param pileupElement The pileup element to update
     * @param refBase       The reference base at this locus
     */
    public synchronized void updateDataForPileupElement(final PileupElement pileupElement, final byte refBase) {
        final int offset = pileupElement.getOffset();
        final ReadCovariates readCovariates = covariateKeySetFrom(pileupElement.getRead());

        tempQualArray[EventType.BASE_SUBSTITUTION.index] = pileupElement.getQual();
        tempErrorArray[EventType.BASE_SUBSTITUTION.index] = !BaseUtils.basesAreEqual(pileupElement.getBase(), refBase);
        tempQualArray[EventType.BASE_INSERTION.index] = pileupElement.getBaseInsertionQual();
        tempErrorArray[EventType.BASE_INSERTION.index] = (pileupElement.getRead().getReadNegativeStrandFlag()) ? pileupElement.isAfterInsertion() : pileupElement.isBeforeInsertion();
        tempQualArray[EventType.BASE_DELETION.index] = pileupElement.getBaseDeletionQual();
        tempErrorArray[EventType.BASE_DELETION.index] = (pileupElement.getRead().getReadNegativeStrandFlag()) ? pileupElement.isAfterDeletedBase() : pileupElement.isBeforeDeletedBase();

        for (final EventType eventType : EventType.values()) {
            final int[] keys = readCovariates.getKeySet(offset, eventType);
            final int eventIndex = eventType.index;
            final byte qual = tempQualArray[eventIndex];
            final boolean isError = tempErrorArray[eventIndex];

            final NestedIntegerArray<RecalDatum> rgRecalTable = recalibrationTables.getTable(RecalibrationTables.TableType.READ_GROUP_TABLE);
            final RecalDatum rgPreviousDatum = rgRecalTable.get(keys[0], eventIndex);
            final RecalDatum rgThisDatum = createDatumObject(qual, isError);
            if (rgPreviousDatum == null)                                                                                // key doesn't exist yet in the map so make a new bucket and add it
                rgRecalTable.put(rgThisDatum, keys[0], eventIndex);
            else
                rgPreviousDatum.combine(rgThisDatum);

            final NestedIntegerArray<RecalDatum> qualRecalTable = recalibrationTables.getTable(RecalibrationTables.TableType.QUALITY_SCORE_TABLE);
            final RecalDatum qualPreviousDatum = qualRecalTable.get(keys[0], keys[1], eventIndex);
            if (qualPreviousDatum == null)
                qualRecalTable.put(createDatumObject(qual, isError), keys[0], keys[1], eventIndex);
            else
                qualPreviousDatum.increment(isError);

            for (int i = 2; i < covariates.length; i++) {
                if (keys[i] < 0)
                    continue;
                final NestedIntegerArray<RecalDatum> covRecalTable = recalibrationTables.getTable(i);
                final RecalDatum covPreviousDatum = covRecalTable.get(keys[0], keys[1], keys[i], eventIndex);
                if (covPreviousDatum == null)
                    covRecalTable.put(createDatumObject(qual, isError), keys[0], keys[1], keys[i], eventIndex);
                else
                    covPreviousDatum.increment(isError);
            }
        }
    }
}
