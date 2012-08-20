package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import java.util.HashMap;
import java.util.Map;

/**
 * An object that keeps track of the base counts as well as the sum of the base, insertion and deletion qualities of each base.
 *
 * @author Mauricio Carneiro
 * @since 6/15/12
 */
public class BaseAndQualsCounts extends BaseCounts {
    private final Map<BaseIndex, Long> sumInsertionQuals;
    private final Map<BaseIndex, Long> sumDeletionQuals;

    public BaseAndQualsCounts() {
        super();
        this.sumInsertionQuals = new HashMap<BaseIndex, Long>();
        this.sumDeletionQuals  = new HashMap<BaseIndex, Long>();
        for (BaseIndex i : BaseIndex.values()) {
            sumInsertionQuals.put(i, 0L);
            sumDeletionQuals.put(i, 0L);
        }
    }

    public void incr(byte base, byte baseQual, byte insQual, byte delQual) {
        super.incr(base, baseQual);
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) {                                                                                                // do not allow Ns
            sumInsertionQuals.put(i, sumInsertionQuals.get(i) + insQual);
            sumDeletionQuals.put(i, sumDeletionQuals.get(i) + delQual);
        }
    }

    public void decr(byte base, byte baseQual, byte insQual, byte delQual) {
        super.decr(base, baseQual);
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) {                                                                                                // do not allow Ns
            sumInsertionQuals.put(i, sumInsertionQuals.get(i) - insQual);
            sumDeletionQuals.put(i, sumDeletionQuals.get(i) - delQual);
        }
    }

    public byte averageInsertionQualsOfMostCommonBase() {
        return getGenericAverageQualOfMostCommonBase(sumInsertionQuals);
    }

    public byte averageDeletionQualsOfMostCommonBase() {
        return getGenericAverageQualOfMostCommonBase(sumDeletionQuals);
    }

    private byte getGenericAverageQualOfMostCommonBase(Map<BaseIndex, Long> sumQuals) {
        BaseIndex base = BaseIndex.byteToBase(baseWithMostCounts());
        return (byte) (sumQuals.get(base) / getCount(base));
    }
}
