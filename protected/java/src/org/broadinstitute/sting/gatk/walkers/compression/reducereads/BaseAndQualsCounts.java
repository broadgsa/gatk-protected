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

    public void incr(final byte base, final byte baseQual, final byte insQual, final byte delQual) {
        super.incr(base, baseQual);
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) {                                                                                                // do not allow Ns
            sumInsertionQuals.put(i, sumInsertionQuals.get(i) + insQual);
            sumDeletionQuals.put(i, sumDeletionQuals.get(i) + delQual);
        }
    }

    public void decr(final byte base, final byte baseQual, final byte insQual, final byte delQual) {
        super.decr(base, baseQual);
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) {                                                                                                // do not allow Ns
            sumInsertionQuals.put(i, sumInsertionQuals.get(i) - insQual);
            sumDeletionQuals.put(i, sumDeletionQuals.get(i) - delQual);
        }
    }

    public byte averageInsertionQualsOfBase(final BaseIndex base) {
        return getGenericAverageQualOfBase(base, sumInsertionQuals);
    }

    public byte averageDeletionQualsOfBase(final BaseIndex base) {
        return getGenericAverageQualOfBase(base, sumDeletionQuals);
    }

    private byte getGenericAverageQualOfBase(final BaseIndex base, final Map<BaseIndex, Long> sumQuals) {
        return (byte) (sumQuals.get(base) / getCount(base));
    }
}
