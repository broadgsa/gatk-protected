package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

/**
 * An object that keeps track of the base counts as well as the sum of the base, insertion and deletion qualities of each base.
 *
 * @author Mauricio Carneiro
 * @since 6/15/12
 */
public class BaseAndQualsCounts extends BaseCounts {
    private final long[] sumInsertionQuals;
    private final long[] sumDeletionQuals;

    public BaseAndQualsCounts() {
        super();
        this.sumInsertionQuals = new long[BaseIndex.values().length];
        this.sumDeletionQuals  = new long[BaseIndex.values().length];
        for (final BaseIndex i : BaseIndex.values()) {
            sumInsertionQuals[i.index] = 0L;
            sumDeletionQuals[i.index] = 0L;
        }
    }

    public void incr(final byte base, final byte baseQual, final byte insQual, final byte delQual) {
        final BaseIndex i = BaseIndex.byteToBase(base);
        super.incr(i, baseQual);
        sumInsertionQuals[i.index] += insQual;
        sumDeletionQuals[i.index] += delQual;
    }

    public void decr(final byte base, final byte baseQual, final byte insQual, final byte delQual) {
        final BaseIndex i = BaseIndex.byteToBase(base);
        super.decr(i, baseQual);
        sumInsertionQuals[i.index] -= insQual;
        sumDeletionQuals[i.index] -= delQual;
    }

    public byte averageInsertionQualsOfBase(final BaseIndex base) {
        return getGenericAverageQualOfBase(base, sumInsertionQuals);
    }

    public byte averageDeletionQualsOfBase(final BaseIndex base) {
        return getGenericAverageQualOfBase(base, sumDeletionQuals);
    }

    private byte getGenericAverageQualOfBase(final BaseIndex base, final long[] sumQuals) {
        return (byte) (sumQuals[base.index] / countOfBase(base));
    }
}
