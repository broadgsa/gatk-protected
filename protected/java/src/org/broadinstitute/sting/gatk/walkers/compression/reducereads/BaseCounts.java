package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;


/**
 * An object to keep track of the number of occurrences of each base and it's quality.
 *
 * User: depristo
 * Date: 4/8/11
 * Time: 2:55 PM
 */

 public class BaseCounts {
    public final static BaseIndex MAX_BASE_INDEX_WITH_NO_COUNTS = BaseIndex.N;
    public final static byte MAX_BASE_WITH_NO_COUNTS = MAX_BASE_INDEX_WITH_NO_COUNTS.getByte();

    private final int[] counts;       // keeps track of the base counts
    private final long[] sumQuals;    // keeps track of the quals of each base
    private int totalCount = 0;       // keeps track of total count since this is requested so often

    public BaseCounts() {
        counts = new int[BaseIndex.values().length];
        sumQuals = new long[BaseIndex.values().length];
        for (final BaseIndex i : BaseIndex.values()) {
            counts[i.index] = 0;
            sumQuals[i.index] = 0L;
        }
    }

    public static BaseCounts createWithCounts(int[] countsACGT) {
        BaseCounts baseCounts = new BaseCounts();
        baseCounts.counts[BaseIndex.A.index] = countsACGT[0];
        baseCounts.counts[BaseIndex.C.index] = countsACGT[1];
        baseCounts.counts[BaseIndex.G.index] = countsACGT[2];
        baseCounts.counts[BaseIndex.T.index] = countsACGT[3];
        baseCounts.totalCount = countsACGT[0] + countsACGT[1] + countsACGT[2] + countsACGT[3];
        return baseCounts;
    }

    @Requires("other != null")
    public void add(final BaseCounts other) {
        for (final BaseIndex i : BaseIndex.values()) {
            final int otherCount = other.counts[i.index];
            counts[i.index] += otherCount;
            totalCount += otherCount;
        }
    }

    @Requires("other != null")
    public void sub(final BaseCounts other) {
        for (final BaseIndex i : BaseIndex.values()) {
            final int otherCount = other.counts[i.index];
            counts[i.index] -= otherCount;
            totalCount -= otherCount;
        }
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(final byte base) {
        final BaseIndex i = BaseIndex.byteToBase(base);
        counts[i.index]++;
        totalCount++;
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(final BaseIndex base, final byte qual) {
        counts[base.index]++;
        totalCount++;
        sumQuals[base.index] += qual;
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) - 1")
    public void decr(final byte base) {
        final BaseIndex i = BaseIndex.byteToBase(base);
        counts[i.index]--;
        totalCount--;
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) - 1")
    public void decr(final BaseIndex base, final byte qual) {
        counts[base.index]--;
        totalCount--;
        sumQuals[base.index] -= qual;
    }

    @Ensures("result >= 0")
    public long getSumQuals(final byte base) {
        return getSumQuals(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public long getSumQuals(final BaseIndex base) {
        return sumQuals[base.index];
    }

    @Ensures("result >= 0")
    public byte averageQuals(final byte base) {
        return (byte) (getSumQuals(base) / countOfBase(base));
    }

    @Ensures("result >= 0")
    public byte averageQuals(final BaseIndex base) {
        return (byte) (getSumQuals(base) / countOfBase(base));
    }

    @Ensures("result >= 0")
    public int countOfBase(final byte base) {
        return countOfBase(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public int countOfBase(final BaseIndex base) {
        return counts[base.index];
    }

    @Ensures("result >= 0")
    public long sumQualsOfBase(final BaseIndex base) {
        return sumQuals[base.index];
    }

    @Ensures("result >= 0")
    public byte averageQualsOfBase(final BaseIndex base) {
        return (byte) (sumQualsOfBase(base) / countOfBase(base));
    }


    @Ensures("result >= 0")
    public int totalCount() {
        return totalCount;
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     *
     * @param base     base
     * @return the proportion of this base over all other bases
     */
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportion(final byte base) {
        return baseCountProportion(BaseIndex.byteToBase(base));
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     *
     * @param baseIndex    base
     * @return the proportion of this base over all other bases
     */
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportion(final BaseIndex baseIndex) {
        return (totalCount == 0) ? 0.0 : (double)counts[baseIndex.index] / (double)totalCount;
    }

    @Ensures("result != null")
    public String toString() {
        StringBuilder b = new StringBuilder();
        for (final BaseIndex i : BaseIndex.values()) {
            b.append(i.toString()).append("=").append(counts[i.index]).append(",");
        }
        return b.toString();
    }

    public byte baseWithMostCounts() {
        return baseIndexWithMostCounts().getByte();
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostCounts() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (counts[i.index] > counts[maxI.index])
                maxI = i;
        }
        return maxI;
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostCountsWithoutIndels() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (i.isNucleotide() && counts[i.index] > counts[maxI.index])
                maxI = i;
        }
        return maxI;
    }

    private boolean hasHigherCount(final BaseIndex targetIndex, final BaseIndex testIndex) {
        final int targetCount = counts[targetIndex.index];
        final int testCount = counts[testIndex.index];
        return  ( targetCount > testCount || (targetCount == testCount && sumQuals[targetIndex.index] > sumQuals[testIndex.index]) );
    }

    public byte baseWithMostProbability() {
        return baseIndexWithMostProbability().getByte();
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostProbability() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (sumQuals[i.index] > sumQuals[maxI.index])
                maxI = i;
        }
        return (sumQuals[maxI.index] > 0L ? maxI : baseIndexWithMostCounts());
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostProbabilityWithoutIndels() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (final BaseIndex i : BaseIndex.values()) {
            if (i.isNucleotide() && sumQuals[i.index] > sumQuals[maxI.index])
                maxI = i;
        }
        return (sumQuals[maxI.index] > 0L ? maxI : baseIndexWithMostCountsWithoutIndels());
    }

    @Ensures("result >=0")
    public int totalCountWithoutIndels() {
        return totalCount - counts[BaseIndex.D.index] - counts[BaseIndex.I.index];
    }

    /**
     * Calculates the proportional count of a base compared to all other bases except indels (I and D)
     *
     * @param base      base
     * @return the proportion of this base over all other bases except indels
     */
    @Requires("index.isNucleotide()")
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportionWithoutIndels(final BaseIndex base) {
        final int total = totalCountWithoutIndels();
        return (total == 0) ? 0.0 : (double)counts[base.index] / (double)total;
    }

    public int[] countsArray() {
        return counts.clone();
    }
}
