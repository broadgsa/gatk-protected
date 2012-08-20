package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.EnumMap;
import java.util.Map;

/**
 * An object to keep track of the number of occurences of each base and it's quality.
 *
 * User: depristo
 * Date: 4/8/11
 * Time: 2:55 PM
 */

 public class BaseCounts {
    public final static BaseIndex MAX_BASE_INDEX_WITH_NO_COUNTS = BaseIndex.N;
    public final static byte MAX_BASE_WITH_NO_COUNTS = MAX_BASE_INDEX_WITH_NO_COUNTS.getByte();

    private final Map<BaseIndex, Integer> counts;   // keeps track of the base counts
    private final Map<BaseIndex, Long> sumQuals;    // keeps track of the quals of each base

    public BaseCounts() {
        counts = new EnumMap<BaseIndex, Integer>(BaseIndex.class);
        sumQuals = new EnumMap<BaseIndex, Long>(BaseIndex.class);
        for (BaseIndex i : BaseIndex.values()) {
            counts.put(i, 0);
            sumQuals.put(i, 0L);
        }
    }

    public static BaseCounts createWithCounts(int[] countsACGT) {
        BaseCounts baseCounts = new BaseCounts();
        baseCounts.counts.put(BaseIndex.A, countsACGT[0]);
        baseCounts.counts.put(BaseIndex.C, countsACGT[1]);
        baseCounts.counts.put(BaseIndex.G, countsACGT[2]);
        baseCounts.counts.put(BaseIndex.T, countsACGT[3]);
        return baseCounts;
    }

    @Requires("other != null")
    public void add(BaseCounts other) {
        for (BaseIndex i : BaseIndex.values())
            counts.put(i, counts.get(i) + other.counts.get(i));
    }

    @Requires("other != null")
    public void sub(BaseCounts other) {
        for (BaseIndex i : BaseIndex.values())
            counts.put(i, counts.get(i) - other.counts.get(i));
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(byte base) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) // no Ns
            counts.put(i, counts.get(i) + 1);
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(byte base, byte qual) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) { // no Ns
            counts.put(i, counts.get(i) + 1);
            sumQuals.put(i, sumQuals.get(i) + qual);
        }
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) - 1")
    public void decr(byte base) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) // no Ns
            counts.put(i, counts.get(i) - 1);
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) - 1")
    public void decr(byte base, byte qual) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if (i != null) { // no Ns
            counts.put(i, counts.get(i) - 1);
            sumQuals.put(i, sumQuals.get(i) - qual);
        }
    }



    @Ensures("result >= 0")
    public int getCount(byte base) {
        return getCount(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public int getCount(BaseIndex base) {
        return counts.get(base);
    }

    @Ensures("result >= 0")
    public long getSumQuals(byte base) {
        return getSumQuals(BaseIndex.byteToBase(base));
    }

    @Ensures("result >= 0")
    public long getSumQuals(BaseIndex base) {
        return sumQuals.get(base);
    }

    @Ensures("result >= 0")
    public byte averageQuals(byte base) {
        return (byte) (getSumQuals(base) / getCount(base));
    }

    @Ensures("result >= 0")
    public byte averageQuals(BaseIndex base) {
        return (byte) (getSumQuals(base) / getCount(base));
    }

    public byte baseWithMostCounts() {
        return baseIndexWithMostCounts().getByte();
    }

    @Ensures("result >= 0")
    public int countOfMostCommonBase() {
        return counts.get(baseIndexWithMostCounts());
    }

    @Ensures("result >= 0")
    public long sumQualsOfMostCommonBase() {
        return sumQuals.get(baseIndexWithMostCounts());
    }

    @Ensures("result >= 0")
    public byte averageQualsOfMostCommonBase() {
        return (byte) (sumQualsOfMostCommonBase() / countOfMostCommonBase());
    }


    @Ensures("result >= 0")
    public int totalCount() {
        int sum = 0;
        for (int c : counts.values())
            sum += c;

        return sum;
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     *
     * @param base
     * @return the proportion of this base over all other bases
     */
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportion(byte base) {
        return (double) counts.get(BaseIndex.byteToBase(base)) / totalCount();
    }

    /**
     * Given a base , it returns the proportional count of this base compared to all other bases
     *
     * @param baseIndex
     * @return the proportion of this base over all other bases
     */
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportion(BaseIndex baseIndex) {
        int total = totalCount();
        if (total == 0)
            return 0.0;
        return (double) counts.get(baseIndex) / totalCount();
    }


    @Ensures("result != null")
    public String toString() {
        StringBuilder b = new StringBuilder();
        for (Map.Entry<BaseIndex, Integer> elt : counts.entrySet()) {
            b.append(elt.toString()).append("=").append(elt.getValue()).append(",");
        }
        return b.toString();
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostCounts() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (BaseIndex i : counts.keySet())
            if (counts.get(i) > counts.get(maxI))
                maxI = i;
        return maxI;
    }

    @Ensures("result != null")
    public BaseIndex baseIndexWithMostCountsWithoutIndels() {
        BaseIndex mostCounts = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for (BaseIndex index : counts.keySet())
            if (index.isNucleotide() && counts.get(index) > counts.get(mostCounts))
                mostCounts = index;
        return mostCounts;
    }

    @Ensures("result >=0")
    public int totalCountWithoutIndels() {
        int sum = 0;
        for (BaseIndex index : counts.keySet())
            if (index.isNucleotide())
                sum += counts.get(index);
        return sum;
    }

    /**
     * Calculates the proportional count of a base compared to all other bases except indels (I and D)
     *
     * @param index
     * @return the proportion of this base over all other bases except indels
     */
    @Requires("index.isNucleotide()")
    @Ensures({"result >=0.0", "result<= 1.0"})
    public double baseCountProportionWithoutIndels(BaseIndex index) {
        int total = totalCountWithoutIndels();
        if (total == 0)
            return 0.0;
        return (double) counts.get(index) / totalCountWithoutIndels();
    }
}
