package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.walkers.indels.HaplotypeIndelErrorModel;
import org.broadinstitute.sting.utils.MathUtils;

public class PoolGenotypePriors implements GenotypePriors {
    private final double[] flatPriors;
    private final double heterozygosity;
    private final int samplesPerPool;
    private double[] priors = null;

    /**
     * Create a new DiploidGenotypePriors object with flat priors for each diploid genotype
     */
    public PoolGenotypePriors(double heterozygosity, int samplesPerPool) {
        flatPriors = new double[2*samplesPerPool+1];
        for (int k=0; k <flatPriors.length; k++)
            flatPriors[k] = Math.log10(heterozygosity);
        priors = flatPriors.clone();
        this.samplesPerPool = samplesPerPool;
        
        this.heterozygosity = heterozygosity;
    }


    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors;
    }

    public double getHeterozygosity() { return heterozygosity; }
    public int getNSamplesPerPool() { return samplesPerPool; }

    public boolean validate(boolean throwException) {
        try {

            for (int i=0; i < priors.length; i++ ) {
                if ( ! MathUtils.wellFormedDouble(priors[i]) || ! MathUtils.isNegativeOrZero(priors[i]) ) {
                    String bad = String.format("Prior %f is badly formed %b", priors[i], MathUtils.isNegativeOrZero(priors[i]));
                    throw new IllegalStateException(String.format("At %d: %s", i, bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

}

