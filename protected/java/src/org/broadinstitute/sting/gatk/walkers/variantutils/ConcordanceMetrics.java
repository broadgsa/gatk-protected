package org.broadinstitute.sting.gatk.walkers.variantutils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.util.*;

/**
 * A class for tabulating and evaluating a callset-by-callset genotype concordance table
 * */
public class ConcordanceMetrics {

    private Map<String,GenotypeConcordanceTable> perSampleGenotypeConcordance;
    private GenotypeConcordanceTable overallGenotypeConcordance;
    private SiteConcordanceTable overallSiteConcordance;

    public ConcordanceMetrics(VCFHeader evaluate, VCFHeader truth) {
        HashSet<String> overlappingSamples = new HashSet<String>(evaluate.getGenotypeSamples());
        overlappingSamples.retainAll(truth.getGenotypeSamples());
        perSampleGenotypeConcordance = new HashMap<String, GenotypeConcordanceTable>(overlappingSamples.size());
        for ( String sample : overlappingSamples ) {
            perSampleGenotypeConcordance.put(sample,new GenotypeConcordanceTable());
        }
        overallGenotypeConcordance = new GenotypeConcordanceTable();
        overallSiteConcordance = new SiteConcordanceTable();
    }

    public GenotypeConcordanceTable getOverallGenotypeConcordance() {
        return overallGenotypeConcordance;
    }

    public SiteConcordanceTable getOverallSiteConcordance() {
        return overallSiteConcordance;
    }

    public GenotypeConcordanceTable getGenotypeConcordance(String sample) {
        GenotypeConcordanceTable table = perSampleGenotypeConcordance.get(sample);
        if ( table == null )
            throw new ReviewedStingException("Attempted to request the concordance table for a sample on which it was not calculated");
        return table;
    }

    public Map<String,GenotypeConcordanceTable> getPerSampleGenotypeConcordance() {
        return Collections.unmodifiableMap(perSampleGenotypeConcordance);
    }

    public Map<String,Double> getPerSampleNRD() {
        Map<String,Double> nrd = new HashMap<String,Double>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            nrd.put(sampleTable.getKey(),calculateNRD(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(nrd);
    }

    public Double getOverallNRD() {
        return calculateNRD(overallGenotypeConcordance);
    }

    public Map<String,Double> getPerSampleNRS() {
        Map<String,Double> nrs = new HashMap<String,Double>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            nrs.put(sampleTable.getKey(),calculateNRS(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(nrs);
    }

    public Double getOverallNRS() {
        return calculateNRS(overallGenotypeConcordance);
    }

    @Requires({"eval != null","truth != null"})
    public void update(VariantContext eval, VariantContext truth) {
        overallSiteConcordance.update(eval,truth);
        Set<String> alleleTruth = new HashSet<String>(8);
        alleleTruth.add(truth.getReference().getBaseString());
        for ( Allele a : truth.getAlternateAlleles() ) {
            alleleTruth.add(a.getBaseString());
        }
        for ( String sample : perSampleGenotypeConcordance.keySet() ) {
            Genotype evalGenotype = eval.getGenotype(sample);
            Genotype truthGenotype = truth.getGenotype(sample);
            perSampleGenotypeConcordance.get(sample).update(evalGenotype,truthGenotype,alleleTruth);
            overallGenotypeConcordance.update(evalGenotype,truthGenotype,alleleTruth);
        }
    }

    private static double calculateNRD(GenotypeConcordanceTable table) {
        return calculateNRD(table.getTable());
    }

    private static double calculateNRD(int[][] concordanceCounts) {
        int correct = 0;
        int total = 0;
        correct += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HET.ordinal()];
        correct += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += correct;
        total += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HET.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HOM_REF.ordinal()];
        total += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HOM_REF.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HET.ordinal()];
        // NRD is by definition incorrec/total = 1.0-correct/total
        return 1.0 - ( (double) correct)/( (double) total);
    }

    private static double calculateNRS(GenotypeConcordanceTable table) {
        return calculateNRS(table.getTable());
    }

    private static double calculateNRS(int[][] concordanceCounts) {
        long confirmedVariant = 0;
        long unconfirmedVariant = 0;
        for ( GenotypeType truthState : Arrays.asList(GenotypeType.HET,GenotypeType.HOM_VAR) ) {
            for ( GenotypeType evalState : GenotypeType.values() ) {
                if ( evalState == GenotypeType.MIXED )
                    continue;
                if ( evalState.equals(GenotypeType.HET) || evalState.equals(GenotypeType.HOM_VAR) )
                    confirmedVariant += concordanceCounts[evalState.ordinal()][truthState.ordinal()];
                else
                    unconfirmedVariant += concordanceCounts[evalState.ordinal()][truthState.ordinal()];
            }
        }

        return ( (double) confirmedVariant ) / ( (double) ( confirmedVariant + unconfirmedVariant ) );
    }


    class GenotypeConcordanceTable {

        private int[][] genotypeCounts;
        private int nMismatchingAlt;

        public GenotypeConcordanceTable() {
            genotypeCounts = new int[GenotypeType.values().length][GenotypeType.values().length];
            nMismatchingAlt = 0;
        }

        @Requires({"eval!=null","truth != null","truthAlleles != null"})
        public void update(Genotype eval, Genotype truth, Set<String> truthAlleles) {
            // this is slow but correct
            boolean matchingAlt = true;
            if ( eval.isCalled() && truth.isCalled() ) {
                // by default, no-calls "match" between alleles, so if
                // one or both sites are no-call or unavailable, the alt alleles match
                // otherwise, check explicitly: if the eval has an allele that's not ref, no-call, or present in truth
                // the alt allele is mismatching - regardless of whether the genotype is correct.
                for ( Allele evalAllele : eval.getAlleles() ) {
                    matchingAlt &= truthAlleles.contains(evalAllele.getBaseString());
                }
            }

            if ( matchingAlt ) {
                genotypeCounts[eval.getType().ordinal()][truth.getType().ordinal()]++;
            } else {
                nMismatchingAlt++;
            }
        }

        public int[][] getTable() {
            return genotypeCounts;
        }

        public int getnMismatchingAlt() {
            return nMismatchingAlt;
        }

        public int getnEvalGenotypes(GenotypeType type) {
            int nGeno = 0;
            for ( GenotypeType comptype : GenotypeType.values() )
                nGeno += genotypeCounts[type.ordinal()][comptype.ordinal()];
            return nGeno;
        }

        public int getnCalledEvalGenotypes() {
            int nGeno = 0;
            for ( GenotypeType evalType : Arrays.asList(GenotypeType.HOM_REF,GenotypeType.HOM_VAR,GenotypeType.HET) ) {
                nGeno += getnEvalGenotypes(evalType);
            }

            return nGeno + nMismatchingAlt;
        }

        public int getnCompGenotypes(GenotypeType type) {
            int nGeno = 0;
            for ( GenotypeType evaltype : GenotypeType.values() )
                nGeno += genotypeCounts[evaltype.ordinal()][type.ordinal()];
            return nGeno;
        }

        public int getnCalledCompGenotypes() {
            int nGeno = 0;
            for ( GenotypeType compType : Arrays.asList(GenotypeType.HOM_REF,GenotypeType.HOM_VAR,GenotypeType.HET) ) {
                nGeno += getnCompGenotypes(compType);
            }
            return nGeno;
        }

        public int get(GenotypeType evalType, GenotypeType compType) {
            return genotypeCounts[evalType.ordinal()][compType.ordinal()];
        }
    }

    class SiteConcordanceTable {

        private int[] siteConcordance;

        public SiteConcordanceTable() {
            siteConcordance = new int[SiteConcordanceType.values().length];
        }

        public void update(VariantContext evalVC, VariantContext truthVC) {
            SiteConcordanceType matchType = getMatchType(evalVC,truthVC);
            siteConcordance[matchType.ordinal()]++;
        }

        @Requires({"evalVC != null","truthVC != null"})
        private SiteConcordanceType getMatchType(VariantContext evalVC, VariantContext truthVC) {
            if ( evalVC.isMonomorphicInSamples() )
                return SiteConcordanceType.TRUTH_ONLY;
            if ( truthVC.isMonomorphicInSamples() )
                return SiteConcordanceType.EVAL_ONLY;

            boolean evalSusbsetTruth = VariantContextUtils.allelesAreSubset(evalVC,truthVC);
            boolean truthSubsetEval = VariantContextUtils.allelesAreSubset(truthVC,evalVC);

            if ( evalSusbsetTruth && truthSubsetEval )
                return SiteConcordanceType.ALLELES_MATCH;
            else if ( evalSusbsetTruth )
                return SiteConcordanceType.EVAL_SUBSET_TRUTH;
            else if ( truthSubsetEval )
                return SiteConcordanceType.EVAL_SUPERSET_TRUTH;

            return SiteConcordanceType.ALLELES_DO_NOT_MATCH;
        }

        public int[] getSiteConcordance() {
            return siteConcordance;
        }

        public int get(SiteConcordanceType type) {
            return getSiteConcordance()[type.ordinal()];
        }
    }

    enum SiteConcordanceType {
        ALLELES_MATCH,
        EVAL_SUPERSET_TRUTH,
        EVAL_SUBSET_TRUTH,
        ALLELES_DO_NOT_MATCH,
        EVAL_ONLY,
        TRUTH_ONLY
    }
}
