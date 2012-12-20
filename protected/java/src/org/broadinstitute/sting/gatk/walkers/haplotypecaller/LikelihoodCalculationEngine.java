/*
 * Copyright (c) 2011 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.pairhmm.*;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;

public class LikelihoodCalculationEngine {

    private static final double LOG_ONE_HALF = -Math.log10(2.0);
    private final byte constantGCP;
    private final boolean DEBUG;
    private final PairHMM pairHMM;

    public LikelihoodCalculationEngine( final byte constantGCP, final boolean debug, final PairHMM.HMM_IMPLEMENTATION hmmType ) {

        switch (hmmType) {
            case EXACT:
                pairHMM = new ExactPairHMM();
                break;
            case ORIGINAL:
                pairHMM = new OriginalPairHMM();
                break;
            case CACHING:
                pairHMM = new CachingPairHMM();
                break;
            case LOGLESS_CACHING:
                pairHMM = new LoglessCachingPairHMM();
                break;
            default:
                throw new UserException.BadArgumentValue("pairHMM", "Specified pairHMM implementation is unrecognized or incompatible with the HaplotypeCaller. Acceptable options are ORIGINAL, EXACT, CACHING, and LOGLESS_CACHING.");
        }

        this.constantGCP = constantGCP;
        DEBUG = debug;
    }

    public Map<String, PerReadAlleleLikelihoodMap> computeReadLikelihoods( final ArrayList<Haplotype> haplotypes, final HashMap<String, ArrayList<GATKSAMRecord>> perSampleReadList ) {

        final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = new HashMap<String, PerReadAlleleLikelihoodMap>();
        int X_METRIC_LENGTH = 0;
        for( final Map.Entry<String, ArrayList<GATKSAMRecord>> sample : perSampleReadList.entrySet() ) {
            for( final GATKSAMRecord read : sample.getValue() ) {
                final int readLength = read.getReadLength();
                if( readLength > X_METRIC_LENGTH ) { X_METRIC_LENGTH = readLength; }
            }
        }
        int Y_METRIC_LENGTH = 0;
        for( final Haplotype h : haplotypes ) {
            final int haplotypeLength = h.getBases().length;
            if( haplotypeLength > Y_METRIC_LENGTH ) { Y_METRIC_LENGTH = haplotypeLength; }
        }

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        X_METRIC_LENGTH += 2;
        Y_METRIC_LENGTH += 2;

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        pairHMM.initialize(X_METRIC_LENGTH, Y_METRIC_LENGTH);

        // for each sample's reads
        for( final Map.Entry<String, ArrayList<GATKSAMRecord>> sampleEntry : perSampleReadList.entrySet() ) {
            //if( DEBUG ) { System.out.println("Evaluating sample " + sample + " with " + perSampleReadList.get( sample ).size() + " passing reads"); }
            // evaluate the likelihood of the reads given those haplotypes
            stratifiedReadMap.put(sampleEntry.getKey(), computeReadLikelihoods(haplotypes, sampleEntry.getValue()));
        }
        return stratifiedReadMap;
    }

    private PerReadAlleleLikelihoodMap computeReadLikelihoods( final ArrayList<Haplotype> haplotypes, final ArrayList<GATKSAMRecord> reads) {

        final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = PerReadAlleleLikelihoodMap.getBestAvailablePerReadAlleleLikelihoodMap();
        final int numHaplotypes = haplotypes.size();
        for( final GATKSAMRecord read : reads ) {
            final byte[] overallGCP = new byte[read.getReadLength()];
            Arrays.fill( overallGCP, constantGCP ); // Is there a way to derive empirical estimates for this from the data?
            Haplotype previousHaplotypeSeen = null;
            final byte[] readQuals = read.getBaseQualities();
            final byte[] readInsQuals = read.getBaseInsertionQualities();
            final byte[] readDelQuals = read.getBaseDeletionQualities();
            for( int kkk = 0; kkk < readQuals.length; kkk++ ) {
                readQuals[kkk] = ( readQuals[kkk] > (byte) read.getMappingQuality() ? (byte) read.getMappingQuality() : readQuals[kkk] ); // cap base quality by mapping quality
                //readQuals[kkk] = ( readQuals[kkk] > readInsQuals[kkk] ? readInsQuals[kkk] : readQuals[kkk] ); // cap base quality by base insertion quality, needs to be evaluated
                //readQuals[kkk] = ( readQuals[kkk] > readDelQuals[kkk] ? readDelQuals[kkk] : readQuals[kkk] ); // cap base quality by base deletion quality, needs to be evaluated
                readQuals[kkk] = ( readQuals[kkk] < (byte) 18 ? QualityUtils.MIN_USABLE_Q_SCORE : readQuals[kkk] );
            }

            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                final Haplotype haplotype = haplotypes.get(jjj);
                final int haplotypeStart = ( previousHaplotypeSeen == null ? 0 : computeFirstDifferingPosition(haplotype.getBases(), previousHaplotypeSeen.getBases()) );
                previousHaplotypeSeen = haplotype;

                perReadAlleleLikelihoodMap.add(read, Allele.create(haplotype.getBases()),
                        pairHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotype.getBases(), read.getReadBases(),
                                readQuals, readInsQuals, readDelQuals, overallGCP, haplotypeStart, jjj == 0));
            }
        }
        return perReadAlleleLikelihoodMap;
    }

    private static int computeFirstDifferingPosition( final byte[] b1, final byte[] b2 ) {
        for( int iii = 0; iii < b1.length && iii < b2.length; iii++ ) {
            if( b1[iii] != b2[iii] ) {
                return iii;
            }
        }
        return Math.min(b1.length, b2.length);
    }

    @Requires({"alleleOrdering.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == alleleOrdering.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final String sample,
                                                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap,
                                                                 final List<Allele> alleleOrdering ) {
        final TreeSet<String> sampleSet = new TreeSet<String>();
        sampleSet.add(sample);
        return computeDiploidHaplotypeLikelihoods(sampleSet, stratifiedReadMap, alleleOrdering);
    }

    @Requires({"alleleOrdering.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == alleleOrdering.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final Set<String> samples,
                                                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap,
                                                                 final List<Allele> alleleOrdering ) {

        final int numHaplotypes = alleleOrdering.size();
        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            Arrays.fill(haplotypeLikelihoodMatrix[iii], Double.NEGATIVE_INFINITY);
        }

        // compute the diploid haplotype likelihoods
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            final Allele iii_allele = alleleOrdering.get(iii);
            for( int jjj = 0; jjj <= iii; jjj++ ) {
                final Allele jjj_allele = alleleOrdering.get(jjj);
                double haplotypeLikelihood = 0.0;
                for( final String sample : samples ) {
                    for( final Map.Entry<GATKSAMRecord, Map<Allele,Double>> entry : stratifiedReadMap.get(sample).getLikelihoodReadMap().entrySet() ) {
                        // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                        // First term is approximated by Jacobian log with table lookup.
                        haplotypeLikelihood += ReadUtils.getMeanRepresentativeReadCount( entry.getKey() ) *
                                ( MathUtils.approximateLog10SumLog10(entry.getValue().get(iii_allele), entry.getValue().get(jjj_allele)) + LOG_ONE_HALF );
                    }
                }
                haplotypeLikelihoodMatrix[iii][jjj] = haplotypeLikelihood;
            }
        }

        // normalize the diploid likelihoods matrix
        return normalizeDiploidLikelihoodMatrixFromLog10( haplotypeLikelihoodMatrix );
    }

    @Requires({"likelihoodMatrix.length == likelihoodMatrix[0].length"})
    @Ensures({"result.length == result[0].length", "result.length == likelihoodMatrix.length"})
    protected static double[][] normalizeDiploidLikelihoodMatrixFromLog10( final double[][] likelihoodMatrix ) {
        final int numHaplotypes = likelihoodMatrix.length;
        double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int index = 0;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ){
                genotypeLikelihoods[index++] = likelihoodMatrix[iii][jjj];
            }
        }
        genotypeLikelihoods = MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
        index = 0;
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ){
                likelihoodMatrix[iii][jjj] = genotypeLikelihoods[index++];
            }
        }
        return likelihoodMatrix;
    }

    /*
    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.size() <= haplotypes.size()"})
    public ArrayList<Haplotype> selectBestHaplotypes( final ArrayList<Haplotype> haplotypes ) {

        // BUGBUG: This function needs a lot of work. Need to use 4-gamete test or Tajima's D to decide to break up events into separate pieces for genotyping

        final int numHaplotypes = haplotypes.size();
        final Set<String> sampleKeySet = haplotypes.get(0).getSampleKeySet(); // BUGBUG: assume all haplotypes saw the same samples
        final ArrayList<Integer> bestHaplotypesIndexList = new ArrayList<Integer>();
        bestHaplotypesIndexList.add(0); // always start with the reference haplotype
        final double[][][] haplotypeLikelihoodMatrix = new double[sampleKeySet.size()][numHaplotypes][numHaplotypes];

        int sampleCount = 0;
        for( final String sample : sampleKeySet ) {
            haplotypeLikelihoodMatrix[sampleCount++] = computeDiploidHaplotypeLikelihoods( haplotypes, sample );
        }

        int hap1 = 0;
        int hap2 = 0;
        int chosenSample = 0;
        //double bestElement = Double.NEGATIVE_INFINITY;
        final int maxChosenHaplotypes = Math.min( 15, sampleKeySet.size() * 2 + 1 );
        while( bestHaplotypesIndexList.size() < maxChosenHaplotypes ) {
            double maxElement = Double.NEGATIVE_INFINITY;
            for( int kkk = 0; kkk < sampleCount; kkk++ ) {
                for( int iii = 0; iii < numHaplotypes; iii++ ) {
                    for( int jjj = 0; jjj <= iii; jjj++ ) {
                        if( haplotypeLikelihoodMatrix[kkk][iii][jjj] > maxElement ) {
                            maxElement = haplotypeLikelihoodMatrix[kkk][iii][jjj];
                            hap1 = iii;
                            hap2 = jjj;
                            chosenSample = kkk;
                        }
                    }
                }
            }
            if( maxElement == Double.NEGATIVE_INFINITY ) { break; }

            if( !bestHaplotypesIndexList.contains(hap1) ) { bestHaplotypesIndexList.add(hap1); }
            if( !bestHaplotypesIndexList.contains(hap2) ) { bestHaplotypesIndexList.add(hap2); }

            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    haplotypeLikelihoodMatrix[chosenSample][iii][jjj] = Double.NEGATIVE_INFINITY;
                }
            }
        }

        if( DEBUG ) { System.out.println("Chose " + (bestHaplotypesIndexList.size() - 1) + " alternate haplotypes to genotype in all samples."); }

        final ArrayList<Haplotype> bestHaplotypes = new ArrayList<Haplotype>();
        for( final int hIndex : bestHaplotypesIndexList ) {
            bestHaplotypes.add( haplotypes.get(hIndex) );
        }
        return bestHaplotypes;
    }
    */

    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.size() <= haplotypes.size()"})
    public ArrayList<Haplotype> selectBestHaplotypes( final ArrayList<Haplotype> haplotypes, final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap ) {

        final int numHaplotypes = haplotypes.size();
        final Set<String> sampleKeySet = stratifiedReadMap.keySet();
        final ArrayList<Integer> bestHaplotypesIndexList = new ArrayList<Integer>();
        bestHaplotypesIndexList.add( findReferenceIndex(haplotypes) ); // always start with the reference haplotype
        final List<Allele> haplotypesAsAlleles = new ArrayList<Allele>();
        for( final Haplotype h : haplotypes ) { haplotypesAsAlleles.add(Allele.create(h.getBases())); }

        final double[][] haplotypeLikelihoodMatrix = computeDiploidHaplotypeLikelihoods( sampleKeySet, stratifiedReadMap, haplotypesAsAlleles ); // all samples pooled together

        int hap1 = 0;
        int hap2 = 0;
        //double bestElement = Double.NEGATIVE_INFINITY;
        final int maxChosenHaplotypes = Math.min( 13, sampleKeySet.size() * 2 + 1 );
        while( bestHaplotypesIndexList.size() < maxChosenHaplotypes ) {
            double maxElement = Double.NEGATIVE_INFINITY;
            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    if( haplotypeLikelihoodMatrix[iii][jjj] > maxElement ) {
                        maxElement = haplotypeLikelihoodMatrix[iii][jjj];
                        hap1 = iii;
                        hap2 = jjj;
                    }
                }
            }
            if( maxElement == Double.NEGATIVE_INFINITY ) { break; }
            if( DEBUG ) { System.out.println("Chose haplotypes " + hap1 + " and " + hap2 + " with diploid likelihood = " + haplotypeLikelihoodMatrix[hap1][hap2]); }
            haplotypeLikelihoodMatrix[hap1][hap2] = Double.NEGATIVE_INFINITY;

            if( !bestHaplotypesIndexList.contains(hap1) ) { bestHaplotypesIndexList.add(hap1); }
            if( !bestHaplotypesIndexList.contains(hap2) ) { bestHaplotypesIndexList.add(hap2); }
        }

        if( DEBUG ) { System.out.println("Chose " + (bestHaplotypesIndexList.size() - 1) + " alternate haplotypes to genotype in all samples."); }

        final ArrayList<Haplotype> bestHaplotypes = new ArrayList<Haplotype>();
        for( final int hIndex : bestHaplotypesIndexList ) {
            bestHaplotypes.add( haplotypes.get(hIndex) );
        }
        return bestHaplotypes;
    }

    public static int findReferenceIndex( final List<Haplotype> haplotypes ) {
        for( final Haplotype h : haplotypes ) {
            if( h.isReference() ) { return haplotypes.indexOf(h); }
        }
        throw new ReviewedStingException( "No reference haplotype found in the list of haplotypes!" );
    }
}