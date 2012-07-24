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
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public class LikelihoodCalculationEngine {

    private static final double LOG_ONE_HALF = -Math.log10(2.0);
    private static final double BEST_LIKELIHOOD_THRESHOLD = 0.1;
    private final byte constantGCP;
    private final boolean DEBUG;
    private final PairHMM pairHMM;

    public LikelihoodCalculationEngine( final byte constantGCP, final boolean debug, final boolean noBanded ) {
        pairHMM = new PairHMM( noBanded );
        this.constantGCP = constantGCP;
        DEBUG = debug;
    }

    public void computeReadLikelihoods( final ArrayList<Haplotype> haplotypes, final HashMap<String, ArrayList<GATKSAMRecord>> perSampleReadList ) {
        final int numHaplotypes = haplotypes.size();

        int X_METRIC_LENGTH = 0;
        for( final String sample : perSampleReadList.keySet() ) {
            for( final GATKSAMRecord read : perSampleReadList.get(sample) ) {
                final int readLength = read.getReadLength();
                if( readLength > X_METRIC_LENGTH ) { X_METRIC_LENGTH = readLength; }
            }
        }
        int Y_METRIC_LENGTH = 0;
        for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
            final int haplotypeLength = haplotypes.get(jjj).getBases().length;
            if( haplotypeLength > Y_METRIC_LENGTH ) { Y_METRIC_LENGTH = haplotypeLength; }
        }

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        X_METRIC_LENGTH += 2;
        Y_METRIC_LENGTH += 2;

        // initial arrays to hold the probabilities of being in the match, insertion and deletion cases
        final double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        PairHMM.initializeArrays(matchMetricArray, XMetricArray, YMetricArray, X_METRIC_LENGTH);

        // for each sample's reads
        for( final String sample : perSampleReadList.keySet() ) {
            //if( DEBUG ) { System.out.println("Evaluating sample " + sample + " with " + perSampleReadList.get( sample ).size() + " passing reads"); }
            // evaluate the likelihood of the reads given those haplotypes
            computeReadLikelihoods( haplotypes, perSampleReadList.get(sample), sample, matchMetricArray, XMetricArray, YMetricArray );
        }
    }

    private void computeReadLikelihoods( final ArrayList<Haplotype> haplotypes, final ArrayList<GATKSAMRecord> reads, final String sample,
                                         final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        final int numHaplotypes = haplotypes.size();
        final int numReads = reads.size();
        final double[][] readLikelihoods = new double[numHaplotypes][numReads];
        for( int iii = 0; iii < numReads; iii++ ) {
            final GATKSAMRecord read = reads.get(iii);

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
                readQuals[kkk] = ( readQuals[kkk] < (byte) 17 ? QualityUtils.MIN_USABLE_Q_SCORE : readQuals[kkk] );
            }

            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                final Haplotype haplotype = haplotypes.get(jjj);
                final int haplotypeStart = ( previousHaplotypeSeen == null ? 0 : computeFirstDifferingPosition(haplotype.getBases(), previousHaplotypeSeen.getBases()) );
                previousHaplotypeSeen = haplotype;

                readLikelihoods[jjj][iii] = pairHMM.computeReadLikelihoodGivenHaplotype(haplotype.getBases(), read.getReadBases(),
                        readQuals, readInsQuals, readDelQuals, overallGCP,
                        haplotypeStart, matchMetricArray, XMetricArray, YMetricArray);
            }
        }
        for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
            haplotypes.get(jjj).addReadLikelihoods( sample, readLikelihoods[jjj] );
        }
    }

    private static int computeFirstDifferingPosition( final byte[] b1, final byte[] b2 ) {
        for( int iii = 0; iii < b1.length && iii < b2.length; iii++ ){
            if( b1[iii] != b2[iii] ) {
                return iii;
            }
        }
        return b1.length;
    }

    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == haplotypes.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final ArrayList<Haplotype> haplotypes, final String sample ) {
        // set up the default 1-to-1 haplotype mapping object, BUGBUG: target for future optimization?
        final ArrayList<ArrayList<Haplotype>> haplotypeMapping = new ArrayList<ArrayList<Haplotype>>();
        for( final Haplotype h : haplotypes ) {
            final ArrayList<Haplotype> list = new ArrayList<Haplotype>();
            list.add(h);
            haplotypeMapping.add(list);
        }
        return computeDiploidHaplotypeLikelihoods( sample, haplotypeMapping );
    }
    
    @Requires({"haplotypeMapping.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == haplotypeMapping.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final String sample, final ArrayList<ArrayList<Haplotype>> haplotypeMapping ) {

        final int numHaplotypes = haplotypeMapping.size();
        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            Arrays.fill(haplotypeLikelihoodMatrix[iii], Double.NEGATIVE_INFINITY);
        }

        // compute the diploid haplotype likelihoods
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ) {                
                for( final Haplotype iii_mapped : haplotypeMapping.get(iii) ) {
                    final double[] readLikelihoods_iii = iii_mapped.getReadLikelihoods(sample);
                    for( final Haplotype jjj_mapped : haplotypeMapping.get(jjj) ) {
                        final double[] readLikelihoods_jjj = jjj_mapped.getReadLikelihoods(sample);
                        double haplotypeLikelihood = 0.0;
                        for( int kkk = 0; kkk < readLikelihoods_iii.length; kkk++ ) {
                            // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                            // First term is approximated by Jacobian log with table lookup.
                            haplotypeLikelihood += MathUtils.approximateLog10SumLog10(readLikelihoods_iii[kkk], readLikelihoods_jjj[kkk]) + LOG_ONE_HALF;
                        }
                        haplotypeLikelihoodMatrix[iii][jjj] = Math.max(haplotypeLikelihoodMatrix[iii][jjj], haplotypeLikelihood); // MathUtils.approximateLog10SumLog10(haplotypeLikelihoodMatrix[iii][jjj], haplotypeLikelihood); // BUGBUG: max or sum?
                    }
                }       
            }
        }

        // normalize the diploid likelihoods matrix
        return normalizeDiploidLikelihoodMatrixFromLog10( haplotypeLikelihoodMatrix );        
    }

    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.length == result[0].length", "result.length == haplotypes.size()"})
    public static double[][] computeDiploidHaplotypeLikelihoods( final ArrayList<Haplotype> haplotypes, final Set<String> samples ) {
        // set up the default 1-to-1 haplotype mapping object, BUGBUG: target for future optimization?
        final ArrayList<ArrayList<Haplotype>> haplotypeMapping = new ArrayList<ArrayList<Haplotype>>();
        for( final Haplotype h : haplotypes ) {
            final ArrayList<Haplotype> list = new ArrayList<Haplotype>();
            list.add(h);
            haplotypeMapping.add(list);
        }

        final int numHaplotypes = haplotypeMapping.size();
        final double[][] haplotypeLikelihoodMatrix = new double[numHaplotypes][numHaplotypes];
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            Arrays.fill(haplotypeLikelihoodMatrix[iii], Double.NEGATIVE_INFINITY);
        }

        // compute the diploid haplotype likelihoods
        for( int iii = 0; iii < numHaplotypes; iii++ ) {
            for( int jjj = 0; jjj <= iii; jjj++ ) {
                for( final Haplotype iii_mapped : haplotypeMapping.get(iii) ) {
                    for( final Haplotype jjj_mapped : haplotypeMapping.get(jjj) ) {
                        double haplotypeLikelihood = 0.0;
                        for( final String sample : samples ) {
                            final double[] readLikelihoods_iii = iii_mapped.getReadLikelihoods(sample);
                            final double[] readLikelihoods_jjj = jjj_mapped.getReadLikelihoods(sample);
                            for( int kkk = 0; kkk < readLikelihoods_iii.length; kkk++ ) {
                                // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                                // First term is approximated by Jacobian log with table lookup.
                                haplotypeLikelihood += MathUtils.approximateLog10SumLog10(readLikelihoods_iii[kkk], readLikelihoods_jjj[kkk]) + LOG_ONE_HALF;
                            }
                        }
                        haplotypeLikelihoodMatrix[iii][jjj] = Math.max(haplotypeLikelihoodMatrix[iii][jjj], haplotypeLikelihood); // MathUtils.approximateLog10SumLog10(haplotypeLikelihoodMatrix[iii][jjj], haplotypeLikelihood); // BUGBUG: max or sum?
                    }
                }
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
    public ArrayList<Haplotype> selectBestHaplotypes( final ArrayList<Haplotype> haplotypes ) {

        final int numHaplotypes = haplotypes.size();
        final Set<String> sampleKeySet = haplotypes.get(0).getSampleKeySet(); // BUGBUG: assume all haplotypes saw the same samples
        final ArrayList<Integer> bestHaplotypesIndexList = new ArrayList<Integer>();
        bestHaplotypesIndexList.add(0); // always start with the reference haplotype
        final double[][] haplotypeLikelihoodMatrix = computeDiploidHaplotypeLikelihoods( haplotypes, sampleKeySet ); // all samples pooled together

        int hap1 = 0;
        int hap2 = 0;
        //double bestElement = Double.NEGATIVE_INFINITY;
        final int maxChosenHaplotypes = Math.min( 8, sampleKeySet.size() * 2 + 1 );
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

    public static Map<String, Map<Allele, List<GATKSAMRecord>>> partitionReadsBasedOnLikelihoods( final GenomeLocParser parser, final HashMap<String, ArrayList<GATKSAMRecord>> perSampleReadList, final HashMap<String, ArrayList<GATKSAMRecord>> perSampleFilteredReadList, final Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>> call) {
        final Map<String, Map<Allele, List<GATKSAMRecord>>> returnMap = new HashMap<String, Map<Allele, List<GATKSAMRecord>>>();
        final GenomeLoc callLoc = parser.createGenomeLoc(call.getFirst());
        for( final String sample : perSampleReadList.keySet() ) {
            final Map<Allele, List<GATKSAMRecord>> alleleReadMap = new HashMap<Allele, List<GATKSAMRecord>>();
            final ArrayList<GATKSAMRecord> readsForThisSample = perSampleReadList.get(sample);
            for( int iii = 0; iii < readsForThisSample.size(); iii++ ) {
                final GATKSAMRecord read = readsForThisSample.get(iii); // BUGBUG: assumes read order in this list and haplotype likelihood list are the same!
                // only count the read if it overlaps the event, otherwise it is not added to the output read list at all
                if( callLoc.overlapsP(parser.createGenomeLoc(read)) ) {
                    final double likelihoods[] = new double[call.getFirst().getAlleles().size()];
                    int count = 0;
                    for( final Allele a : call.getFirst().getAlleles() ) { // find the allele with the highest haplotype likelihood
                        double maxLikelihood = Double.NEGATIVE_INFINITY;
                        for( final Haplotype h : call.getSecond().get(a) ) { // use the max likelihood from all the haplotypes which mapped to this allele (achieved via the haplotype mapper object)
                            final double likelihood = h.getReadLikelihoods(sample)[iii];
                            if( likelihood > maxLikelihood ) {
                                maxLikelihood = likelihood;
                            }
                        }
                        likelihoods[count++] = maxLikelihood;
                    }
                    final int bestAllele = MathUtils.maxElementIndex(likelihoods);
                    final double bestLikelihood = likelihoods[bestAllele];
                    Allele allele = Allele.NO_CALL;
                    boolean isInformativeRead = false;
                    for( final double likelihood : likelihoods ) {
                        if( bestLikelihood - likelihood > BEST_LIKELIHOOD_THRESHOLD ) {
                            isInformativeRead = true;
                            break;
                        }
                    }
                    // uninformative reads get the no call Allele
                    if( isInformativeRead ) {
                        allele = call.getFirst().getAlleles().get(bestAllele);
                    }
                    List<GATKSAMRecord> readList = alleleReadMap.get(allele);
                    if( readList == null ) {
                        readList = new ArrayList<GATKSAMRecord>();
                        alleleReadMap.put(allele, readList);
                    }
                    readList.add(read);
                }
            }
            // add all filtered reads to the NO_CALL list because they weren't given any likelihoods
            List<GATKSAMRecord> readList = alleleReadMap.get(Allele.NO_CALL);
            if( readList == null ) {
                readList = new ArrayList<GATKSAMRecord>();
                alleleReadMap.put(Allele.NO_CALL, readList);
            }
            for( final GATKSAMRecord read : perSampleFilteredReadList.get(sample) ) {
                // only count the read if it overlaps the event, otherwise it is not added to the output read list at all
                if( callLoc.overlapsP(parser.createGenomeLoc(read)) ) {
                    readList.add(read);
                }
            }
            returnMap.put(sample, alleleReadMap);
        }
        return returnMap;
    }
}