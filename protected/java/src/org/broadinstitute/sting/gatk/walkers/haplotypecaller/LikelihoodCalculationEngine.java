/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pairhmm.*;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.variant.variantcontext.Allele;

import java.util.*;

public class LikelihoodCalculationEngine {

    private static final double LOG_ONE_HALF = -Math.log10(2.0);
    private final byte constantGCP;
    private final boolean DEBUG;
    private final PairHMM pairHMM;

    public LikelihoodCalculationEngine( final byte constantGCP, final boolean debug, final PairHMM.HMM_IMPLEMENTATION hmmType ) {

        switch (hmmType) {
            case EXACT:
                pairHMM = new Log10PairHMM(true);
                break;
            case ORIGINAL:
                pairHMM = new Log10PairHMM(false);
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

    public Map<String, PerReadAlleleLikelihoodMap> computeReadLikelihoods( final List<Haplotype> haplotypes, final Map<String, List<GATKSAMRecord>> perSampleReadList ) {

        final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = new HashMap<String, PerReadAlleleLikelihoodMap>();
        int X_METRIC_LENGTH = 0;
        for( final Map.Entry<String, List<GATKSAMRecord>> sample : perSampleReadList.entrySet() ) {
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
        for( final Map.Entry<String, List<GATKSAMRecord>> sampleEntry : perSampleReadList.entrySet() ) {
            //if( DEBUG ) { System.out.println("Evaluating sample " + sample + " with " + perSampleReadList.get( sample ).size() + " passing reads"); }
            // evaluate the likelihood of the reads given those haplotypes
            stratifiedReadMap.put(sampleEntry.getKey(), computeReadLikelihoods(haplotypes, sampleEntry.getValue()));
        }
        return stratifiedReadMap;
    }

    private PerReadAlleleLikelihoodMap computeReadLikelihoods( final List<Haplotype> haplotypes, final List<GATKSAMRecord> reads) {
        // first, a little set up to get copies of the Haplotypes that are Alleles (more efficient than creating them each time)
        final int numHaplotypes = haplotypes.size();
        final Map<Haplotype, Allele> alleleVersions = new HashMap<Haplotype, Allele>(numHaplotypes);
        for ( final Haplotype haplotype : haplotypes ) {
            alleleVersions.put(haplotype, Allele.create(haplotype, true));
        }

        final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        for( final GATKSAMRecord read : reads ) {
            final byte[] overallGCP = new byte[read.getReadLength()];
            Arrays.fill( overallGCP, constantGCP ); // Is there a way to derive empirical estimates for this from the data?
            Haplotype previousHaplotypeSeen = null;
            final byte[] readQuals = read.getBaseQualities();
            final byte[] readInsQuals = read.getBaseInsertionQualities();
            final byte[] readDelQuals = read.getBaseDeletionQualities();
            for( int kkk = 0; kkk < readQuals.length; kkk++ ) {
                readQuals[kkk] = (byte) Math.min( 0xff & readQuals[kkk], read.getMappingQuality()); // cap base quality by mapping quality, as in UG
                //readQuals[kkk] = ( readQuals[kkk] > readInsQuals[kkk] ? readInsQuals[kkk] : readQuals[kkk] ); // cap base quality by base insertion quality, needs to be evaluated
                //readQuals[kkk] = ( readQuals[kkk] > readDelQuals[kkk] ? readDelQuals[kkk] : readQuals[kkk] ); // cap base quality by base deletion quality, needs to be evaluated
                // TODO -- why is Q18 hard-coded here???
                readQuals[kkk] = ( readQuals[kkk] < (byte) 18 ? QualityUtils.MIN_USABLE_Q_SCORE : readQuals[kkk] );
            }

            for( int jjj = 0; jjj < numHaplotypes; jjj++ ) {
                final Haplotype haplotype = haplotypes.get(jjj);

                final int haplotypeStart = ( previousHaplotypeSeen == null ? 0 : PairHMM.findFirstPositionWhereHaplotypesDiffer(haplotype.getBases(), previousHaplotypeSeen.getBases()) );
                previousHaplotypeSeen = haplotype;

                perReadAlleleLikelihoodMap.add(read, alleleVersions.get(haplotype),
                        pairHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotype.getBases(), read.getReadBases(),
                                readQuals, readInsQuals, readDelQuals, overallGCP, haplotypeStart, jjj == 0));
            }
        }
        return perReadAlleleLikelihoodMap;
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

    @Requires({"haplotypes.size() > 0"})
    @Ensures({"result.size() <= haplotypes.size()"})
    public List<Haplotype> selectBestHaplotypes( final List<Haplotype> haplotypes, final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap, final int maxNumHaplotypesInPopulation ) {

        final int numHaplotypes = haplotypes.size();
        final Set<String> sampleKeySet = stratifiedReadMap.keySet();
        final List<Integer> bestHaplotypesIndexList = new ArrayList<Integer>();
        bestHaplotypesIndexList.add( findReferenceIndex(haplotypes) ); // always start with the reference haplotype
        final List<Allele> haplotypesAsAlleles = new ArrayList<Allele>();
        for( final Haplotype h : haplotypes ) { haplotypesAsAlleles.add(Allele.create(h, true)); }

        final double[][] haplotypeLikelihoodMatrix = computeDiploidHaplotypeLikelihoods( sampleKeySet, stratifiedReadMap, haplotypesAsAlleles ); // all samples pooled together

        int hap1 = 0;
        int hap2 = 0;
        //double bestElement = Double.NEGATIVE_INFINITY;
        final int maxChosenHaplotypes = Math.min( maxNumHaplotypesInPopulation, sampleKeySet.size() * 2 + 1 );
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

        final List<Haplotype> bestHaplotypes = new ArrayList<Haplotype>();
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