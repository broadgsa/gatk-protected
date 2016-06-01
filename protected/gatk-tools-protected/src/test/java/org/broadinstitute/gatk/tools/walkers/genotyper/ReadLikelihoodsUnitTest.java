/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.genotyper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.genotyper.*;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test code for {@link ReadLikelihoods}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadLikelihoodsUnitTest
{
    private static final double EPSILON = 1e-6;
    private static final int ODD_READ_START = 101;
    private static final int EVEN_READ_START = 1;

    @Test(dataProvider = "dataSets")
    public void testInstantiationAndQuery(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> result = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);

        Assert.assertEquals(result.sampleCount(), samples.length);
        Assert.assertEquals(result.alleleCount(), alleles.length);


        testSampleQueries(samples, reads, result);
        testAlleleQueries(alleles, result);
        testLikelihoodMatrixQueries(samples, result, null);
    }

    @Test(dataProvider = "dataSets")
    public void testLikelihoodFillingAndQuery(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> result = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] likelihoods = fillWithRandomLikelihoods(samples, alleles, result);
        testLikelihoodMatrixQueries(samples, result, likelihoods);
    }

    private double[][][] fillWithRandomLikelihoods(final String[] samples, final Allele[] alleles, final ReadLikelihoods<Allele> result) {
        final Random rnd = Utils.getRandomGenerator();
        final double[][][] likelihoods = new double[samples.length][alleles.length][];
        for (int s = 0; s < likelihoods.length; s++) {
            final ReadLikelihoods.Matrix<Allele> sampleLikelihoods = result.sampleMatrix(s);
            for (int a = 0; a < likelihoods[s].length; a++) {
                likelihoods[s][a] = new double[result.sampleReadCount(s)];
                for (int r = 0; r < likelihoods[s][a].length; r++)
                    sampleLikelihoods.set(a,r,likelihoods[s][a][r] = -Math.abs(rnd.nextGaussian()));
            }
        }
        return likelihoods;
    }

    @Test(dataProvider = "dataSets")
    public void testBestAlleles(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        fillWithRandomLikelihoods(samples,alleles,original);
        final int alleleCount = alleles.length;
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            final ReadLikelihoods.Matrix<Allele> sampleMatrix = original.sampleMatrix(s);
            final double[] bestLkArray = new double[sampleReadCount];
            final int[] bestIndexArray = new int[sampleReadCount];
            final double[] confidenceArray = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                int bestAlleleIndex = -1;
                double bestAlleleLk = Double.NEGATIVE_INFINITY;
                double secondBestAlleleLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < alleleCount; a++) {
                    final double lk = sampleMatrix.get(a,r);
                    if (lk > bestAlleleLk) {
                        secondBestAlleleLk = bestAlleleLk;
                        bestAlleleLk = lk;
                        bestAlleleIndex = a;
                    } else if (lk > secondBestAlleleLk) {
                        secondBestAlleleLk = lk;
                    }
                }
                bestLkArray[r] = bestAlleleLk;
                confidenceArray[r] = bestAlleleLk - secondBestAlleleLk;
                bestIndexArray[r] = bestAlleleIndex;
            }
            final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = original.bestAlleles();
            for (final ReadLikelihoods<Allele>.BestAllele bestAllele : bestAlleles) {
                final int readIndex = original.readIndex(s,bestAllele.read);
                if (readIndex == -1) continue;
                Assert.assertEquals(bestLkArray[readIndex],bestAllele.likelihood);
                Assert.assertEquals(bestAllele.allele,alleles[bestIndexArray[readIndex]]);
                Assert.assertEquals(bestAllele.confidence,confidenceArray[readIndex],EPSILON);
            }
        }
    }

    @Test(dataProvider = "dataSets")
    public void testBestAlleleMap(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        fillWithRandomLikelihoods(samples,alleles,original);
        final Map<Allele,List<GATKSAMRecord>> expected = new HashMap<>(alleles.length);
        for (final Allele allele : alleles)
            expected.put(allele,new ArrayList<GATKSAMRecord>());

        final int alleleCount = alleles.length;
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            final ReadLikelihoods.Matrix<Allele> sampleMatrix = original.sampleMatrix(s);
            for (int r = 0; r < sampleReadCount; r++) {
                int bestAlleleIndex = -1;
                double bestAlleleLk = Double.NEGATIVE_INFINITY;
                double secondBestAlleleLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < alleleCount; a++) {
                    final double lk = sampleMatrix.get(a,r);
                    if (lk > bestAlleleLk) {
                        secondBestAlleleLk = bestAlleleLk;
                        bestAlleleLk = lk;
                        bestAlleleIndex = a;
                    } else if (lk > secondBestAlleleLk) {
                        secondBestAlleleLk = lk;
                    }
                }
                if ((bestAlleleLk - secondBestAlleleLk) > ReadLikelihoods.BestAllele.INFORMATIVE_THRESHOLD)
                    expected.get(alleles[bestAlleleIndex]).add(sampleMatrix.readAt(r));
            }
        }

        final Map<Allele,List<GATKSAMRecord>> actual = original.readsByBestAlleleMap();

        Assert.assertEquals(actual.size(),alleles.length);
        for (final Allele allele : alleles) {
            final List<GATKSAMRecord> expectedList = expected.get(allele);
            final List<GATKSAMRecord> actualList = actual.get(allele);
            final Set<GATKSAMRecord> expectedSet = new HashSet<>(expectedList);
            final Set<GATKSAMRecord> actualSet = new HashSet<>(actualList);
            Assert.assertEquals(actualSet,expectedSet);
        }
    }

    @Test(dataProvider = "dataSets")
    public void testFilterPoorlyModeledReads(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);

        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int r = 0; r < sampleReadCount; r++) {
                if ((r & 1) == 0) continue;
                for (int a = 0; a < alleles.length; a++)
                    original.sampleMatrix(s).set(a,r,-10000);
            }
        }

        final ReadLikelihoods<Allele> result = original.clone();
        result.filterPoorlyModeledReads(2.0);

        for (int s = 0; s < samples.length; s++) {
            final int oldSampleReadCount = original.sampleReadCount(s);
            final int newSampleReadCount = result.sampleReadCount(s);
            Assert.assertEquals(newSampleReadCount,(oldSampleReadCount + 1) / 2);
            final ReadLikelihoods.Matrix<Allele> newSampleMatrix = result.sampleMatrix(s);
            final ReadLikelihoods.Matrix<Allele> oldSampleMatrix = original.sampleMatrix(s);
            for (int r = 0 ; r < newSampleReadCount; r++) {
                Assert.assertEquals(original.readIndex(s, result.sampleReads(s).get(r)), r * 2);
                for (int a = 0; a < alleles.length; a++) {
                    Assert.assertEquals(newSampleMatrix.get(a,r),oldSampleMatrix.get(a,r*2));
                }
            }
        }
    }

    @Test(dataProvider = "dataSets")
    public void testFilterReadsToOverlap(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final GenomeLoc evenReadOverlap = locParser.createGenomeLoc(SAM_HEADER.getSequenceDictionary().getSequences().get(0).getSequenceName(),EVEN_READ_START ,EVEN_READ_START );
        fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result = original.clone();
        result.filterToOnlyOverlappingUnclippedReads(evenReadOverlap);
        final double[][][] newLikelihoods = new double[samples.length][alleles.length][];
        for (int s = 0; s < samples.length ; s++)
            for (int a = 0; a < alleles.length; a++) {
                newLikelihoods[s][a] = new double[(original.sampleReadCount(s) + 1) / 2];
                final ReadLikelihoods.Matrix<Allele> sampleMatrix = original.sampleMatrix(s);
                for (int r = 0; r < newLikelihoods[s][a].length; r++) {
                    Assert.assertEquals(result.readIndex(s,sampleMatrix.readAt(r << 1)),r);
                    newLikelihoods[s][a][r] = sampleMatrix.get(a, r << 1);
                }
            }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    @Test(dataProvider = "marginalizationDataSets")
    public void testMarginalizationWithOverlap(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads, final Map<Allele,List<Allele>> newToOldAlleleMapping) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final GenomeLoc evenReadOverlap = locParser.createGenomeLoc(SAM_HEADER.getSequenceDictionary().getSequences().get(0).getSequenceName(),EVEN_READ_START ,EVEN_READ_START );
        fillWithRandomLikelihoods(samples, alleles, original);
        final ReadLikelihoods<Allele> marginalized = original.marginalize(newToOldAlleleMapping,evenReadOverlap);
        Assert.assertNotNull(marginalized);
        Assert.assertEquals(newToOldAlleleMapping.size(),marginalized.alleleCount());
        for (int a = 0; a < marginalized.alleleCount(); a++) {
            final List<Allele> oldAlleles = newToOldAlleleMapping.get(marginalized.alleleAt(a));
            Assert.assertNotNull(oldAlleles);
            for (int s = 0; s < samples.length; s++) {
                final ReadLikelihoods.Matrix<Allele> oldSmapleLikelihoods = original.sampleMatrix(s);
                final ReadLikelihoods.Matrix<Allele> sampleLikelihoods = marginalized.sampleMatrix(s);
                final int sampleReadCount = sampleLikelihoods.readCount();
                final int oldSampleReadCount = oldSmapleLikelihoods.readCount();
                Assert.assertEquals(sampleReadCount,(oldSampleReadCount + 1) / 2);
                for (int r = 0; r < sampleReadCount; r++) {
                    double oldBestLk = Double.NEGATIVE_INFINITY;
                    for (final Allele oldAllele : oldAlleles) {
                        oldBestLk = Math.max(oldSmapleLikelihoods.get(original.alleleIndex(oldAllele),r << 1), oldBestLk);
                    }
                    Assert.assertEquals(sampleLikelihoods.get(a,r),oldBestLk);
                }
            }
        }
    }

    @Test(dataProvider = "marginalizationDataSets")
    public void testMarginalization(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads, final Map<Allele,List<Allele>> newToOldAlleleMapping) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        fillWithRandomLikelihoods(samples, alleles, original);
        final ReadLikelihoods<Allele> marginalized = original.marginalize(newToOldAlleleMapping);
        Assert.assertNotNull(marginalized);
        Assert.assertEquals(newToOldAlleleMapping.size(),marginalized.alleleCount());
        for (int a = 0; a < marginalized.alleleCount(); a++) {
            final List<Allele> oldAlleles = newToOldAlleleMapping.get(marginalized.alleleAt(a));
            Assert.assertNotNull(oldAlleles);
            for (int s = 0; s < samples.length; s++) {
                final ReadLikelihoods.Matrix<Allele> oldSmapleLikelihoods = original.sampleMatrix(s);
                final ReadLikelihoods.Matrix<Allele> sampleLikelihoods = marginalized.sampleMatrix(s);
                final int sampleReadCount = sampleLikelihoods.readCount();
                final int oldSampleReadCount = oldSmapleLikelihoods.readCount();
                Assert.assertEquals(oldSampleReadCount,sampleReadCount);
                for (int r = 0; r < sampleReadCount; r++) {
                    double oldBestLk = Double.NEGATIVE_INFINITY;
                    for (final Allele oldAllele : oldAlleles) {
                        oldBestLk = Math.max(oldSmapleLikelihoods.get(original.alleleIndex(oldAllele),r), oldBestLk);
                    }
                    Assert.assertEquals(sampleLikelihoods.get(a,r),oldBestLk);
                }
            }
        }
    }

    @Test(dataProvider = "dataSets")
    public void testNormalizeBestToZero(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result= original.clone();
        result.normalizeLikelihoods(true, Double.NEGATIVE_INFINITY);
        testAlleleQueries(alleles,result);
        final int alleleCount = alleles.length;
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][alleles.length][];
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int a = 0; a < alleleCount; a++)
                newLikelihoods[s][a] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestLk = originalLikelihoods[s][0][r];
                for (int a = 1; a < alleleCount; a++) {
                    bestLk = Math.max(bestLk,originalLikelihoods[s][a][r]);
                }
                for (int a = 0; a < alleleCount; a++) {
                    newLikelihoods[s][a][r] = originalLikelihoods[s][a][r] - bestLk;
                }
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    @Test(dataProvider = "dataSets")
    public void testNormalizeCapWorstLK(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result= original.clone();
        result.normalizeLikelihoods(false, - 0.001);
        testAlleleQueries(alleles,result);
        final int alleleCount = alleles.length;
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][alleles.length][];
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int a = 0; a < alleleCount; a++)
                newLikelihoods[s][a] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestAltLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < alleleCount; a++) {
                    if (alleles[a].isReference())
                        continue;
                    bestAltLk = Math.max(bestAltLk,originalLikelihoods[s][a][r]);
                }
                if (bestAltLk == Double.NEGATIVE_INFINITY)
                    for (int a = 0; a < alleleCount; a++) {
                        newLikelihoods[s][a][r] = originalLikelihoods[s][a][r];
                    }
                else
                    for (int a = 0; a < alleleCount; a++) {
                        newLikelihoods[s][a][r] = Math.max(originalLikelihoods[s][a][r],bestAltLk - 0.001);
                    }
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    @Test(dataProvider = "dataSets")
    public void testNormalizeCapWorstLKAndBestToZero(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result= original.clone();
        result.normalizeLikelihoods(true, - 0.001);
        testAlleleQueries(alleles,result);
        final int alleleCount = alleles.length;
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][alleles.length][];
        for (int s = 0; s < samples.length; s++) {
            final int sampleReadCount = original.sampleReadCount(s);
            for (int a = 0; a < alleleCount; a++)
                newLikelihoods[s][a] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestAltLk = Double.NEGATIVE_INFINITY;
                double bestLk = Double.NEGATIVE_INFINITY;
                for (int a = 0; a < alleleCount; a++) {
                    bestLk = Math.max(bestLk,originalLikelihoods[s][a][r]);
                    if (alleles[a].isReference())
                        continue;
                    bestAltLk = Math.max(bestAltLk,originalLikelihoods[s][a][r]);
                }
                if (bestAltLk == Double.NEGATIVE_INFINITY)
                    for (int a = 0; a < alleleCount; a++) {
                        newLikelihoods[s][a][r] = originalLikelihoods[s][a][r] - bestLk;
                    }
                else
                    for (int a = 0; a < alleleCount; a++) {
                        newLikelihoods[s][a][r] = Math.max(originalLikelihoods[s][a][r],bestAltLk - 0.001) - bestLk;
                    }
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }


    @Test(dataProvider = "dataSets")
    public void testAddMissingAlleles(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result = original.clone();

        // If all the alleles pass are present in the read-likelihoods collection there is no change.
        result.addMissingAlleles(result.alleles(),Double.NEGATIVE_INFINITY);
        testLikelihoodMatrixQueries(samples,result,originalLikelihoods);

        // If the allele list passed is empty there is no effect.
        result.addMissingAlleles(Collections.EMPTY_LIST,Double.NEGATIVE_INFINITY);
        testLikelihoodMatrixQueries(samples,result,originalLikelihoods);

        final Allele newOne;
        final Allele newTwo;
        final Allele newThree;

        // We add a single missing.
        result.addMissingAlleles(Arrays.asList(newOne = Allele.create("ACCCCCAAAATTTAAAGGG".getBytes(),false)),-12345.6);
        Assert.assertEquals(result.alleleCount(), original.alleleCount() + 1);

        // We add too more amongst exisisting alleles:
        result.addMissingAlleles(Arrays.asList(newTwo = Allele.create("ATATATTATATTAATATT".getBytes(), false),result.alleleAt(1),
                result.alleleAt(0),newThree = Allele.create("TGTGTGTATTG".getBytes(),false),Allele.create("ACCCCCAAAATTTAAAGGG".getBytes(),false)),-6.54321);

        Assert.assertEquals(original.alleleCount()+3,result.alleleCount());

        final List<Allele> expectedAlleles = new ArrayList<>(original.alleles());
        expectedAlleles.add(newOne); expectedAlleles.add(newTwo); expectedAlleles.add(newThree);

        Assert.assertEquals(result.alleles(),expectedAlleles);

        final double[][][] newLikelihoods = new double[originalLikelihoods.length][][];
        for (int s = 0; s < samples.length; s++) {
            newLikelihoods[s] = Arrays.copyOf(originalLikelihoods[s],originalLikelihoods[s].length + 3);
            final int sampleReadCount = original.sampleReadCount(s);
            final int originalAlleleCount = originalLikelihoods[s].length;
            newLikelihoods[s][originalAlleleCount] = new double[sampleReadCount];
            Arrays.fill(newLikelihoods[s][originalAlleleCount],-12345.6);
            newLikelihoods[s][originalAlleleCount+1] = new double[sampleReadCount];
            Arrays.fill(newLikelihoods[s][originalAlleleCount+1],-6.54321);
            newLikelihoods[s][originalAlleleCount+2] = new double[sampleReadCount];
            Arrays.fill(newLikelihoods[s][originalAlleleCount+2],-6.54321);
        }
            testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }


    @Test(dataProvider = "dataSets")
    public void testAddNonRefAllele(final String[] samples, final Allele[] alleles, final Map<String,List<GATKSAMRecord>> reads) {
        final ReadLikelihoods<Allele> original = new ReadLikelihoods<>(new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), reads);
        final double[][][] originalLikelihoods = fillWithRandomLikelihoods(samples,alleles,original);
        final ReadLikelihoods<Allele> result = original.clone();
        result.addNonReferenceAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        Assert.assertEquals(result.alleleCount(),original.alleleCount() + 1);
        Assert.assertEquals(result.alleleIndex(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE),result.alleleCount() - 1);
        final double[][][] newLikelihoods = new double[originalLikelihoods.length][][];
        for (int s = 0; s < samples.length; s++) {
            newLikelihoods[s] = Arrays.copyOf(originalLikelihoods[s],originalLikelihoods[s].length + 1);
            final int sampleReadCount = original.sampleReadCount(s);
            final int ordinaryAlleleCount = originalLikelihoods[s].length;
            newLikelihoods[s][ordinaryAlleleCount] = new double[sampleReadCount];
            for (int r = 0; r < sampleReadCount; r++) {
                double bestLk = newLikelihoods[s][0][r];
                double secondBestLk = Double.NEGATIVE_INFINITY;
                for (int a = 1; a < ordinaryAlleleCount; a++) {
                    final double lk = originalLikelihoods[s][a][r];
                    if (lk > bestLk) {
                        secondBestLk = bestLk;
                        bestLk = lk;
                    } else if (lk > secondBestLk) {
                        secondBestLk = lk;
                    }
                }
                final double expectedNonRefLk = Double.isInfinite(secondBestLk) ? bestLk : secondBestLk;
                newLikelihoods[s][ordinaryAlleleCount][r] = expectedNonRefLk;
            }
        }
        testLikelihoodMatrixQueries(samples,result,newLikelihoods);
    }

    private void testLikelihoodMatrixQueries(String[] samples, ReadLikelihoods<Allele> result, final double[][][] likelihoods) {
        for (final String sample : samples) {
            final int sampleIndex = result.sampleIndex(sample);
            final int sampleReadCount = result.sampleReadCount(sampleIndex);
            final int alleleCount = result.alleleCount();
            Assert.assertEquals(result.alleleCount(), alleleCount);
            for (int a = 0; a < alleleCount; a++) {
                Assert.assertEquals(result.sampleReadCount(sampleIndex),sampleReadCount);
                for (int r = 0; r < sampleReadCount; r++)
                    Assert.assertEquals(result.sampleMatrix(sampleIndex).get(a,r),
                            likelihoods == null ? 0.0 : likelihoods[sampleIndex][a][r], EPSILON);
            }
        }
    }

    private void testAlleleQueries(Allele[] alleles, ReadLikelihoods<Allele> result) {
        final Set<Integer> alleleIndices = new HashSet<>();
        for (final Allele allele : alleles) {
            final int alleleIndex = result.alleleIndex(allele);
            Assert.assertTrue(alleleIndex >= 0);
            Assert.assertFalse(alleleIndices.contains(alleleIndex));
            alleleIndices.add(alleleIndex);
            Assert.assertSame(allele,alleles[alleleIndex]);
        }
    }

    private void testSampleQueries(String[] samples, Map<String, List<GATKSAMRecord>> reads, ReadLikelihoods<Allele> result) {
        final Set<Integer> sampleIds = new HashSet<>(samples.length);
        for (final String sample : samples) {
            final int sampleIndex = result.sampleIndex(sample);
            Assert.assertTrue(sampleIndex >= 0);
            Assert.assertFalse(sampleIds.contains(sampleIndex));
            sampleIds.add(sampleIndex);

            final List<GATKSAMRecord> sampleReads = result.sampleReads(sampleIndex);
            final Set<GATKSAMRecord> sampleReadsSet = new HashSet<>(sampleReads);
            final List<GATKSAMRecord> expectedSampleReadArray = reads.get(sample);
            final Set<GATKSAMRecord> expectedSampleReadsSet = new HashSet<>(expectedSampleReadArray);
            Assert.assertEquals(sampleReadsSet,expectedSampleReadsSet);

            final int sampleReadCount = sampleReads.size();
            for (int r = 0; r < sampleReadCount; r++) {
                Assert.assertSame(sampleReads.get(r), expectedSampleReadArray.get(r));
                final int readIndex = result.readIndex(sampleIndex, sampleReads.get(r));
                Assert.assertEquals(readIndex,r);
            }
        }
    }

    private String[][] SAMPLE_SETS = new String[][] {
            {"A","B","C"},
            {"A"},
            {"C","A","D","E","Salsa","Gazpacho"},
    };

    private Allele[][] ALLELE_SETS = new Allele[][] {
            {Allele.create("A",true), Allele.create("T"), Allele.create("C")},
            {Allele.create("A",true)},
            {Allele.create("ATTTA"), Allele.create("A",true)},
            {Allele.create("A"), Allele.create("AT",true)},
                    {Allele.create("A",false), Allele.create("AT",false)},
    };

    @DataProvider(name="marginalizationDataSets")
    public Object[][] marginalizationDataSets() {
        try {
            final Random rnd = Utils.getRandomGenerator();
            final Object[][] result = new Object[SAMPLE_SETS.length * ALLELE_SETS.length * ALLELE_SETS.length][];
            int nextIndex = 0;
            for (int s = 0; s < SAMPLE_SETS.length; s++) {
                for (int a = 0; a < ALLELE_SETS.length; a++) {
                    for (int b = 0; b < ALLELE_SETS.length; b++) {
                        if (ALLELE_SETS[b].length < ALLELE_SETS[a].length)
                            result[nextIndex++] = new Object[]{SAMPLE_SETS[s], ALLELE_SETS[a],
                                    dataSetReads(SAMPLE_SETS[s], rnd), randomAlleleMap(ALLELE_SETS[a], ALLELE_SETS[b])
                            };
                    }
                }
            }
            return Arrays.copyOf(result,nextIndex);
        }catch (final Throwable e) {
            throw new RuntimeException(e);
        }
    }

    private Map<Allele,List<Allele>> randomAlleleMap(final Allele[] fromAlleles, final Allele[] toAlleles) {
        final Map<Allele,List<Allele>> result = new HashMap<>(toAlleles.length);
        for (final Allele toAllele : toAlleles )
            result.put(toAllele,new ArrayList<Allele>(fromAlleles.length));
        final ArrayList<Allele> remaining = new ArrayList<>(Arrays.asList(fromAlleles));
        int nextToIndex = 0;
        final Random rnd = Utils.getRandomGenerator();
        for (int i = 0; i < fromAlleles.length; i++) {
            final int fromAlleleIndex = rnd.nextInt(remaining.size());
            result.get(toAlleles[nextToIndex]).add(remaining.remove(fromAlleleIndex));
            nextToIndex = (nextToIndex + 1) % toAlleles.length;
        }
        return result;
    }


    @DataProvider(name="dataSets")
    public Object[][] dataSets() {
        try {
            final Random rnd = Utils.getRandomGenerator();
            final Object[][] result = new Object[SAMPLE_SETS.length * ALLELE_SETS.length][];
            int nextIndex = 0;
            for (int s = 0; s < SAMPLE_SETS.length; s++)
                for (int a = 0; a < ALLELE_SETS.length; a++) {
                    result[nextIndex++] = new Object[]{SAMPLE_SETS[s], ALLELE_SETS[a],
                            dataSetReads(SAMPLE_SETS[s], rnd)
                    };
                }
            return result;
        }catch (final Throwable e) {
            throw new RuntimeException(e);
        }
    }

    private Map<String,List<GATKSAMRecord>> dataSetReads(final String[] samples,
                                                         final Random rnd) {
        final Map<String,List<GATKSAMRecord>> result = new HashMap<>(samples.length);
        for (final String sample : samples) {
            final int readCount = rnd.nextInt(100);
            final List<GATKSAMRecord> reads = new ArrayList<>(readCount);
            for (int r = 0; r < readCount; r++) {
                final int alignmentStart = (r & 1) == 0 ? EVEN_READ_START : ODD_READ_START;
                reads.add(ArtificialSAMUtils.createArtificialRead(SAM_HEADER,
                        "RRR" + sample + "00" + r, 0, alignmentStart ,"AAAAA".getBytes(), new byte[] {30,30,30,30,30}, "5M"));
            }
            result.put(sample,reads);
        }
        return result;
    }

    @Test(dataProvider="readCountsAndAlleleCountDataSkippingNoAlleleAndWithReference")
    public void testInstantiationAndBasicQueries(final int[] readCounts, final int alleleCount, final boolean hasReference) {
        final SampleList sampleList = sampleList(readCounts);

        final AlleleList<Allele> alleleList = alleleList(alleleCount,hasReference);
        final Map<String,List<GATKSAMRecord>> sampleToReads = ReadLikelihoodsUnitTester.sampleToReads(sampleList, readCounts);
        final ReadLikelihoods<Allele> subject = new ReadLikelihoods<>(sampleList,alleleList,sampleToReads);

        AlleleListUnitTester.assertAlleleList(subject, AlleleListUtils.asList(alleleList));
        SampleListUnitTester.assertSampleList(subject,SampleListUtils.asList(sampleList));

        if (hasReference) {
            final int referenceIndex = AlleleListUtils.indexOfReference(alleleList);
            Assert.assertTrue(referenceIndex >= 0);
            Assert.assertEquals(AlleleListUtils.indexOfReference(alleleList),referenceIndex);
        } else {
            Assert.assertEquals(AlleleListUtils.indexOfReference(subject), -1);
        }

        testLikelihoodMatrixQueries(alleleList, sampleList, sampleToReads, subject);
        testAlleleQueries(alleleList, subject);
        testSampleQueries(sampleList, sampleToReads, subject);
    }

    @Test(dataProvider="readCountsAndAlleleCountDataSkippingNoLikelihoodsOrNoAlleleAndWithReference")
    public void testLikelihoodWriting(final int[] readCounts, final int alleleCount, final boolean hasReference) {
        final SampleList sampleList = sampleList(readCounts);

        final AlleleList<Allele> alleleList = alleleList(alleleCount,hasReference);
        final Map<String,List<GATKSAMRecord>> sampleToReads = ReadLikelihoodsUnitTester.sampleToReads(sampleList,readCounts);
        final ReadLikelihoods<Allele> subject = new ReadLikelihoods<>(sampleList,alleleList,sampleToReads);

        final int sampleCount = readCounts.length;
        int totalLikelihoodsSet = 0;
        int expectedLikelihoodsSet = 0;
        for (int s = 0; s < sampleCount; s++) {
            expectedLikelihoodsSet += readCounts[s] * alleleCount;
            final ReadLikelihoods.Matrix<Allele> matrix = subject.sampleMatrix(s);
            final int readCount = matrix.readCount();
            for (int a = 0; a < alleleCount; a++)
                for (int r = 0; r < readCount; r++)  {
                    final double likelihood = testLikelihood(s, a, r);
                    Assert.assertNotEquals(likelihood,0); //Paranoia
                    totalLikelihoodsSet++;
                    matrix.set(a,r,likelihood);
                    Assert.assertEquals(matrix.get(a, r),likelihood);
                }

        }
        Assert.assertEquals(totalLikelihoodsSet,expectedLikelihoodsSet);
    }

    @Test(dependsOnMethods={"testLikelihoodWriting","testInstantiationAndBasicQueries"},
            dataProvider="readCountsAndAlleleCountDataSkippingNoAlleleAndWithReference")
    public void testMapConversion(final int[] readCounts, final int alleleCount, final boolean hasReference) {
        final SampleList sampleList = sampleList(readCounts);

        final AlleleList<Allele> alleleList = alleleList(alleleCount,hasReference);
        final Map<String,List<GATKSAMRecord>> sampleToReads = ReadLikelihoodsUnitTester.sampleToReads(sampleList,readCounts);

        final Set<Allele> alleleWithLikelihoodsSet = new HashSet<>();
        final Set<GATKSAMRecord> readsWithLikelihoodsSet = new HashSet<>();
        final Map<String,PerReadAlleleLikelihoodMap> map = new HashMap<>(sampleList.sampleCount());
        final int sampleCount = sampleList.sampleCount();
        for (int s = 0; s < sampleCount; s++) {
            final String sample = sampleList.sampleAt(s);
            final PerReadAlleleLikelihoodMap perSampleMap = new PerReadAlleleLikelihoodMap();
            final List<GATKSAMRecord> reads = sampleToReads.get(sample);
            for (int a = 0; a < alleleCount; a++)
                for (int r = 0; r < reads.size(); r++) {
                    perSampleMap.add(reads.get(r), alleleList.alleleAt(a), testLikelihood(s, a, r));
                    alleleWithLikelihoodsSet.add(alleleList.alleleAt(a));
                    readsWithLikelihoodsSet.add(reads.get(r));
                }
            map.put(sample,perSampleMap);

        }

        ReadLikelihoods<Allele> subject = ReadLikelihoods.fromPerAlleleReadLikelihoodsMap(map);

        for (int s = 0; s < sampleCount; s++) {
            final String sample = sampleList.sampleAt(s);
            final int sIndex = subject.sampleIndex(sample);
            Assert.assertTrue(sIndex >= 0);
            Assert.assertTrue(sIndex < sampleCount);
            final int sampleReadCount = sampleToReads.get(sample).size();
            final ReadLikelihoods.Matrix<Allele> sampleLikelihoods = subject.sampleMatrix(sIndex);
            for (int a = 0; a < alleleCount; a++) {
                final Allele allele = alleleList.alleleAt(a);
                final int aIndex = subject.alleleIndex(allele);
                Assert.assertEquals(aIndex >= 0,alleleWithLikelihoodsSet.contains(allele));
                Assert.assertTrue(aIndex < alleleCount);
                if (aIndex == -1) continue;
                for (int r = 0; r < sampleReadCount; r++) {
                    final GATKSAMRecord read = sampleToReads.get(sample).get(r);
                    final int rIndex = subject.readIndex(sIndex,read);
                    final int rIndex2 = sampleLikelihoods.readIndex(read);
                    Assert.assertEquals(rIndex,rIndex2);
                    Assert.assertEquals(rIndex >= 0,readsWithLikelihoodsSet.contains(read));
                    Assert.assertTrue(rIndex < sampleReadCount);
                    if (rIndex == -1)
                        continue;
                    final double likelihood = sampleLikelihoods.get(aIndex,rIndex);
                    Assert.assertEquals(likelihood,testLikelihood(s,a,r));
                }
            }
        }
    }

    private double testLikelihood(final int sampleIndex, final int alleleIndex, final int readIndex) {
        return - Math.abs(31 * (sampleIndex + 1) + 101 * alleleIndex + 1009 * readIndex);
    }


    private final Random rnd = Utils.getRandomGenerator();

    private void testLikelihoodMatrixQueries(final AlleleList<Allele> alleles, final SampleList samples,
                                             final Map<String,List<GATKSAMRecord>> sampleToReads, ReadLikelihoods<Allele> result) {
        for (final String sample : SampleListUtils.asList(samples)) {
            final int sampleIndex = result.sampleIndex(sample);
            final ReadLikelihoods.Matrix<Allele> likelihoodMatrix = result.sampleMatrix(sampleIndex);
            final int sampleReadCount = sampleToReads.get(sample).size();
            final List<GATKSAMRecord> reads = sampleToReads.get(sample);
            Assert.assertEquals(likelihoodMatrix.alleleCount(), alleles.alleleCount());
            Assert.assertEquals(likelihoodMatrix.readCount(), sampleReadCount);
            for (int a = 0; a < likelihoodMatrix.alleleCount(); a++) {
                Assert.assertEquals(likelihoodMatrix.alleleAt(a),alleles.alleleAt(a));
                for (int r = 0; r < sampleReadCount; r++) {
                    Assert.assertEquals(likelihoodMatrix.readAt(r),reads.get(r));
                    Assert.assertEquals(likelihoodMatrix.get(a, r), 0.0);
                }
            }
        }
    }

    private void testAlleleQueries(final AlleleList<Allele> alleles, ReadLikelihoods<Allele> result) {
        final Set<Integer> alleleIndices = new HashSet<>();
        for (final Allele allele : AlleleListUtils.asList(alleles)) {
            final int alleleIndex = result.alleleIndex(allele);
            Assert.assertTrue(alleleIndex >= 0);
            Assert.assertFalse(alleleIndices.contains(alleleIndex));
            alleleIndices.add(alleleIndex);
            Assert.assertSame(allele,alleles.alleleAt(alleleIndex));
        }
    }

    private void testSampleQueries(final SampleList samples, Map<String, List<GATKSAMRecord>> reads,
                                   final ReadLikelihoods<Allele> result) {
        final Set<Integer> sampleIds = new HashSet<>(samples.sampleCount());
        for (final String sample : SampleListUtils.asList(samples)) {
            final int sampleIndex = result.sampleIndex(sample);
            Assert.assertTrue(sampleIndex >= 0);
            Assert.assertFalse(sampleIds.contains(sampleIndex));
            sampleIds.add(sampleIndex);

            final List<GATKSAMRecord> sampleReads = result.sampleReads(sampleIndex);
            final Set<GATKSAMRecord> sampleReadsSet = new HashSet<>(sampleReads);
            final List<GATKSAMRecord> expectedSampleReadArray = reads.get(sample);
            final Set<GATKSAMRecord> expectedSampleReadsSet = new HashSet<>(expectedSampleReadArray);
            Assert.assertEquals(sampleReadsSet,expectedSampleReadsSet);

            final int sampleReadCount = sampleReads.size();
            for (int r = 0; r < sampleReadCount; r++) {
                Assert.assertSame(sampleReads.get(r), expectedSampleReadArray.get(r));
                final int readIndex = result.readIndex(sampleIndex, sampleReads.get(r));
                Assert.assertEquals(readIndex,r);
            }
        }
    }

    private AlleleList<Allele> alleleList(final int alleleCount, final boolean hasReference) {
        final Allele[] alleles = AlleleListUnitTester.generateRandomAlleles(alleleCount,100);
        if (hasReference) {
            final int referenceIndex = rnd.nextInt(alleleCount);
            alleles[referenceIndex] = Allele.create(alleles[referenceIndex].getBases(),true);
        }
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
        if (alleleList.alleleCount() != alleles.length)
            throw new SkipException("repeated alleles, should be infrequent");
        return alleleList;
    }

    private SAMFileHeader SAM_HEADER = ArtificialSAMUtils.createArtificialSamHeader(10, 0, 1000);
    final GenomeLocParser locParser = new GenomeLocParser(SAM_HEADER.getSequenceDictionary());


    private int[][] READ_COUNTS = new int[][] {
            {},
            { 100 },
            { 0 },
            { 0, 0, 0 },
            { 1, 0, 1 },
            { 100, 10 , 100},
            { 1000, 10, 100, 20, 23 }
    };

    private int[] ALLELE_COUNTS = new int[] { 0, 1, 2, 3, 10, 20 };

    @DataProvider(name="readCountsAndAlleleCountData")
    public Object[][] readCountsAndAlleleCountData() {
        final Object[][] result = new Object[READ_COUNTS.length * ALLELE_COUNTS.length * 2][];
        int index = 0;
        for (final int[] readCounts : READ_COUNTS)
            for (final int alleleCount : ALLELE_COUNTS) {
                result[index++] = new Object[]{ readCounts, alleleCount, false};
                result[index++] = new Object[]{ readCounts, alleleCount, true};
            }
        return result;
    }

    @DataProvider(name="readCountsAndAlleleCountDataSkippingNoAlleleAndWithReference")
    public Object[][] readCountsAndAlleleCountDataSkippingNoAlleleAndWithReference() {
        final Object[][] raw = readCountsAndAlleleCountData();
        final List<Object[]> result = new ArrayList<>(raw.length);
        for (final Object[] paramSet : raw)
            if (!paramSet[2].equals(true) || !paramSet[1].equals(0))
                result.add(paramSet);
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider(name="readCountsAndAlleleCountDataSkippingNoLikelihoodsOrNoAlleleAndWithReference")
    public Object[][] readCountsAndAlleleCountDataSkippingNoLikelihoodsOrNoAlleleAndWithReference() {
        final Object[][] raw = readCountsAndAlleleCountDataSkippingNoAlleleAndWithReference();
        final List<Object[]> result = new ArrayList<>(raw.length);
        for (final Object[] paramSet : raw) {
            final int[] readCounts = (int[]) paramSet[0];
            final long totalReadCount = MathUtils.sum(readCounts);
            if (totalReadCount > 0)
                result.add(paramSet);
        }
        return result.toArray(new Object[result.size()][]);
    }

    private SampleList sampleList(final int[] readCounts) {
        final List<String> samples = new ArrayList<>(readCounts.length);
        for (int i = 0; i < readCounts.length; i++)
            samples.add("SAMPLE_" + i);
        return new IndexedSampleList(samples);
    }

}
