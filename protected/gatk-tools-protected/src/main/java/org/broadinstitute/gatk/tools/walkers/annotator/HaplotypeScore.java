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

package org.broadinstitute.gatk.tools.walkers.annotator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.io.Serializable;
import java.util.*;

/**
 * Consistency of the site with strictly two segregating haplotypes
 *
 * <p>For diploid organisms, barring chromosomal abnormalities, we expect that any given sample has no more than 2 segregating haplotypes at a given site. If there is evidence for more
 * than 2 segregating haplotypes, the read data should be considered suspect and the evidence artifactual. Higher scores are indicative of regions with bad alignments, typically leading to artifactual SNP and indel calls.</p>
 *
 * <h3>Caveats</h3>
 * <p>HaplotypeCaller does not output this annotation because it already evaluates haplotype segregation internally. This annotation is only informative (and available) for variants called by Unified Genotyper.</p>
 */
public class HaplotypeScore extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
    private final static Logger logger = Logger.getLogger(HaplotypeScore.class);
    private boolean walkerIdentityCheckWarningLogged = false;

    private final static boolean DEBUG = false;
    private final static int MIN_CONTEXT_WING_SIZE = 10;
    private final static int MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER = 50;
    private final static char REGEXP_WILDCARD = '.';

    @Override
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        // Can only call from UnifiedGenotyper
        if ( !(walker instanceof UnifiedGenotyper) ) {
            if ( !walkerIdentityCheckWarningLogged ) {
                if ( walker != null )
                    logger.warn("Annotation will not be calculated, must be called from UnifiedGenotyper, not " + walker.getClass().getName());
                else
                    logger.warn("Annotation will not be calculated, must be called from UnifiedGenotyper");
                walkerIdentityCheckWarningLogged = true;
            }
            return null;
        }

        if (vc.isSNP() && stratifiedContexts != null)
            return annotatePileup(ref, stratifiedContexts, vc);
        else
            return null;
    }

    private Map<String, Object> annotatePileup(final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc) {

        if (stratifiedContexts.isEmpty()) // empty means that call was made by someone else and we have no data here
            return null;

        final AlignmentContext context = AlignmentContextUtils.joinContexts(stratifiedContexts.values());

        final int contextWingSize = Math.min((ref.getWindow().size() - 1) / 2, MIN_CONTEXT_WING_SIZE);
        final int contextSize = contextWingSize * 2 + 1;

        final int locus = ref.getLocus().getStart() + (ref.getLocus().getStop() - ref.getLocus().getStart()) / 2;

        final ReadBackedPileup pileup = context.getBasePileup();

        // Compute all haplotypes consistent with the current read pileup
        final List<Haplotype> haplotypes = computeHaplotypes(pileup, contextSize, locus, vc);

        final MathUtils.RunningAverage scoreRA = new MathUtils.RunningAverage();
        if (haplotypes != null) {
            for (final Genotype genotype : vc.getGenotypes()) {
                final AlignmentContext thisContext = stratifiedContexts.get(genotype.getSampleName());
                if (thisContext != null) {
                    final ReadBackedPileup thisPileup = thisContext.getBasePileup();
                    scoreRA.add(scoreReadsAgainstHaplotypes(haplotypes, thisPileup, contextSize, locus)); // Taking the simple average of all sample's score since the score can be negative and the RMS doesn't make sense
                }
            }
        }

        // annotate the score in the info field
        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%.4f", scoreRA.mean()));
        return map;
    }

    private static class HaplotypeComparator implements Comparator<Haplotype>, Serializable {

        public int compare(Haplotype a, Haplotype b) {
            if (a.getQualitySum() < b.getQualitySum())
                return 1;
            if (a.getQualitySum() > b.getQualitySum()) {
                return -1;
            }
            return 0;
        }
    }

    private List<Haplotype> computeHaplotypes(final ReadBackedPileup pileup, final int contextSize, final int locus, final VariantContext vc) {
        // Compute all possible haplotypes consistent with current pileup

        int haplotypesToCompute = vc.getAlternateAlleles().size() + 1;

        final PriorityQueue<Haplotype> candidateHaplotypeQueue = new PriorityQueue<>(100, new HaplotypeComparator());
        final PriorityQueue<Haplotype> consensusHaplotypeQueue = new PriorityQueue<>(MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER, new HaplotypeComparator());

        for (final PileupElement p : pileup) {
            final Haplotype haplotypeFromRead = getHaplotypeFromRead(p, contextSize, locus);
            if ( haplotypeFromRead != null )
                candidateHaplotypeQueue.add(haplotypeFromRead);
        }

        // Now that priority queue has been built with all reads at context, we need to merge and find possible segregating haplotypes
        Haplotype elem;
        while ((elem = candidateHaplotypeQueue.poll()) != null) {
            boolean foundHaplotypeMatch = false;
            Haplotype lastCheckedHaplotype = null;
            for (final Haplotype haplotypeFromList : consensusHaplotypeQueue) {
                final Haplotype consensusHaplotype = getConsensusHaplotype(elem, haplotypeFromList);
                if (consensusHaplotype != null) {
                    foundHaplotypeMatch = true;
                    if (consensusHaplotype.getQualitySum() > haplotypeFromList.getQualitySum()) {
                        consensusHaplotypeQueue.remove(haplotypeFromList);
                        consensusHaplotypeQueue.add(consensusHaplotype);
                    }
                    break;
                } else {
                    lastCheckedHaplotype = haplotypeFromList;
                }
            }

            if (!foundHaplotypeMatch && consensusHaplotypeQueue.size() < MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER) {
                consensusHaplotypeQueue.add(elem);
            } else if (!foundHaplotypeMatch && lastCheckedHaplotype != null && elem.getQualitySum() > lastCheckedHaplotype.getQualitySum()) {
                consensusHaplotypeQueue.remove(lastCheckedHaplotype);
                consensusHaplotypeQueue.add(elem);
            }
        }

        // Now retrieve the N most popular haplotypes
        if (consensusHaplotypeQueue.size() > 0) {
            // The consensus haplotypes are in a quality-ordered priority queue, so the best haplotypes are just the ones at the front of the queue
            final Haplotype haplotype1 = consensusHaplotypeQueue.poll();

            List<Haplotype> hlist = new ArrayList<>();
            hlist.add(new Haplotype(haplotype1.getBases(), 60));

            for (int k = 1; k < haplotypesToCompute; k++) {
                Haplotype haplotype2 = consensusHaplotypeQueue.poll();
                if (haplotype2 == null) {
                    haplotype2 = haplotype1;
                } // Sometimes only the reference haplotype can be found
                hlist.add(new Haplotype(haplotype2.getBases(), 20));
            }
            return hlist;
        } else
            return null;
    }

    /**
     * Return a haplotype object constructed from the read or null if read's cigar is null
     *
     * @param p                pileup element representing the read
     * @param contextSize      the context size to use
     * @param locus            the position
     * @return possibly null Haplotype object constructed from the read
     */
    private Haplotype getHaplotypeFromRead(final PileupElement p, final int contextSize, final int locus) {
        final GATKSAMRecord read = p.getRead();
        if ( read.getCigar() == null )
            return null;

        final byte[] haplotypeBases = new byte[contextSize];
        Arrays.fill(haplotypeBases, (byte) REGEXP_WILDCARD);
        final byte[] baseQualities = new byte[contextSize];
        Arrays.fill(baseQualities, (byte)0);

        byte[] readBases = read.getReadBases();
        readBases = AlignmentUtils.readToAlignmentByteArray(read.getCigar(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getBaseQualities();
        readQuals = AlignmentUtils.readToAlignmentByteArray(read.getCigar(), readQuals); // Shift the location of the qual scores based on the Cigar string

        final int readOffsetFromPileup = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), p, read.getAlignmentStart(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1) / 2;

        for (int i = 0; i < contextSize; i++) {
            final int baseOffset = i + baseOffsetStart;
            if (baseOffset < 0) {
                continue;
            }
            if (baseOffset >= readBases.length) {
                break;
            }
            if (readQuals[baseOffset] == PileupElement.DELETION_BASE) {
                readQuals[baseOffset] = PileupElement.DELETION_QUAL;
            }
            if (!BaseUtils.isRegularBase(readBases[baseOffset])) {
                readBases[baseOffset] = (byte) REGEXP_WILDCARD;
                readQuals[baseOffset] = (byte) 0;
            } // N's shouldn't be treated as distinct bases
            readQuals[baseOffset] = (byte) Math.min((int) readQuals[baseOffset], p.getMappingQual());
            if (((int) readQuals[baseOffset]) < 5) {
                readQuals[baseOffset] = (byte) 0;
            } // quals less than 5 are used as codes and don't have actual probabilistic meaning behind them
            haplotypeBases[i] = readBases[baseOffset];
            baseQualities[i] = readQuals[baseOffset];
        }

        return new Haplotype(haplotypeBases, baseQualities);
    }

    private Haplotype getConsensusHaplotype(final Haplotype haplotypeA, final Haplotype haplotypeB) {
        final byte[] a = haplotypeA.getBases();
        final byte[] b = haplotypeB.getBases();

        if (a.length != b.length) {
            throw new ReviewedGATKException("Haplotypes a and b must be of same length");
        }

        byte chA, chB;
        final byte wc = (byte) REGEXP_WILDCARD;

        final int length = a.length;
        final byte[] consensusChars = new byte[length];
        final int[] consensusQuals = new int[length];

        final int[] qualsA = haplotypeA.getQuals();
        final int[] qualsB = haplotypeB.getQuals();

        for (int i = 0; i < length; i++) {
            chA = a[i];
            chB = b[i];

            if ((chA != chB) && (chA != wc) && (chB != wc))
                return null;

            if ((chA == wc) && (chB == wc)) {
                consensusChars[i] = wc;
                consensusQuals[i] = 0;
            } else if ((chA == wc)) {
                consensusChars[i] = chB;
                consensusQuals[i] = qualsB[i];
            } else if ((chB == wc)) {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i];
            } else {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i] + qualsB[i];
            }
        }

        return new Haplotype(consensusChars, consensusQuals);
    }

    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(final List<Haplotype> haplotypes, final ReadBackedPileup pileup, final int contextSize, final int locus) {
        if (DEBUG) System.out.printf("HAP1: %s%n", haplotypes.get(0));
        if (DEBUG) System.out.printf("HAP2: %s%n", haplotypes.get(1));

        final ArrayList<double[]> haplotypeScores = new ArrayList<>();
        for (final PileupElement p : pileup) {
            // Score all the reads in the pileup, even the filtered ones
            final double[] scores = new double[haplotypes.size()];
            for (int i = 0; i < haplotypes.size(); i++) {
                final Haplotype haplotype = haplotypes.get(i);
                final double score = scoreReadAgainstHaplotype(p, contextSize, haplotype, locus);
                scores[i] = score;
                if (DEBUG) {
                    System.out.printf("  vs. haplotype %d = %f%n", i, score);
                }
            }
            haplotypeScores.add(scores);
        }

        double overallScore = 0.0;
        for (final double[] readHaplotypeScores : haplotypeScores) {
            overallScore += MathUtils.arrayMin(readHaplotypeScores);
        }

        return overallScore;
    }

    private double scoreReadAgainstHaplotype(final PileupElement p, final int contextSize, final Haplotype haplotype, final int locus) {
        double expected = 0.0;
        double mismatches = 0.0;

        final GATKSAMRecord read = p.getRead();
        if ( read.getCigar() == null )
            return 0.0;

        // What's the expected mismatch rate under the model that this read is actually sampled from
        // this haplotype?  Let's assume the consensus base c is a random choice one of A, C, G, or T, and that
        // the observed base is actually from a c with an error rate e.  Since e is the rate at which we'd
        // see a miscalled c, the expected mismatch rate is really e.  So the expected number of mismatches
        // is just sum_i e_i for i from 1..n for n sites
        //
        // Now, what's the probabilistic sum of mismatches?  Suppose that the base b is equal to c.  Well, it could
        // actually be a miscall in a matching direction, which would happen at a e / 3 rate.  If b != c, then
        // the chance that it is actually a mismatch is 1 - e, since any of the other 3 options would be a mismatch.
        // so the probability-weighted mismatch rate is sum_i ( matched ? e_i / 3 : 1 - e_i ) for i = 1 ... n
        final byte[] haplotypeBases = haplotype.getBases();
        byte[] readBases = read.getReadBases();

        readBases = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getBaseQualities();
        readQuals = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readQuals); // Shift the location of the qual scores based on the Cigar string
        int readOffsetFromPileup = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, read.getAlignmentStart(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1) / 2;

        for (int i = 0; i < contextSize; i++) {
            final int baseOffset = i + baseOffsetStart;
            if (baseOffset < 0) {
                continue;
            }
            if (baseOffset >= readBases.length) {
                break;
            }

            final byte haplotypeBase = haplotypeBases[i];
            final byte readBase = readBases[baseOffset];

            final boolean matched = (readBase == haplotypeBase || haplotypeBase == (byte) REGEXP_WILDCARD);
            byte qual = readQuals[baseOffset];
            if (qual == PileupElement.DELETION_BASE) {
                qual = PileupElement.DELETION_QUAL;
            } // calcAlignmentByteArrayOffset fills the readQuals array with DELETION_BASE at deletions
            qual = (byte) Math.min((int) qual, p.getMappingQual());
            if (((int) qual) >= 5) { // quals less than 5 are used as codes and don't have actual probabilistic meaning behind them
                final double e = QualityUtils.qualToErrorProb(qual);
                expected += e;
                mismatches += matched ? e : 1.0 - e / 3.0;
            }

            // a more sophisticated calculation would include the reference quality, but it's nice to actually penalize
            // the mismatching of poorly determined regions of the consensus
        }

        return mismatches - expected;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.HAPLOTYPE_SCORE_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    private static class Haplotype  {
        private final byte[] bases;
        private final int[] quals;
        private int qualitySum = -1;

        public Haplotype( final byte[] bases, final int[] quals ) {
            this.bases = bases;
            this.quals = quals;
        }

        public Haplotype( final byte[] bases, final int qual ) {
            this.bases = bases;
            quals = new int[bases.length];
            Arrays.fill(quals, qual);
        }

        public Haplotype( final byte[] bases, final byte[] quals ) {
            this.bases = bases;
            this.quals = new int[quals.length];
            for ( int i = 0 ; i < quals.length; i++ )
                this.quals[i] = (int)quals[i];
        }

        public double getQualitySum() {
            if ( qualitySum == -1 ) {
                qualitySum = 0;
                for ( final int qual : quals ) {
                    qualitySum += qual;
                }
            }
            return qualitySum;
        }

        public int[] getQuals() {
            return quals.clone();
        }

        public byte[] getBases() {
            return bases.clone();
        }
    }
}
