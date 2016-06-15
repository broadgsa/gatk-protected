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

package org.broadinstitute.gatk.tools.walkers.cancer;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.cancer.m2.MuTect2;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Detect clustering of variants near the ends of reads
 *
 * <p> This annotation detects clustering of evidence for a somatic variant near the ends of reads. To turn on the annotation and the accompanying filter (clustered_read_position), add --enable_clustered_read_position_filter flag in the commandline.
 *
 *
 * <h3>Statistical notes</h3>
 * <p> ClusteredReadPosition produces four INFO field annotations. At a given somatic variant site, MEDIAN_LEFT_OFFSET is the median of the number of bases from the left end of the tumor read to the variant. MEDIAN_RIGHT_OFFSET is similar, but counts from the right end of the read. MAD_LEFT_OFFSET and MAD_RIGHT_OFFSET measure the median absolute deviations. The median gives us the offset of a representative read, while the median absolute deviation captures the spread. We filter a variant if MEDIAN_LEFT_OFFSET <= 10 and MAD_LEFT_OFFSET <= 3, or if MEDIAN_RIGHT_OFFSET <= 10 and MAD_RIGHT_OFFSET <= 3.
 *
 *
 * <h3>Caveat</h3>
 * <p> ClusteredReadPosition is available with MuTect2 only </p>
 *
 * <h3>RelatedAnnotation</h3>
 * <li><b><a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php">ReadPosRankSum</a></b> is a similar annotation designed for germline variants.
 *
 */
public class ClusteredReadPosition extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {
    private final static Logger logger = Logger.getLogger(ClusteredReadPosition.class);
    private String tumorSampleName = null;

    @Override
    public List<String> getKeyNames() { return Arrays.asList(
            GATKVCFConstants.MEDIAN_LEFT_OFFSET_KEY,
            GATKVCFConstants.MEDIAN_RIGHT_OFFSET_KEY,
            GATKVCFConstants.MAD_MEDIAN_LEFT_OFFSET_KEY,
            GATKVCFConstants.MAD_MEDIAN_RIGHT_OFFSET_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        List<VCFInfoHeaderLine> descriptions = new ArrayList<>();
        for (final String infoFieldKey : getKeyNames()){
            descriptions.add(GATKVCFHeaderLines.getInfoLine(infoFieldKey));
        }
        return descriptions;

        // the following causes a cryptic class not found error, similar to the one in computeReadPositionStats
        // return getKeyNames().stream().map(GATKVCFHeaderLines::getInfoLine).collect(Collectors.toList());
    }

    @Override
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        // TODO: might make sense to move this code to SomaticGenoypingEngine.
        // FIXME: checking walker is mutect2 is not ideal...moving this annotation to SomaticGenoypingEngine will solve it

        // populate tumorSampleName the first time we call this method. skip afterwards.
        if (tumorSampleName == null){
            if (walker instanceof MuTect2) {
                tumorSampleName = ((MuTect2) walker).getTumorSampleName();
            } else {
                throw new IllegalStateException("ClusteredReadPosition: walker is not MuTect2");
            }
        }

        // we skip multi-allelic sites
        if (vc.getAlternateAlleles().size() > 1){
            return null;
        }

        final Map<String, Object> result = new HashMap<>();

        if ( stratifiedPerReadAlleleLikelihoodMap != null ) {
            final PerReadAlleleLikelihoodMap likelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(tumorSampleName);
            if ( likelihoodMap != null && !likelihoodMap.isEmpty() ) {
                final Optional<MedianStatistics> readPositionStatsOption = computeReadPositionStats(vc, likelihoodMap);
                if (readPositionStatsOption.isPresent()){
                    MedianStatistics readPositionStats = readPositionStatsOption.get();
                    result.put(GATKVCFConstants.MEDIAN_LEFT_OFFSET_KEY, readPositionStats.getLeftMedian());
                    result.put(GATKVCFConstants.MEDIAN_RIGHT_OFFSET_KEY, readPositionStats.getRightMedian());
                    result.put(GATKVCFConstants.MAD_MEDIAN_LEFT_OFFSET_KEY, readPositionStats.getLeftMAD());
                    result.put(GATKVCFConstants.MAD_MEDIAN_RIGHT_OFFSET_KEY, readPositionStats.getRightMAD());
                } else {
                    return null;
                }
            }
        }

        return result;
    }

    /**
     *
     * @param vc
     * @param pralm
     * @return median of left and right offsets and their median absolute deviations. does not return null.
     */
    private Optional<MedianStatistics> computeReadPositionStats(final VariantContext vc,
                                                                final PerReadAlleleLikelihoodMap pralm) {
        final int variantStartPosition = vc.getStart();
        final List<Integer> tumorLeftOffsets = new ArrayList<>();
        final List<Integer> tumorRightOffsets = new ArrayList<>();
        for ( final Map.Entry<GATKSAMRecord, Map<Allele,Double>> readAlleleLikelihood : pralm.getLikelihoodReadMap().entrySet() ) {
            final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(readAlleleLikelihood.getValue());
            final GATKSAMRecord read = readAlleleLikelihood.getKey();
            if ( mostLikelyAllele.getMostLikelyAllele().isReference() || ! mostLikelyAllele.isInformative() || ! isUsableRead(read)) {
                continue;
            }

            final Pair<OptionalInt, OptionalInt> offsetPair = getVariantPositionInRead(read, variantStartPosition);
            final OptionalInt variantPositionInReadFromLeft = offsetPair.getFirst();
            final OptionalInt variantPositionInReadFromRight = offsetPair.getSecond();

            // suffices to check only the left offset because the right offset depends on it
            if ( variantPositionInReadFromLeft.isPresent() ) {
                tumorLeftOffsets.add(variantPositionInReadFromLeft.getAsInt());
                tumorRightOffsets.add(variantPositionInReadFromRight.getAsInt());
            }
        }

        if (tumorLeftOffsets.isEmpty() || tumorRightOffsets.isEmpty()) {
            // This condition seems to arise when the reads as aligned in the bam (as represented by PRALM) do not contain the alt read found by HaplotypeCaller
            logger.warn("At Position " + vc.getContig() + ": " + vc.getStart() + " , the left or right offset list is empty");
            return Optional.empty();
        }

        // The following (mapToDouble() in particular) causes ClusteredReadPosition to be not added to ClassMap
        // leftMedian = median.evaluate(tumorLeftOffsets.stream().mapToDouble( x -> x ).toArray());
        // rightMedian = median.evaluate(tumorRightOffsets.stream().mapToDouble( x -> x).toArray());

        // until we understand why mapToDouble() causes the above error, have to compute medians in two steps
        // first use a for loop to manually cast integer to doubles, then call median :: evaluate
        double[] tumorLeftOffsetsDouble = new double[tumorLeftOffsets.size()];
        double[] tumorRightOffsetsDouble = new double[tumorRightOffsets.size()];
        for (int i = 0; i < tumorLeftOffsets.size(); i++){
            tumorLeftOffsetsDouble[i] = (double) tumorLeftOffsets.get(i);
            tumorRightOffsetsDouble[i] = (double) tumorRightOffsets.get(i);
        }

        Median median = new Median();
        double leftMedian = median.evaluate(tumorLeftOffsetsDouble);
        double rightMedian = median.evaluate(tumorRightOffsetsDouble);
        double leftMAD = calculateMAD(tumorLeftOffsets, leftMedian);
        double rightMAD = calculateMAD(tumorRightOffsets, rightMedian);

        return( Optional.of(new MedianStatistics(leftMedian, rightMedian, leftMAD, rightMAD) ) );
    }

    private static class MedianStatistics {
        private double leftMedian;
        private double rightMedian;
        private double leftMAD;
        private double rightMAD;

        public MedianStatistics(double leftMedian, double rightMedian, double leftMAD, double rightMAD) {
            this.leftMedian = leftMedian;
            this.rightMedian = rightMedian;
            this.leftMAD = leftMAD;
            this.rightMAD = rightMAD;
        }

        public double getLeftMedian() {
            return leftMedian;
        }

        public double getRightMedian() {
            return rightMedian;
        }

        public double getLeftMAD() {
            return leftMAD;
        }

        public double getRightMAD() {
            return rightMAD;
        }
    }


    /**
     Examples below show how we compute the position of the variant with respect to the left and right end of the reads.
     Note that a variant may be SNP, deletion, or insertion, and we are counting the number of bases from the left/right end of the read to that variant.
     We first compute the left offset. Then, right offset = read length - left offset.
     This means that if there is an insertion between the either end of a read and the variant, we count the inserted bases. Conversely, we do not count the deleted bases between the end of a read and a variant.
     We count soft-clipped bases.

     example 1 : SNP

     right offset: 9 8 7 6 5 4 3 2 1 0
     ref:          _ _ _ _ _ _ _ _ _ _
     read:         _ _ _ _ x _ _ _ _ _
     left offset:  0 1 2 3 4 5 6 7 8 9

     left-offset = 4. right offset = 5.
     read.getReadLength() = 10. numReadBasesToVariant = 5.

     example 2: deletion

     We count from the left end of the read to the last non-deleted base i.e. the first deleted base is not counted.
     From the right end, we count bases to the *end* of the deletion.

     right offset: 9 8 7 6 5 4 3 2 1 0
     ref:          _ _ _ _ _ _ _ _ _ _
     read:         _ _ _ _|- - - -|_ _
     left offset:  0 1 2 3 4 5 6 7 8 9

     left-offset = 3. right-offset = 2.
     read.getReadLength() = 6. numReadBasesToVariant = 4

     example 3: insertion

     For insertions, we count from the left to the first inserted base. From the right, we count all the way to the first inserted base.
     In the future, we may modify this; it might be desirable to count from the right to the *last* inserted base.

     right offset: 9 8 7 6 5 4 3 2 1 0
     ref:          _ _ _ _     _ _ _ _
     read:         _ _ _ I I I _ _ _ _
     left offset:  0 1 2 3 4 5 6 7 8 9

     left-offset = 3. right offset = 6
     read.getReadLength() = 10. numReadBasesToVariant = 4.

    */

    /**
     * The function assumes that read contains the variant allele.
     *
     * @param read
     * @param variantStartPosition the location of the variant in the reference
     * @return
     */

    protected Pair<OptionalInt, OptionalInt> getVariantPositionInRead(final GATKSAMRecord read, final int variantStartPosition) {
        final Pair<Integer, Boolean> refPositionAndDeletionFlag = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), variantStartPosition, true);
        // the +1 is needed there because getReadCoordinateForReferenceCoordinate() returns the number of read bases from the left end to the variant - 1
        int numReadBasesFromLeftEndToVariant = refPositionAndDeletionFlag.getFirst() + 1;

        // we don't take advantage of fallsInsideOrJustBeforeDeletionOrSkippedRegion flag now, but we might want to, so I will leave it here in comments.
        // boolean fallsInsideOrJustBeforeDeletionOrSkippedRegion = refPositionAndDeletionFlag.getSecond();

        if ( numReadBasesFromLeftEndToVariant == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            return new Pair(OptionalInt.empty(), OptionalInt.empty());
        } else {
            int leftOffset = numReadBasesFromLeftEndToVariant - 1;
            int rightOffset = read.getReadLength() - numReadBasesFromLeftEndToVariant;
            return new Pair(OptionalInt.of(leftOffset), OptionalInt.of(rightOffset));
        }
    }

    /**
     * Can the read be used in comparative tests between ref / alt bases?
     *
     * @param read   the read to consider
     * @return false if MQ is either 0 or unavailable. true otherwise.
     */
    private boolean isUsableRead(final GATKSAMRecord read) {
        return( read.getMappingQuality() != 0 || read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
    }

    /**
     *
     * @param offsets a list of integers
     * @param median median of the list offsets.
     * @return median absolute deviation (median of the list of deviations from the median)
     */
    private double calculateMAD(final List<Integer> offsets, final double median) {
        // This code is concise but somehow leads to ClusteredReadPosition class being removed from ClassMap.
        // mapToDouble() seems to be the trigger
        // return new Median().evaluate(offsets.stream().mapToDouble(x -> Math.abs(x - median)).toArray());

        double[] medianAbsoluteDeviations = new double[offsets.size()];
        for (int i = 0; i < offsets.size(); i++){
            medianAbsoluteDeviations[i] = Math.abs(offsets.get(i) - median);
        }

        return new Median().evaluate(medianAbsoluteDeviations);
    }
}
