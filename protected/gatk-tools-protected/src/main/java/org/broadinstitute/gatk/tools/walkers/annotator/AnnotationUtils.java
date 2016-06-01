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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.gatk.tools.walkers.variantutils.CombineGVCFs;
import org.broadinstitute.gatk.tools.walkers.variantutils.GenotypeGVCFs;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.gatk.tools.walkers.cancer.m2.MuTect2;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class AnnotationUtils {

    public static final String ANNOTATION_HC_WARN_MSG = " annotation will not be calculated, must be called from HaplotypeCaller or MuTect2";
    public static final int WARNINGS_LOGGED_SIZE = 3;

    /**
     * Helper function to convert a List of Doubles to a comma-separated String
     * @param valueList the ArrayList with Double data
     * @return a comma-separated String
     */
    public static String encodeValueList( final List<Double> valueList, final String precisionFormat ) {
        List<String> outputList = new ArrayList<>();
        for (Double d : valueList) {
            outputList.add(String.format(precisionFormat, d));
        }
        return StringUtils.join(outputList, ",");
    }

    /**
     * Helper function to convert a List of Strings to a comma-separated String
     * @param stringList the ArrayList with String data
     * @return a comma-separated String
     */
    public static String encodeStringList( final List<String> stringList) {
        return StringUtils.join(stringList, ",");
    }

    /**
     * Checks if the walker is compatible with allele-specific annotations
     */
     public static boolean walkerSupportsAlleleSpecificAnnotations(final AnnotatorCompatible walker) {
         return ((walker instanceof HaplotypeCaller) || (walker instanceof CombineGVCFs) || (walker instanceof GenotypeGVCFs));
     }

    /**
     * Checks if the walker should get raw annotation data
     */
    public static boolean walkerRequiresRawData(final AnnotatorCompatible walker) {
        return ((walker instanceof HaplotypeCaller && ((HaplotypeCaller) walker).emitReferenceConfidence()) || walker instanceof CombineGVCFs);
    }

    /**
     * Checks if the input data is appropriate
     *
     * @param annotation the input genotype annotation key name(s)
     * @param walker input walker
     * @param map input map for each read, holds underlying alleles represented by an aligned read, and corresponding relative likelihood.
     * @param g input genotype
     * @param warningsLogged array that enforces the warning is logged once for each caller
     * @param logger logger specific for each caller
     *
     * @return true if the walker is a HaplotypeCaller, the likelihood map is non-null and the genotype is non-null and called, false otherwise
     * @throws IllegalArgumentException if annotation, walker, g, warningsLogged, or logger are null.
     * @throws ReviewedGATKException if the size of warningsLogged is less than 3.
     */
    public static boolean isAppropriateInput(final String annotation, final AnnotatorCompatible walker, final PerReadAlleleLikelihoodMap map, final Genotype g, final boolean[] warningsLogged, final Logger logger) {

        if ( annotation == null ){
            throw new IllegalArgumentException("The input annotation cannot be null");
        }

        if ( walker == null ) {
            throw new IllegalArgumentException("The input walker cannot be null");
        }

        if ( g == null ) {
            throw new IllegalArgumentException("The input genotype cannot be null");
        }

        if ( warningsLogged == null ){
            throw new IllegalArgumentException("The input warnings logged cannot be null");
        }

        if ( logger == null ){
            throw new IllegalArgumentException("The input logger cannot be null");
        }

        if ( warningsLogged.length < WARNINGS_LOGGED_SIZE ){
            throw new ReviewedGATKException("Warnings logged array must have at least " + WARNINGS_LOGGED_SIZE + " elements, but has " + warningsLogged.length);
        }

        if ( !(walker instanceof HaplotypeCaller) && !(walker instanceof MuTect2)) {
            if ( !warningsLogged[0] ) {
                logger.warn(annotation + ANNOTATION_HC_WARN_MSG + ", not " + walker.getClass().getSimpleName());
                warningsLogged[0] = true;
            }
            return false;
        }

        if ( map == null ){
            if ( !warningsLogged[1] ) {
                logger.warn("Annotation will not be calculated, can only be used with likelihood based annotations in the HaplotypeCaller");
                warningsLogged[1] = true;
            }
            return false;
        }

        if ( !g.isCalled() ){
            if ( !warningsLogged[2] ) {
                logger.warn("Annotation will not be calculated, genotype is not called");
                warningsLogged[2] = true;
            }
            return false;
        }

        return true;
    }


    //this method is intended to reconcile uniquified sample names
    // it comes into play when calling this annotation from GenotypeGVCFs with --uniquifySamples because founderIds
    // is derived from the sampleDB, which comes from the input sample names, but vc will have uniquified (i.e. different)
    // sample names. Without this check, the founderIds won't be found in the vc and the annotation won't be calculated.
    protected static Set<String> validateFounderIDs(final Set<String> founderIds, final VariantContext vc) {
        Set<String> vcSamples = new HashSet<>();
        Set<String> returnIDs = founderIds;
        vcSamples.addAll(vc.getSampleNames());
        if (!vcSamples.isEmpty()) {
            if (founderIds != null) {
                vcSamples.removeAll(founderIds);
                if (vcSamples.equals(vc.getSampleNames()))
                    returnIDs = vc.getSampleNames();
            }
        }
        return returnIDs;
    }

    /**
     * Get the position of a variant within a read with respect to the closer end, accounting for hard clipped bases and low quality ends
     * Used by ReadPosRankSum annotations
     *
     * @param read  a read containing the variant
     * @param initialReadPosition   the position based on the modified, post-hard-clipped CIGAR
     * @return read position
     */
    public static int getFinalVariantReadPosition(final GATKSAMRecord read, final int initialReadPosition) {
        final int numAlignedBases = getNumAlignedBases(read);

        int readPos = initialReadPosition;
        //TODO: this doesn't work for the middle-right position if we index from zero
        if (initialReadPosition > numAlignedBases / 2) {
            readPos = numAlignedBases - (initialReadPosition + 1);
        }
        return readPos;

    }

    /**
     *
     * @param read  a read containing the variant
     * @return  the number of hard clipped and low qual bases at the read start (where start is the leftmost end w.r.t. the reference)
     */
    public static int getNumClippedBasesAtStart(final SAMRecord read) {
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        final CigarElement first = c.getCigarElement(0);

        int numStartClippedBases = 0;
        if (first.getOperator() == CigarOperator.H) {
            numStartClippedBases = first.getLength();
        }
        final byte[] unclippedReadBases = read.getReadBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        //TODO: this code may not even get used because HaplotypeCaller already hard clips low quality tails
        for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numStartClippedBases++;
            else
                break;

        }

        return numStartClippedBases;
    }


    /**
     *
     * @param read  a read containing the variant
     * @return  number of non-hard clipped, aligned bases (excluding low quality bases at either end)
     */
    //TODO: this is bizarre -- this code counts hard clips, but then subtracts them from the read length, which already doesn't count hard clips
    public static int getNumAlignedBases(final GATKSAMRecord read) {
        return read.getReadLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }


    /**
     *
     * @param read  a read containing the variant
     * @return  number of hard clipped and low qual bases at the read end (where end is right end w.r.t. the reference)
     */
    public static int getNumClippedBasesAtEnd(final GATKSAMRecord read) {
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        final byte[] unclippedReadBases = read.getReadBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        //TODO: this code may not even get used because HaplotypeCaller already hard clips low quality tails
        for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }

        return numEndClippedBases;
    }
}
