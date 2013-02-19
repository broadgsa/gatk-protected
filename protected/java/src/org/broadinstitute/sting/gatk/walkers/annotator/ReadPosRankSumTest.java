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

package org.broadinstitute.sting.gatk.walkers.annotator;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.gatk.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.variant.variantcontext.Allele;

import java.util.*;

/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele; if the alternate allele is only seen near the ends of reads this is indicative of error).
 * Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
 */
public class ReadPosRankSumTest extends RankSumTest implements StandardAnnotation {

    public List<String> getKeyNames() {
        return Arrays.asList("ReadPosRankSum");
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("ReadPosRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"));
    }

    protected void fillQualsFromPileup(final List<Allele> allAlleles,
                                       final int refLoc,
                                       final ReadBackedPileup pileup,
                                       final PerReadAlleleLikelihoodMap alleleLikelihoodMap,
                                       final List<Double> refQuals, final List<Double> altQuals) {

        if (alleleLikelihoodMap == null) {
            // use old UG SNP-based version if we don't have per-read allele likelihoods
            for ( final PileupElement p : pileup ) {
                if ( isUsableBase(p) && p.getRead().getCigar() != null ) {
                    int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, 0, 0);

                    readPos = getFinalReadPosition(p.getRead(),readPos);

                    if ( allAlleles.get(0).equals(Allele.create(p.getBase(), true)) ) {
                        refQuals.add((double)readPos);
                    } else if ( allAlleles.contains(Allele.create(p.getBase()))) {
                        altQuals.add((double)readPos);
                    }
                }
            }
            return;
        }

        for (Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : alleleLikelihoodMap.getLikelihoodReadMap().entrySet()) {
            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if (a.isNoCall())
                continue; // read is non-informative

            final GATKSAMRecord read = el.getKey();
            if ( read.getSoftStart() + read.getCigar().getReadLength() <= refLoc ) { // make sure the read actually covers the requested ref loc
                continue;
            }
            final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate( read.getSoftStart(), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true );
            if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || read.getCigar() == null )
                continue;
            int readPos = AlignmentUtils.calcAlignmentByteArrayOffset( read.getCigar(), offset, false, 0, 0 );
            final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
            if (readPos > numAlignedBases / 2)
                readPos = numAlignedBases - (readPos + 1);

            if (a.isReference())
                refQuals.add((double)readPos);
            else if (allAlleles.contains(a))
                altQuals.add((double)readPos);
        }
    }

    int getFinalReadPosition(GATKSAMRecord read, int initialReadPosition) {
        final int numAlignedBases = getNumAlignedBases(read);

        int readPos = initialReadPosition;
        if (initialReadPosition > numAlignedBases / 2) {
            readPos = numAlignedBases - (initialReadPosition + 1);
        }
        return readPos;

    }
    int getNumClippedBasesAtStart(SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        final CigarElement first = c.getCigarElement(0);

        int numStartClippedBases = 0;
        if (first.getOperator() == CigarOperator.H) {
            numStartClippedBases = first.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numStartClippedBases++;
            else
                break;

        }

        return numStartClippedBases;
    }

    int getNumAlignedBases(SAMRecord read) {
        return read.getReadLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }

    int getNumClippedBasesAtEnd(SAMRecord read) {
        // compute total number of clipped bases (soft or hard clipped)
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        byte[] unclippedReadBases = read.getReadBases();
        byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < PairHMMIndelErrorModel.BASE_QUAL_THRESHOLD)
                numEndClippedBases++;
            else
                break;
        }


        return numEndClippedBases;
    }

    int getOffsetFromClippedReadStart(SAMRecord read, int offset) {
        return offset - getNumClippedBasesAtStart(read);
    }
}
