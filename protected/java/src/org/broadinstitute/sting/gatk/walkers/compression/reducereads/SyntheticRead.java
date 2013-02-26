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

package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Running Consensus is a read that is compressed as a sliding window travels over the reads
 * and keeps track of all the bases that are outside of variant regions.
 *
 * Consensus reads have qual fields that correspond to the number of reads that had the base
 * and passed the minimum quality threshold.
 *
 * The mapping quality of a consensus read is the average RMS of the mapping qualities of all reads
 * that compose the consensus
 *
 * @author Mauricio Carneiro
 * @since 8/26/11
 */
public class SyntheticRead {
    // Rather than storing a separate list for each attribute in SingleBaseInfo, store one list to reduce
    // memory footprint.
    // TODO: better name
    private static class SingleBaseInfo {
        byte baseIndexOrdinal; // enum BaseIndex.ordinal
        byte count;
        byte qual;
        byte insertionQual;
        byte deletionQual;

        SingleBaseInfo(byte baseIndexOrdinal, byte count, byte qual, byte insertionQual, byte deletionQual) {
            this.baseIndexOrdinal = baseIndexOrdinal;
            this.count = count;
            this.qual = qual;
            this.insertionQual = insertionQual;
            this.deletionQual = deletionQual;
        }
    }

    // This class is merely sharing of code for convertVariableGivenBases().
    private abstract class SingleBaseInfoIterator implements Iterator<Byte> {
        final Iterator<SingleBaseInfo> it;

        SingleBaseInfoIterator() {
            this.it = basesCountsQuals.iterator();
        }

        public boolean hasNext() {
            return it.hasNext();
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
    
        
    // Map from ordinal to enum value.
    private static final BaseIndex[] BaseIndexByOrdinal = new BaseIndex[BaseIndex.values().length];
    static 
    {
        for (final BaseIndex baseIndex : BaseIndex.values()) {
            BaseIndexByOrdinal[baseIndex.ordinal()] = baseIndex;
        }
    }
    
    
    private final List<SingleBaseInfo> basesCountsQuals;
    private double mappingQuality;                                                                                      // the average of the rms of the mapping qualities of all the reads that contributed to this consensus
    private String readTag;

    // Information to produce a GATKSAMRecord
    private SAMFileHeader header;
    private GATKSAMReadGroupRecord readGroupRecord;
    private String contig;
    private int contigIndex;
    private String readName;
    private int refStart;
    private boolean hasIndelQualities = false;
    private boolean isNegativeStrand = false;

    /**
     * Full initialization of the running consensus if you have all the information and are ready to
     * start adding to the running consensus.
     *
     * @param header          GATKSAMRecord file header
     * @param readGroupRecord Read Group for the GATKSAMRecord
     * @param contig          the read's contig name
     * @param contigIndex     the read's contig index
     * @param readName        the read's name
     * @param refStart        the alignment start (reference based)
     * @param readTag         the reduce reads tag for the synthetic read
     */
    public SyntheticRead(SAMFileHeader header, GATKSAMReadGroupRecord readGroupRecord, String contig, int contigIndex, String readName, int refStart, String readTag, boolean hasIndelQualities, boolean isNegativeRead) {
        final int initialCapacity = 10000;
        basesCountsQuals = new ArrayList<SingleBaseInfo>(initialCapacity);
        mappingQuality = 0.0;

        this.readTag = readTag;
        this.header = header;
        this.readGroupRecord = readGroupRecord;
        this.contig = contig;
        this.contigIndex = contigIndex;
        this.readName = readName;
        this.refStart = refStart;
        this.hasIndelQualities = hasIndelQualities;
        this.isNegativeStrand = isNegativeRead;
    }

    public SyntheticRead(List<BaseIndex> bases, List<Byte> counts, List<Byte> quals, List<Byte> insertionQuals, List<Byte> deletionQuals, double mappingQuality, String readTag, SAMFileHeader header, GATKSAMReadGroupRecord readGroupRecord, String contig, int contigIndex, String readName, int refStart, boolean hasIndelQualities, boolean isNegativeRead) {
        basesCountsQuals = new ArrayList<SingleBaseInfo>(bases.size());
        for (int i = 0; i < bases.size(); ++i) {
            basesCountsQuals.add(new SingleBaseInfo(bases.get(i).getOrdinalByte(), counts.get(i), quals.get(i), insertionQuals.get(i), deletionQuals.get(i)));
        }
        this.mappingQuality = mappingQuality;
        this.readTag = readTag;
        this.header = header;
        this.readGroupRecord = readGroupRecord;
        this.contig = contig;
        this.contigIndex = contigIndex;
        this.readName = readName;
        this.refStart = refStart;
        this.hasIndelQualities = hasIndelQualities;
        this.isNegativeStrand = isNegativeRead;
    }

    /**
     * Easy access to keep adding to a running consensus that has already been
     * initialized with the correct read name and refStart
     *
     * @param base   the base to add
     * @param count  number of reads with this base
     */
    @Requires("count <= Byte.MAX_VALUE")
    public void add(BaseIndex base, byte count, byte qual, byte insQual, byte delQual, double mappingQuality) {
        basesCountsQuals.add(new SingleBaseInfo(base.getOrdinalByte(), count, qual, insQual, delQual));
        this.mappingQuality += mappingQuality;
    }

    public BaseIndex getBase(final int readCoordinate) {
        return BaseIndexByOrdinal[basesCountsQuals.get(readCoordinate).baseIndexOrdinal];
    }

    public int getRefStart() {
        return refStart;
    }

    /**
     * Creates a GATKSAMRecord of the synthetic read. Will return null if the read is invalid.
     *
     * Invalid reads are :
     *   - exclusively composed of deletions
     *
     * @return a GATKSAMRecord or null
     */
    public GATKSAMRecord close () {
        if (isAllDeletions())
            return null;

        GATKSAMRecord read = new GATKSAMRecord(header);
        read.setReferenceName(contig);
        read.setReferenceIndex(contigIndex);
        read.setReadPairedFlag(false);
        read.setReadUnmappedFlag(false);
        read.setReadNegativeStrandFlag(isNegativeStrand);
        read.setCigar(buildCigar());                                        // the alignment start may change while building the cigar (leading deletions)
        read.setAlignmentStart(refStart);
        read.setReadName(readName);
        read.setBaseQualities(convertBaseQualities(), EventType.BASE_SUBSTITUTION);
        read.setReadBases(convertReadBases());
        read.setMappingQuality((int) Math.ceil(mappingQuality / basesCountsQuals.size()));
        read.setReadGroup(readGroupRecord);
        read.setAttribute(readTag, convertBaseCounts());

        if (hasIndelQualities) {
            read.setBaseQualities(convertInsertionQualities(), EventType.BASE_INSERTION);
            read.setBaseQualities(convertDeletionQualities(), EventType.BASE_DELETION);
        }

        return read;
    }

    /**
     * Checks if the synthetic read is composed exclusively of deletions
     *
     * @return true if it is, false if it isn't.
     */
    private boolean isAllDeletions() {
        for (SingleBaseInfo b : basesCountsQuals)
            if (b.baseIndexOrdinal != BaseIndex.D.getOrdinalByte())
                return false;
        return true;
    }

    public int size () {
        return basesCountsQuals.size();
    }

    private byte [] convertBaseQualities() {
        return convertVariableGivenBases(new SingleBaseInfoIterator() {
                public Byte next() {
                    return it.next().qual;
                }
            });
    }

    private byte [] convertInsertionQualities() {
        return convertVariableGivenBases(new SingleBaseInfoIterator() {
                public Byte next() {
                    return it.next().insertionQual;
                }
            });
    }

    private byte [] convertDeletionQualities() {
        return convertVariableGivenBases(new SingleBaseInfoIterator() {
                public Byte next() {
                    return it.next().deletionQual;
                }
            });
    }

    protected byte [] convertBaseCounts() {
        byte[] countsArray = convertVariableGivenBases(new SingleBaseInfoIterator() {
                public Byte next() {
                    return it.next().count;
                }
            });

        if (countsArray.length == 0)
            throw new ReviewedStingException("Reduced read has counts array of length 0");

        byte[] compressedCountsArray = new byte [countsArray.length];
        compressedCountsArray[0] = countsArray[0];
        for (int i = 1; i < countsArray.length; i++)
            compressedCountsArray[i] = (byte) MathUtils.bound(countsArray[i] - compressedCountsArray[0], Byte.MIN_VALUE, Byte.MAX_VALUE);

        return compressedCountsArray;
    }

    private byte [] convertReadBases() {
        byte [] readArray = new byte[getReadLengthWithNoDeletions()];
        int i = 0;
        for (final SingleBaseInfo singleBaseInfo : basesCountsQuals) {
            final BaseIndex baseIndex = BaseIndexByOrdinal[singleBaseInfo.baseIndexOrdinal];
            if (baseIndex != BaseIndex.D)
                readArray[i++] = baseIndex.getByte();
        }
        
        return readArray;
    }

    /**
     * Builds the cigar string for the synthetic read
     *
     * Warning: if the synthetic read has leading deletions, it will shift the refStart (alignment start) of the read.
     *
     * @return the cigar string for the synthetic read
     */
    private Cigar buildCigar() {
        LinkedList<CigarElement> cigarElements = new LinkedList<CigarElement>();
        CigarOperator cigarOperator = null;
        int length = 0;
        for (final SingleBaseInfo singleBaseInfo : basesCountsQuals) {
            final BaseIndex b = BaseIndexByOrdinal[singleBaseInfo.baseIndexOrdinal];
            CigarOperator op;
            switch (b) {
                case D:
                    op = CigarOperator.DELETION;
                    break;
                case I:
                    throw new ReviewedStingException("Trying to create an insertion in a synthetic read. This operation is currently unsupported.");
                default:
                    op = CigarOperator.MATCH_OR_MISMATCH;
                    break;
            }
            if (cigarOperator == null) {
                if (op == CigarOperator.D)                                  // read cannot start with a deletion
                    refStart++;                                             // if it does, we need to move the reference start forward
                else
                    cigarOperator = op;
            }
            else if (cigarOperator != op) {                                 // if this is a new operator, we need to close the previous one
                cigarElements.add(new CigarElement(length, cigarOperator)); // close previous operator
                cigarOperator = op;
                length = 0;
            }

            if (cigarOperator != null)                                      // only increment the length of the cigar element if we really added it to the read (no leading deletions)
                length++;
        }
        if (length > 0 && cigarOperator != CigarOperator.D)                 // read cannot end with a deletion
            cigarElements.add(new CigarElement(length, cigarOperator));     // add the last cigar element

        return new Cigar(cigarElements);
    }

    /**
     * Shared functionality for all conversion utilities
     *
     * @param variableIterator the list to convert
     * @return a converted variable given the bases and skipping deletions
     */

    private byte [] convertVariableGivenBases (Iterator<Byte> variableIterator) {
        byte [] variableArray = new byte[getReadLengthWithNoDeletions()];
        int i = 0;
        for (final SingleBaseInfo singleBaseInfo : basesCountsQuals) {
            byte count = variableIterator.next();
            if (singleBaseInfo.baseIndexOrdinal != BaseIndex.D.getOrdinalByte())
                variableArray[i++] = count;
        }
        return variableArray;

    }

    /**
     * Shared functionality for all conversion utilities
     *
     * @return the length of the read with no deletions
     */
    private int getReadLengthWithNoDeletions() {
        int readLength = basesCountsQuals.size();
        for (final SingleBaseInfo singleBaseInfo : basesCountsQuals)
            if (singleBaseInfo.baseIndexOrdinal == BaseIndex.D.getOrdinalByte())
                readLength--;
        return readLength;
    }


}
