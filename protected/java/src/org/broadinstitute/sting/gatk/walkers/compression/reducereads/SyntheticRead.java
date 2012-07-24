package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.gatk.walkers.bqsr.EventType;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
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
    private List<BaseIndex> bases;
    private List<Byte> counts;
    private List<Byte> quals;
    private List<Byte> insertionQuals;
    private List<Byte> deletionQuals;
    private double mappingQuality;                                                                                      // the average of the rms of the mapping qualities of all the reads that contributed to this consensus
    private String readTag;

    // Information to produce a GATKSAMRecord
    private SAMFileHeader header;
    private GATKSAMReadGroupRecord readGroupRecord;
    private String contig;
    private int contigIndex;
    private String readName;
    private Integer refStart;
    private boolean hasIndelQualities = false;

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
    public SyntheticRead(SAMFileHeader header, GATKSAMReadGroupRecord readGroupRecord, String contig, int contigIndex, String readName, Integer refStart, String readTag, boolean hasIndelQualities) {
        final int initialCapacity = 10000;
        bases = new ArrayList<BaseIndex>(initialCapacity);
        counts = new ArrayList<Byte>(initialCapacity);
        quals = new ArrayList<Byte>(initialCapacity);
        insertionQuals = new ArrayList<Byte>(initialCapacity);
        deletionQuals = new ArrayList<Byte>(initialCapacity);
        mappingQuality = 0.0;

        this.readTag = readTag;
        this.header = header;
        this.readGroupRecord = readGroupRecord;
        this.contig = contig;
        this.contigIndex = contigIndex;
        this.readName = readName;
        this.refStart = refStart;
        this.hasIndelQualities = hasIndelQualities;
    }

    public SyntheticRead(List<BaseIndex> bases, List<Byte> counts, List<Byte> quals, List<Byte> insertionQuals, List<Byte> deletionQuals, double mappingQuality, String readTag, SAMFileHeader header, GATKSAMReadGroupRecord readGroupRecord, String contig, int contigIndex, String readName, Integer refStart, boolean hasIndelQualities) {
        this.bases = bases;
        this.counts = counts;
        this.quals = quals;
        this.insertionQuals = insertionQuals;
        this.deletionQuals = deletionQuals;
        this.mappingQuality = mappingQuality;
        this.readTag = readTag;
        this.header = header;
        this.readGroupRecord = readGroupRecord;
        this.contig = contig;
        this.contigIndex = contigIndex;
        this.readName = readName;
        this.refStart = refStart;
        this.hasIndelQualities = hasIndelQualities;
    }

    /**
     * Easy access to keep adding to a running consensus that has already been
     * initialized with the correct read name and refStart
     *
     * @param base   the base to add
     * @param count  number of reads with this base
     */
    @Requires("count < Byte.MAX_VALUE")
    public void add(BaseIndex base, byte count, byte qual, byte insQual, byte delQual, double mappingQuality) {
        counts.add(count);
        bases.add(base);
        quals.add(qual);
        insertionQuals.add(insQual);
        deletionQuals.add(delQual);
        this.mappingQuality += mappingQuality;
    }

    public BaseIndex getBase(int readCoordinate) {
        return bases.get(readCoordinate);
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
        read.setCigar(buildCigar());                                        // the alignment start may change while building the cigar (leading deletions)
        read.setAlignmentStart(refStart);
        read.setReadName(readName);
        read.setBaseQualities(convertBaseQualities(), EventType.BASE_SUBSTITUTION);
        read.setReadBases(convertReadBases());
        read.setMappingQuality((int) Math.ceil(mappingQuality / bases.size()));
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
        for (BaseIndex b : bases)
            if (b != BaseIndex.D)
                return false;
        return true;
    }

    public int size () {
        return bases.size();
    }

    private byte [] convertBaseQualities() {
        return convertVariableGivenBases(bases, quals);
    }

    private byte [] convertInsertionQualities() {
        return convertVariableGivenBases(bases, insertionQuals);
    }

    private byte [] convertDeletionQualities() {
        return convertVariableGivenBases(bases, deletionQuals);
    }

    protected byte [] convertBaseCounts() {
        byte[] countsArray = convertVariableGivenBases(bases, counts);

        if (countsArray.length == 0)
            throw new ReviewedStingException("Reduced read has counts array of length 0");

        byte[] compressedCountsArray = new byte [countsArray.length];
        compressedCountsArray[0] = countsArray[0];
        for (int i = 1; i < countsArray.length; i++)
            compressedCountsArray[i] = (byte) MathUtils.bound(countsArray[i] - compressedCountsArray[0], Byte.MIN_VALUE, Byte.MAX_VALUE);

        return compressedCountsArray;
    }

    private byte [] convertReadBases() {
        byte [] readArray = new byte[getReadLengthWithNoDeletions(bases)];
        int i = 0;
        for (BaseIndex baseIndex : bases)
            if (baseIndex != BaseIndex.D)
                readArray[i++] = baseIndex.getByte();

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
        for (BaseIndex b : bases) {
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
     * @param bases    the read bases
     * @param variable the list to convert
     * @return a converted variable given the bases and skipping deletions
     */

    private static byte [] convertVariableGivenBases (List<BaseIndex> bases, List<Byte> variable) {
        byte [] variableArray = new byte[getReadLengthWithNoDeletions(bases)];
        int i = 0;
        Iterator<Byte> variableIterator = variable.iterator();
        for (BaseIndex baseIndex : bases) {
            byte count = variableIterator.next();
            if (baseIndex != BaseIndex.D)
                variableArray[i++] = count;
        }
        return variableArray;

    }

    /**
     * Shared functionality for all conversion utilities
     *
     * @param bases  the read bases
     * @return the length of the read with no deletions
     */
    private static int getReadLengthWithNoDeletions(List<BaseIndex> bases) {
        int readLength = bases.size();
        for (BaseIndex baseIndex : bases)
            if (baseIndex == BaseIndex.D)
                readLength--;
        return readLength;
    }


}
