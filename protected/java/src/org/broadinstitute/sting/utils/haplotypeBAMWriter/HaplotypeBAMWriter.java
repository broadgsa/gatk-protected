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

package org.broadinstitute.sting.utils.haplotypeBAMWriter;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs.Path;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * A BAMWriter that aligns reads to haplotypes and emits their best alignments to a BAM file
 *
 * User: depristo
 * Date: 2/22/13
 * Time: 2:59 PM
 */
public abstract class HaplotypeBAMWriter {
    /**
     * Allows us to write out unique names for our synthetic haplotype reads
     */
    private long uniqueNameCounter = 1;

    protected final static String READ_GROUP_ID = "ArtificialHaplotype";
    protected final static String HAPLOTYPE_TAG = "HC";

    final SAMFileWriter bamWriter;
    final SAMFileHeader bamHeader;

    /**
     * Possible modes for writing haplotypes to BAMs
     */
    public static enum Type {
        /**
         * A mode that's for method developers.  Writes out all of the possible
         * haplotypes considered, as well as reads aligned to each
         */
        ALL_POSSIBLE_HAPLOTYPES,

        /**
         * A mode for users.  Writes out the reads aligned only to the called
         * haplotypes.  Useful to understand why the caller is calling what it is
         */
        CALLED_HAPLOTYPES
    }

    /**
     * Create a new HaplotypeBAMWriter of type writing SAMRecords to writer
     *
     * @param type the type of the writer we want to create
     * @param stingSAMWriter the destination, must not be null
     * @param header the header of the input BAMs used to make calls, must not be null
     * @return a new HaplotypeBAMWriter
     */
    public static HaplotypeBAMWriter create(final Type type, final StingSAMFileWriter stingSAMWriter, final SAMFileHeader header) {
        if ( header == null ) throw new IllegalArgumentException("header cannot be null");
        if ( stingSAMWriter == null ) throw new IllegalArgumentException("writer cannot be null");
        if ( type == null ) throw new IllegalArgumentException("type cannot be null");

        // prepare the bam header
        final SAMFileHeader bamHeader = new SAMFileHeader();
        bamHeader.setSequenceDictionary(header.getSequenceDictionary());
        bamHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        // include the original read groups plus a new artificial one for the haplotypes
        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>(header.getReadGroups());
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(READ_GROUP_ID);
        rg.setSample("HC");
        rg.setSequencingCenter("BI");
        readGroups.add(rg);
        bamHeader.setReadGroups(readGroups);

        // TODO -- this will be a performance problem at high-scale
        stingSAMWriter.setPresorted(false);
        stingSAMWriter.writeHeader(bamHeader);
        return create(type, stingSAMWriter);
    }

    /**
     * Create a new HaplotypeBAMWriter of type writing SAMRecords to writer
     *
     * Note that writer must have its presorted bit set to false, as reads
     * may come in out of order during writing
     *
     * @param type the type of the writer we want to create
     * @param writer the destination, must not be null
     * @return a new HaplotypeBAMWriter
     */
    public static HaplotypeBAMWriter create(final Type type, final SAMFileWriter writer) {
        if ( writer == null ) throw new IllegalArgumentException("writer cannot be null");
        if ( type == null ) throw new IllegalArgumentException("type cannot be null");

        switch ( type ) {
            case ALL_POSSIBLE_HAPLOTYPES: return new AllHaplotypeBAMWriter(writer);
            case CALLED_HAPLOTYPES: return new CalledHaplotypeBAMWriter(writer);
            default: throw new IllegalArgumentException("Unknown type " + type);
        }
    }

    /**
     * Create a new HaplotypeBAMWriter writing its output to bamWriter
     *
     * Assumes that the header has been fully initialized with a single
     * read group READ_GROUP_ID
     *
     * @param bamWriter our output destination
     */
    protected HaplotypeBAMWriter(SAMFileWriter bamWriter) {
        this.bamWriter = bamWriter;
        this.bamHeader = bamWriter.getFileHeader();
    }

    /**
     * Write out a BAM representing for the haplotype caller at this site
     *
     * @param haplotypes a list of all possible haplotypes at this loc
     * @param paddedReferenceLoc the span of the based reference here
     * @param bestHaplotypes a list of the best (a subset of all) haplotypes that actually went forward into genotyping
     * @param calledHaplotypes a list of the haplotypes at where actually called as non-reference
     * @param stratifiedReadMap a map from sample -> likelihoods for each read for each of the best haplotypes
     */
    public abstract void writeReadsAlignedToHaplotypes(final List<Haplotype> haplotypes,
                                                       final GenomeLoc paddedReferenceLoc,
                                                       final List<Haplotype> bestHaplotypes,
                                                       final Set<Haplotype> calledHaplotypes,
                                                       final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap);

    /**
     * Write out read aligned to haplotype to the BAM file
     *
     * Aligns reads the haplotype, and then projects this alignment of read -> hap onto the reference
     * via the alignment of haplotype (via its getCigar) method.
     *
     * @param originalRead the read we want to write aligned to the reference genome
     * @param haplotype the haplotype that the read should be aligned to, before aligning to the reference
     * @param referenceStart the start of the reference that haplotype is aligned to.  Provides global coordinate frame.
     * @param isInformative true if the read is differentially informative for one of the haplotypes
     */
    protected void writeReadAgainstHaplotype(final GATKSAMRecord originalRead,
                                             final Haplotype haplotype,
                                             final int referenceStart,
                                             final boolean isInformative) {
        final GATKSAMRecord alignedToRef = createReadAlignedToRef(originalRead, haplotype, referenceStart, isInformative);
        if ( alignedToRef != null )
            bamWriter.addAlignment(alignedToRef);
    }

    /**
     * Aligns reads the haplotype, and then projects this alignment of read -> hap onto the reference
     * via the alignment of haplotype (via its getCigar) method.
     *
     * @param originalRead the read we want to write aligned to the reference genome
     * @param haplotype the haplotype that the read should be aligned to, before aligning to the reference
     * @param referenceStart the start of the reference that haplotype is aligned to.  Provides global coordinate frame.
     * @param isInformative true if the read is differentially informative for one of the haplotypes
     * @return a GATKSAMRecord aligned to reference, or null if no meaningful alignment is possible
     */
    protected GATKSAMRecord createReadAlignedToRef(final GATKSAMRecord originalRead,
                                                   final Haplotype haplotype,
                                                   final int referenceStart,
                                                   final boolean isInformative) {
        if ( originalRead == null ) throw new IllegalArgumentException("originalRead cannot be null");
        if ( haplotype == null ) throw new IllegalArgumentException("haplotype cannot be null");
        if ( haplotype.getCigar() == null ) throw new IllegalArgumentException("Haplotype cigar not set " + haplotype);
        if ( referenceStart < 1 ) throw new IllegalArgumentException("reference start much be >= 1 but got " + referenceStart);

        try {
            // compute the smith-waterman alignment of read -> haplotype
            final SWPairwiseAlignment swPairwiseAlignment = new SWPairwiseAlignment(haplotype.getBases(), originalRead.getReadBases(), Path.NEW_SW_PARAMETERS);
            //swPairwiseAlignment.printAlignment(haplotype.getBases(), originalRead.getReadBases());
            if ( swPairwiseAlignment.getAlignmentStart2wrt1() == -1 )
                // sw can fail (reasons not clear) so if it happens just don't write the read
                return null;
            final Cigar swCigar = AlignmentUtils.consolidateCigar(swPairwiseAlignment.getCigar());

            // since we're modifying the read we need to clone it
            final GATKSAMRecord read = (GATKSAMRecord)originalRead.clone();

            addHaplotypeTag(read, haplotype);

            // uninformative reads are set to zero mapping quality to enhance visualization
            if ( !isInformative )
                read.setMappingQuality(0);

            // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
            final Cigar extendedHaplotypeCigar = haplotype.getConsolidatedPaddedCigar(1000);
            final int readStartOnHaplotype = AlignmentUtils.calcFirstBaseMatchingReferenceInCigar(extendedHaplotypeCigar, swPairwiseAlignment.getAlignmentStart2wrt1());
            final int readStartOnReference = referenceStart + haplotype.getAlignmentStartHapwrtRef() + readStartOnHaplotype;
            read.setAlignmentStart(readStartOnReference);

            // compute the read -> ref alignment by mapping read -> hap -> ref from the
            // SW of read -> hap mapped through the given by hap -> ref
            final Cigar haplotypeToRef = AlignmentUtils.trimCigarByBases(extendedHaplotypeCigar, swPairwiseAlignment.getAlignmentStart2wrt1(), extendedHaplotypeCigar.getReadLength() - 1);
            final Cigar readToRefCigarRaw = AlignmentUtils.applyCigarToCigar(swCigar, haplotypeToRef);
            final Cigar readToRefCigarClean = AlignmentUtils.cleanUpCigar(readToRefCigarRaw);
            final Cigar readToRefCigar = AlignmentUtils.leftAlignIndel(readToRefCigarClean, haplotype.getBases(),
                    originalRead.getReadBases(), swPairwiseAlignment.getAlignmentStart2wrt1(), 0, true);

            read.setCigar(readToRefCigar);

            if ( readToRefCigar.getReadLength() != read.getReadLength() )
                throw new IllegalStateException("Cigar " + readToRefCigar + " with read length " + readToRefCigar.getReadLength()
                        + " != read length " + read.getReadLength() + " for read " + read.format() + "\nhapToRef " + haplotypeToRef + " length " + haplotypeToRef.getReadLength() + "/" + haplotypeToRef.getReferenceLength()
                        + "\nreadToHap " + swCigar + " length " + swCigar.getReadLength() + "/" + swCigar.getReferenceLength());

            return read;
        } catch ( CloneNotSupportedException e ) {
            throw new IllegalStateException("GATKSAMRecords should support clone but this one does not " + originalRead);
        }
    }

    /**
     * Add a haplotype tag to the read based on haplotype
     *
     * @param read the read to add the tag to
     * @param haplotype the haplotype that gives rises to read
     */
    private void addHaplotypeTag(final GATKSAMRecord read, final Haplotype haplotype) {
        // add a tag to the read that indicates which haplotype it best aligned to.  It's a uniquish integer
        read.setAttribute(HAPLOTYPE_TAG, haplotype.hashCode());
    }

    /**
     * Write out haplotypes as reads to the BAM, marking specifically those that are among the best haplotypes
     *
     * @param haplotypes a collection of haplotypes to write to the BAM
     * @param bestHaplotypes a subset of haplotypes that contains those that are best "either good or called"
     * @param paddedReferenceLoc the genome loc of the padded reference
     */
    protected void writeHaplotypesAsReads(final Collection<Haplotype> haplotypes,
                                          final Set<Haplotype> bestHaplotypes,
                                          final GenomeLoc paddedReferenceLoc) {
        for ( final Haplotype haplotype : haplotypes )
            writeHaplotype(haplotype, paddedReferenceLoc, bestHaplotypes.contains(haplotype));
    }

    /**
     * Write out a representation of this haplotype as a read
     *
     * @param haplotype a haplotype to write out.  Cannot be null
     * @param paddedRefLoc the reference location.  Cannot be null
     * @param isAmongBestHaplotypes true if among the best haplotypes, false if it was just one possible but not so good
     */
    private void writeHaplotype(final Haplotype haplotype,
                                final GenomeLoc paddedRefLoc,
                                final boolean isAmongBestHaplotypes) {
        final GATKSAMRecord record = new GATKSAMRecord(bamHeader);
        record.setReadBases(haplotype.getBases());
        record.setAlignmentStart(paddedRefLoc.getStart() + haplotype.getAlignmentStartHapwrtRef());
        record.setBaseQualities(Utils.dupBytes((byte) '!', haplotype.getBases().length));
        record.setCigar(AlignmentUtils.consolidateCigar(haplotype.getCigar()));
        record.setMappingQuality(isAmongBestHaplotypes ? 60 : 0);
        record.setReadName("HC" + uniqueNameCounter++);
        addHaplotypeTag(record, haplotype);
        record.setReadUnmappedFlag(false);
        record.setReferenceIndex(paddedRefLoc.getContigIndex());
        record.setAttribute(SAMTag.RG.toString(), READ_GROUP_ID);
        record.setFlags(16);
        bamWriter.addAlignment(record);
    }
}