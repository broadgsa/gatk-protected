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
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.*;


public class ArtificialReadPileupTestProvider {
    final String refBases = "ACAGAGCTGACCCTCCCTCCCCTCTCCCAGTGCAACAGCACGGGCGGCGACTGCTTTTACCGAGGCTACACGTCAGGCGTGGCGGCTGTCCAGGACTGGTACCACTTCCACTATGTGGATCTCTGCTGAGGACCAGGAAAGCCAGCACCCGCAGAGACTCTTCCCCAGTGCTCCATACGATCACCATTCTCTGCAGAAGGTCAGACGTCACTGGTGGCCCCCCAGCCTCCTCAGCAGGGAAGGATACTGTCCCGCAGATGAGATGAGCGAGAGCCGCCAGACCCACGTGACGCTGCACGACATCGACCCTCAGGCCTTGGACCAGCTGGTGCAGTTTGCCTACACGGCTGAGATTGTGGTGGGCGAGGGC";
    final int contigStart = 1;
    final int contigStop = refBases.length();
    final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, contigStop - contigStart + 1);
//    final GATKSAMReadGroupRecord artificialGATKRG = new GATKSAMReadGroupRecord("synthetic");
    final String artificialContig = "chr1";
  //  final int artificialContigIndex = 0;
    final String artificialReadName = "synth";
    final int artificialRefStart = 1;
    final int artificialMappingQuality = 60;
    Map<String, GATKSAMReadGroupRecord> sample2RG = new HashMap<String, GATKSAMReadGroupRecord>();
    List<SAMReadGroupRecord> sampleRGs;
    List<String> sampleNames = new ArrayList<String>();
    private String sampleName(int i) { return sampleNames.get(i); }
    private GATKSAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }
    public final int locStart = 105; // start position where we desire artificial variant
    private final int readLength = 10; // desired read length in pileup
    public final int readOffset = 4;
    private final int readStart = locStart - readOffset;
    public final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    public final GenomeLoc loc = genomeLocParser.createGenomeLoc(artificialContig,locStart,locStart);
    public final GenomeLoc window = genomeLocParser.createGenomeLoc(artificialContig,locStart-100,locStart+100);
    public final String windowBases = refBases.substring(locStart-100-1,locStart+100);
    public final ReferenceContext referenceContext = new ReferenceContext(genomeLocParser,loc,window,windowBases.getBytes());

    byte BASE_QUAL = 50;

    public ArtificialReadPileupTestProvider(final int numSamples, final String SAMPLE_PREFIX) {
        sampleRGs = new ArrayList<SAMReadGroupRecord>();

        for ( int i = 0; i < numSamples; i++ ) {
            sampleNames.add(String.format("%s%04d", SAMPLE_PREFIX, i));
            GATKSAMReadGroupRecord rg = createRG(sampleName(i));
            sampleRGs.add(rg);
            sample2RG.put(sampleName(i), rg);
        }

    }

    public ArtificialReadPileupTestProvider(final int numSamples, final String SAMPLE_PREFIX, final byte q) {
        this(numSamples,SAMPLE_PREFIX);
        BASE_QUAL = q;
    }
    public List<String> getSampleNames() {
        return sampleNames;
    }
    public byte getRefByte() {
        return referenceContext.getBase();
    }

    public ReferenceContext getReferenceContext()   { return referenceContext;}
    public GenomeLocParser getGenomeLocParser()     { return genomeLocParser; }

    public Map<String,AlignmentContext> getAlignmentContextFromAlleles(int eventLength, String altBases, int[] numReadsPerAllele) {
        return getAlignmentContextFromAlleles(eventLength, altBases, numReadsPerAllele, false, BASE_QUAL);
    }
    public Map<String,AlignmentContext> getAlignmentContextFromAlleles(final int eventLength,
                                                                       final String altBases,
                                                                       final int[] numReadsPerAllele,
                                                                       final boolean addBaseErrors,
                                                                       final int phredScaledBaseErrorRate) {
        final String refChar = new String(new byte[]{referenceContext.getBase()});

        String refAllele, altAllele;
        if (eventLength == 0)  {
            // SNP case
            refAllele = refChar;
            altAllele = altBases.substring(0,1);

        } else if (eventLength>0){
            // insertion
            refAllele = refChar;
            altAllele = refChar+altBases/*.substring(0,eventLength)*/;
        }
        else {
            // deletion
            refAllele = new String(referenceContext.getForwardBases()).substring(0,Math.abs(eventLength)+1);
            altAllele = refChar;
        }

        Map<String,AlignmentContext> contexts = new HashMap<String,AlignmentContext>();

        for (String sample: sampleNames) {
            AlignmentContext context = new AlignmentContext(loc, generateRBPForVariant(loc, refAllele, altAllele, altBases, numReadsPerAllele, sample, addBaseErrors, phredScaledBaseErrorRate));
            contexts.put(sample,context);

        }

        return contexts;
    }

    private GATKSAMReadGroupRecord createRG(String name) {
        GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(name);
        rg.setPlatform("ILLUMINA");
        rg.setSample(name);
        return rg;
    }

    private ReadBackedPileup generateRBPForVariant( GenomeLoc loc, String refAllele, String altAllele, String altBases,
                                                    int[] numReadsPerAllele, String sample, boolean addErrors, int phredScaledErrorRate) {
        List<PileupElement> pileupElements = new ArrayList<PileupElement>();
        final int refAlleleLength = refAllele.length();

        pileupElements.addAll(createPileupElements(refAllele, loc, numReadsPerAllele[0], sample, readStart, altBases, addErrors, phredScaledErrorRate, refAlleleLength, true));
        pileupElements.addAll(createPileupElements(altAllele, loc, numReadsPerAllele[1], sample, readStart, altBases, addErrors, phredScaledErrorRate, refAlleleLength, false));
        return new ReadBackedPileupImpl(loc,pileupElements);
    }

    private List<PileupElement> createPileupElements(String allele, GenomeLoc loc, int numReadsPerAllele, String sample, int readStart, String altBases, boolean addErrors, int phredScaledErrorRate, int refAlleleLength, boolean isReference) {

        int alleleLength = allele.length();
        List<PileupElement> pileupElements = new ArrayList<PileupElement>();

        int readCounter = 0;
        for ( int d = 0; d < numReadsPerAllele; d++ ) {
            byte[] readBases = trueHaplotype(allele, refAlleleLength, readLength);
            if (addErrors)
                addBaseErrors(readBases, phredScaledErrorRate);

            byte[] readQuals = new byte[readBases.length];
            Arrays.fill(readQuals, (byte)phredScaledErrorRate);

            GATKSAMRecord read = new GATKSAMRecord(header);
            read.setBaseQualities(readQuals);
            read.setReadBases(readBases);
            read.setReadName(artificialReadName+readCounter++);

            boolean isBeforeDeletion = alleleLength<refAlleleLength;
            boolean isBeforeInsertion = alleleLength>refAlleleLength;

            int eventLength = alleleLength - refAlleleLength;
            if (isReference)
                read.setCigarString(readBases.length + "M");
            else {
                if (isBeforeDeletion || isBeforeInsertion)
                    read.setCigarString((readOffset+1)+"M"+ Math.abs(eventLength) + (isBeforeDeletion?"D":"I") +
                            (readBases.length-readOffset)+"M");
                else // SNP case
                    read.setCigarString(readBases.length+"M");
            }

            read.setReadPairedFlag(false);
            read.setAlignmentStart(readStart);
            read.setMappingQuality(artificialMappingQuality);
            read.setReferenceName(loc.getContig());
            read.setReadNegativeStrandFlag(false);
            read.setReadGroup(sampleRG(sample));

            pileupElements.add(LocusIteratorByState.createPileupForReadAndOffset(read, readOffset));
        }

        return pileupElements;
    }

    /**
     * Create haplotype with desired allele and reference context
     * @param allele                             Desired allele string
     * @param refAlleleLength                    Length of reference allele.
     * @param desiredLength                      Desired haplotype length
     * @return                                   String with haplotype formed by (prefix)+allele bases + postfix
     */
    private byte[] trueHaplotype(final String allele, final int refAlleleLength, final int desiredLength) {
        // create haplotype based on a particular allele
        final int startIdx= locStart - readOffset-1;

        final String prefix = refBases.substring(startIdx, locStart-1);
        final String postfix = refBases.substring(locStart+refAlleleLength-1,startIdx + desiredLength);

        return (prefix+allele+postfix).getBytes();
    }

    private void addBaseErrors(final byte[] readBases, final int phredScaledErrorRate) {
        double errorProbability = QualityUtils.qualToErrorProb((byte)phredScaledErrorRate);

        for (int k=0; k < readBases.length; k++) {
            if (Utils.getRandomGenerator().nextDouble() < errorProbability) {
                // random offset
                int offset = BaseUtils.simpleBaseToBaseIndex(readBases[k]);          //0..3
                offset += (Utils.getRandomGenerator().nextInt(3)+1);  // adds 1,2 or 3
                offset %= 4;
                readBases[k] = BaseUtils.baseIndexToSimpleBase(offset);

            }

        }

    }
}
