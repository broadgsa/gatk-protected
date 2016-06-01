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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import htsjdk.samtools.GATKBin;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math.distribution.ExponentialDistribution;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
* Mock-up active region data used in testing.
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
public class ActiveRegionTestDataSet {

    private final byte[] referenceBytes;
    protected String reference;
    protected String[] haplotypeCigars;
    protected List<String> haplotypeStrings;
    protected String[] readCigars;
    protected byte[] bq;
    protected byte[] dq;
    protected byte[] iq;
    protected int kmerSize;
    private List<Haplotype> haplotypeList;
    private List<GATKSAMRecord> readList;
    private AssemblyResultSet assemblyResultSet;
    private Map<String,GATKSAMRecord> readBySequence;
    private String stringRepresentation;
    private List<List<Civar.ElementOffset>> readEventOffsetList;
    private GenomeLocParser genomeLocParser;

    /** Create a new active region data test set */
    public ActiveRegionTestDataSet(final int kmerSize, final String reference, final String[] haplotypes,
                                   final String[] readCigars, final byte[] bq, final byte[] dq, final byte[] iq) {
        this.reference = reference;
        this.referenceBytes = reference.getBytes();
        this.haplotypeCigars = haplotypes;
        this.readCigars = readCigars;
        this.bq = bq;
        this.dq = dq;
        this.iq = iq;
        this.kmerSize = kmerSize;
        this.genomeLocParser = new GenomeLocParser(ArtificialSAMUtils.createArtificialSamHeader(1,1,reference.length()).getSequenceDictionary());
    }

    public String getReference() {
        return reference;
    }

    public String toString() {
        if (stringRepresentation == null)
            return super.toString();
        else return stringRepresentation;
    }

    public AssemblyResultSet assemblyResultSet() {
        if (assemblyResultSet == null) {
            final ReadThreadingGraph rtg = new ReadThreadingGraph(kmerSize);
            rtg.addSequence("anonymous", this.getReference().getBytes(), true);
            for (final String haplotype : this.haplotypesStrings()) {
                rtg.addSequence("anonymous", haplotype.getBytes(), false);
            }
            rtg.buildGraphIfNecessary();
            if (rtg.hasCycles())
                throw new RuntimeException("there is cycles in the reference with kmer size " + kmerSize + ". Don't use this size for the benchmark or change the reference");

            List<Haplotype> haplotypeList = this.haplotypeList();

            assemblyResultSet = new AssemblyResultSet();
            final AssemblyResult ar = new AssemblyResult((haplotypeList.size() > 1 ?
                    AssemblyResult.Status.ASSEMBLED_SOME_VARIATION : AssemblyResult.Status.JUST_ASSEMBLED_REFERENCE),rtg.convertToSequenceGraph());
            ar.setThreadingGraph(rtg);

            for (final Haplotype h : haplotypeList)
                assemblyResultSet.add(h, ar);
        }
        return assemblyResultSet;
    }

    public List<String> haplotypesStrings() {
        if (haplotypeStrings != null) {
            return haplotypeStrings;
        }
        final List<String> result = new ArrayList<>(haplotypeCigars.length);
        String reference = this.reference;
        for (final String cigar : haplotypeCigars) {
            if (cigar.matches("^Civar:.*$")) {
                stringRepresentation = cigar.substring(6);
                result.addAll(expandAllCombinations(cigar.substring(6),reference));
            } else if (cigar.matches("^.*\\d+.*$")) {
                result.add(applyCigar(reference, cigar,0,true));
            } else {
                result.add(cigar);
            }
        }
        haplotypeStrings = result;
        return result;
    }

    private List<String> expandAllCombinations(final String cigarString, final String reference) {
        final Civar civar = Civar.fromCharSequence(cigarString);
        final List<Civar> unrolledCivars = civar.optionalizeAll().unroll();
        List<String> result = new ArrayList<>(unrolledCivars.size());
        for (final Civar c : unrolledCivars) {
            result.add(c.applyTo(reference));
        }
        return result;
    }

    private List<Haplotype> expandAllHaplotypeCombinations(final String civarString, final String reference) {
        final Civar civar = Civar.fromCharSequence(civarString);
        final List<Civar> unrolledCivars = civar.optionalizeAll().unroll();
        List<Haplotype> result = new ArrayList<>(unrolledCivars.size());
        for (final Civar c : unrolledCivars) {
            final String baseString = c.applyTo(reference);
            final Haplotype haplotype = new Haplotype(baseString.getBytes(),baseString.equals(reference));
            haplotype.setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,reference.length()));
            try {
            haplotype.setCigar(c.toCigar(reference.length()));
            } catch (final RuntimeException ex) {
                c.applyTo(reference);
                c.toCigar(reference.length());
                throw new RuntimeException("" + c + " " + ex.getMessage(),ex);
            }
            result.add(haplotype);
        }
        return result;
    }


    public List<Haplotype> haplotypeList() {
        if (haplotypeList == null) {

          final List<Haplotype> result = new ArrayList<>(haplotypeCigars.length);
          final String reference = this.reference;
          for (final String cigar : haplotypeCigars) {
              if (cigar.matches("^Civar:.*$")) {
                  stringRepresentation = cigar.substring(6);
                  result.addAll(expandAllHaplotypeCombinations(cigar.substring(6), reference));
              } else if (cigar.matches("^.*\\d+.*$")) {
                  result.add(cigarToHaplotype(reference, cigar, 0, true));
              } else {
                  final Haplotype h = new Haplotype(cigar.getBytes());
                  h.setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,reference.length()));
                  result.add(h);
              }
          }
          haplotypeList = result;
        }
        return haplotypeList;
    }


    protected SAMSequenceDictionary artificialSAMSequenceDictionary() {
        return new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("00",reference.length())));
    }

    protected SAMFileHeader artificialSAMFileHeader() {
        return ArtificialSAMUtils.createArtificialSamHeader(artificialSAMSequenceDictionary());
    }

    public List<GATKSAMRecord> readList() {
        if (readList == null) {
            final SAMFileHeader header = artificialSAMFileHeader();
            readList = new ArrayList<>(readCigars.length);
            final List<String> haplotypes = haplotypesStrings();
            int count = 0;
            for (final String descr : readCigars) {
                String sequence;
                if (descr.matches("^\\d+:\\d+:.+$")) {
                    final String[] parts = descr.split(":");
                    int allele = Integer.valueOf(parts[0]);
                    int offset = Integer.valueOf(parts[1]);
                    final String cigar = parts[2];
                    final String base = allele == 0 ? reference : haplotypes.get(allele - 1);
                    sequence = applyCigar(base, cigar, offset, false);
                    final GATKSAMRecord samRecord = ArtificialSAMUtils.createArtificialRead(header, "read_" + count, 0, 1, sequence.getBytes(), Arrays.copyOf(bq, sequence.length()));
                    readList.add(new MyGATKSAMRecord(samRecord));
                } else if (descr.matches("^\\*:\\d+:\\d+$")) {
                    int readCount = Integer.valueOf(descr.split(":")[1]);
                    int readLength = Integer.valueOf(descr.split(":")[2]);
                    readList.addAll(generateSamRecords(haplotypes, readCount, readLength, header, count));
                } else {
                    sequence = descr;
                    final GATKSAMRecord samRecord = ArtificialSAMUtils.createArtificialRead(header, "read_" + count, 0, 1, sequence.getBytes(), Arrays.copyOf(bq, sequence.length()));
                    readList.add(new MyGATKSAMRecord(samRecord));
                }
                count = readList.size();
            }
        }
        return readList;
    }

    public List<List<Civar.ElementOffset>> readEventOffsetList() {
        if (haplotypeCigars.length != 1 || !haplotypeCigars[0].startsWith("Civar:"))
            throw new UnsupportedOperationException();
        if (readEventOffsetList == null) {
            final Civar civar = Civar.fromCharSequence(haplotypeCigars[0].substring(6));
            final List<Civar> unrolledCivars = civar.optionalizeAll().unroll();

            readEventOffsetList = new ArrayList<>(readCigars.length);
            int count = 0;
            for (final String descr : readCigars) {
                if (descr.matches("^\\d+:\\d+:.+$")) {
                    throw new UnsupportedOperationException();
                } else if (descr.matches("^\\*:\\d+:\\d+$")) {
                    int readCount = Integer.valueOf(descr.split(":")[1]);
                    int readLength = Integer.valueOf(descr.split(":")[2]);
                    readEventOffsetList.addAll(generateElementOffsetRecords(haplotypesStrings(), unrolledCivars, readCount, readLength, count));
                } else {
                    throw new UnsupportedOperationException();
                }
                count = readEventOffsetList.size();
            }
            readEventOffsetList = Collections.unmodifiableList(readEventOffsetList);
        }
        return readEventOffsetList;
    }




    @SuppressWarnings("unused")
    public String cigarToSequence(final String cigar) {
            String reference = this.reference;
            return applyCigar(reference, cigar,0,true);
    }

    @SuppressWarnings("unused")
    public GATKSAMRecord readFromString(final String readSequence) {
        if (readBySequence == null) {
            final List<GATKSAMRecord> readList = readList();
            readBySequence = new HashMap<>(readList.size());
            for (final GATKSAMRecord r : readList)
                readBySequence.put(r.getReadString(),r);
        }
        return readBySequence.get(readSequence);
    }

    public List<Civar> unrolledCivars() {
        if (haplotypeCigars.length != 1 || !haplotypeCigars[0].startsWith("Civar:"))
            throw new UnsupportedOperationException();
        final Civar civar = Civar.fromCharSequence(haplotypeCigars[0].substring(6));
        return civar.optionalizeAll().unroll();
    }

    public void introduceErrors(final Random rnd) {
        final List<GATKSAMRecord> reads = readList();
        final ArrayList<GATKSAMRecord> result = new ArrayList<>(reads.size());
        for (final GATKSAMRecord read : reads) {
            result.add(new MyGATKSAMRecord(read,rnd));
        }
        readList = result;
    }

    private class MyGATKSAMRecord extends GATKSAMRecord {
            protected MyGATKSAMRecord(final GATKSAMRecord r) {
                super(r);
                this.setMappingQuality(100);
                GATKBin.setReadIndexingBin(this, -1);
            }

        ExponentialDistribution indelLengthDist = MathUtils.exponentialDistribution(1.0 / 0.9);

        public MyGATKSAMRecord(final GATKSAMRecord r, final Random rnd) {
            super(r);
            this.setMappingQuality(100);
            // setting read indexing bin last

            final byte[] bases = new byte[r.getReadBases().length];

            final byte[] readBases = r.getReadBases();
            final byte[] bq = r.getBaseQualities();
            final byte[] iq = r.getBaseInsertionQualities();
            final byte[] dq = r.getBaseDeletionQualities();
            int refOffset = r.getAlignmentStart() - 1;
            int readOffset = 0;
            for (int i = 0; i < r.getReadBases().length;) {
                double p = rnd.nextDouble();
                double iqp = QualityUtils.qualToErrorProb(iq[i]);
                if (p < iqp) { // insertion
                    final int length = Math.min(generateIndelLength(rnd),r.getReadBases().length - i);
                    final int refStart = rnd.nextInt(reference.length() - length);
                    System.arraycopy(referenceBytes,refStart,bases,i,length);
                    i += length;
                    continue;
                }
                p -= iqp;
                double dqp = QualityUtils.qualToErrorProb(dq[i]);
                if (p < dqp) {
                    final int length = generateIndelLength(rnd);
                    refOffset += length;
                    refOffset = refOffset % referenceBytes.length;
                    readOffset += length;
                    continue;
                }
                p -= dqp;
                double bqp = QualityUtils.qualToErrorProb(bq[i]);
                byte b = readOffset < readBases.length ? readBases[readOffset] : referenceBytes[refOffset];
                byte nb;
                if (p < bqp) {
                   switch (b) {
                       case 'A': nb = 'C'; break;
                       case 'T': nb = 'A'; break;
                       case 'C': nb = 'G'; break;
                       case 'G': nb = 'B'; break;
                       default: nb = 'A';
                   }
                } else
                    nb = b;

                bases[i++] = nb;
                refOffset++;
                refOffset = refOffset % referenceBytes.length;
                readOffset++;
            }
            this.setReadBases(bases);
            this.setBaseQualities(r.getBaseQualities());
            this.setReadName(r.getReadName());

            GATKBin.setReadIndexingBin(this, -1);
        }

        private int generateIndelLength(final Random rnd) {
            final int length;
            try {
                length = (int) Math.round(indelLengthDist.inverseCumulativeProbability(rnd.nextDouble()) + 1);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
            return length;
        }

        @Override
            public byte[] getBaseDeletionQualities() {
                return Arrays.copyOf(dq,getReadLength());
            }

            @Override
            public byte[] getBaseInsertionQualities() {
                return Arrays.copyOf(iq,getReadLength());
            }

            @Override
            public int getMappingQuality() {
                return 100;
            }

            @Override
            public int hashCode() {
                return getReadName().hashCode();
            }

            @Override
            public boolean equals(Object o) {
                if (o instanceof GATKSAMRecord) {
                    return getReadName().equals(((GATKSAMRecord)o).getReadName());
                } else {
                    return false;
                }
            }

            public String toString() {
                return super.toString() + " " + this.getReadString();
            }
    }


    public List<String> readStrings() {
        final List<String> result = new ArrayList<>(readCigars.length);
        final List<String> haplotypes = haplotypesStrings();
        for (final String descr : readCigars) {
            String sequence;
            if (descr.matches("^\\d+:\\d+:.+$")) {
                final String[] parts = descr.split(":");
                int allele = Integer.valueOf(parts[0]);
                int offset = Integer.valueOf(parts[1]);
                final String cigar = parts[2];
                final String base = allele == 0 ? reference : haplotypes.get(allele - 1);
                sequence = applyCigar(base, cigar, offset, false);
                result.add(sequence);
            } else if (descr.matches("\\*:^\\d+:\\d+")) {
                int readCount = Integer.valueOf(descr.split(":")[1]);
                int readLength = Integer.valueOf(descr.split(":")[2]);
                result.addAll(generateReads(haplotypes, readCount, readLength));
            } else {
                sequence = descr;
                result.add(sequence);
            }
        }
        return result;
    }

    private List<String> generateReads(final List<String> haplotypes, final int readCount, final int readLength) {
        final List<String> result = new ArrayList<>(readCount);
        for (int i = 0; i < readCount; i++) {
            int hi = i % haplotypes.size();
            final String h = haplotypes.get(hi);
            int offset = i % h.length() - readLength;
            result.add(h.substring(offset,offset + readLength));
        }
        return result;
    }

    private List<MyGATKSAMRecord> generateSamRecords(final List<String> haplotypes, final int readCount, final int readLength, final SAMFileHeader header, final int idStart) {
        int id = idStart;
        final List<MyGATKSAMRecord> result = new ArrayList<>(readCount);
        for (int i = 0; i < readCount; i++) {
            int hi = i % haplotypes.size();
            final String h = haplotypes.get(hi);
            int offset = h.length() <= readLength ? 0 : i % (h.length() - readLength);
            int to = Math.min(h.length(),offset + readLength);
            byte[] bases = h.substring(offset,to).getBytes();
            byte[] quals = Arrays.copyOf(bq,to - offset);
            final GATKSAMRecord samRecord = ArtificialSAMUtils.createArtificialRead(header,"read_" + id++,0,offset + 1,bases, quals);
            result.add(new MyGATKSAMRecord(samRecord));
        }
        return result;
    }


    private List<List<Civar.ElementOffset>> generateElementOffsetRecords(final List<String> haplotypes, final List<Civar> unrolledCivars, final int readCount, final int readLength, final int count) {

        final List<List<Civar.ElementOffset>> result = new ArrayList<>(readCount);
        for (int i = 0; i < readCount; i++) {
            int hi = i % unrolledCivars.size();
            final Civar c = unrolledCivars.get(hi);
            final String h = haplotypes.get(hi);
            int offset = h.length() <= readLength ? 0 : i % (h.length() - readLength);
            int to = Math.min(h.length(),offset + readLength);
            result.add(c.eventOffsets(reference,offset,to));
        }
        return result;
    }

    private static final Pattern cigarPattern = Pattern.compile("(\\d+)([=A-Z])");


    private Haplotype cigarToHaplotype(final String reference, final String cigar, final int offset, final boolean global) {
        final String sequence = applyCigar(reference,cigar,offset,global);
        final Haplotype haplotype = new Haplotype(sequence.getBytes(),reference.equals(sequence));
        haplotype.setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,reference.length()));
        haplotype.setCigar(Civar.fromCharSequence(cigar).toCigar(reference.length()));
        return haplotype;
    }

    private String applyCigar(final String reference, final String cigar, final int offset, final boolean global) {
        final Matcher pm = cigarPattern.matcher(cigar);
        StringBuffer sb = new StringBuffer();
        int index = offset;
        while (pm.find()) {
            int length = Integer.valueOf(pm.group(1));
            char operator = pm.group(2).charAt(0);
            switch (operator) {
                case '=' :
                    try {
                      sb.append(reference.substring(index, index + length));
                    } catch (Exception e) {
                      throw new RuntimeException(" " + index + " " + (index + length) + " " + reference.length() + " " + cigar,e);
                    }
                    index += length; break;
                case 'D' :
                    index += length; break;
                case 'I' :
                    String insert = cigar.substring(pm.end(),pm.end() + length).toUpperCase();
                    sb.append(insert); break;
                case 'V' :
                    sb.append(transversionV(reference.charAt(index))); index++; break;
                case 'W' :
                        sb.append(transversionW(reference.charAt(index))); index++; break;
                case 'T' :
                    sb.append(transition(reference.charAt(index))); index++; break;
                default:
                    throw new UnsupportedOperationException("cigar operator " + operator + " not supported.");
            }
        }
        if (global && index != reference.length()) {
            throw new RuntimeException(" haplotype cigar does not explain reference length (" + index + " != " + reference.length() + ") on cigar " + cigar);
        } else if (index > reference.length()) {
            throw new RuntimeException(" index beyond end ");
        }
        return sb.toString();
    }

    protected int kmerSize() {
        return kmerSize;
    }

    private char transversionV(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'C';
            case 'G': return 'T';
            case 'C': return 'A';
            case 'T': return 'G';
            default:
                return c;
        }

    }

    private char transversionW(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'T';
            case 'G': return 'C';
            case 'T': return 'A';
            case 'C': return 'G';
            default:
                return c;
        }

    }

    private char transition(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'G';
            case 'G': return 'A';
            case 'T': return 'C';
            case 'C': return 'T';
            default:
                return c;
        }

    }
}
