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

package org.broadinstitute.gatk.utils.haplotypeBAMWriter;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class HaplotypeBAMWriterUnitTest extends BaseTest {
    private final static boolean DEBUG = false;
    final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

    private GATKSAMRecord makeRead(final String baseString) {
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        final byte[] bases = baseString.getBytes();
        read.setReadBases(bases.clone());
        read.setBaseQualities(Utils.dupBytes((byte)30, read.getReadLength()));
        return read;
    }

    private Haplotype makeHaplotype(final String bases, final String cigar) {
        final Haplotype hap = new Haplotype(bases.getBytes());
        hap.setCigar(TextCigarCodec.decode(cigar));
        return hap;
    }

    private static class MockDestination extends ReadDestination {
        private final static SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();

        private MockDestination() {
            super(header, "foo");
        }

        @Override
        public void add(GATKSAMRecord read) {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    @Test
    public void testCreate() throws Exception {
        final MockDestination writer = new MockDestination();
        Assert.assertTrue(HaplotypeBAMWriter.create(HaplotypeBAMWriter.Type.CALLED_HAPLOTYPES, writer) instanceof CalledHaplotypeBAMWriter);
        Assert.assertTrue(HaplotypeBAMWriter.create(HaplotypeBAMWriter.Type.ALL_POSSIBLE_HAPLOTYPES, writer) instanceof AllHaplotypeBAMWriter);
    }


    //////////////////////////////////////////
    // Test HaplotypeBAMWriter.createReadAlignedToRef() //
    //////////////////////////////////////////

    @DataProvider(name = "ReadAlignedToRefData")
    public Object[][] makeReadAlignedToRefData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String hapBases = "ACTGAAGGTTCC";
        final Haplotype allM = makeHaplotype(hapBases, hapBases.length() + "M");

        // make sure we get back a cigar of the right length
        for ( int i = -1; i < hapBases.length(); i++ ) {
            final GATKSAMRecord read = makeRead(hapBases);
            if ( i != -1 ) read.getReadBases()[i] = (byte)'A';
            tests.add(new Object[]{read, allM, 10, 10, allM.getCigar().toString()});
        }

        // make sure insertions at the front are correctly handled
        for ( int padFront = 1; padFront < 10; padFront++ ) {
            final GATKSAMRecord read = makeRead(Utils.dupString("N", padFront) + hapBases);
            tests.add(new Object[]{read, allM, 10, 10, padFront + "I" + allM.getCigar().toString()});
        }

        // make sure insertions at the back are correctly handled
        for ( int padBack = 1; padBack < 10; padBack++ ) {
            final GATKSAMRecord read = makeRead(hapBases + Utils.dupString("N", padBack));
            tests.add(new Object[]{read, allM, 10, 10, allM.getCigar().toString() + padBack + "I"});
        }

        // make sure refStart and hapStart are respected
        for ( int refStart = 1; refStart < 10; refStart++ ) {
            for ( int hapStart = refStart; hapStart < 10 + refStart; hapStart++ ) {
                final Haplotype hap = new Haplotype(allM.getBases());
                hap.setCigar(allM.getCigar());
                hap.setAlignmentStartHapwrtRef(hapStart);

                final GATKSAMRecord read = makeRead(new String(hap.getBases()));
                tests.add(new Object[]{read, hap, refStart, refStart + hapStart, allM.getCigar().toString()});
            }
        }

        // example case of bad alignment because SW doesn't necessarily left-align indels
        {
            final String hap = "ACTGTGGGTTCCTCTTATTTTATTTCTACATCAATGTTCATATTTAACTTATTATTTTATCTTATTTTTAAATTTCTTTTATGTTGAGCCTTGATGAAAGCCATAGGTTCTCTCATATAATTGTATGTGTATGTATGTATATGTACATAATATATACATATATGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATATATACATATATG";
            final String hapCigar = "399M";
            final String readBases = "ATGTACATAATATATACATATATGTATATGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATAT";
            final GATKSAMRecord read = makeRead(readBases);
            final int refStart = 10130100;
            final int hapStart = 500;
            final String badCigar = "31M6D211M";
            final String goodCigar = "28M6D214M";
            final Haplotype badHap = new Haplotype(hap.getBytes());
            badHap.setCigar(TextCigarCodec.decode(hapCigar));
            badHap.setAlignmentStartHapwrtRef(hapStart);

            final int expectedPos = 10130740;
            tests.add(new Object[]{read, badHap, refStart, expectedPos, goodCigar});
        }

        return tests.toArray(new Object[][]{});
    }


    @Test(dataProvider = "ReadAlignedToRefData", enabled = true)
    public void testReadAlignedToRef(final GATKSAMRecord read, final Haplotype haplotype, final int refStart, final int expectedReadStart, final String expectedReadCigar) throws Exception {
        final HaplotypeBAMWriter writer = new CalledHaplotypeBAMWriter(new MockDestination());
        final GATKSAMRecord originalReadCopy = (GATKSAMRecord)read.clone();

        if ( expectedReadCigar == null ) {
            Assert.assertNull(AlignmentUtils.createReadAlignedToRef(read, haplotype, haplotype, refStart, true));
        } else {
            final Cigar expectedCigar = TextCigarCodec.decode(expectedReadCigar);
            final GATKSAMRecord alignedRead = AlignmentUtils.createReadAlignedToRef(read, haplotype, haplotype, refStart, true);

            Assert.assertEquals(alignedRead.getReadName(), originalReadCopy.getReadName());
            Assert.assertEquals(alignedRead.getAlignmentStart(), expectedReadStart);
            Assert.assertEquals(alignedRead.getReadBases(), originalReadCopy.getReadBases());
            Assert.assertEquals(alignedRead.getBaseQualities(), originalReadCopy.getBaseQualities());
            Assert.assertEquals(alignedRead.getAlignmentStart(), expectedReadStart);
            Assert.assertEquals(alignedRead.getCigar(), expectedCigar);
            Assert.assertNotNull(alignedRead.getAttribute("HC"));
        }

        Assert.assertEquals(read, originalReadCopy, "createReadAlignedToRef seems be modifying the original read!");
    }

    private static class Mutation implements Comparable<Mutation> {
        int pos, len;
        CigarOperator operator;

        private Mutation(int pos, int len, CigarOperator operator) {
            this.pos = pos;
            this.len = len;
            this.operator = operator;
        }
        public int getNMismatches() { return len; }

        @Override
        public int compareTo(Mutation o) {
            return Integer.valueOf(pos).compareTo(o.pos);
        }

        private String apply(final String seq) {
            switch ( operator ) {
                case M:
                    final byte[] bases = seq.getBytes();
                    if ( pos < seq.length() )
                        bases[pos] = (byte)(bases[pos] == 'A' ? 'C' : 'A');
                    return new String(bases);
                case I: {
                    final String prefix = seq.substring(0, pos);
                    final String postfix = seq.substring(pos, seq.length());
                    return prefix + "GTCAGTTA".substring(0, len) + postfix;
                } case D: {
                    final String prefix = seq.substring(0, pos);
                    final String postfix = seq.substring(pos + len, seq.length());
                    return prefix + postfix;
                }default:
                    throw new IllegalStateException("Unexpected operator " + operator);
            }
        }
    }

    private static class MutatedSequence {
        int numMismatches;
        String seq;

        private MutatedSequence(int numMismatches, String seq) {
            this.numMismatches = numMismatches;
            this.seq = seq;
        }
    }

    private MutatedSequence mutateSequence(final String hapIn, final List<Mutation> mutations) {
        Collections.sort(mutations);
        int mismatches = 0;
        String hap = hapIn;
        for ( final Mutation mut : mutations ) {
            hap = mut.apply(hap);
            mismatches += mut.getNMismatches();
        }
        return new MutatedSequence(mismatches, hap);
    }

    @DataProvider(name = "ComplexReadAlignedToRef")
    public Object[][] makeComplexReadAlignedToRef() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final List<Mutation> allMutations = Arrays.asList(
                new Mutation(1, 1, CigarOperator.M),
                new Mutation(2, 1, CigarOperator.M),
                new Mutation(3, 1, CigarOperator.I),
                new Mutation(7, 1, CigarOperator.D)
        );

        int i = 0;
        final String referenceBases  = "ACTGACTGACTG";
        final String paddedReference = "NNNN" + referenceBases + "NNNN";
        for ( final List<Mutation> mutations : Utils.makePermutations(allMutations, 3, false) ) {
            final MutatedSequence hap = mutateSequence(referenceBases, mutations);
            final Haplotype haplotype = new Haplotype(hap.seq.getBytes());
            final SWPairwiseAlignment align = new SWPairwiseAlignment(paddedReference.getBytes(), hap.seq.getBytes());
            haplotype.setAlignmentStartHapwrtRef(align.getAlignmentStart2wrt1());
            haplotype.setCigar(align.getCigar());

            for ( final List<Mutation> readMutations : Utils.makePermutations(allMutations, 3, false) ) {
                final MutatedSequence readBases = mutateSequence(hap.seq, readMutations);
                final GATKSAMRecord read = makeRead(readBases.seq);
                tests.add(new Object[]{i++, read, paddedReference, haplotype, hap.numMismatches + readBases.numMismatches});
            }
        }

        // for convenient testing of a single failing case
        //tests.add(new Object[]{makeRead("ACCGGGACTGACTG"), reference, makeHaplotype("AAAGGACTGACTG", "1M1I11M"), 2});

        return tests.toArray(new Object[][]{});
    }


    @Test(dataProvider = "ComplexReadAlignedToRef", enabled = !DEBUG)
    public void testReadAlignedToRefComplexAlignment(final int testIndex, final GATKSAMRecord read, final String reference, final Haplotype haplotype, final int expectedMaxMismatches) throws Exception {
        final HaplotypeBAMWriter writer = new CalledHaplotypeBAMWriter(new MockDestination());
        final GATKSAMRecord alignedRead = AlignmentUtils.createReadAlignedToRef(read, haplotype, new Haplotype(reference.getBytes(),true), 1, true);
        if ( alignedRead != null ) {
            final int mismatches = AlignmentUtils.getMismatchCount(alignedRead, reference.getBytes(), alignedRead.getAlignmentStart() - 1).numMismatches;
            Assert.assertTrue(mismatches <= expectedMaxMismatches,
                    "Alignment of read to ref looks broken.  Expected at most " + expectedMaxMismatches + " but saw " + mismatches
                            + " for readBases " + new String(read.getReadBases()) + " with cigar " + read.getCigar() + " reference " + reference + " haplotype "
                            + haplotype + " with cigar " + haplotype.getCigar() + " aligned read cigar " + alignedRead.getCigarString() + " @ " + alignedRead.getAlignmentStart());
        }
    }
}
