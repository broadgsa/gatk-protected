/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;


public class ConstrainedMateFixingManagerUnitTest extends BaseTest {

    private static SAMFileHeader header;
    private static GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(3, 1, 10000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @Test
    public void testSecondaryAlignmentsDoNotInterfere() {
        final List<GATKSAMRecord> properReads = ArtificialSAMUtils.createPair(header, "foo", 1, 10, 30, true, false);
        final GATKSAMRecord read1 = properReads.get(0);
        read1.setAlignmentStart(8); // move the read
        read1.setFlags(99);   // first in proper pair, mate negative strand

        final GATKSAMRecord read2Primary = properReads.get(1);
        read2Primary.setFlags(147);   // second in pair, mate unmapped, not primary alignment

        Assert.assertEquals(read1.getInferredInsertSize(), 21);

        final GATKSAMRecord read2NonPrimary = new GATKSAMRecord(read2Primary);
        read2NonPrimary.setFlags(393);   // second in proper pair, on reverse strand

        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, genomeLocParser, 1000, 1000, 1000);
        manager.addRead(read1, true, false);
        manager.addRead(read2NonPrimary, false, false);
        manager.addRead(read2Primary, false, false);

        Assert.assertEquals(manager.getNReadsInQueue(), 3);

        for ( final SAMRecord read : manager.getReadsInQueueForTesting() ) {
            if ( read.getFirstOfPairFlag() ) {
                Assert.assertEquals(read.getFlags(), 99);
                Assert.assertEquals(read.getInferredInsertSize(), 23);
            } else if ( read.getNotPrimaryAlignmentFlag() ) {
                Assert.assertEquals(read.getFlags(), 393);
                Assert.assertEquals(read.getInferredInsertSize(), -21);
            } else {
                Assert.assertEquals(read.getFlags(), 147);
                Assert.assertEquals(read.getInferredInsertSize(), -23);
            }
        }
    }

    @Test
    public void testSecondaryAlignmentsDoNotCauseAccidentalRemovalOfMate() {
        final List<GATKSAMRecord> properReads = ArtificialSAMUtils.createPair(header, "foo", 1, 530, 1594, true, false);
        final GATKSAMRecord read1 = properReads.get(0);
        read1.setFlags(99);   // first in proper pair, mate negative strand

        final GATKSAMRecord read2Primary = properReads.get(1);
        read2Primary.setFlags(147);   // second in pair, mate unmapped, not primary alignment
        read2Primary.setAlignmentStart(1596); // move the read

        final GATKSAMRecord read2NonPrimary = new GATKSAMRecord(read2Primary);
        read2NonPrimary.setReadName("foo");
        read2NonPrimary.setFlags(393);   // second in proper pair, on reverse strand
        read2NonPrimary.setAlignmentStart(451);
        read2NonPrimary.setMateAlignmentStart(451);

        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, genomeLocParser, 10000, 200, 10000);
        manager.addRead(read2NonPrimary, false, false);
        manager.addRead(read1, false, false);

        for ( int i = 0; i < ConstrainedMateFixingManager.EMIT_FREQUENCY; i++ )
            manager.addRead(ArtificialSAMUtils.createArtificialRead(header, "foo" + i, 0, 1500, 10), false, false);

        Assert.assertTrue(manager.forMateMatching.containsKey("foo"));
    }

    @Test
    public void testSupplementaryAlignmentsDoNotCauseBadMateFixing() {
        final List<GATKSAMRecord> properReads = ArtificialSAMUtils.createPair(header, "foo", 1, 1000, 2000, true, false);
        final GATKSAMRecord read1 = properReads.get(0);
        read1.setFlags(99);   // first in pair, negative strand

        final GATKSAMRecord read2 = properReads.get(1);
        read2.setFlags(161);   // second in pair, mate negative strand

        final GATKSAMRecord read2Supp = new GATKSAMRecord(read2);
        read2Supp.setReadName("foo");
        read2Supp.setFlags(2209);   // second in pair, mate negative strand, supplementary
        read2Supp.setAlignmentStart(100);
        read2Supp.setMateAlignmentStart(1000);

        final DummyWriter writer = new DummyWriter();
        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(writer, genomeLocParser, 10000, 200, 10000);
        manager.addRead(read2Supp, false, false);
        manager.addRead(read1, false, false);
        manager.addRead(read2, false, false);
        manager.close(); // "write" the reads to our dummy writer

        // check to make sure that none of the mate locations were changed, which is the problem brought to us by a user
        for ( final SAMRecord read : writer.reads ) {
            final int start = read.getAlignmentStart();
            switch (start) {
                case 100:
                    Assert.assertEquals(read.getMateAlignmentStart(), 1000);
                    break;
                case 1000:
                    Assert.assertEquals(read.getMateAlignmentStart(), 2000);
                    break;
                case 2000:
                    Assert.assertEquals(read.getMateAlignmentStart(), 1000);
                    break;
                default:
                    Assert.assertTrue(false, "We saw a read located at the wrong position");
            }
        }
    }

    private class DummyWriter implements SAMFileWriter {

        public List<SAMRecord> reads;

        public DummyWriter() { reads = new ArrayList<>(10); }

        public void addAlignment(final SAMRecord alignment) { reads.add(alignment);}

        public SAMFileHeader getFileHeader() { return null; }

        public void setProgressLogger(final ProgressLoggerInterface progress) {}

        public void close() {}
    }
}