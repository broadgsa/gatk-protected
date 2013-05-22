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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.UnvalidatingGenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.collections.PrimitivePair;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class LocalAssemblyEngineUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;
    private IndexedFastaSequenceFile seq;
    private SAMFileHeader header;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
        header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
    }

    private enum Assembler {DEBRUIJN_ASSEMBLER, READ_THREADING_ASSEMBLER}
    private LocalAssemblyEngine createAssembler(final Assembler type) {
        switch ( type ) {
            case DEBRUIJN_ASSEMBLER: return new DeBruijnAssembler();
            case READ_THREADING_ASSEMBLER: return new ReadThreadingAssembler();
            default: throw new IllegalStateException("Unexpected " + type);
        }
    }

    @DataProvider(name = "AssembleIntervalsData")
    public Object[][] makeAssembleIntervalsData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String contig = "20";
        final int start = 10000000;
        final int end   = 10100000;
        final int windowSize = 100;
        final int stepSize = 200;
        final int nReadsToUse = 5;

        for ( final Assembler assembler : Assembler.values() ) {
            for ( int startI = start; startI < end; startI += stepSize) {
                final int endI = startI + windowSize;
                final GenomeLoc refLoc = genomeLocParser.createGenomeLoc(contig, startI, endI);
                tests.add(new Object[]{assembler, refLoc, nReadsToUse});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "AssembleIntervalsWithVariantData")
    public Object[][] makeAssembleIntervalsWithVariantData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String contig = "20";
        final int start = 10000000;
        final int end   = 10001000;
        final int windowSize = 100;
        final int stepSize = 200;
        final int variantStepSize = 1;
        final int nReadsToUse = 5;

        for ( final Assembler assembler : Assembler.values() ) {
            for ( int startI = start; startI < end; startI += stepSize) {
                final int endI = startI + windowSize;
                final GenomeLoc refLoc = genomeLocParser.createGenomeLoc(contig, startI, endI);
                for ( int variantStart = windowSize / 2 - 10; variantStart < windowSize / 2 + 10; variantStart += variantStepSize ) {
                    tests.add(new Object[]{assembler, refLoc, nReadsToUse, variantStart});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AssembleIntervalsData")
    public void testAssembleRef(final Assembler assembler, final GenomeLoc loc, final int nReadsToUse) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();

        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        for ( int i = 0; i < nReadsToUse; i++ ) {
            final byte[] bases = refBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, refBases.length);
            final String cigar = refBases.length + "M";
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, loc.getContig(), loc.getContigIndex(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        // TODO -- generalize to all assemblers
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertEquals(haplotypes, Collections.singletonList(refHaplotype));
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndSNP(final Assembler assembler, final GenomeLoc loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        final Allele refBase = Allele.create(refBases[variantSite], true);
        final Allele altBase = Allele.create((byte)(refBase.getBases()[0] == 'A' ? 'C' : 'A'), false);
        final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite, Arrays.asList(refBase, altBase));
        testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndDeletion(final Assembler assembler, final GenomeLoc loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        for ( int deletionLength = 1; deletionLength < 10; deletionLength++ ) {
            final Allele refBase = Allele.create(new String(refBases).substring(variantSite, variantSite + deletionLength + 1), true);
            final Allele altBase = Allele.create(refBase.getBases()[0], false);
            final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite + deletionLength, Arrays.asList(refBase, altBase));
            testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
        }
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndInsertion(final Assembler assembler, final GenomeLoc loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        for ( int insertionLength = 1; insertionLength < 10; insertionLength++ ) {
            final Allele refBase = Allele.create(refBases[variantSite], false);
            final Allele altBase = Allele.create(new String(refBases).substring(variantSite, variantSite + insertionLength + 1), true);
            final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite + insertionLength, Arrays.asList(refBase, altBase));
            testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
        }
    }

    private void testAssemblyWithVariant(final Assembler assembler, final byte[] refBases, final GenomeLoc loc, final int nReadsToUse, final VariantContext site) {
        final String preRef = new String(refBases).substring(0, site.getStart());
        final String postRef = new String(refBases).substring(site.getEnd() + 1, refBases.length);
        final byte[] altBases = (preRef + site.getAlternateAllele(0).getBaseString() + postRef).getBytes();

//        logger.warn("ref " + new String(refBases));
//        logger.warn("alt " + new String(altBases));

        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        for ( int i = 0; i < nReadsToUse; i++ ) {
            final byte[] bases = altBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, altBases.length);
            final String cigar = altBases.length + "M";
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, loc.getContig(), loc.getContigIndex(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Haplotype altHaplotype = new Haplotype(altBases, false);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertEquals(haplotypes, Arrays.asList(refHaplotype, altHaplotype));
    }


    private List<Haplotype> assemble(final Assembler assembler, final byte[] refBases, final GenomeLoc loc, final List<GATKSAMRecord> reads) {
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);

        final ActiveRegion activeRegion = new ActiveRegion(loc, null, true, genomeLocParser, 0);
        activeRegion.addAll(reads);
        final LocalAssemblyEngine engine = createAssembler(assembler);
//        logger.warn("Assembling " + activeRegion + " with " + engine);
        return engine.runLocalAssembly(activeRegion, refHaplotype, refBases, loc, Collections.<VariantContext>emptyList());
    }

    @DataProvider(name = "SimpleAssemblyTestData")
    public Object[][] makeSimpleAssemblyTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String contig  = "20";
        final int start      = 10000000;
        final int windowSize = 200;
        final int end        = start + windowSize;

        final Map<Assembler, Integer> edgeExcludesByAssembler = new EnumMap<>(Assembler.class);
        edgeExcludesByAssembler.put(Assembler.DEBRUIJN_ASSEMBLER, 26);
        edgeExcludesByAssembler.put(Assembler.READ_THREADING_ASSEMBLER, 25); // TODO -- decrease to zero when the edge calling problem is fixed

        final String ref = new String(seq.getSubsequenceAt(contig, start, end).getBases());
        final GenomeLoc refLoc = genomeLocParser.createGenomeLoc(contig, start, end);

        for ( final Assembler assembler : Assembler.values() ) {
            final int excludeVariantsWithXbp = edgeExcludesByAssembler.get(assembler);
            for ( int snpPos = 0; snpPos < windowSize; snpPos++) {
                if ( snpPos > excludeVariantsWithXbp && (windowSize - snpPos) >= excludeVariantsWithXbp ) {
                    final byte[] altBases = ref.getBytes();
                    altBases[snpPos] = 'N';
                    final String alt = new String(altBases);
                    tests.add(new Object[]{"SNP at " + snpPos, assembler, refLoc, ref, alt});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SimpleAssemblyTestData")
    public void testSimpleAssembly(final String name, final Assembler assembler, final GenomeLoc loc, final String ref, final String alt) {
        final byte[] refBases = ref.getBytes();
        final byte[] altBases = alt.getBytes();

        final List<GATKSAMRecord> reads = new LinkedList<>();
        for ( int i = 0; i < 20; i++ ) {
            final byte[] bases = altBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, altBases.length);
            final String cigar = altBases.length + "M";
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, loc.getContig(), loc.getContigIndex(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Haplotype altHaplotype = new Haplotype(altBases, false);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertTrue(haplotypes.size() > 0, "Failed to find ref haplotype");
        Assert.assertEquals(haplotypes.get(0), refHaplotype);

        Assert.assertEquals(haplotypes.size(), 2, "Failed to find single alt haplotype");
        Assert.assertEquals(haplotypes.get(1), altHaplotype);
    }
}
