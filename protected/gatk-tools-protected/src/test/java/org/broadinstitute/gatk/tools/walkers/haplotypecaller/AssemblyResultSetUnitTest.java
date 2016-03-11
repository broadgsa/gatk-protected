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

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.RandomDNA;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests for {@link AssemblyResultSet}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class AssemblyResultSetUnitTest extends BaseTest
{
    private GenomeLocParser genomeLocParser;
    private SAMFileHeader header;

    @BeforeClass
    public void init() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }


    @Test
    public void testEmptyResultSet() {
        final AssemblyResultSet subject = new AssemblyResultSet();

        Assert.assertEquals(subject.getHaplotypeList().size(), 0);
        Assert.assertEquals(subject.getHaplotypeCount(),0);
        Assert.assertEquals(subject.getReferenceHaplotype(),null);
        Assert.assertEquals(subject.getFullReferenceWithPadding(),null);
        Assert.assertEquals(subject.getPaddedReferenceLoc(),null);
        Assert.assertEquals(subject.getRegionForGenotyping(),null);
        Assert.assertEquals(subject.getUniqueReadThreadingGraph(10),null);
        Assert.assertFalse(subject.hasMultipleKmerSizes());
    }

    @Test
    public void testAddReferenceHaplotype() {

        final Haplotype ref = new Haplotype("ACGT".getBytes(),true);
        ref.setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,ref.length() + 1 ));
        final AssemblyResultSet subject = new AssemblyResultSet();

        Assert.assertTrue(subject.add(ref));
        Assert.assertFalse(subject.add(ref));

        Assert.assertEquals(subject.getReferenceHaplotype(),ref);
        Assert.assertEquals(subject.getHaplotypeCount(),1);
        Assert.assertEquals(subject.getHaplotypeList().size(),1);
    }

    @Test(dataProvider="assemblyResults")
    public void testAddManyHaplotypes(final java.util.List<AssemblyResult> assemblyResults,
                                      final java.util.List<java.util.List<Haplotype>> haplotypes) {
        final AssemblyResultSet subject = new AssemblyResultSet();
        for (int i = 0; i < haplotypes.size(); i++) {
            final int haplotypeCountBefore = subject.getHaplotypeCount();
            final java.util.List<Haplotype> haplos = haplotypes.get(i);
            final AssemblyResult ar = assemblyResults.get(i);
            for (final Haplotype h : haplos) {
                Assert.assertTrue(subject.add(h, ar));
                Assert.assertFalse(subject.add(h,ar));
                if (h.isReference())
                    Assert.assertEquals(subject.getReferenceHaplotype(),h);
            }
            final int haplotypeCountAfter = subject.getHaplotypeCount();
            Assert.assertEquals(haplos.size(),haplotypeCountAfter - haplotypeCountBefore);
            Assert.assertTrue(subject.getMaximumKmerSize() >= ar.getKmerSize());
            Assert.assertTrue(subject.getMinimumKmerSize() <= ar.getKmerSize());
            Assert.assertEquals(subject.getUniqueReadThreadingGraph(ar.getKmerSize()), ar.getThreadingGraph());
        }
    }

    @Test(dataProvider="trimmingData")
    public void testTrimTo(final Map<Haplotype,AssemblyResult> haplotypesAndResultSets, final ActiveRegion original) {
        final AssemblyResultSet subject = new AssemblyResultSet();
        for (final Map.Entry<Haplotype,AssemblyResult> entry : haplotypesAndResultSets.entrySet())
            subject.add(entry.getKey(),entry.getValue());
        subject.setRegionForGenotyping(original);
        final GenomeLoc originalLocation = original.getExtendedLoc();
        final int length = originalLocation.size();
        final GenomeLoc newLocation = originalLocation.setStop(originalLocation.setStart(originalLocation,originalLocation.getStart() + length / 2),originalLocation.getStop() - length / 2);
        final ActiveRegion newRegion = original.trim(newLocation);

        final Map<Haplotype,Haplotype> originalHaplotypesByTrimmed = new HashMap<>(haplotypesAndResultSets.size());
        for (final Haplotype h : haplotypesAndResultSets.keySet())
            originalHaplotypesByTrimmed.put(h.trim(newRegion.getExtendedLoc()), h);

        final AssemblyResultSet trimmed = subject.trimTo(newRegion);

        Assert.assertFalse(subject.wasTrimmed());
        Assert.assertTrue(trimmed.wasTrimmed());

        for (final Haplotype h : trimmed.getHaplotypeList()) {
            Assert.assertEquals(h.getGenomeLocation(),newLocation);
            Assert.assertEquals(h.getBases().length,newLocation.size());
        }
    }

    @DataProvider(name="trimmingData")
    public Iterator<Object[]> trimmingData() {
        final ActiveRegion activeRegion = new ActiveRegion(genomeLocParser.createGenomeLoc("chr1",1000,1100),genomeLocParser,25);
        final int length = activeRegion.getExtendedLoc().size();
        final RandomDNA rnd = new RandomDNA(13); // keep it prepoducible by fixing the seed to lucky 13.
        final ActiveRegionTestDataSet actd = new ActiveRegionTestDataSet(10,new String(rnd.nextBases(length)),new String[] {
                "Civar:*1T*" }, new String[0], new byte[0], new byte[0], new byte[0]);

        final List<Haplotype> haplotypes = actd.haplotypeList();
        for (final Haplotype h : haplotypes)
            h.setGenomeLocation(activeRegion.getExtendedLoc());

        final ReadThreadingGraph rtg = new ReadThreadingGraph(10);
        for (final Haplotype h : haplotypes)
            rtg.addSequence("seq-" + Math.abs(h.hashCode()), h.getBases(), h.isReference());
        final SeqGraph seqGraph = rtg.convertToSequenceGraph();
        final AssemblyResult ar = new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION,seqGraph);
        ar.setThreadingGraph(rtg);
        final Map<Haplotype,AssemblyResult> result =
                new HashMap<>();
        for (final Haplotype h : haplotypes)
            result.put(h,ar);
        return Collections.singleton(new Object[] {result,activeRegion}).iterator();

    }




    @DataProvider(name="assemblyResults")
    public java.util.Iterator<Object[]> assemblyResults() {
        final int size = THREE_KS_GRAPH_AND_HAPLOTYPES.length * (1 + TEN_KS_GRAPH_AND_HAPLOTYPES.length);
        final Object[][] result = new Object[size][];

        for (int i = 0; i < THREE_KS_GRAPH_AND_HAPLOTYPES.length; i++) {
            final ReadThreadingGraph rtg = new ReadThreadingGraph((String) THREE_KS_GRAPH_AND_HAPLOTYPES[i][0]);
            final AssemblyResult ar = new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION,rtg.convertToSequenceGraph());
            ar.setThreadingGraph(rtg);
            final Object[] haplotypeStrings = (Object[]) THREE_KS_GRAPH_AND_HAPLOTYPES[i][1];
            final Haplotype[] haplotypes = new Haplotype[haplotypeStrings.length];
            for (int j = 0; j < haplotypeStrings.length; j++) {
                haplotypes[j] = new Haplotype(((String)haplotypeStrings[j]).getBytes(),j == 0);
                haplotypes[j].setGenomeLocation(genomeLocParser.createGenomeLoc("chr1",1,haplotypes[j].length() + 1));
            }
            result[i] = new Object[] { Collections.singletonList(ar),Arrays.asList(Arrays.asList(haplotypes))};
            for (int j = 0; j < TEN_KS_GRAPH_AND_HAPLOTYPES.length; j++) {
                final ReadThreadingGraph rtg10 = new ReadThreadingGraph((String) TEN_KS_GRAPH_AND_HAPLOTYPES[j][0]);
                final AssemblyResult ar10 = new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION,rtg10.convertToSequenceGraph());
                ar10.setThreadingGraph(rtg10);
                final Object[] haplotypeStrings10 = (Object[]) TEN_KS_GRAPH_AND_HAPLOTYPES[j][1];
                final Haplotype[] haplotype10 = new Haplotype[haplotypeStrings10.length];
                for (int k = 0; k < haplotypeStrings10.length; k++) {
                    haplotype10[k] = new Haplotype(((String)haplotypeStrings10[k]).getBytes(),false);
                    haplotype10[k].setGenomeLocation(genomeLocParser.createGenomeLoc("chr1", 1, haplotype10[k].length() + 1));
                }

                result[THREE_KS_GRAPH_AND_HAPLOTYPES.length + i * TEN_KS_GRAPH_AND_HAPLOTYPES.length + j] = new Object[] { Arrays.asList(ar,ar10),
                        Arrays.asList( Arrays.asList(haplotypes), Arrays.asList(haplotype10)) };
            }
        }
        return Arrays.asList(result).iterator();
    }


    private static final Object[][] THREE_KS_GRAPH_AND_HAPLOTYPES = new Object[][] {
            {"[ks=3]{REF: ACT}",new Object[] {"ACT"}},
            {"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                    "{ (3) -> A -> G ->          (2) }" +
                    "{  (1) -> A -> G ->  (2) }",new Object[] {"ACTTGA","ACTAGGA","ACTTAGGA"}},
            {"[ks=3]{REF: ACT -> C(1) -> G}{ACT -> C(1) -> G}{ACT -> C(1) -> G}", new Object[] {"ACTCG"}} ,
            {"[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> C -> G -> T }" +
                    "{A(1) -> T -> A(2) }", new Object[] {"ACTAGACGT","ACTATACGT"}}  ,
            {"[ks=3]{REF: ACT -> A -> T(2) -> C -> A -> G -> T -> A -> C -> G -> T -> A(1) -> T}" +
                    "{ ACT -> A -> T(2) -> C -> T -> A -> C -> G -> T -> A(1) -> T}",
                           new Object[] {"ACTATCAGTACGTAT","ACTATCTACGTAT"}} ,
            {"[ks=3]{REF: ACT -> A -> T    -> C -> A -> G -> T -> A -> C -> G -> T -> A    -> T}",
                           new Object[] {"ACTATCAGTACGTAT"}},
            {"[ks=3]{REF: ACT -> A -> T(1) }" +
                    "{ ACT -> A -> T(1) }", new Object[] {"ACTAT"}},
            {"[ks=3]{REF: TTT -> A(1) -> C -> T(2)}{ A(1) -> T(2) } ", new Object[] {"TTTACT","TTTAT"}}
    };

    private static final Object[][] TEN_KS_GRAPH_AND_HAPLOTYPES = new Object[][] {
            {"[ks=10]{ACTAGTAAAT -> A -> T -> A -> A -> T -> A", new Object[] {"ACTAGTAAATATAATA"}},
            {"[ks=10]{ATAGTAATAA(1) -> A -> C -> T -> A(2) -> C}{ (1) -> C -> C -> C -> A(2) -> C}",
                    new Object[] {"ATAGTAATAAACTAC","ATAGTAATAACCCAC"}},

    };

}
