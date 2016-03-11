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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Tests {@link org.broadinstitute.gatk.tools.walkers.variantutils.ReferenceConfidenceVariantContextMerger}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class VariantContextMergerUnitTest  extends BaseTest {
    Allele Aref, T, C, G, Cref, ATC, ATCATC;
    Allele ATCATCT;
    Allele ATref;
    Allele Anoref;
    Allele GT;
    Allele del;

    private GenomeLocParser genomeLocParser;

    @BeforeSuite
    public void setup() throws IOException {
        // alleles
        Aref = Allele.create("A", true);
        Cref = Allele.create("C", true);
        T = Allele.create("T");
        C = Allele.create("C");
        G = Allele.create("G");
        ATC = Allele.create("ATC");
        ATCATC = Allele.create("ATCATC");
        ATCATCT = Allele.create("ATCATCT");
        ATref = Allele.create("AT",true);
        Anoref = Allele.create("A",false);
        del = Allele.SPAN_DEL;
        GT = Allele.create("GT",false);
        genomeLocParser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(new File(hg18Reference)));
    }

    @Test(dataProvider = "referenceConfidenceMergeData")
    public void testReferenceConfidenceMerge(final String testID, final List<VariantContext> toMerge, final GenomeLoc loc,
                                             final boolean returnSiteEvenIfMonomorphic, final boolean uniquifySamples, final VariantContext expectedResult) {
        final VariantContext result = ReferenceConfidenceVariantContextMerger.merge(toMerge, loc, returnSiteEvenIfMonomorphic ? (byte) 'A' : null, true, uniquifySamples, null);
        if ( result == null ) {
            Assert.assertTrue(expectedResult == null);
            return;
        }
        Assert.assertEquals(result.getAlleles(), expectedResult.getAlleles(),testID);
        Assert.assertEquals(result.getNSamples(), expectedResult.getNSamples(),testID);
        for ( final Genotype expectedGenotype : expectedResult.getGenotypes() ) {
            Assert.assertTrue(result.hasGenotype(expectedGenotype.getSampleName()), "Missing " + expectedGenotype.getSampleName());
            // use string comparisons to test equality for now
            Assert.assertEquals(result.getGenotype(expectedGenotype.getSampleName()).toString(), expectedGenotype.toString());
        }
    }

    @Test
    public void testGenerateADWithNewAlleles() {

        final int[] originalAD = new int[] {1,2,0};
        final int[] indexesOfRelevantAlleles = new int[] {0,1,2,2};

        final int[] newAD = ReferenceConfidenceVariantContextMerger.generateAD(originalAD, indexesOfRelevantAlleles);
        Assert.assertEquals(newAD, new int[]{1,2,0,0});
    }


    @Test(expectedExceptions = UserException.class)
    public void testGetIndexesOfRelevantAllelesWithNoALT() {

        final List<Allele> alleles1 = new ArrayList<>(1);
        alleles1.add(Allele.create("A", true));
        final List<Allele> alleles2 = new ArrayList<>(1);
        alleles2.add(Allele.create("A", true));
        GenotypeBuilder builder = new GenotypeBuilder();
        ReferenceConfidenceVariantContextMerger.getIndexesOfRelevantAlleles(alleles1, alleles2, -1, builder.make());
        Assert.fail("We should have thrown an exception because the <ALT> allele was not present");
    }

    @Test(dataProvider = "getIndexesOfRelevantAllelesData")
    public void testGetIndexesOfRelevantAlleles(final int allelesIndex, final List<Allele> allAlleles) {
        final List<Allele> myAlleles = new ArrayList<>(3);

        // always add the reference and <ALT> alleles
        myAlleles.add(allAlleles.get(0));
        myAlleles.add(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        // optionally add another alternate allele
        if ( allelesIndex > 0 )
            myAlleles.add(allAlleles.get(allelesIndex));

        GenotypeBuilder builder = new GenotypeBuilder();

        final int[] indexes = ReferenceConfidenceVariantContextMerger.getIndexesOfRelevantAlleles(myAlleles, allAlleles, -1, builder.make());

        Assert.assertEquals(indexes.length, allAlleles.size());

        for ( int i = 0; i < allAlleles.size(); i++ ) {
            if ( i == 0 )
                Assert.assertEquals(indexes[i], 0);    // ref should always match
            else if ( i == allelesIndex )
                Assert.assertEquals(indexes[i], 2);    // allele
            else
                Assert.assertEquals(indexes[i], 1);    // <ALT>
        }
    }


    @DataProvider(name = "referenceConfidenceMergeData")
    public Object[][] makeReferenceConfidenceMergeData() {
        final List<Object[]> tests = new ArrayList<>();
        final int start = 10;
        final GenomeLoc loc = new UnvalidatingGenomeLoc("20", 0, start, start);
        final VariantContext VCbase = new VariantContextBuilder("test", "20", start, start, Arrays.asList(Aref)).make();
        final VariantContext VCbase2 = new VariantContextBuilder("test2", "20", start, start, Arrays.asList(Aref)).make();
        final VariantContext VCprevBase = new VariantContextBuilder("test", "20", start-1, start-1, Arrays.asList(Aref)).make();

        final int[] standardPLs = new int[]{30, 20, 10, 71, 72, 73};
        final int[] reorderedSecondAllelePLs = new int[]{30, 71, 73, 20, 72, 10};

        final List<Allele> noCalls = new ArrayList<>(2);
        noCalls.add(Allele.NO_CALL);
        noCalls.add(Allele.NO_CALL);

        final List<Allele> A_ALT = Arrays.asList(Aref, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gA_ALT = new GenotypeBuilder("A").PL(new int[]{0, 100, 1000}).alleles(noCalls).make();
        final VariantContext vcA_ALT = new VariantContextBuilder(VCbase).alleles(A_ALT).genotypes(gA_ALT).make();

        final Allele AAref = Allele.create("AA", true);
        final List<Allele> AA_ALT = Arrays.asList(AAref, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gAA_ALT = new GenotypeBuilder("AA").PL(new int[]{0, 80, 800}).alleles(noCalls).make();
        final VariantContext vcAA_ALT = new VariantContextBuilder(VCprevBase).alleles(AA_ALT).genotypes(gAA_ALT).make();

        final List<Allele> A_C = Arrays.asList(Aref, C);
        final Genotype gA_C = new GenotypeBuilder("A_C").PL(new int[]{30, 20, 10}).alleles(noCalls).make();
        final List<Allele> A_C_ALT = Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gA_C_ALT = new GenotypeBuilder("A_C").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcA_C = new VariantContextBuilder(VCbase2).alleles(A_C_ALT).genotypes(gA_C).make();
        final VariantContext vcA_C_ALT = new VariantContextBuilder(VCbase).alleles(A_C_ALT).genotypes(gA_C_ALT).make();

        final List<Allele> A_G_ALT = Arrays.asList(Aref, G, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gA_G_ALT = new GenotypeBuilder("A_G").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcA_G_ALT = new VariantContextBuilder(VCbase).alleles(A_G_ALT).genotypes(gA_G_ALT).make();

        final List<Allele> A_C_G = Arrays.asList(Aref, C, G);
        final Genotype gA_C_G = new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30}).alleles(noCalls).make();
        final List<Allele> A_C_G_ALT = Arrays.asList(Aref, C, G, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gA_C_G_ALT = new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30, 71, 72, 73, 74}).alleles(noCalls).make();
        final VariantContext vcA_C_G = new VariantContextBuilder(VCbase2).alleles(A_C_G_ALT).genotypes(gA_C_G).make();
        final VariantContext vcA_C_G_ALT = new VariantContextBuilder(VCbase).alleles(A_C_G_ALT).genotypes(gA_C_G_ALT).make();

        final List<Allele> A_ATC_ALT = Arrays.asList(Aref, ATC, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gA_ATC_ALT = new GenotypeBuilder("A_ATC").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcA_ATC_ALT = new VariantContextBuilder(VCbase).alleles(A_ATC_ALT).genotypes(gA_ATC_ALT).make();

        final Allele A = Allele.create("A", false);
        final List<Allele> AA_A_ALT = Arrays.asList(AAref, A, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final Genotype gAA_A_ALT = new GenotypeBuilder("AA_A").PL(standardPLs).alleles(noCalls).make();
        final VariantContext vcAA_A_ALT = new VariantContextBuilder(VCprevBase).alleles(AA_A_ALT).genotypes(gAA_A_ALT).make();
        final List<Allele> A_C_del = Arrays.asList(Aref, C, del);

        // first test the case of a single record
        tests.add(new Object[]{"test00",Arrays.asList(vcA_C_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(gA_C).make()});

        // now, test pairs:
        // a SNP with another SNP
        tests.add(new Object[]{"test01",Arrays.asList(vcA_C_ALT, vcA_G_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(gA_C_ALT, new GenotypeBuilder("A_G").PL(reorderedSecondAllelePLs).alleles(noCalls).make()).make()});
        // a SNP with an indel
        tests.add(new Object[]{"test02",Arrays.asList(vcA_C_ALT, vcA_ATC_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(Arrays.asList(Aref, C, ATC)).genotypes(gA_C_ALT, new GenotypeBuilder("A_ATC").PL(reorderedSecondAllelePLs).alleles(noCalls).make()).make()});
        // a SNP with 2 SNPs
        tests.add(new Object[]{"test03",Arrays.asList(vcA_C_ALT, vcA_C_G_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(gA_C_ALT, gA_C_G).make()});
        // a SNP with a ref record
        tests.add(new Object[]{"test04",Arrays.asList(vcA_C_ALT, vcA_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(gA_C, gA_ALT).make()});

        // spanning records:
        // a SNP with a spanning ref record
        tests.add(new Object[]{"test05",Arrays.asList(vcA_C_ALT, vcAA_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(gA_C, gAA_ALT).make()});
        // a SNP with a spanning deletion
        tests.add(new Object[]{"test06",Arrays.asList(vcA_C_ALT, vcAA_A_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(A_C_del).genotypes(new GenotypeBuilder("A_C").PL(new int[]{30, 20, 10, 71, 72, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("AA_A").PL(new int[]{30, 71, 73, 20, 72, 10}).alleles(noCalls).make()).make()});

        // combination of all
        tests.add(new Object[]{"test07",Arrays.asList(vcA_C_ALT, vcA_G_ALT, vcA_ATC_ALT, vcA_C_G_ALT, vcA_ALT, vcAA_ALT, vcAA_A_ALT),
                loc, false, false,
                new VariantContextBuilder(VCbase).alleles(Arrays.asList(Aref, C, G, ATC, del)).genotypes(new GenotypeBuilder("A_C").PL(new int[]{30, 20, 10, 71, 72, 73, 71, 72, 73, 73, 71, 72, 73, 73, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_G").PL(new int[]{30, 71, 73, 20, 72, 10, 71, 73, 72, 73, 71, 73, 72, 73, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_ATC").PL(new int[]{30, 71, 73, 71, 73, 73, 20, 72, 72, 10, 71, 73, 73, 72, 73}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C_G").PL(new int[]{40, 20, 30, 20, 10, 30, 71, 72, 73, 74, 71, 72, 73, 74, 74}).alleles(noCalls).make(),
                        new GenotypeBuilder("A").PL(new int[]{0, 100, 1000, 100, 1000, 1000, 100, 1000, 1000, 1000, 100, 1000, 1000, 1000, 1000}).alleles(noCalls).make(),
                        new GenotypeBuilder("AA").PL(new int[]{0, 80, 800, 80, 800, 800, 80, 800, 800, 800, 80, 800, 800, 800, 800}).alleles(noCalls).make(),
                        new GenotypeBuilder("AA_A").PL(new int[]{30, 71, 73, 71, 73, 73, 71, 73, 73, 73, 20, 72, 72, 72, 10}).alleles(noCalls).make()).make()});

        // just spanning ref contexts, trying both instances where we want/do not want ref-only contexts
        tests.add(new Object[]{"test08",Arrays.asList(vcAA_ALT),

                loc, false, false,
                null});
        tests.add(new Object[]{"test09", Arrays.asList(vcAA_ALT),
                loc, true, false,
                new VariantContextBuilder(VCbase).alleles(Arrays.asList(Allele.create("A", true))).genotypes(new GenotypeBuilder("AA").PL(new int[]{0}).alleles(noCalls).make()).make()});

        // test uniquification of sample names
        tests.add(new Object[]{"test10",Arrays.asList(vcA_C, vcA_C_ALT), loc, false, true,
                new VariantContextBuilder(VCbase).alleles(A_C).genotypes(
                        new GenotypeBuilder("A_C.test2").PL(new int[]{30, 20, 10}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C.test").PL(new int[]{30, 20, 10}).alleles(noCalls).make()).make()});

        tests.add(new Object[]{"test11",Arrays.asList(vcA_C_G, vcA_C_G_ALT), loc, false, true,
                new VariantContextBuilder(VCbase).alleles(A_C_G).genotypes(
                        new GenotypeBuilder("A_C_G.test2").PL(new int[]{40, 20, 30, 20, 10, 30}).alleles(noCalls).make(),
                        new GenotypeBuilder("A_C_G.test").PL(new int[]{40, 20, 30, 20, 10, 30}).alleles(noCalls).make()).make()});

        final Object[][] result = tests.toArray(new Object[][]{});
        return result;
    }
    @DataProvider(name = "getIndexesOfRelevantAllelesData")
    public Object[][] makeGetIndexesOfRelevantAllelesData() {
        final int totalAlleles = 5;
        final List<Allele> alleles = new ArrayList<>(totalAlleles);
        alleles.add(Allele.create("A", true));
        for ( int i = 1; i < totalAlleles; i++ )
            alleles.add(Allele.create(Utils.dupString('A', i + 1), false));

        final List<Object[]> tests = new ArrayList<>();

        for ( int alleleIndex = 0; alleleIndex < totalAlleles; alleleIndex++ ) {
            tests.add(new Object[]{alleleIndex, alleles});
        }

        return tests.toArray(new Object[][]{});
    }
}
