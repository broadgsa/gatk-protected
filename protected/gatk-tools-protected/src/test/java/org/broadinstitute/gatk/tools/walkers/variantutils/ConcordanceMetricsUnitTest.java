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

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.StringBufferInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ConcordanceMetricsUnitTest extends BaseTest {

    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);
    }
    public static String HEADER_BASE = "##fileformat=VCFv4.0\n" +
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n" +
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    public static String TEST_1_HEADER = HEADER_BASE + "test1_sample1\ttest1_sample2\ttest1_sample3\n";
    public static String TEST_2_HEADER = HEADER_BASE + "test2_sample1\ttest2_sample2\n";
    public static String TEST_3_HEADER_1 = HEADER_BASE + "test3_sample1\ttest3_sample2\ttest3_sample3\ttest3_sample4\ttest3_sample5\n";
    public static String TEST_3_HEADER_2 = HEADER_BASE + "test3_sample6\ttest3_sample7\ttest3_sample8\ttest3_sample9\ttest3_sample10\n";
    public static String TEST_3_HEADER_3 = HEADER_BASE + "test3_sample3\ttest3_sample6\ttest3_sample7\ttest3_sample8\ttest3_sample9\ttest3_sample10\n";


    private Pair<VariantContext,VariantContext> getData1() {

        Allele reference_A = Allele.create(BaseUtils.Base.A.base,true);
        Allele alt_C = Allele.create(BaseUtils.Base.C.base);

        Genotype sam_1_1_eval = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_eval = GenotypeBuilder.create("test1_sample2", Arrays.asList(reference_A,alt_C));
        Genotype sam_1_3_eval = GenotypeBuilder.create("test1_sample3", Arrays.asList(reference_A,alt_C));

        Genotype sam_1_1_truth = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_truth = GenotypeBuilder.create("test1_sample2", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_3_truth = GenotypeBuilder.create("test1_sample3", Arrays.asList(alt_C,alt_C));

        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 3, 3);
        VariantContextBuilder eval_1_builder = new VariantContextBuilder();
        VariantContextBuilder truth_1_builder = new VariantContextBuilder();

        eval_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        eval_1_builder.genotypes(Arrays.asList(sam_1_1_eval,sam_1_2_eval,sam_1_3_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_1_1_truth,sam_1_2_truth,sam_1_3_truth));

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());

        Pair<VariantContext,VariantContext> testData = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        return testData;
    }

    @Test(enabled=true)
    public void testSimpleComparison() {
        Pair<VariantContext,VariantContext> data = getData1();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);
        Assert.assertEquals(eval.getGenotype("test1_sample2").getType().ordinal(), 2);
        Assert.assertEquals(truth.getGenotype("test1_sample2").getType().ordinal(),1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getnMismatchingAlt(),0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[2][1],1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][3],1);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getTable()[1][1],1);
    }

    private Pair<VariantContext,VariantContext> getData2() {

        Allele reference_A = Allele.create(BaseUtils.Base.A.base,true);
        Allele alt_C = Allele.create(BaseUtils.Base.C.base);
        Allele alt_T = Allele.create(BaseUtils.Base.T.base);

        Genotype sam_1_1_eval = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_eval = GenotypeBuilder.create("test1_sample2", Arrays.asList(reference_A,alt_T));
        Genotype sam_1_3_eval = GenotypeBuilder.create("test1_sample3", Arrays.asList(reference_A,alt_C));

        Genotype sam_1_1_truth = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_truth = GenotypeBuilder.create("test1_sample2", Arrays.asList(reference_A,alt_C));
        Genotype sam_1_3_truth = GenotypeBuilder.create("test1_sample3", Arrays.asList(alt_C,alt_C));

        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 3, 3);
        VariantContextBuilder eval_1_builder = new VariantContextBuilder();
        VariantContextBuilder truth_1_builder = new VariantContextBuilder();

        eval_1_builder.alleles(Arrays.asList(reference_A,alt_C,alt_T));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        eval_1_builder.genotypes(Arrays.asList(sam_1_1_eval,sam_1_2_eval,sam_1_3_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_1_1_truth,sam_1_2_truth,sam_1_3_truth));

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());

        Pair<VariantContext,VariantContext> testData = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        return testData;
    }

    @Test(enabled=true)
    public void testMismatchingAlleleInAlleleSubset() {
        Pair<VariantContext,VariantContext> data = getData2();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);
        Assert.assertEquals(eval.getGenotype("test1_sample2").getType().ordinal(), 2);
        Assert.assertEquals(truth.getGenotype("test1_sample2").getType().ordinal(),2);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getnMismatchingAlt(),1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][3],1);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getTable()[1][1],1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.EVAL_SUPERSET_TRUTH.ordinal()],1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.ALLELES_DO_NOT_MATCH.ordinal()],0);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.ALLELES_MATCH.ordinal()],0);

        // now flip them around

        eval = data.getSecond();
        truth = data.getFirst();
        codec = new VCFCodec();
        evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);
        Assert.assertEquals(eval.getGenotype("test1_sample2").getType().ordinal(), 2);
        Assert.assertEquals(truth.getGenotype("test1_sample2").getType().ordinal(),2);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getnMismatchingAlt(),1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[1][2],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[1][2],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[3][2],1);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getTable()[1][1],1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.EVAL_SUPERSET_TRUTH.ordinal()],0);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.EVAL_SUBSET_TRUTH.ordinal()],1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.ALLELES_DO_NOT_MATCH.ordinal()],0);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.ALLELES_MATCH.ordinal()],0);
    }

    private Pair<VariantContext,VariantContext> getData3() {

        Allele reference_ACT = Allele.create(new byte[]{BaseUtils.Base.A.base,BaseUtils.Base.C.base,BaseUtils.Base.T.base},true);
        Allele alt_AC = Allele.create(new byte[]{BaseUtils.Base.A.base,BaseUtils.Base.C.base});
        Allele alt_A = Allele.create(BaseUtils.Base.A.base);
        Allele alt_ATT = Allele.create(new byte[]{BaseUtils.Base.A.base,BaseUtils.Base.T.base,BaseUtils.Base.T.base});

        Genotype sam_1_1_eval = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_ACT,alt_ATT));
        Genotype sam_1_2_eval = GenotypeBuilder.create("test1_sample2", Arrays.asList(alt_A,alt_A));
        Genotype sam_1_3_eval = GenotypeBuilder.create("test1_sample3", Arrays.asList(reference_ACT,alt_A));

        Genotype sam_1_1_truth = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_ACT,alt_AC));
        Genotype sam_1_2_truth = GenotypeBuilder.create("test1_sample2", Arrays.asList(alt_A,alt_A));
        Genotype sam_1_3_truth = GenotypeBuilder.create("test1_sample3", Arrays.asList(reference_ACT,alt_A));

        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 3, 5);
        VariantContextBuilder eval_1_builder = new VariantContextBuilder();
        VariantContextBuilder truth_1_builder = new VariantContextBuilder();

        eval_1_builder.alleles(Arrays.asList(reference_ACT,alt_ATT,alt_A));
        truth_1_builder.alleles(Arrays.asList(reference_ACT,alt_AC,alt_A));
        eval_1_builder.genotypes(Arrays.asList(sam_1_1_eval,sam_1_2_eval,sam_1_3_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_1_1_truth,sam_1_2_truth,sam_1_3_truth));

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());

        Pair<VariantContext,VariantContext> testData = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        return testData;
    }

    @Test(enabled=true)
    public void testComplex() {
        Pair<VariantContext,VariantContext> data = getData3();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample1").getnMismatchingAlt(),1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[3][3],1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[1][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][2],1);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getTable()[3][3],1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.EVAL_SUPERSET_TRUTH.ordinal()],0);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.ALLELES_DO_NOT_MATCH.ordinal()],1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().getSiteConcordance()[ConcordanceMetrics.SiteConcordanceType.ALLELES_MATCH.ordinal()],0);
    }

    private Pair<VariantContext,VariantContext> getData4() {

        Allele reference_A = Allele.create(BaseUtils.Base.A.base,true);
        Allele alt_C = Allele.create(BaseUtils.Base.C.base);
        Allele alt_T = Allele.create(BaseUtils.Base.T.base);

        Genotype sam_1_1_eval = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_eval = GenotypeBuilder.create("test1_sample2", Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
        Genotype sam_1_3_eval = GenotypeBuilder.create("test1_sample3", Arrays.asList(reference_A,alt_C));

        Genotype sam_1_1_truth = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_truth = GenotypeBuilder.create("test1_sample2", Arrays.asList(reference_A,alt_C));
        Genotype sam_1_3_truth = GenotypeBuilder.create("test1_sample3", Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));

        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 3, 3);
        VariantContextBuilder eval_1_builder = new VariantContextBuilder();
        VariantContextBuilder truth_1_builder = new VariantContextBuilder();

        eval_1_builder.alleles(Arrays.asList(reference_A,alt_C,alt_T));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        eval_1_builder.genotypes(Arrays.asList(sam_1_1_eval,sam_1_2_eval,sam_1_3_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_1_1_truth,sam_1_2_truth,sam_1_3_truth));

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());

        Pair<VariantContext,VariantContext> testData = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        return testData;
    }

    @Test(enabled=true)
    public void testNoCalls() {
        Pair<VariantContext,VariantContext> data = getData4();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getnMismatchingAlt(),0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[0][2],1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][3],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][0],1);
    }

    private Pair<VariantContext,VariantContext> getData5() {

        Allele reference_A = Allele.create(BaseUtils.Base.A.base,true);
        Allele alt_C = Allele.create(BaseUtils.Base.C.base);
        Allele alt_T = Allele.create(BaseUtils.Base.T.base);

        Genotype sam_1_1_eval = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_eval = GenotypeBuilder.create("test1_sample2", new ArrayList<Allele>(0));
        Genotype sam_1_3_eval = GenotypeBuilder.create("test1_sample3", Arrays.asList(reference_A,alt_C));

        Genotype sam_1_1_truth = GenotypeBuilder.create("test1_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_1_2_truth = GenotypeBuilder.create("test1_sample2", Arrays.asList(reference_A,alt_C));
        Genotype sam_1_3_truth = GenotypeBuilder.create("test1_sample3", new ArrayList<Allele>(0));

        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 3, 3);
        VariantContextBuilder eval_1_builder = new VariantContextBuilder();
        VariantContextBuilder truth_1_builder = new VariantContextBuilder();

        eval_1_builder.alleles(Arrays.asList(reference_A,alt_C,alt_T));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        eval_1_builder.genotypes(Arrays.asList(sam_1_1_eval,sam_1_2_eval,sam_1_3_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_1_1_truth,sam_1_2_truth,sam_1_3_truth));

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());

        Pair<VariantContext,VariantContext> testData = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        return testData;
    }

    @Test(enabled=true)
    public void testMissing() {
        Pair<VariantContext,VariantContext> data = getData5();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);
        Assert.assertTrue(eval.getGenotype("test1_sample2").getType().equals(GenotypeType.UNAVAILABLE));
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getnMismatchingAlt(),0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[0][2],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[4][2],1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][3],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][0],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][4],1);
    }

    private List<Pair<VariantContext,VariantContext>> getData6() {

        Allele reference_A = Allele.create(BaseUtils.Base.A.base,true);
        Allele alt_C = Allele.create(BaseUtils.Base.C.base);


        // site 1 -
        //  sample 1: hom-ref/hom-ref
        //  sample 2: het/hom-ref

        Genotype sam_2_1_1_eval = GenotypeBuilder.create("test2_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_2_2_1_eval = GenotypeBuilder.create("test2_sample2", Arrays.asList(reference_A,alt_C));

        Genotype sam_2_1_1_truth = GenotypeBuilder.create("test2_sample1", Arrays.asList(reference_A,reference_A));
        Genotype sam_2_2_1_truth = GenotypeBuilder.create("test2_sample2", Arrays.asList(reference_A,reference_A));

        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 3, 3);
        VariantContextBuilder eval_1_builder = new VariantContextBuilder();
        VariantContextBuilder truth_1_builder = new VariantContextBuilder();

        eval_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        eval_1_builder.genotypes(Arrays.asList(sam_2_1_1_eval,sam_2_2_1_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_2_1_1_truth,sam_2_2_1_truth));

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());

        Pair<VariantContext,VariantContext> testDataSite1 = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        reference_A = Allele.create(BaseUtils.Base.A.base,true);
        Allele alt_T = Allele.create(BaseUtils.Base.T.base);

        // site 2 -
        //  sample 1: no-call/hom-ref
        //  sample 2: hom-var/hom-var

        Genotype sam_2_1_2_eval = GenotypeBuilder.create("test2_sample1",Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
        Genotype sam_2_2_2_eval = GenotypeBuilder.create("test2_sample2",Arrays.asList(alt_T,alt_T));
        Genotype sam_2_1_2_truth = GenotypeBuilder.create("test2_sample1",Arrays.asList(reference_A,reference_A));
        Genotype sam_2_2_2_truth = GenotypeBuilder.create("test2_sample2",Arrays.asList(alt_T,alt_T));

        loc = genomeLocParser.createGenomeLoc("chr1", 4, 4);
        eval_1_builder = new VariantContextBuilder();
        truth_1_builder = new VariantContextBuilder();

        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        eval_1_builder.alleles(Arrays.asList(reference_A,alt_T));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_T));
        eval_1_builder.genotypes(Arrays.asList(sam_2_1_2_eval,sam_2_2_2_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_2_1_2_truth,sam_2_2_2_truth));

        Pair<VariantContext,VariantContext> testDataSite2 = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        Allele alt_G = Allele.create(BaseUtils.Base.G.base);

        // site 3 -
        //  sample 1: alleles do not match
        //  sample 2: het/het
        Genotype sam_2_1_3_eval = GenotypeBuilder.create("test2_sample1",Arrays.asList(alt_G,alt_T));
        Genotype sam_2_2_3_eval = GenotypeBuilder.create("test2_sample2",Arrays.asList(reference_A,alt_T));
        Genotype sam_2_1_3_truth = GenotypeBuilder.create("test2_sample1",Arrays.asList(alt_T,alt_T));
        Genotype sam_2_2_3_truth = GenotypeBuilder.create("test2_sample2",Arrays.asList(reference_A,alt_T));

        loc = genomeLocParser.createGenomeLoc("chr1",5,5);
        eval_1_builder = new VariantContextBuilder();
        truth_1_builder = new VariantContextBuilder();
        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        eval_1_builder.alleles(Arrays.asList(reference_A,alt_T,alt_G));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_T));
        eval_1_builder.genotypes(Arrays.asList(sam_2_1_3_eval,sam_2_2_3_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_2_1_3_truth,sam_2_2_3_truth));

        Pair<VariantContext,VariantContext> testDataSite3 = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        // site 4 -
        //  sample 1: unavailable/het
        //  sample 2: unavailable/ref
        Genotype sam_2_1_4_eval = GenotypeBuilder.create("test2_sample1",new ArrayList<Allele>(0));
        Genotype sam_2_2_4_eval = GenotypeBuilder.create("test2_sample2",new ArrayList<Allele>(0));
        Genotype sam_2_1_4_truth = GenotypeBuilder.create("test2_sample1",Arrays.asList(reference_A,alt_T));
        Genotype sam_2_2_4_truth = GenotypeBuilder.create("test2_sample2",Arrays.asList(reference_A,reference_A));

        loc = genomeLocParser.createGenomeLoc("chr1",6,6);
        eval_1_builder = new VariantContextBuilder();
        truth_1_builder = new VariantContextBuilder();
        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        eval_1_builder.alleles(Arrays.asList(reference_A,alt_T));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_T));
        eval_1_builder.genotypes(Arrays.asList(sam_2_1_4_eval,sam_2_2_4_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_2_1_4_truth,sam_2_2_4_truth));

        Pair<VariantContext,VariantContext> testDataSite4 = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        // site 5 -
        //  sample 1: hom-var/no-call
        //  sample 2: het/het
        Genotype sam_2_1_5_eval = GenotypeBuilder.create("test2_sample1",Arrays.asList(alt_C,alt_C));
        Genotype sam_2_2_5_eval = GenotypeBuilder.create("test2_sample2",Arrays.asList(reference_A,alt_C));
        Genotype sam_2_1_5_truth = GenotypeBuilder.create("test2_sample1",Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
        Genotype sam_2_2_5_truth = GenotypeBuilder.create("test2_sample2",Arrays.asList(reference_A,alt_C));

        loc = genomeLocParser.createGenomeLoc("chr1",7,7);
        eval_1_builder = new VariantContextBuilder();
        truth_1_builder = new VariantContextBuilder();
        eval_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        truth_1_builder.loc(loc.getContig(),loc.getStart(),loc.getStop());
        eval_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        truth_1_builder.alleles(Arrays.asList(reference_A,alt_C));
        eval_1_builder.genotypes(Arrays.asList(sam_2_1_5_eval,sam_2_2_5_eval));
        truth_1_builder.genotypes(Arrays.asList(sam_2_1_5_truth,sam_2_2_5_truth));

        Pair<VariantContext,VariantContext> testDataSite5 = new Pair<VariantContext, VariantContext>(eval_1_builder.make(),truth_1_builder.make());

        return Arrays.asList(testDataSite1,testDataSite2,testDataSite3,testDataSite4,testDataSite5);
    }

    @Test(enabled=true)
    public void testMultiSite() {
        int[][] sample1_expected = new int[GenotypeType.values().length][GenotypeType.values().length];
        int[][] sample2_expected = new int[GenotypeType.values().length][GenotypeType.values().length];
        // order: no-call,ref,het,hom-var,unavailable,mixed
        sample1_expected[0] = new int[]{0,1,0,0,0,0};
        sample2_expected[0] = new int[]{0,0,0,0,0,0};
        sample1_expected[1] = new int[]{0,1,0,0,0,0};
        sample2_expected[1] = new int[]{0,0,0,0,0,0};
        sample1_expected[2] = new int[]{0,0,0,0,0,0};
        sample2_expected[2] = new int[]{0,1,2,0,0,0};
        sample1_expected[3] = new int[]{1,0,0,0,0,0};
        sample2_expected[3] = new int[]{0,0,0,1,0,0};
        sample1_expected[4] = new int[]{0,0,1,0,0,0};
        sample2_expected[4] = new int[]{0,1,0,0,0,0};

        List<Pair<VariantContext,VariantContext>> data = getData6();

        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_2_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_2_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);

        for ( Pair<VariantContext,VariantContext> contextPair : data ) {
            VariantContext eval = contextPair.getFirst();
            VariantContext comp = contextPair.getSecond();
            logger.warn(eval.toString());
            logger.warn(comp.toString());
            Assert.assertTrue(eval != null);
            Assert.assertTrue(comp != null);
            Assert.assertTrue(eval.getGenotype("test2_sample1") != null);
            Assert.assertTrue(comp.getGenotype("test2_sample1") != null);
            Assert.assertTrue(eval.getGenotype("test2_sample2") != null);
            Assert.assertTrue(comp.getGenotype("test2_sample2") != null);
            metrics.update(eval,comp);
        }

        int[][] sample1_observed = metrics.getGenotypeConcordance("test2_sample1").getTable();
        int[][] sample2_observed = metrics.getGenotypeConcordance("test2_sample2").getTable();
        for ( GenotypeType eType : GenotypeType.values() ) {
            for ( GenotypeType cType : GenotypeType.values() ) {
                Assert.assertEquals(sample1_expected[eType.ordinal()][cType.ordinal()],sample1_observed[eType.ordinal()][cType.ordinal()]);
                Assert.assertEquals(sample2_expected[eType.ordinal()][cType.ordinal()],sample2_observed[eType.ordinal()][cType.ordinal()]);
            }
        }
    }

    @Test(enabled=true)
    public void testNRD_testNRS_testMargins() {
        Pair<VariantContext,VariantContext> data = getData3();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        int[][] table = metrics.getOverallGenotypeConcordance().getTable();
        // set up the table
        table[0] = new int[] {30, 12, 7, 5, 6, 0};
        table[1] = new int[] {10, 100, 5, 1, 7, 1};
        table[2] = new int[] {5, 7, 150, 3, 3, 1};
        table[3] = new int[] {3, 2, 6, 50, 1, 0};
        table[4] = new int[] {10, 6, 3, 3, 2, 0};
        table[5] = new int[] {12, 0, 34, 20, 10, 0};
        double EXPEC_NRS = 0.8969957;
        double EXPEC_NRD = 0.1071429;
        double EXPEC_OGC = 0.92592592; // (100+150+50)/(100+5+1+150+7+3+50+2+6)
        Assert.assertEquals(EXPEC_NRS,metrics.getOverallNRS(),1e-7);
        Assert.assertEquals(EXPEC_NRD,metrics.getOverallNRD(),1e-7);
        Assert.assertEquals(EXPEC_OGC,metrics.getOverallOGC(),1e-7);
        int EXPEC_EVAL_REF = 124;
        int EXPEC_EVAL_HET = 169;
        int EXPEC_EVAL_VAR = 62;
        int EXPEC_COMP_REF = 127;
        int EXPEC_COMP_HET = 205;
        int EXPEC_COMP_VAR = 82;
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getnEvalGenotypes(GenotypeType.HOM_REF),EXPEC_EVAL_REF);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getnEvalGenotypes(GenotypeType.HET),EXPEC_EVAL_HET);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getnEvalGenotypes(GenotypeType.HOM_VAR),EXPEC_EVAL_VAR);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getnCompGenotypes(GenotypeType.HOM_REF),EXPEC_COMP_REF);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getnCompGenotypes(GenotypeType.HET),EXPEC_COMP_HET);
        Assert.assertEquals(metrics.getOverallGenotypeConcordance().getnCompGenotypes(GenotypeType.HOM_VAR),EXPEC_COMP_VAR);
    }

    @Test(enabled=true)
    public void testRobustness() {
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_3_HEADER_1))));
        VCFHeader disjointCompHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_3_HEADER_2))));
        VCFHeader overlapCompHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_3_HEADER_3))));
        ConcordanceMetrics disjointMetrics = new ConcordanceMetrics(evalHeader,disjointCompHeader,null);
        ConcordanceMetrics overlapMetrics = new ConcordanceMetrics(evalHeader,overlapCompHeader,null);

        // test what happens if you put in disjoint sets and start making requests
        Assert.assertEquals(0,disjointMetrics.getPerSampleGenotypeConcordance().size());
        String msg = "No Exception Thrown";
        try {
            disjointMetrics.getGenotypeConcordance("test3_sample4");
        } catch ( Exception e) {
            msg = e.getMessage();
        }
        Assert.assertEquals("Attempted to request the concordance table for sample test3_sample4 on which it was not calculated",msg);

        // test that the overlapping sample is in the overlapping table (basically do this without throwing an exception)
        overlapMetrics.getGenotypeConcordance("test3_sample3");

        String msg2 = "No Exception Thrown";
        try {
            disjointMetrics.getGenotypeConcordance("test3_sample4");
        } catch ( Exception e) {
            msg2 = e.getMessage();
        }
        Assert.assertEquals("Attempted to request the concordance table for sample test3_sample4 on which it was not calculated",msg2);

        // test what happens if you try to calculate NRS and NRD on an empty table
        Assert.assertEquals(disjointMetrics.getOverallNRD(), 1.0, 1e-16);
        Assert.assertEquals(disjointMetrics.getOverallNRS(), 0.0, 1e-16);
    }

    public List<Pair<VariantContext,VariantContext>> getData7() {

        Allele ref1 = Allele.create(BaseUtils.Base.T.base,true);
        Allele alt1 = Allele.create(BaseUtils.Base.C.base);
        Allele alt2 = Allele.create(BaseUtils.Base.G.base);
        Allele alt3 = Allele.create(BaseUtils.Base.A.base);

        GenomeLoc loc1 = genomeLocParser.createGenomeLoc("chr1",1,1);
        VariantContextBuilder site1Eval = new VariantContextBuilder();
        VariantContextBuilder site1Comp = new VariantContextBuilder();


        // site 1: eval superset comp
        site1Eval.loc(loc1.getContig(),loc1.getStart(),loc1.getStop());
        site1Comp.loc(loc1.getContig(),loc1.getStart(),loc1.getStop());
        site1Eval.alleles(Arrays.asList(ref1,alt1,alt2));
        site1Comp.alleles(Arrays.asList(ref1,alt2));
        site1Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt2)));
        site1Comp.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt2)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt2)));

        // site 2: eval subset comp
        GenomeLoc loc2 = genomeLocParser.createGenomeLoc("chr1",2,2);
        VariantContextBuilder site2Eval = new VariantContextBuilder();
        VariantContextBuilder site2Comp = new VariantContextBuilder();
        site2Eval.loc(loc2.getContig(),loc2.getStart(),loc2.getStop());
        site2Comp.loc(loc2.getContig(),loc2.getStart(),loc2.getStop());
        site2Eval.alleles(Arrays.asList(ref1,alt1));
        site2Comp.alleles(Arrays.asList(ref1,alt1,alt3));
        site2Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt1)));
        site2Comp.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt3)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt1)));

        // site 3: eval only
        GenomeLoc loc3 = genomeLocParser.createGenomeLoc("chr1",3,3);
        VariantContextBuilder site3Eval = new VariantContextBuilder();
        VariantContextBuilder site3Comp = new VariantContextBuilder();
        site3Eval.loc(loc3.getContig(),loc3.getStart(),loc3.getStop());
        site3Comp.loc(loc3.getContig(),loc3.getStart(),loc3.getStop());
        site3Eval.alleles(Arrays.asList(ref1,alt1));
        site3Comp.alleles(Arrays.asList(ref1,alt1));
        site3Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt1)));
        site3Comp.genotypes(GenotypeBuilder.create("test2_sample1",new ArrayList<Allele>(0)),GenotypeBuilder.create("test2_sample2",new ArrayList<Allele>(0)));

        // site 4: comp only - monomorphic
        GenomeLoc loc4 = genomeLocParser.createGenomeLoc("chr1",4,4);
        VariantContextBuilder site4Eval = new VariantContextBuilder();
        VariantContextBuilder site4Comp = new VariantContextBuilder();
        site4Eval.loc(loc4.getContig(),loc4.getStart(),loc4.getStop());
        site4Comp.loc(loc4.getContig(),loc4.getStart(),loc4.getStop());
        site4Eval.alleles(Arrays.asList(ref1,alt1));
        site4Comp.alleles(Arrays.asList(ref1,alt1));
        site4Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,ref1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,ref1)));
        site4Comp.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt1)));

        // site 5: overlapping
        GenomeLoc loc5 = genomeLocParser.createGenomeLoc("chr1",5,5);
        VariantContextBuilder site5Eval = new VariantContextBuilder();
        VariantContextBuilder site5Comp = new VariantContextBuilder();
        site5Eval.loc(loc5.getContig(),loc5.getStart(),loc5.getStop());
        site5Comp.loc(loc5.getContig(),loc5.getStart(),loc5.getStop());
        site5Eval.alleles(Arrays.asList(ref1,alt1,alt3));
        site5Comp.alleles(Arrays.asList(ref1,alt1,alt3));
        site5Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(alt1,alt3)));
        site5Comp.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(alt1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(alt3,alt3)));

        // site 6: some non-matching alts
        GenomeLoc loc6 = genomeLocParser.createGenomeLoc("chr1",6,6);
        VariantContextBuilder site6Eval = new VariantContextBuilder();
        VariantContextBuilder site6Comp = new VariantContextBuilder();
        site6Eval.loc(loc6.getContig(),loc6.getStart(),loc6.getStop());
        site6Comp.loc(loc6.getContig(),loc6.getStart(),loc6.getStop());
        site6Eval.alleles(Arrays.asList(ref1,alt1,alt2));
        site6Comp.alleles(Arrays.asList(ref1,alt1,alt3));
        site6Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt2)));
        site6Comp.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt3)));

        // site 7: matching with no-calls
        GenomeLoc loc7 = genomeLocParser.createGenomeLoc("chr1",7,7);
        VariantContextBuilder site7Eval = new VariantContextBuilder();
        VariantContextBuilder site7Comp = new VariantContextBuilder();
        site7Eval.loc(loc7.getContig(),loc7.getStart(),loc7.getStop());
        site7Comp.loc(loc7.getContig(),loc7.getStart(),loc7.getStop());
        site7Eval.alleles(Arrays.asList(ref1,alt1));
        site7Comp.alleles(Arrays.asList(ref1,alt1));
        site7Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(Allele.NO_CALL,Allele.NO_CALL)));
        site7Comp.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(ref1,alt1)),GenotypeBuilder.create("test2_sample2",Arrays.asList(ref1,alt1)));

        Pair<VariantContext,VariantContext> site1 = new Pair<VariantContext, VariantContext>(site1Eval.make(),site1Comp.make());
        Pair<VariantContext,VariantContext> site2 = new Pair<VariantContext, VariantContext>(site2Eval.make(),site2Comp.make());
        Pair<VariantContext,VariantContext> site3 = new Pair<VariantContext, VariantContext>(site3Eval.make(),site3Comp.make());
        Pair<VariantContext,VariantContext> site4 = new Pair<VariantContext, VariantContext>(site4Eval.make(),site4Comp.make());
        Pair<VariantContext,VariantContext> site5 = new Pair<VariantContext, VariantContext>(site5Eval.make(),site5Comp.make());
        Pair<VariantContext,VariantContext> site6 = new Pair<VariantContext, VariantContext>(site6Eval.make(),site6Comp.make());
        Pair<VariantContext,VariantContext> site7 = new Pair<VariantContext, VariantContext>(site7Eval.make(),site7Comp.make());

        return Arrays.asList(site1,site2,site3,site4,site5,site6,site7);
    }

    @Test(enabled = true)
    public void testSites() {
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_2_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_2_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);

        List<Pair<VariantContext,VariantContext>> data = getData7();

        int idx = 0;
        int[] expecNotMatch = new int[]{0,0,0,0,0,1,1};
        for ( Pair<VariantContext,VariantContext> varPair : data ) {
            metrics.update(varPair.getFirst(),varPair.getSecond());
            Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.ALLELES_DO_NOT_MATCH),expecNotMatch[idx]);
            logger.info(idx);
            idx++;
        }

        Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.ALLELES_DO_NOT_MATCH),1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.ALLELES_MATCH),2);
        Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.EVAL_ONLY),1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.TRUTH_ONLY),1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.EVAL_SUBSET_TRUTH),1);
        Assert.assertEquals(metrics.getOverallSiteConcordance().get(ConcordanceMetrics.SiteConcordanceType.EVAL_SUPERSET_TRUTH),1);

    }

    private Pair<VariantContext,VariantContext> getMonoallelicData() {

        final Allele ref = Allele.create(BaseUtils.Base.T.base,true);
        final Allele alt = Allele.create(BaseUtils.Base.C.base);

        //Site in eval is monoallelic, both samples are HOM_REF
        //sample1 in comp is HOM_VAR, sample2 is NO_CALL
        //None of these should trigger mismatching alleles
        final GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1",1,1);
        final VariantContextBuilder site1Comp = new VariantContextBuilder();
        final VariantContextBuilder site1Eval = new VariantContextBuilder();
        site1Comp.loc(loc.getContig(), loc.getStart(), loc.getStop());
        site1Eval.loc(loc.getContig(), loc.getStart(), loc.getStop());
        site1Comp.alleles(Arrays.asList(ref));
        site1Eval.alleles(Arrays.asList(ref, alt));
        site1Comp.genotypes(GenotypeBuilder.create("test2_sample1", Arrays.asList(ref, ref)),
                GenotypeBuilder.create("test2_sample2", Arrays.asList(ref, ref)));
        site1Eval.genotypes(GenotypeBuilder.create("test2_sample1",Arrays.asList(alt,alt)),
                GenotypeBuilder.create("test2_sample2",Arrays.asList(Allele.NO_CALL,Allele.NO_CALL)));

        return new Pair<>(site1Eval.make(), site1Comp.make());
    }

    @Test
    public void testMonoallelicSite() {
        final Pair<VariantContext,VariantContext> data = getMonoallelicData();
        final VariantContext eval = data.getFirst();
        final VariantContext truth = data.getSecond();
        final VCFCodec codec = new VCFCodec();
        final VCFHeader evalHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_2_HEADER))));
        final VCFHeader compHeader = (VCFHeader)codec.readActualHeader(codec.makeSourceFromStream(new PositionalBufferedStream(new StringBufferInputStream(TEST_2_HEADER))));
        final ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader,null);
        metrics.update(eval,truth);


        Assert.assertEquals(metrics.getGenotypeConcordance("test2_sample1").getnMismatchingAlt(),0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test2_sample2").getnMismatchingAlt(),0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test2_sample1").getTable()[3][1],1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test2_sample2").getTable()[0][1],1);
    }
}