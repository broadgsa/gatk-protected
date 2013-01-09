package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.variantutils.ConcordanceMetrics;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.GenotypeType;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.testng.annotations.Test;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.vcf.*;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.StringBufferInputStream;
import java.util.ArrayList;
import java.util.Set;
import java.util.Arrays;
import java.util.List;
import net.sf.picard.reference.ReferenceSequenceFile;

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


    private Pair<VariantContext,VariantContext> getData1() {

        Allele reference_A = Allele.create(BaseUtils.A,true);
        Allele alt_C = Allele.create(BaseUtils.C);

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
        VCFHeader evalHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader);
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

        Allele reference_A = Allele.create(BaseUtils.A,true);
        Allele alt_C = Allele.create(BaseUtils.C);
        Allele alt_T = Allele.create(BaseUtils.T);

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
    public void testMismatchingAllele() {
        Pair<VariantContext,VariantContext> data = getData2();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader);
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
    }

    private Pair<VariantContext,VariantContext> getData3() {

        Allele reference_ACT = Allele.create(new byte[]{BaseUtils.A,BaseUtils.C,BaseUtils.T},true);
        Allele alt_AC = Allele.create(new byte[]{BaseUtils.A,BaseUtils.C});
        Allele alt_A = Allele.create(BaseUtils.A);
        Allele alt_ATT = Allele.create(new byte[]{BaseUtils.A,BaseUtils.T,BaseUtils.T});

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
        VCFHeader evalHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader);
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

        Allele reference_A = Allele.create(BaseUtils.A,true);
        Allele alt_C = Allele.create(BaseUtils.C);
        Allele alt_T = Allele.create(BaseUtils.T);

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
        VCFHeader evalHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader);
        metrics.update(eval,truth);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getnMismatchingAlt(),0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample2").getTable()[0][2],1);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][1],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][3],0);
        Assert.assertEquals(metrics.getGenotypeConcordance("test1_sample3").getTable()[2][0],1);
    }

    private Pair<VariantContext,VariantContext> getData5() {

        Allele reference_A = Allele.create(BaseUtils.A,true);
        Allele alt_C = Allele.create(BaseUtils.C);
        Allele alt_T = Allele.create(BaseUtils.T);

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
        VCFHeader evalHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader);
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

    @Test(enabled=true)
    public void testNRD_testNRS() {
        Pair<VariantContext,VariantContext> data = getData3();
        VariantContext eval = data.getFirst();
        VariantContext truth = data.getSecond();
        VCFCodec codec = new VCFCodec();
        VCFHeader evalHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        VCFHeader compHeader = (VCFHeader)codec.readHeader(new AsciiLineReader(new PositionalBufferedStream(new StringBufferInputStream(TEST_1_HEADER))));
        ConcordanceMetrics metrics = new ConcordanceMetrics(evalHeader,compHeader);
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
        Assert.assertEquals(EXPEC_NRS,metrics.getOverallNRS(),1e-7);
        Assert.assertEquals(EXPEC_NRD,metrics.getOverallNRD(),1e-7);
    }
}