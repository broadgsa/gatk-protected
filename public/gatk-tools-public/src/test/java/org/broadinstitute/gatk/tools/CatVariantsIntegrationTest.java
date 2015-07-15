/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.MD5DB;
import org.broadinstitute.gatk.utils.MD5Mismatch;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.runtime.ProcessController;
import org.broadinstitute.gatk.utils.runtime.ProcessSettings;
import org.broadinstitute.gatk.utils.runtime.RuntimeUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class CatVariantsIntegrationTest {
    private final MD5DB md5db = new MD5DB();
    private final File CatVariantsDir = new File(BaseTest.privateTestDir, "CatVariants");
    private final File CatVariantsVcf1 = new File(CatVariantsDir, "CatVariantsTest1.vcf");
    private final File CatVariantsVcf2 = new File(CatVariantsDir, "CatVariantsTest2.vcf");
    private final File CatVariantsBcf1 = new File(CatVariantsDir, "CatVariantsTest1.bcf");
    private final File CatVariantsBcf2 = new File(CatVariantsDir, "CatVariantsTest2.bcf");

    private class CatVariantsTestProvider extends BaseTest.TestDataProvider {
        private final File file1;
        private final File file2;
        public final File outputFile;
        public final String md5;

        private CatVariantsTestProvider(final File file1, final File file2, final File outputFile, final String md5) {
            super(CatVariantsTestProvider.class);

            this.file1 = file1;
            this.file2 = file2;
            this.outputFile = outputFile;
            this.md5 = md5;
        }

        public final String getCmdLine() {
            return String.format("java -cp \"%s\" %s -R %s -V %s -V %s -out %s",
                    StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                    CatVariants.class.getCanonicalName(), BaseTest.b37KGReference, file1, file2, outputFile);
        }

        public String toString() {
            return String.format("CatVariantsTestProvider %s + %s -> %s", file1.getName(), file2.getName(), outputFile.getName());
        }
    }

    @DataProvider(name = "ExtensionsTest")
    public Object[][] makeExtensionsTestProvider() {
        final File catVariantsTempList1 = BaseTest.createTempListFile("CatVariantsTest1", CatVariantsVcf1.getAbsolutePath());
        final File catVariantsTempList2 = BaseTest.createTempListFile("CatVariantsTest2", CatVariantsVcf2.getAbsolutePath());

        new CatVariantsTestProvider(CatVariantsVcf1, CatVariantsVcf2, BaseTest.createTempFile("CatVariantsTest", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider(CatVariantsBcf1, CatVariantsBcf2, BaseTest.createTempFile("CatVariantsTest", ".bcf"), "6a57fcbbf3cae490896d13a288670d83");

        for (String extension1 : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS) {
            for (String extension2 : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS) {
                final File file1 = new File(CatVariantsDir, "CatVariantsTest1.vcf" + extension1);
                final File file2 = new File(CatVariantsDir, "CatVariantsTest2.vcf" + extension2);
                new CatVariantsTestProvider(file1, file2, BaseTest.createTempFile("CatVariantsTest.", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
                new CatVariantsTestProvider(file1, file2, BaseTest.createTempFile("CatVariantsTest.", ".bcf"), "6a57fcbbf3cae490896d13a288670d83");
                new CatVariantsTestProvider(file1, file2, BaseTest.createTempFile("CatVariantsTest.", ".vcf" + extension1), "33f728ac5c70ce2994f3619a27f47088");
            }
            new CatVariantsTestProvider(CatVariantsVcf1, CatVariantsVcf2, BaseTest.createTempFile("CatVariantsTest.", ".vcf" + extension1), "33f728ac5c70ce2994f3619a27f47088");
            new CatVariantsTestProvider(CatVariantsBcf1, CatVariantsBcf2, BaseTest.createTempFile("CatVariantsTest.", ".vcf" + extension1), "f1a55575f59707f80b8c17e2591fbf53");
        }

        //Test list parsing functionality
        new CatVariantsTestProvider(catVariantsTempList1, CatVariantsVcf2, BaseTest.createTempFile("CatVariantsTest.", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider(CatVariantsVcf1, catVariantsTempList2, BaseTest.createTempFile("CatVariantsTest.", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");
        new CatVariantsTestProvider(catVariantsTempList1, catVariantsTempList2, BaseTest.createTempFile("CatVariantsTest.", ".vcf"), "d0d81eb7fd3905256c4ac7c0fc480094");

        return CatVariantsTestProvider.getTests(CatVariantsTestProvider.class);
    }

    @Test(dataProvider = "ExtensionsTest")
    public void testExtensions(final CatVariantsTestProvider cfg) throws IOException {

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(Utils.escapeExpressions(cfg.getCmdLine()));
        pc.execAndCheck(ps);

        MD5DB.MD5Match result = md5db.testFileMD5("testExtensions", "CatVariantsTestProvider", cfg.outputFile, cfg.md5, false);
        if(result.failed) {
            final MD5Mismatch failure = new MD5Mismatch(result.actualMD5, result.expectedMD5, result.diffEngineOutput);
            Assert.fail(failure.toString());
        }
    }

    @DataProvider(name = "MismatchedExtensionsTest")
    public Object[][] makeMismatchedExtensionsTestProvider() {
        return new Object[][]{
                {".vcf", ".vcf.gz"},
                {".vcf.gz", ".vcf"},
                {".bcf", ".vcf.gz"},
                {".vcf.gz", ".bcf"},
                {".vcf", ".bcf"},
                {".bcf", ".vcf"}
        };
    }

    @Test(dataProvider = "MismatchedExtensionsTest", expectedExceptions = IOException.class)
    public void testMismatchedExtensions1(final String extension1, final String extension2) throws IOException {
        String cmdLine = String.format("java -cp \"%s\" %s -R %s -V %s -V %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                new File(CatVariantsDir, "CatVariantsTest1" + extension1),
                new File(CatVariantsDir, "CatVariantsTest2" + extension2),
                BaseTest.createTempFile("CatVariantsTest", ".bcf"));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(Utils.escapeExpressions(cmdLine));
        pc.execAndCheck(ps);
    }

    @Test(dataProvider = "MismatchedExtensionsTest", expectedExceptions = IOException.class)
    public void testMismatchedExtensions2(final String extension1, final String extension2) throws IOException {

        String cmdLine = String.format("java -cp \"%s\" %s -R %s -V %s -V %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                new File(CatVariantsDir, "CatVariantsTest1" + extension1),
                new File(CatVariantsDir, "CatVariantsTest2" + extension2),
                BaseTest.createTempFile("CatVariantsTest", ".vcf"));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(Utils.escapeExpressions(cmdLine));
        pc.execAndCheck(ps);
    }

    @DataProvider(name = "SortingTest")
    public Object[][] makeVariantSortingTestProvider() {
    	// Make all of our variants REF==A and ALT==T
    	Allele Aref = Allele.create("A", true);
    	Allele Talt = Allele.create("T");
    	List<Allele> varAlleles = Arrays.asList(Aref, Talt);
    	
    	// First, let's make a couple of variants
    	VariantContext vc1_10 = new VariantContextBuilder("Source", "1", 10, 10, varAlleles).make();
    	VariantContext vc1_30 = new VariantContextBuilder("Source", "1", 30, 30, varAlleles).make();
    	VariantContext vc2_20 = new VariantContextBuilder("Source", "2", 20, 20, varAlleles).make();
    	
    	File vc1_10_fn = BaseTest.createTempFile("CatVariantsSortingTest1", "vcf");
    	File vc1_30_fn = BaseTest.createTempFile("CatVariantsSortingTest1", "vcf");
    	File vc2_20_fn = BaseTest.createTempFile("CatVariantsSortingTest1", "vcf");
    	
    	// And now, write exactly one variant per file
    	VariantContextWriterBuilder builder = new VariantContextWriterBuilder();

    	VCFHeader stdHeader = VCFStandardHeaderLines.repairStandardHeaderLines(new VCFHeader());
    	
        VariantContextWriter vc1_10_w = builder
            .setOutputFile(vc1_10_fn)
            .build();
        vc1_10_w.writeHeader(stdHeader);
        vc1_10_w.add(vc1_10);
        vc1_10_w.close();
       
        VariantContextWriter vc1_30_w = builder
            .setOutputFile(vc1_30_fn)
            .build();
        vc1_30_w.writeHeader(stdHeader);
        vc1_30_w.add(vc1_30);
        vc1_30_w.close();
       
        VariantContextWriter vc2_20_w = builder
            .setOutputFile(vc2_20_fn)
            .build();
        vc2_20_w.writeHeader(stdHeader);
        vc2_20_w.add(vc2_20);
        vc2_20_w.close();
    	
    	// and return them in an unsorted order to make sure the sorting works as advertised
    	return new Object[][]{ 
    			{Arrays.asList(vc1_30_fn, vc2_20_fn, vc1_10_fn), 
    		     Arrays.asList(vc1_10, vc1_30, vc2_20) } 
    		};        
    }
    
    @Test(dataProvider = "SortingTest")    
    public void testVariantSorting(final List<File> fileInput, final List<VariantContext> sortedVariants) throws IOException {
        String inputCmd = "-V " + StringUtils.join(fileInput, " -V ");
        
        File combinedFile = BaseTest.createTempFile("CatVariantsSortingTestCombo", ".vcf"); 
        	
        String cmdLine = String.format("java -cp \"%s\" %s -R %s %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                inputCmd,
                combinedFile);

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(Utils.escapeExpressions(cmdLine));
        pc.execAndCheck(ps);
	
        // Now, we need to check and ensure that the chromosome and position of 
        // each variant read in the combined file is exactly what we expect in 
        // the sortedVariants input
        VCFFileReader comboReader = new VCFFileReader(combinedFile);
        Iterator<VariantContext> sortedIter = sortedVariants.iterator();
        for(VariantContext vc : comboReader){
        	// we should have the same number of variants in the combined file
        	// as the list passed to us
        	Assert.assertTrue(sortedIter.hasNext());
        	VariantContext nextSortvc = sortedIter.next();
        	
        	// we should get the same contig and start position for each
        	// variant read!
        	Assert.assertEquals(vc.getContig(), nextSortvc.getContig());
        	Assert.assertEquals(vc.getStart(), nextSortvc.getStart());
        }
        
        // finally, we should NOT have any variants left in the sortedVariants list
        Assert.assertFalse(sortedIter.hasNext());
                
    }
    
    //
    //
    // IndexCreator tests
    //
    //

    private class VCFIndexCreatorTest extends BaseTest.TestDataProvider {
        private final GATKVCFIndexType type;
        private final int parameter;

        private VCFIndexCreatorTest(GATKVCFIndexType type, int parameter) {
            super(VCFIndexCreatorTest.class);

            this.type = type;
            this.parameter = parameter;
        }

        public String toString() {
            return String.format("Index Type %s, Index Parameter %s", type, parameter);
        }

        public Index getIndex(final File vcfFile) {
            switch (type) {
                case DYNAMIC_SEEK : return IndexFactory.createDynamicIndex(vcfFile, new VCFCodec(), IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
                case DYNAMIC_SIZE : return IndexFactory.createDynamicIndex(vcfFile, new VCFCodec(), IndexFactory.IndexBalanceApproach.FOR_SIZE);
                case LINEAR : return IndexFactory.createLinearIndex(vcfFile, new VCFCodec(), parameter);
                case INTERVAL : return IndexFactory.createIntervalIndex(vcfFile, new VCFCodec(), parameter);
                default : throw new TestException("Invalid index type");
            }
        }
    }

    @DataProvider(name = "IndexDataProvider")
    public Object[][] indexCreatorData() {
        new VCFIndexCreatorTest(GATKVCFIndexType.DYNAMIC_SEEK, 0);
        new VCFIndexCreatorTest(GATKVCFIndexType.DYNAMIC_SIZE, 0);
        new VCFIndexCreatorTest(GATKVCFIndexType.LINEAR, 100);
        new VCFIndexCreatorTest(GATKVCFIndexType.LINEAR, 10000);
        new VCFIndexCreatorTest(GATKVCFIndexType.INTERVAL, 20);
        new VCFIndexCreatorTest(GATKVCFIndexType.INTERVAL, 2000);

        return BaseTest.TestDataProvider.getTests(VCFIndexCreatorTest.class);
    }

    @Test(dataProvider = "IndexDataProvider")
    public void testCatVariantsVCFIndexCreation(VCFIndexCreatorTest testSpec) throws IOException{

        String cmdLine = String.format("java -cp \"%s\" %s -R %s -V %s -V %s --variant_index_type %s --variant_index_parameter %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                CatVariantsVcf1,
                CatVariantsVcf2,
                testSpec.type,
                testSpec.parameter,
                BaseTest.createTempFile("CatVariantsVCFIndexCreationTest", ".vcf"));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(Utils.escapeExpressions(cmdLine));
        pc.execAndCheck(ps);
    }

    @Test()
    public void testCatVariantsGVCFIndexCreation() throws IOException{

        String cmdLine = String.format("java -cp \"%s\" %s -R %s -V %s -V %s -out %s",
                StringUtils.join(RuntimeUtils.getAbsoluteClassPaths(), File.pathSeparatorChar),
                CatVariants.class.getCanonicalName(),
                BaseTest.b37KGReference,
                CatVariantsVcf1,
                CatVariantsVcf2,
                BaseTest.createTempFile("CatVariantsGVCFIndexCreationTest", "." + GATKVCFUtils.GVCF_EXT));

        ProcessController pc = ProcessController.getThreadLocal();
        ProcessSettings ps = new ProcessSettings(Utils.escapeExpressions(cmdLine));
        pc.execAndCheck(ps);
    }
}