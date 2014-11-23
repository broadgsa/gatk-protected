package org.broadinstitute.gatk.tools.walkers.phasing;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.broadinstitute.gatk.tools.walkers.phasing.PhasingUtils;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Created by ronlevine and philipdexheimer on 11/22/14.
 *
 * TODO - need to exercise the code paths in matchHaplotypeAlleles and reallyMergeIntoMNP
 */
public class PhasingUtilitiesUnitTest extends BaseTest  {

    private final int start = 10;
    private ReferenceSequenceFile referenceFile;
    private Allele Aref;

    @BeforeClass
    public void init() throws FileNotFoundException {
        referenceFile = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        Aref = Allele.create("A", true);
    }

    @Test
    public void TestMatchHaplotypeAlleles() {

        final List<Allele> noCalls = new ArrayList<>(1);
        noCalls.add(Allele.NO_CALL);
        final Genotype gA_ALT = new GenotypeBuilder("A").PL(new int[]{0, 100, 1000}).alleles(noCalls).make();

        // Must match because that same genotype
        PhasingUtils.SameHaplotypeAlleles sameHaplotypeAlleles = PhasingUtils.matchHaplotypeAlleles(gA_ALT, gA_ALT);

        // TODO - make a reference SameHaplotypeAlleles to test against
        Assert.assertEquals(sameHaplotypeAlleles, sameHaplotypeAlleles);
    }

    @Test
    public void TestReallyMergeIntoMNP() {

        final VariantContext VCbase = new VariantContextBuilder("test", "20", start, start, Arrays.asList(Aref)).make();

        // Merging the same variant context
        VariantContext variantContext = PhasingUtils.reallyMergeIntoMNP(VCbase, VCbase, referenceFile);

        // TODO - make a reference VariantContext to test against
        Assert.assertEquals(variantContext, variantContext);
    }
}
