package org.broadinstitute.sting.gatk.datasources.providers;

import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.LocusShard;
import org.broadinstitute.sting.gatk.iterators.GenomeLocusIterator;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

import java.util.Collections;
/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/** Tests for viewing the reference from the perspective of a locus. */

public class LocusReferenceViewUnitTest extends ReferenceViewTemplate {

//
//    /** Multiple-base pair queries should generate exceptions. */
//    @Test(expected = InvalidPositionException.class)
//    public void testSingleBPFailure() {
//        Shard shard = new LocusShard(GenomeLocParser.createGenomeLoc(0, 1, 50));
//
//        ShardDataProvider dataProvider = new ShardDataProvider(shard, null, sequenceFile, null);
//        LocusReferenceView view = new LocusReferenceView(dataProvider);
//
//        view.getReferenceContext(shard.getGenomeLoc()).getBase();
//    }

    @Test
    public void testOverlappingReferenceBases() {
        Shard shard = new LocusShard(Collections.singletonList(GenomeLocParser.createGenomeLoc(0, sequenceFile.getSequence("chrM").length() - 10, sequenceFile.getSequence("chrM").length())));
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, null, shard.getGenomeLocs().get(0), null, sequenceFile, null);
        LocusReferenceView view = new LocusReferenceView(dataProvider);

        char[] results = view.getReferenceBases(GenomeLocParser.createGenomeLoc(0, sequenceFile.getSequence("chrM").length() - 10, sequenceFile.getSequence("chrM").length() + 9));
        Assert.assertEquals(20, results.length);
        for (int x = 0; x < results.length; x++) {
            if (x <= 10) Assert.assertTrue(results[x] != 'X');
            else Assert.assertTrue(results[x] == 'X');
        }
    }


    /** Queries outside the bounds of the shard should result in reference context window trimmed at the shard boundary. */
    @Test
    public void testBoundsFailure() {
        Shard shard = new LocusShard(Collections.singletonList(GenomeLocParser.createGenomeLoc(0, 1, 50)));

        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, null, shard.getGenomeLocs().get(0), null, sequenceFile, null);
        LocusReferenceView view = new LocusReferenceView(dataProvider);

        GenomeLoc locus = GenomeLocParser.createGenomeLoc(0, 50, 51);

        ReferenceContext rc = view.getReferenceContext(locus);
        Assert.assertTrue(rc.getLocus().equals(locus));
        Assert.assertTrue(rc.getWindow().equals(GenomeLocParser.createGenomeLoc(0,50)));
        Assert.assertTrue(rc.getBases().length == 1);
    }


    /**
     * Compares the contents of the fasta and view at a specified location.
     *
     * @param loc
     */
    protected void validateLocation( GenomeLoc loc ) {
        Shard shard = new LocusShard(Collections.singletonList(loc));
        GenomeLocusIterator shardIterator = new GenomeLocusIterator(loc);

        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, null, loc, null, sequenceFile, null);
        LocusReferenceView view = new LocusReferenceView(dataProvider);

        while (shardIterator.hasNext()) {
            GenomeLoc locus = shardIterator.next();

            ReferenceSequence expectedAsSeq = sequenceFile.getSubsequenceAt(locus.getContig(), locus.getStart(), locus.getStop());
            char expected = StringUtil.bytesToString(expectedAsSeq.getBases()).charAt(0);
            char actual = view.getReferenceContext(locus).getBase();

            Assert.assertEquals(String.format("Value of base at position %s in shard %s does not match expected", locus.toString(), shard.getGenomeLocs()),
                    expected,
                    actual);
        }
    }

}