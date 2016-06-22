/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 1:12:35 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Template for testing reference views (ReadReferenceView and LocusReferenceView).
 */

public abstract class ReferenceViewTemplate extends BaseTest {
    /**
     * The fasta, for comparison.
     */
    protected ReferenceSequenceFile sequenceFile = null;
    protected GenomeLocParser genomeLocParser = null;

    //
    // The bulk of sequence retrieval is tested by IndexedFastaSequenceFile, but we'll run a few spot
    // checks here to make sure that data is flowing through the LocusReferenceView.

    /**
     * Initialize the fasta.
     */
    @BeforeClass
    public void initialize() throws FileNotFoundException {
        sequenceFile = new CachingIndexedFastaSequenceFile( new File(hg18Reference) );
        genomeLocParser = new GenomeLocParser(sequenceFile);
    }

    /**
     * Test the initial fasta location.
     */
    @Test
    public void testReferenceStart() {
        validateLocation( genomeLocParser.createGenomeLoc(sequenceFile.getSequenceDictionary().getSequence(0).getSequenceName(),1,25) );
    }

    /**
     * Test the end of a contig.
     */
    @Test
    public void testReferenceEnd() {
        // Test the last 25 bases of the first contig.
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(sequenceFile.getSequenceDictionary().getSequences().size()-1);
        final int contigStart = selectedContig.getSequenceLength() - 24;
        final int contigStop = selectedContig.getSequenceLength();
        validateLocation( genomeLocParser.createGenomeLoc(selectedContig.getSequenceName(),contigStart,contigStop) );
    }

    /**
     * Test the start of the middle contig.
     */
    @Test
    public void testContigStart() {
        // Test the last 25 bases of the first contig.
        int contigPosition = sequenceFile.getSequenceDictionary().getSequences().size()/2;
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(contigPosition);
        validateLocation( genomeLocParser.createGenomeLoc(selectedContig.getSequenceName(),1,25) );
    }


    /**
     * Test the end of the middle contig.
     */
    @Test
    public void testContigEnd() {
        // Test the last 25 bases of the first contig.
        int contigPosition = sequenceFile.getSequenceDictionary().getSequences().size()/2;
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(contigPosition);
        final int contigStart = selectedContig.getSequenceLength() - 24;
        final int contigStop = selectedContig.getSequenceLength();
        validateLocation( genomeLocParser.createGenomeLoc(selectedContig.getSequenceName(),contigStart,contigStop) );
    }

    protected abstract void validateLocation( GenomeLoc loc );
}
