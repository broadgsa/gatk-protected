package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class GATKBAMIndexFromDataSource extends GATKBAMIndex {
    private File sourceFile;
    private SAMFileHeader sourceHeader;
    private BrowseableBAMIndex index;

    public GATKBAMIndexFromDataSource(File sourceFile, SAMFileHeader sourceHeader, BrowseableBAMIndex index) {
        this.sourceFile = sourceFile;
        this.sourceHeader = sourceHeader;
        this.index = index;
    }

    @Override
    public GATKBAMIndexData readReferenceSequence(int referenceSequence) {
        List<SAMSequenceRecord> sequences = sourceHeader.getSequenceDictionary().getSequences();
        if (referenceSequence >= sequences.size())
            throw new ReviewedGATKException("Invalid sequence number " + referenceSequence + " in index file " + sourceFile);


        BinList sourceBins = index.getBinsOverlapping(referenceSequence, 0, sequences.get(referenceSequence).getSequenceLength());

        List<GATKBin> bins = new ArrayList<GATKBin>();
        for (Bin sourceBin : sourceBins) {
            int indexBin = sourceBin.getBinNumber();
            while(indexBin >= bins.size())
                bins.add(null);

            GATKBin bin = new GATKBin(referenceSequence, indexBin);
            List<Chunk> chunks = index.getSpanOverlapping(sourceBin).getChunks();
            List<GATKChunk> gatkChunks = new ArrayList<>(chunks.size());
            for (Chunk chunk : chunks) {
                gatkChunks.add(new GATKChunk(chunk));
            }

            bin.setChunkList(gatkChunks.toArray(new GATKChunk[gatkChunks.size()]));

            bins.set(indexBin, bin);
        }

        // there is no interface to get linear index from HTSJDK
        LinearIndex linearIndex = new LinearIndex(referenceSequence, 0, new long[]{});

        return new GATKBAMIndexData(this,referenceSequence,bins,linearIndex);
    }

    @Override
    public int getLevelSize(int levelNumber) {
        return index.getLevelSize(levelNumber);
    }

    @Override
    public int getLevelForBin(Bin bin) {
        return index.getLevelForBin(bin);
    }

    @Override
    public int getFirstLocusInBin(Bin bin) {
        return index.getFirstLocusInBin(bin);
    }

    @Override
    public int getLastLocusInBin(Bin bin) {
        return index.getLastLocusInBin(bin);
    }

    @Override
    public long getStartOfLastLinearBin() {
        return index.getStartOfLastLinearBin();
    }
}
