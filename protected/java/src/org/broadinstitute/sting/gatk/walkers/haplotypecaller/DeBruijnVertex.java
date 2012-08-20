package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */
// simple node class for storing kmer sequences
public class DeBruijnVertex {

    protected final byte[] sequence;
    public final int kmer;

    public DeBruijnVertex( final byte[] sequence, final int kmer ) {
        this.sequence = sequence;
        this.kmer = kmer;
    }

    @Override
    public boolean equals( Object v ) {
        return v instanceof DeBruijnVertex && Arrays.equals(sequence, ((DeBruijnVertex) v).sequence);
    }

    @Override
    public int hashCode() { // necessary to override here so that graph.containsVertex() works the same way as vertex.equals() as one might expect
        return Arrays.hashCode(sequence);
    }

    public String toString() {
        return new String(sequence);
    }   
    
    public String getSuffixString() {
        return new String( getSuffix() );
    }

    public byte[] getSequence() {
        return sequence;
    }

    public byte[] getSuffix() {
        return Arrays.copyOfRange( sequence, kmer - 1, sequence.length );
    }
}
