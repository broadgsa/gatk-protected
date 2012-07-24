package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

/**
 * Simple byte / base index conversions
 *
 *
 * @author carneiro
 * @since 8/26/11
 */
public enum BaseIndex {
    A ( 'A', 0 ),
    C ( 'C', 1 ),
    G ( 'G', 2 ),
    T ( 'T', 3 ),
    D ( 'D', 4 ),
    I ( 'I', 5 ), // insertion to the right of the base
    N ( 'N', 6 );

    final byte b;
    final int index;

    public byte getByte() { return b; }

    private BaseIndex(char base, int index) {
        this.b = (byte)base;
        this.index = index;
    }

    /**
     * Converts a byte representation of a base to BaseIndex
     *
     * @param base the byte representation of the base
     * @return the BaseIndex representation of the base;
     */
    public static BaseIndex byteToBase(final byte base) {
        switch (base) {
            case 'A':
            case 'a':
                return A;
            case 'C':
            case 'c':
                return C;
            case 'G':
            case 'g':
                return G;
            case 'T':
            case 't':
                return T;
            case 'D':
            case 'd':
            case '-':
                return D;
            case 'I':
            case 'i':
                return I;
            case 'N':
            case 'n':
                return N;
            default: return null;
        }
    }

    /**
     * Definition of a nucleotide for the BaseIndex is anything that has been read as a base
     * by the machine (A,C,G,T), even if it couldn't tell which base it was, but it knows
     * there is a base there (N).
     *
     * @return whether or not it is a nucleotide, given the definition above
     */
    public boolean isNucleotide() {
        return this == A || this == C || this == G || this == T || this == N;
    }

    /**
     * Whether or not this base is an insertion or a deletion
     *
     * @return true for I or D, false otherwise
     */
    public boolean isIndel() {
        return this == D || this == I;
    }
}
