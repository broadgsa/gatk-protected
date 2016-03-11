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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;

import java.lang.reflect.Array;
import java.util.*;

/**
 * Represent a sequence of kmers where any two consecutive kmers overlap in kmer length - 1 elements.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 */
public class KmerSequence implements List<Kmer> {
    private final byte[] sequence;
    private final int start;
    private final int size;
    private final int kmerSize;
    private final int rawLength;

    /**
     * Creates a kmer sequence from a read's sequence.
     *
     * @param read the read to represent as a sequence of kmers.
     * @param kmerSize the kmer size.
     */
    public KmerSequence(final SAMRecord read, final int kmerSize) {
        this(read.getReadBases(), kmerSize);
    }

    /**
     * Creates a kmer sequence from a haplotype's sequence.
     *
     * @param hap the haplotype to represent as a sequence of kmers.
     * @param kmerSize the kmer size.
     */
    @SuppressWarnings("unused")
    public KmerSequence(final Haplotype hap, final int kmerSize) {
        this(hap.getBases(), kmerSize);
    }

    /**
     * Creates a kmer sequence out of a byte sequence.
     *
     * @param sequence the byte array to represent as a kmer sequence.
     * @param kmerSize the kmer size.
     */
    public KmerSequence(final byte[] sequence, final int kmerSize) {
        this(sequence,0,Math.max(0,sequence.length - kmerSize + 1),kmerSize, sequence.length);
    }

    /**
     * Creates a kmer sequence out of a range of a byte array
     *
     * @param sequence the input array.
     * @param start inclusive first position of the array that maps to the first position in the first kmer.
     * @param size number kmers in the output.
     * @param kmerSize kmer length in bases.
     * @param rawLength the of the range in bases.
     */
    protected KmerSequence(final byte[] sequence, final int start, final int size, final int kmerSize, final int rawLength) {
        if (sequence == null) {
            throw new IllegalArgumentException("start must be 0 or greater");
        }
        if (rawLength > sequence.length - start) {
            throw new IllegalArgumentException("the raw sequence length goes beyond the array capacity");
        }
        if (size < 0) {
            throw new IllegalArgumentException("the length cannot be negative");
        }
        if (start < 0) {
            throw new IllegalArgumentException("start must be 0 or greater");
        }
        if (size > 0 && size + kmerSize - 1 > rawLength) {
            throw new IllegalArgumentException(
                    String.format("the kmerSize (%d) + size (%d) - 1 cannot be larger than rawLength (%d)",kmerSize,size,rawLength) );
        }
        this.sequence = sequence;
        this.start = start;
        this.size = size;
        this.kmerSize = kmerSize;
        this.rawLength = rawLength;
    }

    public int kmerSize() {
        return kmerSize;
    }

    public KmerSequence subsequence(final int from, final int to) {
        if (from < 0 || from > to) {
            throw new IllegalArgumentException();
        }
        if (to > size) {
            throw new IllegalArgumentException();
        }
        return new KmerSequence(sequence,this.start + from,to - from,kmerSize,rawLength - from - (size - to));
    }


    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    @Override
    public boolean contains(final Object o) {
        if (o instanceof Kmer) {
            if (o instanceof MyKmer) {
                final MyKmer k = (MyKmer) o;
                if (k.bases == sequence && k.start >= start && k.length == kmerSize && k.start < start + size) {
                    return true;
                }
            }
            final Kmer k = (Kmer) o;
            if (k.length != kmerSize) {
                return false;
            }
            for (int i = 0; i < size; i++) {
                int j;
                for (j = 0; j < kmerSize; j++) {
                    if (sequence[start + i + j] != k.bases[k.start + j]) {
                         break;
                    }
                }
                if (j == kmerSize) {
                    return true;
                }
            }
            return false;
        } else {
            return false;
        }
    }

    @Override
    public Iterator<Kmer> iterator() {
        return new Iterator<Kmer>() {

            private int offset = 0;

            @Override
            public boolean hasNext() {
                return offset < size;
            }

            @Override
            public Kmer next() {
                return new Kmer(sequence,start + offset,kmerSize);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public Object[] toArray() {
        return toArray(new Kmer[size()]);
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T[] toArray(final T[] a) {
        if (a == null) {
            throw new IllegalArgumentException();
        } else if (!a.getClass().getComponentType().isAssignableFrom(Kmer.class)) {
            throw new IllegalArgumentException();
        } else {
            T[] result;
            if (a.length < size) {
                result = (T[]) Array.newInstance(a.getClass().getComponentType(), size);
            } else {
                result = a;
            }
            for (int i = 0; i < size; i++) {
                result[i] = (T) new Kmer(sequence,start + i,kmerSize);
            }
            return result;
        }
    }

    @Override
    public boolean add(final Kmer kmer) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean remove(final Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean containsAll(final Collection<?> c) {
        for (final Object o : c)
            if (!contains(o))
                return false;
        return true;
    }

    @Override
    public boolean addAll(final Collection<? extends Kmer> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(final int index, final Collection<? extends Kmer> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(final Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(final Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Kmer get(final int index) {
        if (index < 0 || index >= size) {
            throw new IllegalArgumentException();
        }
        return new Kmer(sequence,start + index,kmerSize);
    }

    @Override
    public Kmer set(final int index, final Kmer element) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void add(final int index, final Kmer element) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Kmer remove(final int index) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int indexOf(final Object o) {
        if (o instanceof Kmer) {
            final Kmer k = (Kmer) o;
            if (k.length != kmerSize) {
                return -1;
            }
            for (int i = 0; i < size; i++) {
                int j;
                for (j = 0; j < kmerSize; j++) {
                    if (sequence[start + i + j] != k.bases[k.start + j]) {
                        break;
                    }
                }
                if (j == kmerSize) {
                    return i;
                }
            }
            return -1;
        } else {
            return -1;
        }
    }

    @Override
    public int lastIndexOf(final Object o) {
        if (o instanceof Kmer) {
            final Kmer k = (Kmer) o;
            if (k.length != kmerSize) {
                return -1;
            }
            for (int i = size - 1; i >= 0; i--) {
                int j;
                for (j = kmerSize - 1; j >= 0; j--) {
                    if (sequence[start + i + j] != k.bases[k.start + j]) {
                        break;
                    }
                }
                if (j == 0) {
                    return i;
                }
            }
            return -1;
        } else {
            return -1;
        }
    }

    @Override
    public ListIterator<Kmer> listIterator() {
        return new MyListIterator(0);
    }

    @Override
    public ListIterator<Kmer> listIterator(final int index) {
        return new MyListIterator(index);
    }

    @Override
    public List<Kmer> subList(final int fromIndex, final int toIndex) {
        return subsequence(fromIndex,toIndex);
    }

    /**
     * Returns the byte array representation of the kmer sequence.
     * @return never {@code null}.
     */
    public byte[] getBytes() {
        if (start == 0 && rawLength == sequence.length)
            return sequence;
        else
            return Arrays.copyOfRange(sequence, start, rawLength + start);
    }

    /**
     * Internal class that implements the {@link Kmer} more efficiently
     * making reference to the sequence's own byte array.
     */
    protected class MyKmer extends Kmer {

        /**
         * Create a new instance give the offset in the byte array.
         * @param start the start base offset for the kmer.
         */
        public MyKmer(final int start) {
            super(sequence,start,kmerSize);
        }
    }

    /**
     * Iterator implementation of Kmer elements.
     */
    private class MyListIterator implements ListIterator<Kmer> {

        private int i = 0;

        /**
         * Creates a iterator at certain offset in the sequence.
         * @param idx the start position or kmer offset.
         */
        private MyListIterator(final int idx) {
            i = idx;
        }

        @Override
        public boolean hasNext() {
            return i < size;
        }

        @Override
        public Kmer next() {
            return new Kmer(sequence,start + i++,kmerSize);
        }

        @Override
        public boolean hasPrevious() {
            return i > 0;
        }

        @Override
        public Kmer previous() {
            return new Kmer(sequence,start + --i,kmerSize);
        }

        @Override
        public int nextIndex() {
            return i;
        }

        @Override
        public int previousIndex() {
            return i - 1;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void set(final Kmer kmer) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void add(final Kmer kmer) {
            throw new UnsupportedOperationException();
        }

    }

}
