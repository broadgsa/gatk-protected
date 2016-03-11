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

package org.broadinstitute.gatk.utils.collections;

import com.google.java.contract.Requires;

import java.lang.reflect.Array;
import java.util.*;

/**
 * Efficient implementation for a small set of integer primitive values.
 * <p>
 * It includes a increment operation incAll which is convenient when analyzing the read-threading graphs. Nevertheless
 * it can be also be used in general purpose.
 * </p>
 * <p>
 * It does not provide a O(1) look-up of its elements though. These are kept in a sorted array so look up is implemented
 * using a binary search O(log n). Therefore it might not be optimal for problems that require large integer sets.
 * </p>
 * <p>
 * Also note that addition can be costly for large sets unless done in order: O(n).
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CountSet implements Cloneable, Set<Integer> {

    /**
     * The size of the set.
     */
    private int size;

    /**
     * Holds the element of the set within the subrange <code>[0 .. size - 1]</code> in ascending order.
     */
    private int[] elements;

    /**
     * Creates a copy of an existing int-set.
     * @param template the intset to copy values from.
     */
    public CountSet(final CountSet template) {
        elements = template.elements.clone();
        size = template.size;
    }

    /**
     * Creates a new set indicating the expected maximum number of elements it will contain.
     * @param initialCapacity the desired initial capacity of the set.
     * @throws IllegalArgumentException if <code>initialCapacity</code> is negative.
     */
    public CountSet(int initialCapacity) {
        if (initialCapacity < 0)
            throw new IllegalArgumentException();
        elements = new int[initialCapacity];
        size = 0;
    }

    /**
     * Set the set contents to a single integer value.
     * @param value the integer value to set the set to.
     */
    public void setTo(int value) {
        ensureCapacity(1);
        size = 1;
        elements[0] = value;
    }

    /**
     * Set the content of this set to a collection of integers.
     * @param values the new values to be included in the set.
     * @throws NullPointerException if <code>value</code> is <code>null</code>.
     */
    public void setTo(int ... values) {
        ensureCapacity(values.length);
        size = values.length;
        System.arraycopy(values, 0, elements, 0, size);
        Arrays.sort(elements,0,size);
    }

    /**
     * Increase (or decrease) all elements in the set by a number.
     * @param delta the number of add (or substract if negative) to all elements.
     *
     * @return <code>true</code> if the set changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean incAll(final int delta) {
        if (size == 0 || delta == 0)
            return false;
        for (int i = 0; i < size; i++)
            elements[i] += delta;
        return true;
    }

    /**
     * Returns the smallest integer value in the set.
     *
     * @throws NoSuchElementException if the set is empty (thus there is no minimum).
     * @return the smallest integer value in the set.
     */
    public int min() {
        if (size == 0)
            throw new NoSuchElementException("cannot have a min from an empty set");
        return elements[0];
    }

    /**
     * Returns the largest integer value in the set.
     *
     * @throws NoSuchElementException if the set is empty (thus there is no maximum).
     * @return the largest integer value in the set.
     */
    public int max() {
        if (size == 0)
            throw new NoSuchElementException("cannot have a max from an empty set");
        return elements[size - 1];
    }

    /**
     * Adds a range of integer values to the collection.
     *
     * This method avoid the need to explicity indicate all values in that range. Notice that the range is fully inclusive.
     * You can indicate a decrease range (fromValue > toValue).
     *
     * @param fromValue the first value to add in the set (inclusive).
     * @param toValue the last value to add to the set (inclusive).
     * @return <code>true</code> if the set changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean addRange(final int fromValue, final int toValue) {
        final int lowEnd;
        final int highEnd;

        if (fromValue <= toValue) {
            lowEnd = fromValue; highEnd = toValue;
        } else {
            highEnd = fromValue; lowEnd = toValue;
        }

        //TODO to be optimized to add missing sub-ranges in one go:
        boolean result = false;
        for (int i = lowEnd; i <= highEnd; i++)
            result = add(i) | result;
        return result;
    }

    /**
     * Add an integer value to the set.
     * @param value to add to the set.
     * @return <code>true</code> if the set changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean add(final int value) {
        int pos = Arrays.binarySearch(elements,0,size,value);
        if (pos >= 0) return false;
        int insertPos = - pos - 1;
        ensureCapacity(size + 1);
        System.arraycopy(elements, insertPos, elements, insertPos + 1, size - insertPos);
        elements[insertPos] = value;
        size++;
        return true;
    }

    /**
     * Add a arbitrary number of integers to the set.
     *
     * @param values integer to add to the set.
     * @return <code>true</code> if the set changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean addAll(final int ... values) {
        ensureCapacity(size + values.length);
        boolean result = false;
        for (final int v : values)
            result = add(v) | result;
        return result;
    }

    @Override
    public boolean addAll(final Collection<? extends Integer> numbers) {
        ensureCapacity(size + numbers.size());
        boolean result = false;
        for (final Number n : numbers)
            result = add(n.intValue()) | result;
        return result;
    }

    /**
     * Add all values within a range in an integer array.
     *
     * @param source array where the values to add are found.
     * @param fromIndex first position from <code>source</code> to add (inclusive).
     * @param toIndex index after the last position in <code>source</code> to add (thus exclusive).
     * @throws NullPointerException if <code>source</code> is <code>null</code>.
     * @throws NegativeArraySizeException if <code>fromIndex</code> or <code>toIndex</code> are negative.
     * @throws ArrayIndexOutOfBoundsException if <code>fromIndex</code> or <code>toIndex</code> are beyond bounds
     *    allowed <code>[0 .. source.length]</code>.
     * @return <code>true</code> if the set changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean addAll(final int[] source, final int fromIndex, final int toIndex) {
        ensureCapacity(size + source.length);
        boolean result = false;
        for (int i = fromIndex; i < toIndex; i++)
            result = add(source[i]) | result;
        return result;
    }


    /**
     * Add all elements present in a int-set.
     *
     * @param other the other inset.
     *
     * @throws NullPointerException if <code>other</code> is <code>null</code>.
     * @return <code>true</code> if this set changed due to this operation, <code>false</code> otherwise.
     */
    public boolean addAll(final CountSet other) {
        return addAll(other.elements,0,other.size);
    }

    /**
     * Checks whether a integer value is included in the set.
     * @param value the value to check.
     * @return <code>true</code> if <code>value</code> is inside the set, <code>false</code> otherwise.
     */
    public boolean contains(final int value) {
        return Arrays.binarySearch(elements, 0, size, value) >= 0;
    }

    /**
     * Make sure that this int-set has capacity to handle a number of elements.
     * <p/>
     * If the set has already that or greater capacity nothing would be changed.
     *
     * @param capacity the requested capacity.
     */
    private void ensureCapacity(final int capacity) {
        if (elements.length >= capacity) return;
        int newLength = Math.max(elements.length << 1, capacity);
        elements = Arrays.copyOf(elements,newLength);
    }


    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size() == 0;
    }

    @Override
    public boolean contains(final Object o) {
        if (o instanceof Integer) {
            final int i = (Integer)o;
            return contains(i);
        } else
            return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public Iterator<Integer> iterator() {
        return new MyIterator();
    }

    @Override
    public Object[] toArray() {
        final Integer[] result = new Integer[size];
        for (int i = 0; i < size; i++)
            result[i] = elements[i];
        return result;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T[] toArray(final T[] a) {
        if (a == null)
            throw new NullPointerException();

        @SuppressWarnings("unchecked")
        final Class<T> componentClass = (Class) a.getClass().getComponentType();
        if (!componentClass.isAssignableFrom(Integer.class))
            throw new ArrayStoreException();

        @SuppressWarnings("unchecked")
        final T[] dest = (a.length < size) ? (T[]) Array.newInstance(componentClass, size) : a;

        for (int i = 0; i < size; i++)
            dest[i] = (T) (Integer) elements[i];
        return dest;
    }

    /**
     * Copies the content of the set into an integer array. The result can be freely modified by the invoker.
     * @return never <code>null</code> but a zero-length array if the set is empty.
     */
    public int[] toIntArray() {
        return Arrays.copyOfRange(elements,0,size);
    }

    /**
     * Copy the content of the set into an array.
     * @param dest the destination array.
     * @param offset where to store the first element of the set.
     * @throws NullPointerException if <code>dest</code> is <code>null</code>.
     * @throws ArrayIndexOutOfBoundsException if <code>offset</code> is out of range of there is not enough
     * space after <code>offset</code> in the destination array to hold all values in the set.
     */
    public void copyTo(final int[] dest, int offset) {
        if (dest == null)
            throw new NullPointerException();
        if (dest.length < (size + offset))
            throw new ArrayIndexOutOfBoundsException("destination is to short");
        System.arraycopy(elements,0,dest,offset,size);
    }

    /**
     * Copy the content of the set into an array.
     * @param dest the destination array.
     * @throws NullPointerException if <code>dest</code> is <code>null</code>.
     * @throws ArrayIndexOutOfBoundsException if there is not enough
     * space after <code>offset</code> in the destination array to hold all values in the set.
     */
    public void copyTo(final int[] dest) {
        copyTo(dest,0);
    }


    @Override
    public boolean add(final Integer integer) {
        return add((int) integer);
    }

    @Override
    public boolean remove(final Object o) {
        return o instanceof Integer && remove((int)o);
    }

    /**
     * Removes a single integer value for the set.
     * @param i the value to remove.
     * @return <code>true</code> if the set has changed as a result of this invocation, <code>false</code> otherwise.
     */
    public boolean remove(final int i) {
        final int pos = Arrays.binarySearch(elements,0,size,i);
        if (pos < 0)
            return false;
        else {
            removeIndex(pos);
            return true;
        }
    }

    @Override
    public boolean containsAll(final Collection<?> c) {
        for (final Object o : c)
            if (!contains(o))
                return false;
        return true;
    }


    @Override
    public boolean retainAll(final Collection<?> c) {
        if (size == 0)
            return false;
        @SuppressWarnings("all")
        final CountSet retainIndices = new CountSet(c.size() + 2);
        retainIndices.add(-1);
        retainIndices.add(size);
        for (final Object o : c) {
            if (!(o instanceof Integer))
                continue;
            final int pos = Arrays.binarySearch(elements,0,size,(int) o);
            if (pos < 0)
                continue;
            retainIndices.add(pos);
        }
        if (retainIndices.size == 2) {
            size = 0;
            return true;
        } else if (retainIndices.size == size + 2) {
            return false;
        } else {
            for (int idx = retainIndices.size - 1; idx > 0; idx--) {
                final int toIdx = retainIndices.elements[idx];
                final int fromIdx = retainIndices.elements[idx - 1] + 1;
                removeIndices(toIdx,fromIdx);
            }
            return true;
        }
    }

    /**
     * Removes the values found in a range of indexes in {@link #elements}.
     * @param fromIdx first index to remove (inclusive).
     * @param toIdx right after last index to remove (exclusive).
     */
    @Requires("fromIdx >= toIdx & fromIdx >= 0 & toIdx <= size")
    private void removeIndices(final int fromIdx, final int toIdx) {
       System.arraycopy(elements,toIdx,elements,fromIdx,size - toIdx);
       size -= toIdx - fromIdx;
    }

    @Override
    public boolean removeAll(final Collection<?> c) {
        boolean result = false;
        for (final Object o : c)
            result = remove(o) | result;
        return result;
    }

    @Requires("idx >= 0 && idx < size")
    private void removeIndex(int idx) {
        System.arraycopy(elements,idx+1,elements,idx,size - idx - 1);
    }

    @Override
    public void clear() {
        size = 0;
    }

    /**
     * Returns a copy of this set which can be changed without modifying the original one.
     * @return never {@code null}.
     */
    @SuppressWarnings("all")
    public CountSet clone() {
        return new CountSet(this);
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder(2 + size() * 10);
        sb.append('{');
        for (int i = 0; i < size; i++)
            sb.append(elements[i]).append(',');
        sb.replace(sb.length()-1,sb.length(),"}");
        return sb.toString();

    }


    /**
     * Custom iterator class for {@link CountSet IntSets}
     */
    private class MyIterator implements Iterator<Integer> {
        /** What position I am in. */
        private int next = 0;

        @Override
        public boolean hasNext() {
            return next < size;
        }

        @Override
        public Integer next() {
            if (next >= size)
                throw new NoSuchElementException();
            return elements[next];
        }

        @Override
        public void remove() {
            if (next == 0)
                throw new IllegalStateException();
            if (next >= size)
                throw new NoSuchElementException();
            removeIndex(next - 1);
        }
    }


}
