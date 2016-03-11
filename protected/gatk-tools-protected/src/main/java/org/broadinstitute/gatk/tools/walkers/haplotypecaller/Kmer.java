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

import com.google.java.contract.Requires;

import java.util.Arrays;

/**
 * Fast wrapper for byte[] kmers
 *
 * This objects has several important features that make it better than using a raw byte[] for a kmer:
 *
 * -- Can create kmer from a range of a larger byte[], allowing us to avoid Array.copyOfRange
 * -- Fast equals and hashcode methods
 * -- can get actual byte[] of the kmer, even if it's from a larger byte[], and this operation
 *    only does the work of that operation once, updating its internal state
 *
 * User: depristo
 * Date: 4/8/13
 * Time: 7:54 AM
 */
public class Kmer {
    // this values may be updated in the course of interacting with this kmer
    protected byte[] bases;
    protected int start;

    // two constants
    final protected int length;
    final protected int hash;

    /**
     * Create a new kmer using all bases in kmer
     * @param kmer a non-null byte[]
     */
    public Kmer(byte[] kmer) {
        this(kmer, 0, kmer.length);
    }

    /**
     * Create a new kmer based on the string kmer
     *
     * This is not a good method to use for performance
     *
     * @param kmer the bases as a string
     */
    public Kmer(final String kmer) {
        this(kmer.getBytes());
    }

    /**
     * Create a new kmer backed by the bases in bases, spanning start -> start + length
     *
     * Under no circumstances can bases be modified anywhere in the client code.  This does not make a copy
     * of bases for performance reasons
     *
     * @param bases an array of bases
     * @param start the start of the kmer in bases, must be >= 0 and < bases.length
     * @param length the length of the kmer.  Must be >= 0 and start + length < bases.length
     */
    public Kmer(final byte[] bases, final int start, final int length) {
        if ( bases == null ) throw new IllegalArgumentException("bases cannot be null");
        if ( start < 0 ) throw new IllegalArgumentException("start must be >= 0 but got " + start);
        if ( length < 0 ) throw new IllegalArgumentException("length must be >= 0 but got " + length);
        if ( (start + length) > bases.length ) throw new IllegalArgumentException("start + length " + (start + length) + " must be <= bases.length " + bases.length + " but got " + start + " with length " + length);
        this.bases = bases;
        this.start = start;
        this.length = length;
        this.hash = myHashCode(bases, start, length);
    }

    /**
     * Create a new kmer that's a shallow copy of kmer
     * @param kmer the kmer to shallow copy
     */
    public Kmer(final Kmer kmer) {
        this.bases = kmer.bases;
        this.start = kmer.start;
        this.length = kmer.length;
        this.hash = kmer.hash;
    }

    public Kmer(final Kmer kmer, final byte nextChar) {
        final byte[] sequence = new byte[kmer.length];
        System.arraycopy(kmer.bases,kmer.start + 1,sequence,0,kmer.length - 1);
        sequence[kmer.length - 1] = nextChar;
        bases = sequence;
        start = 0;
        length = kmer.length;
        hash = myHashCode(bases,start,length);
    }

    /**
     * Create a derived shallow kmer that starts at newStart and has newLength bases
     * @param newStart the new start of kmer, where 0 means that start of the kmer, 1 means skip the first base
     * @param newLength the new length
     * @return a new kmer based on the data in this kmer.  Does not make a copy, so shares most of the data
     */
    public Kmer subKmer(final int newStart, final int newLength) {
        return new Kmer(bases, start + newStart, newLength);
    }

    /**
     * Get the bases of this kmer.  May create a copy of the bases, depending on how this kmer was constructed.
     *
     * Note that this function is efficient in that if it needs to copy the bases this only occurs once.
     *
     * @return a non-null byte[] containing length() bases of this kmer, regardless of how this kmer was created
     */
    public byte[] bases() {

        if ( start != 0 || bases.length != length ) {
            // update operation.  Rip out the exact byte[] and update start so we don't ever do this again
            bases = Arrays.copyOfRange(bases, start, start + length);
            start = 0;
        }

        return bases;
    }


    /**
     * Copies kmer bytes into a byte array.
     *
     * @param start first position of the kmer to copy
     * @param dest  what array to copy into
     * @param offset what position the first byte to copy should go into the destination array.
     * @param length how many bytes to copy
     *
     * @throws IllegalArgumentException if <code>start</code> is negative or combined with <code>length</code> it goes
     *        beyond the end of the kmer. Also if <code>length</code> is negative.
     * @throws NullPointerException if dest is <code>null</code>
     * @throws ArrayIndexOutOfBoundsException if dest does not have capacity to received the data.
     */
    public void copyTo(final int start, final byte[] dest, final int offset, final int length) {
        if (start + length > this.length) {
            throw new IllegalArgumentException("request goes beyond end of kmer");
        }
        if (length < 0) {
            throw new IllegalArgumentException("requested length cannot be negative");
        }
        System.arraycopy(bases,this.start + start,dest,offset,length);
    }

    /**
     * Copies kmer bytes into a byte array.
     *
     * @param dest  what array to copy into
     * @param offset what position the first byte to copy should go into the destination array.
     *
     * @throws IllegalArgumentException if <code>start</code> is negative or combined with <code>length</code> it goes
     *        beyond the end of the kmer. Also if <code>length</code> is negative.
     * @throws NullPointerException if dest is <code>null</code>
     */
    public void copyTo(final byte[] dest, final int offset) {
        System.arraycopy(bases,start,dest,offset,length);
    }

    /**
     * Backdoor method for fast base peeking: avoids copying like bases() and doesn't modify internal state.
     * Intended to be used for fast computation of neighboring kmers
     * @return                        Reference to complete bases stores in this kmer
     * WARNING: UNSAFE, caller should NEVER modify bases. Speed/safety tradeoff!!
     */
    private byte[] unsafePeekAtBases() {
        return bases;
    }
    /**
     * Get a string representation of the bases of this kmer
     * @return a non-null string
     */
    public String baseString() {
        return new String(bases());
    }

    /**
     * The length of this kmer
     * @return an integer >= 0
     */
    public int length() {
        return length;
    }

    /**
     * Gets a set of differing positions and bases from another k-mer, limiting up to a max distance.
     * For example, if this = "ACATT" and other = "ACGGT":
     * - if maxDistance < 2 then -1 will be returned, since distance between kmers is 2.
     * - If maxDistance >=2, then 2 will be returned, and arrays will be filled as follows:
     * differingIndeces = {2,3}
     * differingBases = {'G','G'}
     * @param other                 Other k-mer to test
     * @param maxDistance           Maximum distance to search. If this and other k-mers are beyond this Hamming distance,
     *                              search is aborted and a null is returned
     * @param differingIndeces      Array with indices of differing bytes in array
     * @param differingBases        Actual differing bases
     * @return                      Set of mappings of form (int->byte), where each elements represents index
     *                              of k-mer array where bases mismatch, and the byte is the base from other kmer.
     *                              If both k-mers differ by more than maxDistance, returns null
     */
    @Requires({"other != null","differingIndeces != null","differingBases != null",
            "differingIndeces.size>=maxDistance","differingBases.size>=maxDistance"})
    public int getDifferingPositions(final Kmer other,
                                                   final int maxDistance,
                                                   final int[] differingIndeces,
                                                   final byte[] differingBases) {


        int dist = 0;
        if (length == other.length()) {
            final byte[] f2 = other.unsafePeekAtBases();
            for (int i=0; i < length; i++)
                if(bases[start+i] != f2[i]) {
                    differingIndeces[dist] = i;
                    differingBases[dist++] = f2[i];
                    if (dist > maxDistance)
                        return -1;
                }

        }
        return dist;
    }

    @Override
    public String toString() {
        return "Kmer{" + new String(bases,start,length) + "}";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || !Kmer.class.isAssignableFrom(o.getClass())) return false;

        final Kmer kmer = (Kmer) o;

        // very fast test.  If hash aren't equal you are done, otherwise compare the bases
        if ( hash != kmer.hash ) return false;
        if ( length != kmer.length ) return false;

        for ( int i = 0; i < length; i++ )
            if ( bases[start + i] != kmer.bases[kmer.start + i] )
                return false;

        return true;
    }

    @Override
    public int hashCode() {
        return hash;
    }

    /**
     * Helper method that computes the hashcode for this kmer based only the bases in
     * a[], starting at start and running length bases
     *
     * @param a a non-null bases array
     * @param start where to start in bases
     * @param length the length of the bases
     * @return a hashcode value appropriate for a[start] -> a[start + length]
     */
    private static int myHashCode(final byte a[], final int start, final int length) {
        if (a == null)
            return 0;

        int result = 1;
        for (int i = 0; i < length; i++)
            result = 31 * result + a[start + i];

        return result;
    }

    public byte base(final int i) {
        return bases[start + i];
    }

    public Kmer shift(final byte nextChar) {
        if (bases.length > start + length && bases[start + length] == nextChar) {
           return new Kmer(bases,start + 1,length);
        } else {
           final byte[] newBases = new byte[length];
           System.arraycopy(bases,start + 1,newBases,0,length - 1);
           newBases[length - 1] = nextChar;
           return new Kmer(newBases,0,length);
        }
    }

    public byte lastBase() {
        return bases[start + length - 1];
    }
}
