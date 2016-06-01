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
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.*;

/**
 * Compact Idosyncratic Variation Alignment Report.
 * <p/>
 * Allows to specify variation of a sequence.
 * <p/>
 */
public class Civar {


    protected List<Element> elements;
    private String string;
    private transient int minimumTemplateSize = -1;
    private transient Boolean expands = null;
    private transient int starCount = -1;
    private Boolean hasEmbeddedCivars = null;
    private Boolean hasOptionalVariation = null;
    private Boolean hasVariation = null;
    private Boolean allVariationIsOptional;


    protected Civar(final List<Element> elements) {
        this.elements = elements;
    }

    public static Civar fromCharSequence(final CharSequence cs, final int from, final int to) {
        return Parser.parse(cs, from, to);
    }

    public static Civar fromCharSequence(final CharSequence cs) {
        return fromCharSequence(cs, 0, cs.length());
    }

    @Override
    public String toString() {
        if (string == null)
            analyzeCivar();
        return string;
    }

    /**
     * Correspond to the minimum template sequence size to which this civar can be applied to.
     */
    public int minimumTemplateSequenceSize() {
        if (minimumTemplateSize < 0)
            analyzeCivar();
        return minimumTemplateSize;
    }

    public boolean expands() {
        if (expands == null)
            analyzeCivar();
        return expands;
    }

    public boolean hasEmbeddedCivars() {
        if (hasEmbeddedCivars == null)
            analyzeCivar();
        return hasEmbeddedCivars;
    }

    public List<Civar> unroll() {
        if (!isUnrolled()) {
            final List<Civar> result = new LinkedList<>();
            unroll(0, new LinkedList<Element>(), "", result);
            return result;
        } else {
            return Collections.singletonList(this);
        }
    }

    public Civar optionalizeAll() {
        if (allElementsAreOptional()) {
            return this;
        }
        final Element[] newElements = new Element[this.elements.size()];
        int next = 0;
        for (final Element e : elements) {
           final Element newElement = e.clone();
           if (newElement.operator() == Operator.EMBEDDED) {
               newElement.embedded = newElement.embedded.optionalizeAll();
           }
           newElement.makeOptional();
           newElements[next++] = newElement;
        }
        return new Civar(Collections.unmodifiableList(Arrays.asList(newElements)));
    }

    private boolean allElementsAreOptional() {
        if (allVariationIsOptional == null)
            analyzeCivar();
        return allVariationIsOptional;
    }


    private void unroll(final int elementIndex, final LinkedList<Element> leftElements, final String leftString, final List<Civar> dest) {
        if (elementIndex == elements.size()) {
            final Civar result = new Civar(Collections.unmodifiableList(new ArrayList<>(leftElements)));
            dest.add(result);
            return;
        }
        final Element currentElement = elements.get(elementIndex);
        if (currentElement.operator() == Operator.EMBEDDED) {
            List<Civar> embeddedUnroll = currentElement.embedded.unroll();
            Civar novar = null;
            for (final Civar c : embeddedUnroll) {
                if (!c.hasVariation()) {
                    novar = c;
                    break;
                }
            }
            if (novar == null && currentElement.isOptional()) {

                embeddedUnroll = new LinkedList<>(embeddedUnroll);
                embeddedUnroll.add(0, currentElement.embedded.novarEquivalent());
            }

            for (final Civar embedded: embeddedUnroll) {
                final Element embeddedElement = new Element(embedded,false);
                leftElements.add(embeddedElement);
                unroll(elementIndex + 1, leftElements, leftString + embeddedElement.toString(), dest);
                leftElements.removeLast();
            }
        } else if (currentElement.isOptional() && currentElement.isVariant()) {
            leftElements.add(currentElement.matchEquivalent());
            unroll(elementIndex + 1, leftElements, leftString + currentElement.matchEquivalent().toString(), dest);
            leftElements.removeLast();
            leftElements.add(currentElement.mandatoryEquivalent());
            unroll(elementIndex + 1, leftElements, leftString + currentElement.mandatoryEquivalent().toString(), dest);
            leftElements.removeLast();
        } else {
            leftElements.add(currentElement);
            unroll(elementIndex + 1, leftElements, leftString + currentElement.toString(), dest);
            leftElements.removeLast();
        }
    }

    private Civar novarEquivalent() {
        int mtss = minimumTemplateSequenceSize();
        int sc = starCount();
        if (mtss > 0) {
            if (sc > 0) {
                final Civar newEmbedded = new Civar(Collections.unmodifiableList(Arrays.asList(
                        new Element(Operator.MATCH,mtss,false,false),new Element(Operator.MATCH,sc,true,false))));
                //newEmbedded.string = newEmbedded.elements.get(0).toString() + newEmbedded.elements.get(1).toString();
                return newEmbedded;
            } else {
                return new Civar(Collections.unmodifiableList(Collections.singletonList(new Element(Operator.MATCH, mtss, false, false))));
            }
        } else if (sc > 0) {
            return new Civar(Collections.unmodifiableList(Collections.singletonList(new Element(Operator.MATCH, sc, true, false))));
        } else {
            return new Civar((List<Element>) (List) Collections.emptyList());
        }
    }

    public String applyTo(final CharSequence seq) {
        return applyTo(seq,0,seq.length());
    }

    public List<ElementOffset> eventOffsets(final CharSequence seq, final int from, final int to) {
        if (!isUnrolled())
            throw new UnsupportedOperationException("you cannot apply an unrolled Civar to a DNA sequence");
        final List<ElementOffset> result = new ArrayList<>(elements().size());
        final CharSequence sequence = seq;

        int sequenceLength = sequence.length();
        int minSeqLen = minimumTemplateSequenceSize();
        if (!expands() && sequenceLength != minSeqLen)
            throw new IllegalArgumentException("the sequence provided does not match this Civar size " + sequence.length() + " != " + minSeqLen);
        if (sequenceLength < minSeqLen)
            throw new IllegalArgumentException("the sequence provided is too small for this Civar " + sequence.length() + " < " + minSeqLen);
        int starCount = starCount();
        int paddingTotal = sequenceLength - minSeqLen;
        int starPadding = starCount == 0 ? 0 : paddingTotal / starCount;
        int excessPadding = starCount == 0 ? paddingTotal : paddingTotal % starCount;
        int nextInSeq = 0;
        int nextElement = 0;

        int outputLength = 0;
        while (nextInSeq < sequenceLength && nextElement < elements().size()) {
            final Element e = elements().get(nextElement++);
            final int outputStart = outputLength;
            final int sequenceStart = nextInSeq;
            int size = e.expands() ? starPadding * e.size() : e.size();
            if (e.expands() && excessPadding != 0) {
                size++;
                excessPadding--;
            }
            switch (e.operator()) {
                case EMBEDDED:
                    throw new IllegalStateException("supposed to be unrolled Civar");
                case DELETION:
                    nextInSeq += size;
                    break;
                case INSERTION:
                    outputLength += size;
                    break;
                default:
                    outputLength += size;
                    nextInSeq += size;
                    break;
            }
            final int outputEnd = outputLength;
            final int sequenceEnd = nextInSeq;
            if (outputEnd > from && outputStart < to) {
                result.add(new ElementOffset(e,Math.max(outputStart - from,0),
                      Math.min(outputEnd - from, to -from),
                      sequenceStart + Math.max(from - outputStart,0),
                      sequenceEnd - Math.max(outputEnd - to,0)));
            }
        }
        if (nextInSeq != sequenceLength) {
            throw new IllegalStateException("probable bug: mismatched sequence and Civar application length " + nextInSeq + " != " + to);
        }
        return result;
    }

    public String applyTo(final CharSequence seq, final int from, final int to) {

        final CharSequence sequence = seq.subSequence(from, to);
        if (!isUnrolled())
            throw new UnsupportedOperationException("you cannot apply an unrolled Civar to a DNA sequence");

        int referenceLength = sequence.length();
        int minSeqLen = minimumTemplateSequenceSize();
        if (!expands() && referenceLength != minSeqLen)
            throw new IllegalArgumentException("the sequence provided does not match this Civar size " + sequence.length() + " != " + minSeqLen);
        if (referenceLength < minSeqLen)
            throw new IllegalArgumentException("the sequence provided is too small for this Civar " + sequence.length() + " < " + minSeqLen);
        int starCount = starCount();
        int paddingTotal = referenceLength - minSeqLen;
        int starPadding = starCount == 0 ? 0 : paddingTotal / starCount;
        int excessPadding = starCount == 0 ? paddingTotal : paddingTotal % starCount;
        int nextInSeq = 0;
        int nextElement = 0;

        final StringBuffer sb = new StringBuffer(sequence.length() * 2);
        while (nextInSeq < to && nextElement < elements.size()) {
            final Element e = elements.get(nextElement++);
            int size = e.expands() ? starPadding * e.size() : e.size();
            if (e.expands() && excessPadding != 0) {
                size++;
                excessPadding--;
            }
            switch (e.operator()) {
                case EMBEDDED:
                    throw new IllegalStateException("supposed to be unrolled Civar");
                case MATCH:
                    sb.append(sequence.subSequence(nextInSeq, nextInSeq += size));
                    break;
                case DELETION:
                    nextInSeq += size;
                    break;
                case INSERTION:
                    final CharSequence xmer = e.xmer().toString().toUpperCase();
                    while (size >= xmer.length()) {
                        sb.append(xmer);
                        size -= xmer.length();
                    }
                    sb.append(xmer.subSequence(0, size));
                    break;
                case TRANSITION:
                    transition(sequence, nextInSeq, sb, size);
                    nextInSeq += size;
                    break;
                case TRANSVERSION:
                    transversion(sequence, nextInSeq, sb, size);
                    nextInSeq += size;
                    break;
                case COMPLEMENT:
                    complement(sequence, nextInSeq, sb, size);
                    nextInSeq += size;
                    break;
            }
        }
        if (nextInSeq != referenceLength) {
            throw new IllegalStateException("probable bug: mismatched sequence and Civar application length " + nextInSeq + " != " + to);
        }
        return sb.toString();
    }


    private static char transversion(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A':
                return 'C';
            case 'G':
                return 'T';
            case 'C':
                return 'A';
            case 'T':
                return 'G';
            default:
                return c;
        }
    }

    private static char complement(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A':
                return 'T';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            default:
                return c;
        }

    }

    private static char transition(final char c) {
        switch (Character.toUpperCase(c)) {
            case 'A':
                return 'G';
            case 'G':
                return 'A';
            case 'T':
                return 'C';
            case 'C':
                return 'T';
            default:
                return c;
        }

    }

    private void transition(final CharSequence charSequence, final int from, final StringBuffer dest, final int length) {
        for (int i = from; i < length + from; i++) {
            dest.append(transition(charSequence.charAt(i)));
        }
    }

    private void transversion(final CharSequence cs, final int from, final StringBuffer dest, final int length) {
        for (int i = from; i < length + from; i++) {
            dest.append(transversion(cs.charAt(i)));
        }
    }

    private void complement(final CharSequence cs, final int from, final StringBuffer dest, final int length) {
        for (int i = from; i < length + from; i++) {
            dest.append(complement(cs.charAt(i)));
        }
    }

    private static List<Element> cle(final List<Element> original) {
        if (original.size() == 0)
            return original;
        if (original.size() == 1) {
            final Element e = original.get(0);
            if (e.operator() == Operator.EMBEDDED) {
                return original;
            } else if (e.size() == 0) {
                return Collections.emptyList();
            } else {
                return original;
            }
        } else {
            final ArrayList<Element> result = new ArrayList<>(original);
            for (int i = 0; i < result.size(); ) {
                final Element e = result.get(i);
                if (e.operator() == Operator.EMBEDDED) {
                    i++;
                } else if (e.size() == 0) {
                    result.remove(i);
                } else if (i == result.size() - 1) {
                    i++;
                } else {
                    final Element next = result.get(i + 1);
                    if (next.operator() != e.operator()) {
                        i++;
                    } else if (e.operator() != Operator.INSERTION) {
                        if (next.expands() == e.expands()) {
                            result.remove(i + 1);
                            result.set(i, new Element(e.operator(), e.size() + next.size(), e.expands(), false, null));
                        } else if (e.expands()) {
                            int j;
                            for (j = i + 1; j < result.size(); j++) {
                                if (result.get(j).operator() == Operator.MATCH && result.get(j).expands()) {
                                    break;
                                }
                            }
                            result.add(j, e);
                            result.remove(i);
                        }
                    } else {         // INSERTION my be fussed if their sizes correspond to the xmers without expansions
                        if (!e.expands() && !next.expands() && e.xmer().length() == e.size() && next.xmer().length() == next.size()) {
                            result.remove(i + 1);
                            result.set(i, new Element(Operator.INSERTION, e.size() + next.size(), false,  false, e.xmer().toString() + next.xmer().toString()));
                        }
                    }
                }
            }
            return result;
        }
    }


    protected void analyzeCivar() {
        int minimumTemplateSize = 0;
        boolean expands = false;
        int starCount = 0;
        boolean hasEmbeddedCivars = false;
        boolean hasOptionalElements = false;
        boolean allElementsAreOptional = true;
        boolean hasVariation = false;
        StringBuffer strBuffer = new StringBuffer(100);
        for (final Element e : elements) {
            strBuffer.append(e.toString());
            if (e.operator() == Operator.EMBEDDED) {
                hasEmbeddedCivars = true;
                if (e.embedded.hasVariation()) {
                    hasVariation = true;
                    allElementsAreOptional &= (e.optional || e.embedded.allElementsAreOptional());
                    hasOptionalElements |= e.optional || e.embedded.hasOptionalVariationElements();
                }
                minimumTemplateSize += e.embedded.minimumTemplateSequenceSize();
                if (e.embedded.expands()) {
                    expands = true;
                    starCount += e.embedded.starCount();
                }
            } else {
                if (e.isVariant()) {
                    hasVariation = true;
                    allElementsAreOptional &= e.optional;
                    hasOptionalElements |= e.optional;
                }
                if (e.expands()) {
                    starCount += e.size();
                    expands = true;
                    continue;
                }
                if (e.operator() == Operator.INSERTION)
                    continue;
                minimumTemplateSize += e.size();
            }
        }
        this.string = strBuffer.toString();
        this.hasVariation = hasVariation;
        this.allVariationIsOptional = allElementsAreOptional;
        this.hasOptionalVariation = hasOptionalElements;
        this.hasEmbeddedCivars = hasEmbeddedCivars;
        this.starCount = starCount;
        this.expands = expands;
        this.minimumTemplateSize = minimumTemplateSize;
    }

    public boolean isUnrolled() {
        return !hasOptionalVariationElements();
    }

    public boolean hasVariation() {
        if (hasVariation == null)
            analyzeCivar();
        return hasVariation;
    }

    private boolean hasOptionalVariationElements() {
        if (hasOptionalVariation == null)
            analyzeCivar();
        return hasOptionalVariation;
    }

    private int starCount() {
        if (starCount == -1)
            analyzeCivar();
        return starCount;
    }

    public List<Element> elements() {
        return elements;
    }


    public enum Operator {

        MATCH('='), INSERTION('I'), DELETION('D'), TRANSITION('T'), EMBEDDED('('),

        /**
         * Transversion that is not to the complement nucleotide
         */
        TRANSVERSION('V'),

        /**
         * Transverison to the complement nucleotide A <-> T or C <-> G
         */
        COMPLEMENT('C'),

        /**
         * Marks
         */
        START('^'), END('$');

        public final char charValue;

        Operator(final char c) {
            charValue = c;
        }

        /**
         * Does this operator represents a SNP.
         *
         * @return <code>true</code> is so, <code>false</code> otherwise.
         */
        public boolean isSnp() {
            switch (this) {
                case TRANSITION:
                case TRANSVERSION:
                case COMPLEMENT:
                    return true;
                default:
                    return false;
            }
        }

        /**
         * Checks whether the operation requires an x-mer.
         *
         * @return <code>true</code> if so, <code>false</code> otherwise
         */
        public boolean requiresXmer() {
            return this == INSERTION;
        }

        public boolean acceptsXmer() {
            return this == INSERTION;
        }

        public static Operator fromChar(final char c) {
            switch (c) {
                case 'I':
                    return INSERTION;
                case 'V':
                    return TRANSVERSION;
                case 'T':
                    return TRANSITION;
                case 'C':
                    return COMPLEMENT;
                case 'M':
                case '=':
                    return MATCH;
                case '?':
                    return START;
                case '$':
                    return END;
                case 'D':
                    return DELETION;
                default:
                    throw new IllegalArgumentException("the chacter " + c + " does not denote a valid Civar operator");
            }
        }
    }


    public static class Element implements Cloneable {

        protected int size;
        protected boolean expands;
        protected Operator o;

        protected CharSequence xmer;
        protected Civar embedded;
        protected boolean optional;


        public Element(final Civar c, boolean optional) {
            this(Operator.EMBEDDED, 1, false, optional);
            embedded = c;
        }

        @Override
        protected Element clone() {
            try {
                return (Element) super.clone();
            } catch (CloneNotSupportedException e) {
                throw new RuntimeException("unexpected exception ", e);
            }
        }


        /**
         * Calculate the length provided the the length of the padding unit and how much
         * exceed padding is still available.
         *
         * @param starPadding
         * @param excessPaddingRemaining
         *
         * @return never {@code null}.
         */
        public int size(final int starPadding, final int excessPaddingRemaining) {
            if (!expands)
                return size;
            else
                return size * (starPadding) + Math.min(excessPaddingRemaining,size);
        }

        /**
         * Returns the size of this element.
         * <p/>
         * If the element is non-explandable this value is the exact size of the indicated operation in term of bases
         * (deleted, inserted, changed, etc...)
         * <p/>
         * If the element is expandable, this is relative size respect with other expandable elements.
         *
         * @return never negative, and typically greater than 0; 0 would be use for possible future symbolic operands.
         */
        public int size() {
            return size;
        }

        /**
         * Checks whether this is an expandable element.
         *
         * @return <code>true</code> if so, <code>false</code> otherwise.
         */
        public boolean expands() {
            return expands;
        }

        public boolean isOptional() {
            return optional;
        }

        /**
         * Returns the appendix sequence for those elements for which it applies.
         *
         * @return never <code>null</code> but perhaps an empty sequence. You must not attempt to modify the returned element.
         */
        public CharSequence xmer() {
            return xmer;
        }


        /**
         * Returns the operator for this element.
         *
         * @return never <code>null</code>.
         */
        public Operator operator() {
            return o;
        }

        protected Element(final Operator o, final int size, final boolean expands, final boolean optional) {
            if (o == null)
                throw new IllegalArgumentException("operator cannot be null");
            if (size < 0)
                throw new IllegalArgumentException("element size cannot be negative");
            this.o = o;
            this.optional = optional;
            this.expands = expands;
            this.size = size;
        }

        /**
         * Constructs an element given its properties.
         *
         * @param o       the operator for this element.
         * @param size    the size of the element
         * @param expands if the element is expandable.
         * @param xmer    the xmer for this element.
         */
        public Element(final Operator o, final int size, final boolean expands, final boolean optional, final String xmer) {
            this(o,size,expands,optional);
            if ((xmer == null || xmer.length() == 0) && o.requiresXmer())
                throw new IllegalArgumentException("operator  " + o + " requires a x-mer");
            if ((xmer != null && !o.acceptsXmer()))
                throw new IllegalArgumentException("operator  " + o + " does not accept a x-mer");
            this.xmer = xmer;
        }

        public String toString() {
            if (operator() == Operator.EMBEDDED) {
                return "(" + embedded.toString() + ")";
            }
            final String sizeString = expands ? (size == 1 ? "*":"" + size) : ("" + size);
            final String sizeAndOperator = sizeString + operator().charValue;
            final String sizeAndOperatorAndXmer = o == Operator.INSERTION ? sizeAndOperator + xmer : sizeAndOperator;
            return optional ? sizeAndOperatorAndXmer + "?" : sizeAndOperatorAndXmer;
        }

        protected void makeOptional() {
            if (o != Operator.MATCH)
                optional = true;
        }

        public boolean isVariant() {
            switch (o) {
                case MATCH: return false;
                default:
                    return true;
            }
        }

        private Element matchEquivalent() {
            switch (o) {
                case MATCH: return this;
                case INSERTION: return new Element(Operator.MATCH,0,false,false);
                case EMBEDDED:
                    int mtss = embedded.minimumTemplateSequenceSize();
                    int sc = embedded.starCount();
                    if (mtss > 0) {
                        if (sc > 0) {
                            final Civar newEmbedded = new Civar(Collections.unmodifiableList(Arrays.asList(
                                    new Element(Operator.MATCH,mtss,false,false),new Element(Operator.MATCH,sc,true,false))));
                            //newEmbedded.string = newEmbedded.elements.get(0).toString() + newEmbedded.elements.get(1).toString();
                            return new Element(newEmbedded,false);
                        } else {
                            return new Element(Operator.MATCH,mtss,false,false);
                        }
                    } else if (sc > 0) {
                        return new Element(Operator.MATCH,sc,true,false);

                    } else {
                        return new Element(Operator.MATCH,0,false,false);
                    }
                default:
                    return new Element(Operator.MATCH,size,expands,false);
            }
        }

        public Element mandatoryEquivalent() {
            final Element result = this.clone();
            result.optional = false;
            return result;
        }

        public int excessPaddingUsed(final int excessPadding) {
            return expands ? Math.min(excessPadding,size) : 0;
        }
    }

    protected static class Parser {

        public static Civar parse(final CharSequence cs, final int from, final int to) {
            if (cs == null)
                throw new NullPointerException("the input char-sequence cannot be null");
            if (from < 0)
                throw new IndexOutOfBoundsException("the from index cannot be negative");
            if (to < from)
                throw new IllegalArgumentException("the to index cannot less than the from index");
            if (to > cs.length())
                throw new IllegalArgumentException("the to index cannot be greater than the end of the sequence");
            if (cs == null)
                throw new IllegalArgumentException("cs cannot be null");

            final String s = cs.subSequence(from, to).toString();
            final LinkedList<Token> tokens = tokenize(s, 0, s.length());
            final LinkedList<Element> elements = elementize(tokens, s);
            final Civar result;
            if (elements.size() == 0) {
                result = new Civar((List<Element>) (List) Collections.emptyList());
            } else if (elements.size() == 1) {
                result = new Civar(Collections.singletonList(elements.getFirst()));
            } else {
                result = new Civar(Collections.unmodifiableList(Arrays.asList(elements.toArray(new Element[elements.size()]))));
            }
            //result.string = s;
            return result;
        }

        @Requires("tokens != null")
        private static LinkedList<Element> elementize(final LinkedList<Token> tokens, final CharSequence cs) {
            LinkedList<Element> elements = new LinkedList<>();
            Stack<Token> stack = new Stack<>();
            while (!tokens.isEmpty()) {
                Token nextToken = tokens.pollFirst();
                stack.push(nextToken);
                switch (nextToken.type) {
                    case OPERATOR:
                        if (nextToken.asOperator() == Operator.INSERTION) break;
                        stack.push(Token.element(popOperation(stack, null,cs)));
                        break;
                    case XMER:
                        stack.push(Token.element(popOperation(stack, stack.pop().asXmer(),cs)));
                        break;
                    case CLOSE_BRACKET:
                        stack.push(Token.element(popEmbedded(stack, cs)));
                        break;
                    case QMARK:
                        stack.pop();
                        final Token t;
                        stack.push(t = Token.element(popOperation(stack,null,cs)));
                        t.asElement().makeOptional();
                        break;
                    default:
                }
            }

            for (int i = 0; i < stack.size(); i++) {
                final Token t = stack.get(i);
                if (t.type == TokenType.ELEMENT) {
                    elements.add(t.asElement());
                } else if (t.type == TokenType.STAR) {
                    elements.add(new Element(Operator.MATCH,1,true, false, null));
                } else if (t.type == TokenType.NUMBER) {
                    if (i < stack.size() - 1) {
                        if (stack.get(i+1).type == TokenType.STAR) {
                            elements.add(new Element(Operator.MATCH,t.asNumber(),true, false, null));
                            i++;
                        } else {
                            elements.add(new Element(Operator.MATCH,t.asNumber(),false, false, null));
                        }
                    } else {
                        elements.add(new Element(Operator.MATCH,t.asNumber(),false, false, null));
                    }
                } else {
                   throw new IllegalArgumentException("Invalid Civar string: " + cs);
                }
            }
            return elements;
        }

        private static Element popEmbedded(final Stack<Token> stack, final CharSequence cs) {
            Token closeBracket = stack.pop();// remove close parentesis.
            LinkedList<Element> embeddedElements = new LinkedList<>();
            while (!stack.isEmpty()) {
                final Token nextToken = stack.pop();
                if (nextToken.type == TokenType.OPEN_BRACKET) {
                    final Civar embeddedCivar = new Civar(Collections.unmodifiableList(embeddedElements));
                    //embeddedCivar.string = cs.subSequence(nextToken.asNumber() + 1, closeBracket.asNumber()).toString();
                    return new Element(embeddedCivar,false);
                } else if (nextToken.type != TokenType.ELEMENT) {
                    throw new IllegalArgumentException("Civar format error");
                } else {
                    embeddedElements.add(0, nextToken.asElement());
                }
            }
            throw new IllegalArgumentException("Civar format error");
        }

        private static Element popOperation(final Stack<Token> stack, final String xmer, final CharSequence cs) {
            if (stack.isEmpty()) {
                throw new IllegalArgumentException("Invalid Civar string: " + cs);
            }
            Token operator = stack.pop();
            if (operator.type == TokenType.ELEMENT) {
                return operator.asElement();
            }
            if (operator.type != TokenType.OPERATOR) {
                throw new IllegalArgumentException("Invalid Civar string:" + operator.type +  " " + cs);
            }
            if (stack.isEmpty()) {
                return new Element(operator.asOperator(), 1, false, false, xmer);
            } else {
                Token sizeOrStar = stack.pop();
                if (sizeOrStar.type == TokenType.STAR) {
                    if (stack.isEmpty()) {
                        return new Element(operator.asOperator(), 1, true, false, xmer);
                    } else if (stack.peek().type == TokenType.NUMBER) {
                        return new Element(operator.asOperator(), stack.pop().asNumber(), true, false, xmer);
                    } else {
                        return new Element(operator.asOperator(), 1, true, false, xmer);
                    }
                } else if (sizeOrStar.type == TokenType.NUMBER) {
                    return new Element(operator.asOperator(), sizeOrStar.asNumber(), false, false, xmer);
                } else {
                    stack.push(sizeOrStar);
                    return new Element(operator.asOperator(), 1, false, false, xmer);
                }
            }
        }


        private static LinkedList<Token> tokenize(final CharSequence cs, final int from, final int to) {
            final LinkedList<Token> tokens = new LinkedList<>();
            int i = from;
            while (i < to) {
                char c = cs.charAt(i++);
                // NUMBER tokens
                if (Character.isDigit(c)) {
                    int num = 0;
                    do {
                        num = num * 10 + (c - '0');
                        if (i == to)
                            break;
                        c = cs.charAt(i++);
                    }
                    while (Character.isDigit(c) || (i-- == 0)); // || (i-- == 0) is a trick to "pushback" the first non digit character.
                    tokens.add(Token.number(num));
                } else if (c == '*') {
                    tokens.add(Token.star());
                } else if (c == '(') {
                    tokens.add(Token.openBracket(i - 1));
                } else if (c == ')') {
                    tokens.add(Token.closeBracket(i - 1));
                } else if (c == '?') {
                    tokens.add(Token.qmark());
                } else if (c == '=') {
                    tokens.add(Token.operator(Operator.MATCH));
                } else if (Character.isLetter(c)) {
                    if (Character.isUpperCase(c)) {
                        tokens.add(Token.operator(Operator.fromChar(c)));
                    } else {
                        int start = i - 1;
                        do {
                            if (i == to)
                                break;
                            c = cs.charAt(i++);
                        } while (Character.isLowerCase(c) || (i-- == 0));
                        tokens.add(Token.xmer(cs.subSequence(start, i)));
                    }
                } else {
                    throw new IllegalArgumentException("the Civar string contains invalid characters starting with '" + c + '"');
                }
            }
            return tokens;
        }


        public Civar parse(final CharSequence cs) {
            if (cs == null)
                throw new IllegalArgumentException("cs cannot be null");
            return parse(cs, 0, cs.length());
        }


    }

    /**
     * Transforms a civar into the equivalent Cigar.
     * @return never {@code null}.
     */
    public Cigar toCigar(final int templateLength) {

        int minSeqLen = minimumTemplateSequenceSize();
        if (!expands() && templateLength != minSeqLen)
            throw new IllegalArgumentException("the sequence provided does not match this Civar size " + templateLength + " != " + minSeqLen);
        if (templateLength < minSeqLen)
            throw new IllegalArgumentException("the sequence provided is too small for this Civar " + templateLength + " < " + minSeqLen);
        int starCount = starCount();
        int paddingTotal = templateLength - minSeqLen;
        int starPadding = starCount == 0 ? 0 : paddingTotal / starCount;
        int excessPadding = starCount == 0 ? paddingTotal : paddingTotal % starCount;

       // We first get the equivalent cigar elements for the elements in the Civar.
       final List<CigarElement> cigarElements = new LinkedList<>();

       for (final Element e : this.elements()) {
            final int size = e.size(starPadding,excessPadding);
            excessPadding -= e.excessPaddingUsed(excessPadding);

            switch (e.operator()) {
                case EMBEDDED:
                    cigarElements.addAll(e.embedded.toCigar(size).getCigarElements());
                    break;
                case MATCH:
                case TRANSITION:
                case COMPLEMENT:
                case TRANSVERSION:
                    cigarElements.add(new CigarElement(size, CigarOperator.M));
                    break;
                case INSERTION:
                    cigarElements.add(new CigarElement(size,CigarOperator.I));
                    break;
                case DELETION:
                    cigarElements.add(new CigarElement(size,CigarOperator.D));
                    break;
                default:
            }
       }

       // No we look for consequitive elements with the same operator and we merge them.
       final ListIterator<CigarElement> it = cigarElements.listIterator();
       while (it.hasNext()) {
            final CigarElement thisElement = it.next();
            if (!it.hasNext())
                continue;
            final CigarElement nextElement = it.next();
            if (thisElement.getOperator() == nextElement.getOperator()) {
                final int nextLength = nextElement.getLength();
                it.remove();
                it.previous();
                it.set(new CigarElement(thisElement.getLength() + nextLength, thisElement.getOperator()));
            } else
                it.previous();
       }
       return new Cigar(cigarElements);
    }



    protected static enum TokenType {
        OPERATOR, XMER, NUMBER, STAR, OPEN_BRACKET, CLOSE_BRACKET, ELEMENT, START, END, QMARK;
    }

    protected static class Token {
        public final TokenType type;
        public final Object content;

        @Requires("type != null && content != null")
        protected Token(final TokenType type, final Object content) {
            this.type = type;
            this.content = content;

        }

        public String toString() {
            switch (this.type) {
                case STAR: return "*";
                case OPEN_BRACKET: return "(";
                case CLOSE_BRACKET: return ")";
                default:
                    return String.valueOf(content);
            }
        }


        public Operator asOperator() {
            return (Operator) content;
        }

        public String asXmer() {
            return ((CharSequence) content).toString();
        }

        public int asNumber() {
            return ((Number) content).intValue();
        }

        public Element asElement() {
            return ((Element) content);
        }

        protected static Token xmer(final CharSequence cs) {
            return new Token(TokenType.XMER, cs);

        }

        protected static Token operator(final Operator o) {
            return new Token(TokenType.OPERATOR, o);
        }

        protected static Token number(final int n) {
            return new Token(TokenType.NUMBER, n);
        }



        protected static Token star() {
            return new Token(TokenType.STAR, null);
        }

        protected static Token openBracket(int offset) {
            return new Token(TokenType.OPEN_BRACKET, offset);
        }

        protected static Token closeBracket(int offset) {
            return new Token(TokenType.CLOSE_BRACKET, offset);
        }

        protected static Token element(final Element e) {
            return new Token(TokenType.ELEMENT, e);
        }


        public static Token qmark() {
            return new Token(TokenType.QMARK, null);
        }
    }


    public static class ElementOffset {
        public final Element element;
        public final int sequenceFrom;
        public final int sequenceTo;
        public final int templateFrom;
        public final int templateTo;


        @Requires("e != null && from >= 0 && to >= from && tFrom >= 0 && tTo >= tFrom")
        protected ElementOffset(final Element e, final int from, final int to, final int tFrom, final int tTo) {
            element = e;
            sequenceFrom = from;
            sequenceTo = to;
            templateFrom = tFrom;
            templateTo = tTo;
        }
    }
}
