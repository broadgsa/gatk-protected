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

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.collections.CountSet;
import org.broadinstitute.gatk.utils.haplotype.EventMap;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.haplotype.HaplotypeSizeAndBaseComparator;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;

/**
 * Collection of read assembly using several kmerSizes.
 *
 * <p>
 *     There could be a different assembly per each kmerSize. In turn, haplotypes are result of one of those
 *     assemblies.
 * </p>
 *
 * <p>
 *     Where there is more than one possible kmerSize that generates a haplotype we consider the smaller one.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.com&gt;
 */
public  class AssemblyResultSet {

    private final Map<Integer,AssemblyResult> assemblyResultByKmerSize;
    private final Set<Haplotype> haplotypes;
    private final Map<Haplotype,AssemblyResult> assemblyResultByHaplotype;
    private ActiveRegion regionForGenotyping;
    private byte[] fullReferenceWithPadding;
    private GenomeLoc paddedReferenceLoc;
    private boolean variationPresent;
    private Haplotype refHaplotype;
    private boolean wasTrimmed = false;
    private final CountSet kmerSizes;
    private TreeSet<VariantContext> variationEvents;
    private boolean debug;
    private static Logger logger = Logger.getLogger(AssemblyResultSet.class);

    /**
     * Constructs a new empty assembly result set.
     */
    public AssemblyResultSet() {
        assemblyResultByKmerSize = new LinkedHashMap<>(4);
        haplotypes = new LinkedHashSet<>(10);
        assemblyResultByHaplotype = new LinkedHashMap<>(10);
        kmerSizes = new CountSet(4);
    }


    /**
     * Change the debug status for this assembly-result-set.
     * @param newValue new value for the debug status.
     */
    void setDebug(final boolean newValue) {
        debug = newValue;
    }

    /**
     * Trims an assembly result set down based on a new set of trimmed haplotypes.
     *
     * @param trimmedActiveRegion the trimmed down active region.
     *
     * @throws NullPointerException if any argument in {@code null} or
     *      if there are {@code null} entries in {@code originalByTrimmedHaplotypes} for trimmed haplotype keys.
     * @throws IllegalArgumentException if there is no reference haplotype amongst the trimmed ones.
     *
     * @return never {@code null}, a new trimmed assembly result set.
     */
    public AssemblyResultSet trimTo(final ActiveRegion trimmedActiveRegion) {

        final Map<Haplotype,Haplotype> originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(trimmedActiveRegion);
        if (refHaplotype == null) throw new IllegalStateException();
        if (trimmedActiveRegion == null) throw new NullPointerException();
        final AssemblyResultSet result = new AssemblyResultSet();

        for (final Haplotype trimmed : originalByTrimmedHaplotypes.keySet()) {
            final Haplotype original = originalByTrimmedHaplotypes.get(trimmed);
            if (original == null)
                throw new NullPointerException("all trimmed haplotypes must have an original one");
            final AssemblyResult as = assemblyResultByHaplotype.get(original);
            if (as == null) result.add(trimmed); else result.add(trimmed, as);
        }

        result.setRegionForGenotyping(trimmedActiveRegion);
        result.setFullReferenceWithPadding(this.fullReferenceWithPadding);
        result.setPaddedReferenceLoc(this.paddedReferenceLoc);
        if (result.refHaplotype == null)
            throw new IllegalStateException("missing reference haplotype in the trimmed set");
        result.wasTrimmed = true;
        return result;
    }

    private Map<Haplotype, Haplotype> calculateOriginalByTrimmedHaplotypes(final ActiveRegion trimmedActiveRegion) {
        if ( debug ) logger.info("Trimming active region " + getRegionForGenotyping() + " with " + getHaplotypeCount() + " haplotypes");

        final List<Haplotype> haplotypeList = getHaplotypeList();

        // trim down the haplotypes
        final Map<Haplotype,Haplotype> originalByTrimmedHaplotypes = new HashMap<>();

        for ( final Haplotype h : haplotypeList ) {
            final Haplotype trimmed = h.trim(trimmedActiveRegion.getExtendedLoc());

            if ( trimmed != null ) {
                if (originalByTrimmedHaplotypes.containsKey(trimmed)) {
                    if (trimmed.isReference()) {
                        originalByTrimmedHaplotypes.remove(trimmed);
                        originalByTrimmedHaplotypes.put(trimmed, h);
                    }
                } else
                    originalByTrimmedHaplotypes.put(trimmed,h);
            } else if (h.isReference())
                throw new IllegalStateException("trimming eliminates the reference haplotype");
            else if ( debug ) {
                logger.info("Throwing out haplotype " + h + " with cigar " + h.getCigar() +
                        " because it starts with or ends with an insertion or deletion when trimmed to " +
                        trimmedActiveRegion.getExtendedLoc());
            }
        }

        // create the final list of trimmed haplotypes
        final List<Haplotype> trimmedHaplotypes = new ArrayList<>(originalByTrimmedHaplotypes.keySet());

        // resort the trimmed haplotypes.
        Collections.sort(trimmedHaplotypes,new HaplotypeSizeAndBaseComparator());
        final Map<Haplotype,Haplotype> sortedOriginalByTrimmedHaplotypes = new LinkedHashMap<>(trimmedHaplotypes.size());
        for (final Haplotype trimmed : trimmedHaplotypes)
            sortedOriginalByTrimmedHaplotypes.put(trimmed,originalByTrimmedHaplotypes.get(trimmed));


        if ( debug ) logger.info("Trimmed region to " + trimmedActiveRegion.getLocation() + " size " +
                trimmedActiveRegion.getLocation().size() + " reduced number of haplotypes from " +
                haplotypeList.size() + " to only " + trimmedHaplotypes.size());
        if ( debug )
            for ( final Haplotype remaining: trimmedHaplotypes )
                logger.info("Remains: " + remaining + " cigar " + remaining.getCigar());
        return sortedOriginalByTrimmedHaplotypes;
    }

    /**
     * Query the reference haplotype in the result set.
     * @return {@code null} if none wasn't yet added, otherwise a reference haplotype.
     */
    public Haplotype getReferenceHaplotype() {
        return refHaplotype;
    }

    /**
     * Checks whether there is any variation present in the assembly result set.
     *
     * <p>
     *     This is equivalent to whether there is more than one haplotype.
     * </p>
     *
     * @return {@code true} if there is variation present, {@code false} otherwise.
     */
    public boolean isVariationPresent() {
        return variationPresent && haplotypes.size() > 1;
    }

    /**
     * Dumps debugging information into a print-writer.
     *
     * @param pw where to dump the information.
     *
     * @throws NullPointerException if {@code pw} is {@code null}.
     */
    public void debugDump(final PrintWriter pw) {
        if (getHaplotypeList().size() == 0) {
            return;
        }
        pw.println("Active Region " + this.regionForGenotyping.getLocation());
        pw.println("Extended Act Region " + this.getRegionForGenotyping().getExtendedLoc());
        pw.println("Ref haplotype coords " + getHaplotypeList().get(0).getGenomeLocation());
        pw.println("Haplotype count " + haplotypes.size());
        final Map<Integer,Integer> kmerSizeToCount = new HashMap<>();

        for (final Map.Entry<Haplotype,AssemblyResult> e : assemblyResultByHaplotype.entrySet()) {
            final AssemblyResult as = e.getValue();
            final int kmerSize = as.getGraph().getKmerSize();
            if (kmerSizeToCount.containsKey(kmerSize)) {
                kmerSizeToCount.put(kmerSize,kmerSizeToCount.get(kmerSize) + 1);
            } else {
                kmerSizeToCount.put(kmerSize,1);
            }
        }
        pw.println("Kmer sizes count " + kmerSizeToCount.entrySet().size() );
        Integer[] kmerSizes = new Integer[kmerSizeToCount.size()];
        kmerSizes = kmerSizeToCount.keySet().toArray(kmerSizes);
        Arrays.sort(kmerSizes);
        pw.println("Kmer sizes values " + Arrays.toString(kmerSizes));
        for (int size : kmerSizes) {
            pw.println("Kmer size " + size + " count " + kmerSizeToCount.get(size));
        }
    }

    /**
     * Adds a haplotype to the result set without indicating a generating assembly result.
     *
     * <p>
     *     It is possible to call this method with the same haplotype several times. In that the second and further
     *     calls won't have any effect (thus returning {@code false}).
     * </p>
     *
     * @param h the haplotype to add to the assembly result set.
     *
     * @throws NullPointerException if {@code h} is {@code null}
     * @throws IllegalArgumentException if {@code h} does not have a genome location.
     *
     * @return {@code true} if the assembly result set has been modified as a result of this call.
     */
    public boolean add(final Haplotype h) {
        if (h == null) throw new NullPointerException("input haplotype cannot be null");
        if (h.getGenomeLocation() == null)
            throw new IllegalArgumentException("the haplotype provided must have a genomic location");
        if (haplotypes.contains(h))
            return false;
        haplotypes.add(h);
        updateReferenceHaplotype(h);
        return true;
    }

    /**
     * Adds simultaneously a haplotype and the generating assembly-result.
     *
     * <p>
     *     Haplotypes and their assembly-result can be added multiple times although just the first call will have
     *     any effect (return value is {@code true}).
     * </p>
     *
     *
     * @param h haplotype to add.
     * @param ar assembly-result that is assumed to have given rise to that haplotype.
     *
     * @throws NullPointerException if {@code h} or {@code ar} is {@code null}.
     * @throws IllegalArgumentException if {@code h} has not defined genome location.
     *
     * @return {@code true} iff this called changes the assembly result set.
     */
    public boolean add(final Haplotype h, final AssemblyResult ar) {
        if (h == null) throw new NullPointerException("input haplotype cannot be null");
        if (ar == null) throw new NullPointerException("input assembly-result cannot be null");
        if (h.getGenomeLocation() == null)
            throw new IllegalArgumentException("the haplotype provided must have a genomic location");

        final boolean assemblyResultAdditionReturn =  add(ar);

        if (haplotypes.contains(h)) {
            final AssemblyResult previousAr = assemblyResultByHaplotype.get(h);
            if (previousAr == null) {
                assemblyResultByHaplotype.put(h, ar);
                return true;
            } else if (!previousAr.equals(ar))
                throw new IllegalStateException("there is already a different assembly result for the input haplotype");
            else
                return assemblyResultAdditionReturn;
        } else {
            haplotypes.add(h);
            assemblyResultByHaplotype.put(h,ar);
            updateReferenceHaplotype(h);
            if (h.isNonReference()) variationPresent = true;
            return true;
        }
    }

    /**
     * Add a assembly-result object.
     *
     * @param ar the assembly result to add.
     *
     * @throws NullPointerException if {@code ar} is {@code null}.
     * @throws IllegalStateException if there is an assembly result with the same kmerSize.
     * @return {@code true} iff this addition changed the assembly result set.
     */
    public boolean add(final AssemblyResult ar) {
        if (ar == null)
            throw new NullPointerException();
        final int kmerSize = ar.getKmerSize();
        if (assemblyResultByKmerSize.containsKey(kmerSize)) {
            if (!assemblyResultByKmerSize.get(kmerSize).equals(ar))
                throw new IllegalStateException("a different assembly result with the same kmerSize was already added");
            return false;
        } else {
            assemblyResultByKmerSize.put(kmerSize, ar);
            kmerSizes.add(kmerSize);
            return true;
        }
    }

    /**
     * Returns the current region for genotyping.
     *
     * @return might be {@code null}.
     */
    public ActiveRegion getRegionForGenotyping() {
        return regionForGenotyping;
    }

    /**
     * Sets the region for genotyping.
     *
     * @param regionForGenotyping the new value.
     */
    public void setRegionForGenotyping(final ActiveRegion regionForGenotyping) {
        this.regionForGenotyping = regionForGenotyping;
    }

    /**
     * Returns the current full reference with padding.
     *
     * @return might be {@code null}.
     */
    public byte[] getFullReferenceWithPadding() {
        return fullReferenceWithPadding;
    }

    /**
     * Sets the full reference with padding base sequence.
     *
     * @param fullReferenceWithPadding the new value.
     */
    public void setFullReferenceWithPadding(final byte[] fullReferenceWithPadding) {
        this.fullReferenceWithPadding = fullReferenceWithPadding;
    }

    /**
     * Returns the padded reference location.
     *
     * @return might be {@code null}
     */
    public GenomeLoc getPaddedReferenceLoc() {
        return paddedReferenceLoc;
    }

    /**
     * Changes the padded reference location.
     * @param paddedReferenceLoc the new value.
     */
    public void setPaddedReferenceLoc(final GenomeLoc paddedReferenceLoc) {
        this.paddedReferenceLoc = paddedReferenceLoc;
    }

    /**
     * Returns the number of haplotypes in the assembly result set.
     * @return {@code 0} or greater.
     */
    public int getHaplotypeCount() {
        return haplotypes.size();
    }

    /**
     * Returns the haplotypes as a list.
     *
     * <p>
     *     The result is unmodifiable.
     * </p>
     *
     * @return never {@code null}, but perhaps a empty list if no haplotype was generated during assembly.
     */
    public List<Haplotype> getHaplotypeList() {
        return Arrays.asList(haplotypes.toArray(new Haplotype[haplotypes.size()]));
    }

    /**
     * Returns the maximum kmerSize available.
     *
     * @throws IllegalStateException if no assembly-result was added to the set, thus there is no kmerSize.
     *
     * @return greater than 0.
     */
    public int getMaximumKmerSize() {
        if (kmerSizes.size() == 0)
            throw new IllegalStateException("there is yet no kmerSize in this assembly result set");
        return kmerSizes.max();
    }

    /**
     * Indicates whether there are more than one kmerSize in the set.
     *
     * @return {@code true} iff there is more than one kmerSize assembly in the set.
     */
    public boolean hasMultipleKmerSizes() {
        return kmerSizes.size() > 1;
    }

    /**
     * Returns the minimum kmerSize available.
     *
     * @throws IllegalStateException if no assembly-result was added to the set, thus there is no kmerSize.
     *
     * @return greater than 0.
     */
    public int getMinimumKmerSize() {
        if (kmerSizes.size() == 0)
            throw new IllegalStateException("there is yet no kmerSize in this assembly result set");
        return kmerSizes.min();
    }

    /**
     * Returns a read-threading graph in the assembly set that has a particular kmerSize.
     *
     * @param kmerSize the requested kmerSize.
     *
     * @return {@code null} if there is no read-threading-graph amongst assembly results with that kmerSize.
     */
    public ReadThreadingGraph getUniqueReadThreadingGraph(final int kmerSize) {
        final AssemblyResult assemblyResult = assemblyResultByKmerSize.get(kmerSize);
        if (assemblyResult == null) return null;
        return assemblyResult.getThreadingGraph();
    }

    /**
     * Checks whether this assembly result set was trimmed.
     *
     * @return {@code true} iff this assembly result set was trimmed.
     */
    public boolean wasTrimmed() {
        return wasTrimmed;
    }

    /**
     * Marks the assembly as not having variation even if it has more than one haplotype.
     */
    public void resetVariationPresent() {
        variationPresent = false;
    }

    /**
     * Dumps debugging information into a logger.
     *
     * @param logger where to dump the information.
     *
     * @throws NullPointerException if {@code logger} is {@code null}.
     */
    public void debugDump(final Logger logger) {
        final StringWriter sw = new StringWriter();
        final PrintWriter pw = new PrintWriter(sw);
        debugDump(pw);
        final String str = sw.toString();
        final String[] lines = str.split("\n");
        for (final String line : lines) {
            if (line.isEmpty()) {
                continue;
            }
            logger.debug(line);
        }
    }

    /**
     * Given whether a new haplotype that has been already added to {@link #haplotypes} collection is the
     * reference haplotype and updates {@link #refHaplotype} accordingly.
     *
     * <p>
     *     This method assumes that the colling code has verified that the haplotype was not already in {@link #haplotypes}
     *     I.e. that it is really a new one. Otherwise it will result in an exception if it happen to be a reference
     *     haplotype and this has already be set. This is the case even if the new haplotypes and the current reference
     *     are equal.
     * </p>
     *
     * @param newHaplotype the new haplotype.
     * @throws NullPointerException if {@code newHaplotype} is {@code null}.
     * @throws IllegalStateException if there is already a reference haplotype.
     */
    private void updateReferenceHaplotype(final Haplotype newHaplotype) {
        if (!newHaplotype.isReference()) return;
        if (refHaplotype == null)
            refHaplotype = newHaplotype;
        else // assumes that we have checked wether the haplotype is already in the collection and so is no need to check equality.
            throw new IllegalStateException("the assembly-result-set already have a reference haplotype that is different");
    }

    /**
     * Returns a sorted set of variant events that best explain the haplotypes found by the assembly
     * across kmerSizes.
     *
     * <p/>
     * The result is sorted incrementally by location.
     *
     * @return never {@code null}, but perhaps an empty collection.
     */
    public TreeSet<VariantContext> getVariationEvents() {
        if (variationEvents == null) {
            final List<Haplotype> haplotypeList = getHaplotypeList();
            EventMap.buildEventMapsForHaplotypes(haplotypeList,fullReferenceWithPadding,paddedReferenceLoc,debug);
            variationEvents = EventMap.getAllVariantContexts(haplotypeList);
        }
        return variationEvents;
    }
}
