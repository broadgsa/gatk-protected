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

import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

/**
 * Merges VariantContexts in a series of haplotypes according to their pairwise LD
 *
 * User: depristo
 * Date: 3/28/13
 * Time: 6:17 PM
 */
public class LDMerger extends MergeVariantsAcrossHaplotypes {
    private final static Logger logger = Logger.getLogger(LDMerger.class);

    private final boolean DEBUG;
    private final int minSamplesToMergeSNPs;
    private final int minSamplesToMergeOtherEvents;

    public LDMerger(boolean DEBUG, int minSamplesToMergeSNPs, int minSamplesToMergeOtherEvents) {
        super();
        this.DEBUG = DEBUG;
        this.minSamplesToMergeSNPs = minSamplesToMergeSNPs;
        this.minSamplesToMergeOtherEvents = minSamplesToMergeOtherEvents;
    }

    protected LDMerger() {
        this(false, 1, 1);
    }

    // TODO -- should be class arguments and static variables in HC
    protected final static int MAX_DISTANCE_BETWEEN_SNPS_TO_MERGE = 6;
    protected final static int MAX_DISTANCE_BETWEEN_OTHER_EVENTS_TO_MERGE = 25;

    /**
     * We require 99% confidence that only the phased haplotypes exist in the population to merge the records
     */
    protected final static double MERGE_EVENTS_PROB_PHASED_THRESHOLD = 0.99;

    /**
     * Merge as many events among the haplotypes as possible based on pairwise LD among variants
     *
     * @param haplotypes a list of haplotypes whose events we want to merge
     * @param readLikelihoods map from sample name -> read likelihoods for each haplotype
     * @param startPosKeySet a set of starting positions of all events among the haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     */
    @Override
    public boolean merge( final List<Haplotype> haplotypes,
                          final ReadLikelihoods<Haplotype> readLikelihoods,
                          final TreeSet<Integer> startPosKeySet,
                          final byte[] ref,
                          final GenomeLoc refLoc ) {
        if ( haplotypes == null ) throw new IllegalArgumentException("haplotypes cannot be null");
        if ( readLikelihoods == null ) throw new IllegalArgumentException("readLikelihoods cannot be null");
        if ( startPosKeySet == null ) throw new IllegalArgumentException("startPosKeySet cannot be null");
        if ( ref == null ) throw new IllegalArgumentException("ref cannot be null");
        if ( refLoc == null ) throw new IllegalArgumentException("refLoc cannot be null");
        if ( refLoc.size() != ref.length ) throw new IllegalArgumentException("refLoc size " + refLoc.size() + " != ref.length " + ref.length + " at " + refLoc);

        if( startPosKeySet.size() <= 1 ) { return false; }

        final int nSamples = readLikelihoods.sampleCount();
        final HaplotypeLDCalculator r2Calculator = new HaplotypeLDCalculator(haplotypes, readLikelihoods);
        boolean somethingWasMerged = false;
        boolean mapWasUpdated = true;
        while( mapWasUpdated ) {
            mapWasUpdated = mergeConsecutiveEventsBasedOnLDOnce(haplotypes, r2Calculator, nSamples, startPosKeySet, ref, refLoc);
            somethingWasMerged |= mapWasUpdated;
        }
        return somethingWasMerged;
    }

    /**
     * Merge the next pair of events, if possible
     *
     * @param haplotypes a list of haplotypes whose events we want to merge
     * @param ldCalculator calculates R^2 for pairs of events on demand
     * @param startPosKeySet a set of starting positions of all events among the haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     * @return true if something was merged, false otherwise
     */
    protected boolean mergeConsecutiveEventsBasedOnLDOnce( final List<Haplotype> haplotypes,
                                                           final HaplotypeLDCalculator ldCalculator,
                                                           final int nSamples,
                                                           final TreeSet<Integer> startPosKeySet,
                                                           final byte[] ref,
                                                           final GenomeLoc refLoc ) {
        // loop over the set of start locations and consider pairs that start near each other
        final Iterator<Integer> iter = startPosKeySet.iterator();
        int thisStart = iter.next();
        while( iter.hasNext() ) {
            final int nextStart = iter.next();
            final LDMergeData toMerge = getPairOfEventsToMerge(haplotypes, thisStart, nextStart);

            if ( toMerge.canBeMerged(nSamples) ) {
                final double pPhased = ldCalculator.computeProbOfBeingPhased(toMerge.firstVC, toMerge.secondVC);

                if( DEBUG ) {
                    logger.info("Found consecutive biallelic events with R^2 = " + String.format("%.4f", pPhased));
                    logger.info("-- " + toMerge.firstVC);
                    logger.info("-- " + toMerge.secondVC);
                }

                if( pPhased > MERGE_EVENTS_PROB_PHASED_THRESHOLD) {
                    final VariantContext mergedVC = createMergedVariantContext(toMerge.firstVC, toMerge.secondVC, ref, refLoc);
                    // if for some reason the merging resulting in a bad allele, mergedVC will be null, and we will just remove first and second
                    replaceVariantContextsInMap(haplotypes, startPosKeySet, mergedVC, toMerge.firstVC, toMerge.secondVC);
                    return true; // break out of tree set iteration since it was just updated, start over from the beginning and keep merging events
                }
            }

            thisStart = nextStart;
        }

        return false;
    }

    /**
     * Info about potential LD merge of two variant contexts
     */
    private class LDMergeData {
        VariantContext firstVC = null, secondVC = null;
        boolean canBeMerged = true;

        /** Tell this object that it cant be merged for some reason */
        public LDMergeData cantBeMerged() {
            canBeMerged = false;
            return this;
        }

        /**
         * Can these two events be merged
         * @param nSamples the number of samples we're considering
         * @return true if we can merge our two variant contexts
         */
        public boolean canBeMerged(final int nSamples) {
            if ( ! canBeMerged || firstVC == null || secondVC == null )
                return false;

            final int distance = secondVC.getStart() - firstVC.getEnd();
            if ( firstVC.isSNP() && secondVC.isSNP() ) {
                return nSamples >= minSamplesToMergeSNPs && distance <= MAX_DISTANCE_BETWEEN_SNPS_TO_MERGE;
            } else {
                return nSamples >= minSamplesToMergeOtherEvents && distance <= MAX_DISTANCE_BETWEEN_OTHER_EVENTS_TO_MERGE;
            }
        }
    }

    /**
     * Get the information about the potential merge of two events starting at thisStart and nextStart
     * @param haplotypes our haplotypes
     * @param thisStart the starting position of the first event to merge
     * @param nextStart the starting position of the next event to merge
     * @return never {@code null}.
     */
    private LDMergeData getPairOfEventsToMerge(final List<Haplotype> haplotypes, final int thisStart, final int nextStart) {
        final LDMergeData mergeData = new LDMergeData();

        for( final Haplotype h : haplotypes ) {
            // only make complex substitutions out of consecutive biallelic sites
            final VariantContext thisHapVC = h.getEventMap().get(thisStart);
            if( thisHapVC != null && !thisHapVC.isSymbolic() ) { // something was found at this location on this haplotype
                if( mergeData.firstVC == null ) {
                    mergeData.firstVC = thisHapVC;
                } else if( !thisHapVC.hasSameAllelesAs( mergeData.firstVC) ) {
                    return mergeData.cantBeMerged();
                }
            }
            final VariantContext nextHapVC = h.getEventMap().get(nextStart);
            if( nextHapVC != null && !nextHapVC.isSymbolic() ) { // something was found at the next location on this haplotype
                if( mergeData.secondVC == null ) {
                    mergeData.secondVC = nextHapVC;
                } else if( !nextHapVC.hasSameAllelesAs( mergeData.secondVC) ) {
                    return mergeData.cantBeMerged();
                }
            }
        }

        // don't try to merge overlapping events
        if ( mergeData.firstVC != null && mergeData.secondVC != null && mergeData.firstVC.getEnd() >= mergeData.secondVC.getStart() )
            return mergeData.cantBeMerged();

        return mergeData;
    }

    // BUGBUG: make this merge function more general
    protected VariantContext createMergedVariantContext( final VariantContext thisVC, final VariantContext nextVC, final byte[] ref, final GenomeLoc refLoc ) {
        final int thisStart = thisVC.getStart();
        final int nextStart = nextVC.getStart();
        byte[] refBases = new byte[]{};
        byte[] altBases = new byte[]{};
        refBases = ArrayUtils.addAll(refBases, thisVC.getReference().getBases());
        altBases = ArrayUtils.addAll(altBases, thisVC.getAlternateAllele(0).getBases());
        int locus;
        for( locus = thisStart + refBases.length; locus < nextStart; locus++ ) {
            final byte refByte = ref[locus - refLoc.getStart()];
            refBases = ArrayUtils.add(refBases, refByte);
            altBases = ArrayUtils.add(altBases, refByte);
        }
        refBases = ArrayUtils.addAll(refBases, ArrayUtils.subarray(nextVC.getReference().getBases(), locus > nextStart ? 1 : 0, nextVC.getReference().getBases().length)); // special case of deletion including the padding base of consecutive indel
        altBases = ArrayUtils.addAll(altBases, nextVC.getAlternateAllele(0).getBases());

        int iii = 0;
        if( refBases.length == altBases.length ) { // insertion + deletion of same length creates an MNP --> trim common prefix bases off the beginning of the allele
            while( iii < refBases.length && refBases[iii] == altBases[iii] ) { iii++; }
            if ( iii == refBases.length ) {
                // we've become a null allele, such as with CA/C + A/AA -> CA/CA => after trimming there's nothing left
                // so return a null variant context so we can eliminate the variants from consideration
                return null;
            }
        }


        final Allele refAllele = Allele.create( ArrayUtils.subarray(refBases, iii, refBases.length), true );
        final Allele altAllele =  Allele.create( ArrayUtils.subarray(altBases, iii, altBases.length), false );
        return new VariantContextBuilder("merged", thisVC.getChr(), thisVC.getStart() + iii, nextVC.getEnd(), Arrays.asList(refAllele, altAllele)).make();
    }

    /**
     * Update the event maps in all haplotypes to replace a replacement of update1 and 2 with replacement
     *
     * @param haplotypes the haplotypes whose event maps we need to update
     * @param startPosKeySet a sorted set of start positions that we must update
     * @param replacement a VariantContext to replace update1 and update2 with.  Can be null, indicating that we just want to remove update1 and update2
     * @param update1 the first VC we want to update
     * @param update2 the second VC we want to update
     */
    private void replaceVariantContextsInMap(final List<Haplotype> haplotypes,
                                             final TreeSet<Integer> startPosKeySet,
                                             final VariantContext replacement,
                                             final VariantContext update1, final VariantContext update2) {
        // remove the old event from the eventMap on every haplotype and the start pos key set, replace with merged event
        for( final Haplotype h : haplotypes ) {
            // if we had both events, add replacement.  In some cases the haplotype may not have both
            // events but they were still merged because the haplotype isn't a particularly informative
            // haplotype in any case.  The order of operations here is important because we are modifying the map
            final boolean shouldAdd = h.getEventMap().containsKey(update1.getStart()) && h.getEventMap().containsKey(update2.getStart());
            h.getEventMap().remove(update1.getStart());
            h.getEventMap().remove(update2.getStart());
            if ( shouldAdd && replacement != null ) {
                h.getEventMap().addVC(replacement, false); // cannot merge we other events at the same position
            }
        }

        startPosKeySet.remove(update1.getStart());
        startPosKeySet.remove(update2.getStart());
        if ( replacement != null ) startPosKeySet.add(replacement.getStart());
    }
}
