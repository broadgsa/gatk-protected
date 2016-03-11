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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.tools.walkers.annotator.AlleleSpecificAnnotationData;
import org.broadinstitute.gatk.tools.walkers.annotator.ReducibleAnnotationData;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Variant context utilities related to merging variant-context instances.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReferenceConfidenceVariantContextMerger {

    private final static Logger logger = Logger.getLogger(ReferenceConfidenceVariantContextMerger.class);

    private static Comparable combineAnnotationValues( final List<Comparable> array ) {
        return MathUtils.median(array); // right now we take the median but other options could be explored
    }

    /**
     * Merges VariantContexts from gVCFs into a single hybrid.
     * Assumes that none of the input records are filtered.
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning (i.e. don't start at the location in loc); if null, we'll return null in this case
     * @param removeNonRefSymbolicAllele if true, remove the <NON_REF> allele from the merged VC
     * @param samplesAreUniquified  if true, sample names have been uniquified
     * @return new VariantContext representing the merge of all VCs or null if it not relevant
     */
    public static VariantContext merge(final List<VariantContext> VCs, final GenomeLoc loc, final Byte refBase, final boolean removeNonRefSymbolicAllele,
                                       final boolean samplesAreUniquified, final VariantAnnotatorEngine annotatorEngine) {
        // this can happen if e.g. you are using a dbSNP file that spans a region with no gVCFs
        if ( VCs == null || VCs.isEmpty() ) {
            return null; }

        // establish the baseline info (sometimes from the first VC)
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();

        // ref allele
        final Allele refAllele = determineReferenceAlleleGivenReferenceBase(VCs, loc, refBase);
        if ( refAllele == null ) {
            return null; }

        // FinalAlleleSet contains the alleles of the new resulting VC
        // Using linked set in order to guarantee a stable order
        final LinkedHashSet<Allele> finalAlleleSet = new LinkedHashSet<>(10);
        // Reference goes first
        finalAlleleSet.add(refAllele);

        final Map<String, Object> attributes = new LinkedHashMap<>();
        final Set<String> rsIDs = new LinkedHashSet<>(1); // most of the time there's one id
        int depth = 0;
        final Map<String, List<ReducibleAnnotationData>> annotationMap = new LinkedHashMap<>();
        final GenotypesContext genotypes = GenotypesContext.create();

        // In this list we hold the mapping of each variant context alleles.
        final List<Pair<VariantContext,List<Allele>>> vcAndNewAllelePairs = new ArrayList<>(VCs.size());
        // Keep track of whether we saw a spanning deletion and a non-spanning event
        boolean sawSpanningDeletion = false;
        boolean sawNonSpanningEvent = false;

        // cycle through and add info from the other VCs
        for ( final VariantContext vc : VCs ) {

            // if this context doesn't start at the current location then it must be a spanning event (deletion or ref block)
            final boolean isSpanningEvent = loc.getStart() != vc.getStart();
            // record whether it's also a spanning deletion/event (we know this because the VariantContext type is no
            // longer "symbolic" but "mixed" because there are real alleles mixed in with the symbolic non-ref allele)
            sawSpanningDeletion |= ( isSpanningEvent && vc.isMixed() ) || vc.getAlternateAlleles().contains(Allele.SPAN_DEL) ||
                    vc.getAlternateAlleles().contains(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED );
            sawNonSpanningEvent |= ( !isSpanningEvent && vc.isMixed() );

            vcAndNewAllelePairs.add(new Pair<>(vc, isSpanningEvent ? replaceWithNoCallsAndDels(vc) : remapAlleles(vc, refAllele, finalAlleleSet)));
        }

        // Add <DEL> and <NON_REF> to the end if at all required in the output.
        if ( sawSpanningDeletion && (sawNonSpanningEvent || !removeNonRefSymbolicAllele) ) finalAlleleSet.add(Allele.SPAN_DEL);
        if (!removeNonRefSymbolicAllele) finalAlleleSet.add(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);

        final List<Allele> allelesList = new ArrayList<>(finalAlleleSet);

        for ( final Pair<VariantContext,List<Allele>> pair : vcAndNewAllelePairs ) {
            final VariantContext vc = pair.getFirst();
            final List<Allele> remappedAlleles = pair.getSecond();

            mergeRefConfidenceGenotypes(genotypes, vc, remappedAlleles, allelesList, samplesAreUniquified);

            // special case DP (add it up) for all events
            if ( vc.hasAttribute(VCFConstants.DEPTH_KEY) ) {
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            } else { // handle the gVCF case from the HaplotypeCaller
                for( final Genotype gt : vc.getGenotypes() ) {
                    depth += (gt.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) ? Integer.parseInt((String)gt.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)) : (gt.hasDP() ? gt.getDP() : 0));
                }
            }

            if ( loc.getStart() != vc.getStart() ) {
                continue;
            }

            // special case ID (just preserve it)
            if ( vc.hasID() ) rsIDs.add(vc.getID());

            // add attributes to annotationMap, store all info field annotations as AlleleSpecificAnnotationData in case they can be parsed that way
            addReferenceConfidenceAttributes(pair, annotationMap);
        }

        //combine the annotations that are reducible and remove them from annotationMap
        Map<String, Object> combinedAnnotations = new HashMap<>();
        if (annotatorEngine != null) {
            combinedAnnotations = annotatorEngine.combineAnnotations(allelesList, annotationMap);
        }
        attributes.putAll(combinedAnnotations);

        // remove stale AC and AF based attributes (including MLEAC and MLEAF lists)
        //these will be recalculated after genotyping
        removeStaleAttributesAfterMerge(annotationMap);

        //annotatorEngine.combineAnnotations removed the successfully combined annotations, so now parse those that are left
        //here we're assuming that things that are left are scalars per sample
        Map<String, List<Comparable>> parsedAnnotationMap = parseRemainingAnnotations(annotationMap);

        // when combining remaining annotations use the median value from all input VCs which had annotations provided
        for ( final Map.Entry<String, List<Comparable>> p : parsedAnnotationMap.entrySet() ) {
            if ( ! p.getValue().isEmpty() ) {
                attributes.put(p.getKey(), combineAnnotationValues(p.getValue()));
            }
        }

        if ( depth > 0 ) {
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));
        }

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        // note that in order to calculate the end position, we need a list of alleles that doesn't include anything symbolic
        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID).alleles(allelesList)
                .chr(loc.getContig()).start(loc.getStart()).computeEndFromAlleles(nonSymbolicAlleles(allelesList), loc.getStart(), loc.getStart())
                .genotypes(genotypes).unfiltered().attributes(new TreeMap<>(attributes)).log10PError(CommonInfo.NO_LOG10_PERROR);  // we will need to re-genotype later

        return builder.make();
    }

    /**
     * parse the annotations that were not identified as reducible annotations and combined by the annotation engine
     * @param annotationMap the map of info field annotation names and the list of their data from the merged VCs
     * @return  info field data parsed as ints or doubles
     */
    private static Map<String, List<Comparable>> parseRemainingAnnotations(final Map<String, List<ReducibleAnnotationData>> annotationMap) {
        final Map<String, List<Comparable>> parsedAnnotations = new HashMap<>();
        for (Map.Entry<String, List<ReducibleAnnotationData>> currentData : annotationMap.entrySet()) {
            List<Comparable> annotationValues = new ArrayList<>();
            for (ReducibleAnnotationData value : currentData.getValue()) {
                try {
                    final String stringValue = value.getRawData();
                    if (stringValue.contains(".")) {
                        annotationValues.add(Double.parseDouble(stringValue));
                    } else if (Character.isDigit(stringValue.charAt(0))){
                        annotationValues.add(Integer.parseInt(stringValue));
                    //TODO: uncomment this to parse dbSNP membership annotation once allele-specific merging for that attribute is added
                    /*} else if (Character.isLetter(stringValue.charAt(0))) {
                        if (stringValue.equalsIgnoreCase("true"))
                            annotationValues.add(true);
                        else if (stringValue.equalsIgnoreCase("false"))
                            annotationValues.add(false);*/
                    }

                } catch (final NumberFormatException e) {
                    logger.warn("WARNING: remaining (non-reducible) annotations are assumed to be ints or doubles or booleans, but " + value.getRawData() + " doesn't parse and will not be annotated in the final VC.");
                }
            }
            parsedAnnotations.put(currentData.getKey(),annotationValues);
        }
        return parsedAnnotations;
    }

    /**
     * @param list  the original alleles list
     * @return a non-null list of non-symbolic alleles
     */
    private static List<Allele> nonSymbolicAlleles(final List<Allele> list) {
        final List<Allele> result = new ArrayList<>(list.size());
        for ( final Allele allele : list ) {
            if ( !allele.isSymbolic() ) {
                result.add(allele);
            }
        }
        return result;
    }

    /**
     * Determines the ref allele given the provided reference base at this position
     *
     * @param VCs     collection of unsorted genomic VCs
     * @param loc     the current location
     * @param refBase the reference allele to use if all contexts in the VC are spanning
     * @return new Allele or null if no reference allele/base is available
     */
    private static Allele determineReferenceAlleleGivenReferenceBase(final List<VariantContext> VCs, final GenomeLoc loc, final Byte refBase) {
        final Allele refAllele = GATKVariantContextUtils.determineReferenceAllele(VCs, loc);
        if ( refAllele == null ) {
            return (refBase == null ? null : Allele.create(refBase, true));
        }
        return refAllele;
    }

    /**
     * Remove the stale attributes from the merged set
     *
     * @param attributes the attribute map
     */
    private static void removeStaleAttributesAfterMerge(final Map<String, List<ReducibleAnnotationData>> attributes) {
        attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
        attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
        attributes.remove(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        attributes.remove(VCFConstants.END_KEY);
    }

    /**
     * Adds attributes to the global map from the new context in a sophisticated manner
     *
     * @param pair                      VariantContext/Allele list pair from which to get attributes
     * @param annotationMap              map of annotations for combining later
     */
    private static void addReferenceConfidenceAttributes(Pair<VariantContext,List<Allele>> pair,
                                                         final Map<String, List<ReducibleAnnotationData>> annotationMap) {
        final Map<String, Object> myAttributes = pair.getFirst().getAttributes(); //these are the info field attributes of the VC in pair
        final List<Allele> sampleAlleles = pair.getSecond();

        for ( final Map.Entry<String, Object> p : myAttributes.entrySet() ) {
            final String key = p.getKey();
            //allele-specific attributes will always be in list form because they've already been parsed per-allele
            //non-allele-specific attributes (DP, etc.) will be a list of length 1
            final List<Object> valueList = pair.getFirst().getAttributeAsList(key);

            // add the existing annotation values to a list for combining later
            List<ReducibleAnnotationData> rawValuesList = annotationMap.get(key);
            if( rawValuesList == null ) {
                rawValuesList = new ArrayList<>();
                annotationMap.put(key, rawValuesList);
            }
            String combinedString = "";
            for(int i=0; i < valueList.size(); i++) {
                if (i > 0)
                    combinedString += ",";
                combinedString += valueList.get(i);
            }
            ReducibleAnnotationData pairData = new AlleleSpecificAnnotationData(sampleAlleles, combinedString);
            rawValuesList.add(pairData);
            annotationMap.put(key, rawValuesList);
        }
    }

    /**
     * This method does a couple of things:
     * <ul><li>
     *     remaps the vc alleles considering the differences between the final reference allele and its own reference,</li>
     * <li>
     *     collects alternative alleles present in variant context and add them to the {@code finalAlleles} set.
     * </li></ul>
     *
     * @param vc           the variant context.
     * @param refAllele    final reference allele.
     * @param finalAlleles where to add the final set of non-ref called alleles.
     * @return never {@code null}
     */
    //TODO as part of a larger refactoring effort {@link #remapAlleles} can be merged with {@link GATKVariantContextUtils#remapAlleles}.
    private static List<Allele> remapAlleles(final VariantContext vc, final Allele refAllele, final LinkedHashSet<Allele> finalAlleles) {

        final Allele vcRef = vc.getReference();
        final byte[] refBases = refAllele.getBases();
        final int extraBaseCount = refBases.length - vcRef.getBases().length;
        if (extraBaseCount < 0) throw new IllegalStateException("the wrong reference was selected");

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());
        result.add(refAllele);

        for (final Allele a : vc.getAlternateAlleles()) {
            if (a.isSymbolic()) {
                result.add(a);
                // we always skip <NON_REF> when adding to finalAlleles; this is done outside if it applies.
                // we also skip <*:DEL> if there isn't a real alternate allele.
                if ( !a.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) && !vc.isSymbolic() )
                    finalAlleles.add(a);
            } else if ( a == Allele.SPAN_DEL ) {
                result.add(a);
                // we skip * if there isn't a real alternate allele.
                if ( !vc.isBiallelic() )
                    finalAlleles.add(a);
            } else if (a.isCalled()) {
                final Allele newAllele;
                if (extraBaseCount > 0) {
                    final byte[] oldBases = a.getBases();
                    final byte[] newBases = Arrays.copyOf(oldBases,oldBases.length + extraBaseCount);
                    System.arraycopy(refBases,refBases.length - extraBaseCount,newBases,oldBases.length,extraBaseCount);
                    newAllele = Allele.create(newBases,false);
                } else
                    newAllele = a;
                result.add(newAllele);
                finalAlleles.add(newAllele);
            } else { // NO_CALL and strange miscellanea
                result.add(a);
            }
        }
        return result;
    }

    /**
     * Replaces any alleles in the VariantContext with NO CALLS or the symbolic deletion allele as appropriate, except for the generic ALT allele
     *
     * @param vc   VariantContext with the alleles to replace
     * @return non-null list of alleles
     */
    private static List<Allele> replaceWithNoCallsAndDels(final VariantContext vc) {
        if ( vc == null ) throw new IllegalArgumentException("VariantContext cannot be null");

        final List<Allele> result = new ArrayList<>(vc.getNAlleles());

        // no-call the reference allele
        result.add(Allele.NO_CALL);

        // handle the alternate alleles
        for ( final Allele allele : vc.getAlternateAlleles() ) {
            final Allele replacement;
            if ( allele.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) )
                replacement = allele;
            else if ( allele.length() < vc.getReference().length() )
                replacement = Allele.SPAN_DEL;
            else
                replacement = Allele.NO_CALL;

            result.add(replacement);
        }
        return result;
    }

    /**
     * Merge into the context a new genotype represented by the given VariantContext for the provided list of target alleles.
     * This method assumes that none of the alleles in the VC overlaps with any of the alleles in the set.
     *
     * @param mergedGenotypes       the genotypes context to add to
     * @param vc                    the Variant Context for the sample
     * @param remappedAlleles       the list of remapped alleles for the sample
     * @param targetAlleles         the list of target alleles
     * @param samplesAreUniquified  true if sample names have been uniquified
     */
    private static void mergeRefConfidenceGenotypes(final GenotypesContext mergedGenotypes,
                                                    final VariantContext vc,
                                                    final List<Allele> remappedAlleles,
                                                    final List<Allele> targetAlleles,
                                                    final boolean samplesAreUniquified) {
        final int maximumPloidy = vc.getMaxPloidy(GATKVariantContextUtils.DEFAULT_PLOIDY);
        // the map is different depending on the ploidy, so in order to keep this method flexible (mixed ploidies)
        // we need to get a map done (lazily inside the loop) for each ploidy, up to the maximum possible.
        final int[][] genotypeIndexMapsByPloidy = new int[maximumPloidy + 1][];
        final int maximumAlleleCount = Math.max(remappedAlleles.size(),targetAlleles.size());
        int[] perSampleIndexesOfRelevantAlleles;

        for (final Genotype g : vc.getGenotypes()) {
            final String name;
            if (samplesAreUniquified)
                name = g.getSampleName() + "." + vc.getSource();
            else
                name = g.getSampleName();
            final int ploidy = g.getPloidy();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy()));
            genotypeBuilder.name(name);
            final boolean hasPL = g.hasPL();
            final boolean hasSAC = g.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY);
            if (hasPL || hasSAC) {
                perSampleIndexesOfRelevantAlleles = getIndexesOfRelevantAlleles(remappedAlleles, targetAlleles, vc.getStart(), g);
                if (g.hasPL()) {
                    // lazy initialization of the genotype index map by ploidy.

                    final int[] genotypeIndexMapByPloidy = genotypeIndexMapsByPloidy[ploidy] == null
                            ? GenotypeLikelihoodCalculators.getInstance(ploidy, maximumAlleleCount).genotypeIndexMap(perSampleIndexesOfRelevantAlleles)
                            : genotypeIndexMapsByPloidy[ploidy];
                    final int[] PLs = generatePL(g, genotypeIndexMapByPloidy);
                    final int[] AD = g.hasAD() ? generateAD(g.getAD(), perSampleIndexesOfRelevantAlleles) : null;
                    genotypeBuilder.PL(PLs).AD(AD);
                }
                if (g.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY)) {
                    final List<Integer> sacIndexesToUse = adaptToSACIndexes(perSampleIndexesOfRelevantAlleles);
                    final int[] SACs = GATKVariantContextUtils.makeNewSACs(g, sacIndexesToUse);
                    genotypeBuilder.attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, SACs);
                }
            }
            mergedGenotypes.add(genotypeBuilder.make());
        }
    }

    /**
      * Adapt the relevant alleles to the SAC indexes
      *
      * @param perSampleIndexesOfRelevantAlleles
      * @return SAC indexes
     */
    private static List<Integer> adaptToSACIndexes(final int[] perSampleIndexesOfRelevantAlleles) {
        if (perSampleIndexesOfRelevantAlleles == null)
            throw new IllegalArgumentException("The per sample index of relevant alleles must not be null");

        final List<Integer> sacIndexesToUse = new ArrayList(2 * perSampleIndexesOfRelevantAlleles.length);

        for (int item : perSampleIndexesOfRelevantAlleles) {
            sacIndexesToUse.add(new Integer(2 * item));
            sacIndexesToUse.add(new Integer(2 * item + 1));
        }

        return sacIndexesToUse;
    }


    /**
     * Composes a new likelihood array given the original genotype and the genotype index map.
     *
     * @param g the original genotype.
     * @param genotypeIndexMapByPloidy genotype index map. The ith element indicates what genotype in {@code g} corresponds
     *                                 to the ith genotype in the return likelihoods array.
     *
     * @throws NullPointerException if {@code g} or {@code genotypeIndexMapByPloidy} is {@code null}, or if {@code g}
     *    does not contain likelihoods.
     * @throws IndexOutOfBoundsException if {@code genotypeIndexMapByPloidy} contain non valid
     *  genotype indices given the likelihood array in {@code g}.
     *
     * @return never {@code null} but an array of exactly {@code genotypeIndexMapByPloidy.length} positions.
     */
    private static int[] generatePL(final Genotype g, final int[] genotypeIndexMapByPloidy) {
        final int[] PLs = new int[genotypeIndexMapByPloidy.length];
        final int[] oldPLs = g.getPL();
        for (int i = 0; i < PLs.length; i++)
            PLs[i] = oldPLs[genotypeIndexMapByPloidy[i]];
        return PLs;
    }

    /**
     * Determines the allele mapping from myAlleles to the targetAlleles, substituting the generic "<ALT>" as appropriate.
     * If the remappedAlleles set does not contain "<ALT>" as an allele, it throws an exception.
     *
     * @param remappedAlleles   the list of alleles to evaluate
     * @param targetAlleles     the target list of alleles
     * @param position          position to output error info
     * @param g                 genotype from which targetAlleles are derived
     * @return non-null array of ints representing indexes
     */
    protected static int[] getIndexesOfRelevantAlleles(final List<Allele> remappedAlleles, final List<Allele> targetAlleles, final int position, final Genotype g) {

        if ( remappedAlleles == null || remappedAlleles.isEmpty()) throw new IllegalArgumentException("The list of input alleles must not be null or empty");
        if ( targetAlleles == null || targetAlleles.isEmpty() ) throw new IllegalArgumentException("The list of target alleles must not be null or empty");

        if ( !remappedAlleles.contains(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) )
            throw new UserException("The list of input alleles must contain " + GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE + " as an allele but that is not the case at position " + position + "; please use the Haplotype Caller with gVCF output to generate appropriate records");

        final int indexOfNonRef = remappedAlleles.indexOf(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        final int[] indexMapping = new int[targetAlleles.size()];

        // the reference likelihoods should always map to each other (even if the alleles don't)
        indexMapping[0] = 0;

        // create the index mapping, using the <NON-REF> allele whenever such a mapping doesn't exist
        for ( int i = 1; i < targetAlleles.size(); i++ ) {
            final Allele targetAllele = targetAlleles.get(i);

            // if there’s more than 1 DEL allele then we need to use the best one
            if ( targetAllele == Allele.SPAN_DEL && g.hasPL() ) {
                final int occurrences = Collections.frequency(remappedAlleles, Allele.SPAN_DEL);
                if ( occurrences > 1 ) {
                    final int indexOfBestDel = indexOfBestDel(remappedAlleles, g.getPL(), g.getPloidy());
                    indexMapping[i] = ( indexOfBestDel == -1 ? indexOfNonRef : indexOfBestDel );
                    continue;
                }
            }

            final int indexOfRemappedAllele = remappedAlleles.indexOf(targetAllele);
            indexMapping[i] = indexOfRemappedAllele == -1 ? indexOfNonRef : indexOfRemappedAllele;
        }

        return indexMapping;
    }

    /**
     * Returns the index of the best spanning deletion allele based on AD counts
     *
     * @param alleles   the list of alleles
     * @param PLs       the list of corresponding PL values
     * @param ploidy    the ploidy of the sample
     * @return the best index or -1 if not found
     */
    private static int indexOfBestDel(final List<Allele> alleles, final int[] PLs, final int ploidy) {
        int bestIndex = -1;
        int bestPL = Integer.MAX_VALUE;

        for ( int i = 0; i < alleles.size(); i++ ) {
            if ( alleles.get(i) == Allele.SPAN_DEL ) {
                final int homAltIndex = findHomIndex(i, ploidy, alleles.size());
                final int PL = PLs[homAltIndex];
                if ( PL < bestPL ) {
                    bestIndex = i;
                    bestPL = PL;
                }
            }
        }

        return bestIndex;
    }

    /**
     * Returns the index of the PL that represents the homozygous genotype of the given i'th allele
     *
     * @param i           the index of the allele with the list of alleles
     * @param ploidy      the ploidy of the sample
     * @param numAlleles  the total number of alleles
     * @return the hom index
     */
    private static int findHomIndex(final int i, final int ploidy, final int numAlleles) {
        // some quick optimizations for the common case
        if ( ploidy == 2 )
            return i + (i * (i + 1) / 2);   // this is straight from the VCF spec on PLs
        if ( ploidy == 1 )
            return i;

        final GenotypeLikelihoodCalculator calculator = GenotypeLikelihoodCalculators.getInstance(ploidy, numAlleles);
        final int[] alleleIndexes = new int[ploidy];
        Arrays.fill(alleleIndexes, i);
        return calculator.allelesToIndex(alleleIndexes);
    }

    /**
     * Generates a new AD array by adding zeros for missing alleles given the set of indexes of the Genotype's current
     * alleles from the original AD.
     *
     * @param originalAD    the original AD to extend
     * @param indexesOfRelevantAlleles the indexes of the original alleles corresponding to the new alleles
     * @return non-null array of new AD values
     */
    protected static int[] generateAD(final int[] originalAD, final int[] indexesOfRelevantAlleles) {
        if ( originalAD == null || indexesOfRelevantAlleles == null ) throw new IllegalArgumentException("The list of input AD values and alleles must not be null");

        final int numADs = indexesOfRelevantAlleles.length;
        final int[] newAD = new int[numADs];

        for ( int i = 0; i < numADs; i++ ) {
            final int oldIndex = indexesOfRelevantAlleles[i];
            if ( oldIndex >= originalAD.length )
                newAD[i] = 0;
            else
                newAD[i] = originalAD[oldIndex];
        }

        return newAD;
    }
}
