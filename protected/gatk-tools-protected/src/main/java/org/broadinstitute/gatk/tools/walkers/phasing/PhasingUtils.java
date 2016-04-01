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

package org.broadinstitute.gatk.tools.walkers.phasing;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.variantcontext.*;

import java.util.*;

/**
 * Utility class for phasing analysis
 */
class PhasingUtils {

    /**
     * Merge variants into a multi-nucleotide polymorphism (MNP)
     *
     * @param genomeLocParser parse the genome locations
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @param referenceFile sequence file containing the reference genome
     * @param alleleMergeRule rule for merging variants
     * @return merged variant or null if the variants are NOT an SNP or MNP, on the same contig, variant location 1 is the same or after the  variant location 2,
     * their genotypes do NOT have the same number of chromosomes, haplotype, number of attributes as chromosomes, are both hetrozygous or do not abide by the merge rule
     */
    static VariantContext mergeIntoMNP(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile, AlleleMergeRule alleleMergeRule) {

        // Check if variants are an SNP or MNP, on the same contig, variant location 1 is not before variant location 2
        if (!mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2))
            return null;

        // Check if variant genotypes have the same number of chromosomes, haplotype, number of attributes as chromosomes, and either genotype is homozygous
        if (!allSamplesAreMergeable(vc1, vc2))
            return null;

        // Check if there's a "point" in merging the VCs (e.g., annotations could be changed)
        if (!alleleMergeRule.allelesShouldBeMerged(vc1, vc2))
            return null;

        return reallyMergeIntoMNP(vc1, vc2, referenceFile);
    }

    /**
     * Find the alleles with the same haplotype
     * assumes alleleSegregationIsKnown
     * TODO - should alleleSegregationIsKnown be called within this method?
     *
     * @param gt1 genotype 1
     * @param gt2 genotype 2
     * @return gt1 and gt2 alleles with the same haplotype
     */
    static SameHaplotypeAlleles matchHaplotypeAlleles(final Genotype gt1, final Genotype gt2) {

        final SameHaplotypeAlleles hapAlleles = new SameHaplotypeAlleles();

        // Get the alleles
        final int numAlleles = gt1.getPloidy();
        final Allele[] site1AllelesArray = gt1.getAlleles().toArray(new Allele[numAlleles]);
        final Allele[] site2AllelesArray = gt2.getAlleles().toArray(new Allele[numAlleles]);

        // locations of the same HP attribute in gt2 to gt2
        final int[] site1ToSite2Inds = new int[numAlleles];

        if (gt1.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY) && gt2.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY)) {
            final String[] hp1 = (String[]) gt1.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY);
            final String[] hp2 = (String[]) gt2.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY);

            // Map of HP attribute to it's array index
            final HashMap<String, Integer> hpNameToSite1Inds = new HashMap<String, Integer>();

            // Hp name to index
            for (int ind1 = 0; ind1 < hp1.length; ++ind1) {
                hpNameToSite1Inds.put(hp1[ind1], ind1);
            }

            // Find the index of the gt2 HP attribute in gt1 HP attribute array
            for (int ind2 = 0; ind2 < hp2.length; ++ind2) {
                final int ind1 = hpNameToSite1Inds.get(hp2[ind2]);

                // attributes are not in the same position in both genotypes
                if (ind2 != ind1)
                    hapAlleles.requiresSwap = true;

                site1ToSite2Inds[ind1] = ind2;
            }
        }
        else { // gt1.isHom() || gt2.isHom() ; so, we trivially merge the corresponding alleles
            for (int ind = 0; ind < site1ToSite2Inds.length; ++ind)
                site1ToSite2Inds[ind] = ind;
        }

        // Get the alleles for gt1 and gt2 with the same haplotype
        for (int ind1 = 0; ind1 < numAlleles; ++ind1) {
            final Allele all1 = site1AllelesArray[ind1];
            final int ind2 = site1ToSite2Inds[ind1];
            final Allele all2 = site2AllelesArray[ind2]; // this is OK, since alleleSegregationIsKnown(gt1, gt2)

            // add the 2 alleles
            hapAlleles.hapAlleles.add(new AlleleOneAndTwo(all1, all2));
        }

        return hapAlleles;
    }

    /**
     * Merge variants into a multi-nucleotide polymorphism (MNP)
     *
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @param referenceFile sequence file containing the reference genome
     * @return variant with the merged MNP
     */
    static VariantContext reallyMergeIntoMNP(VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile) {
        final int startInter = vc1.getEnd() + 1;
        final int endInter = vc2.getStart() - 1;
        byte[] intermediateBases = null;

        // get bases between vc1 and vc2 in the reference sequence
        if (startInter <= endInter) {
            intermediateBases = referenceFile.getSubsequenceAt(vc1.getChr(), startInter, endInter).getBases();
            StringUtil.toUpperCase(intermediateBases);
        }

        // merge the reference bases with vc1 and vc2
        final MergedAllelesData mergeData = new MergedAllelesData(intermediateBases, vc1, vc2); // ensures that the reference allele is added

        final GenotypesContext mergedGenotypes = GenotypesContext.create();
        for (final Genotype gt1 : vc1.getGenotypes()) {
            final Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            // Alleles with the same haplotype
            final SameHaplotypeAlleles hapAlleles = matchHaplotypeAlleles(gt1, gt2);

            boolean isPhased = gt1.isPhased() && gt2.isPhased();
            if (hapAlleles.requiresSwap) // at least one swap of allele order was necessary, so trio-based phasing order cannot be maintained
                isPhased = false;

            final List<Allele> mergedAllelesForSample = new LinkedList<Allele>();
            for (AlleleOneAndTwo all1all2 : hapAlleles.hapAlleles) {
                final Allele mergedAllele = mergeData.ensureMergedAllele(all1all2.all1, all1all2.all2);
                mergedAllelesForSample.add(mergedAllele);
            }

            final double mergedGQ = Math.min(gt1.getLog10PError(), gt2.getLog10PError());

            final Map<String, Object> mergedGtAttribs = new HashMap<String, Object>();

            // get the min read backed phasing quality
            double PQ = Double.MAX_VALUE;
            if (gt1.hasAnyAttribute(VCFConstants.PHASE_QUALITY_KEY)) {
                PQ = Math.min(PQ, (double) gt1.getAnyAttribute(VCFConstants.PHASE_QUALITY_KEY));
            }
            if (gt2.hasAnyAttribute(VCFConstants.PHASE_QUALITY_KEY)) {
                PQ = Math.min(PQ, (double) gt2.getAnyAttribute(VCFConstants.PHASE_QUALITY_KEY));
            }
            if (PQ != Double.MAX_VALUE)
                mergedGtAttribs.put(VCFConstants.PHASE_QUALITY_KEY, PQ);

            if (gt1.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY)) {
                mergedGtAttribs.put(GATKVCFConstants.RBP_HAPLOTYPE_KEY, gt1.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY));
            }
            else if (gt2.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY)) { // gt1 doesn't have, but merged (so gt1 is hom and can take gt2's haplotype names):
                mergedGtAttribs.put(GATKVCFConstants.RBP_HAPLOTYPE_KEY, gt2.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY));
            }

            // make the merged genotype
            final Genotype mergedGt = new GenotypeBuilder(gt1.getSampleName(), mergedAllelesForSample).log10PError(mergedGQ).attributes(mergedGtAttribs).phased(isPhased).make();
            mergedGenotypes.add(mergedGt);
        }

        // get the merged name
        final String mergedName = mergeVariantContextNames(vc1.getSource(), vc2.getSource());
        final double mergedLog10PError = Math.min(vc1.getLog10PError(), vc2.getLog10PError());
        final Set<String> mergedFilters = new HashSet<String>(); // Since vc1 and vc2 were unfiltered, the merged record remains unfiltered
        final Map<String, Object> mergedAttribs = mergeVariantContextAttributes(vc1, vc2);

        // get the merged ID
        final List<String> mergedIDs = new ArrayList<String>();
        if ( vc1.hasID() ) mergedIDs.add(vc1.getID());
        if ( vc2.hasID() ) mergedIDs.add(vc2.getID());
        final String mergedID = mergedIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(VCFConstants.ID_FIELD_SEPARATOR, mergedIDs);

        // make the merged variant context
        final VariantContextBuilder mergedBuilder = new VariantContextBuilder(mergedName, vc1.getChr(), vc1.getStart(), vc2.getEnd(), mergeData.getAllMergedAlleles()).id(mergedID).genotypes(mergedGenotypes).log10PError(mergedLog10PError).filters(mergedFilters).attributes(mergedAttribs);
        VariantContextUtils.calculateChromosomeCounts(mergedBuilder, true);
        return mergedBuilder.make();
    }

    /**
     * Merge variant context names
     *
     * @param name1 variant context 1 name
     * @param name2 variant context 2 name
     * @return merged variant names (name1_name2)
     */
    static String mergeVariantContextNames(String name1, String name2) {
        return name1 + "_" + name2;
    }

    /**
     * Get preset attributes and that are in vc1 or vc2
     * TODO: Will always return an empty map because MERGE_OR_ATTRIBS is empty
     *
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @return merged attributes in vc1 or vc2
     */
    static Map<String, Object> mergeVariantContextAttributes(VariantContext vc1, VariantContext vc2) {
        // Map of attribute name to value
        Map<String, Object> mergedAttribs = new HashMap<String, Object>();

        final List<VariantContext> vcList = new LinkedList<VariantContext>();
        vcList.add(vc1);
        vcList.add(vc2);

        // Attribute of interest
        //String[] MERGE_OR_ATTRIBS = {VCFConstants.DBSNP_KEY};
        final String[] MERGE_OR_ATTRIBS = {};
        for (String orAttrib : MERGE_OR_ATTRIBS) {
            boolean attribVal = false;
            for (VariantContext vc : vcList) {
                // Does the variant have the attribute?
                attribVal = vc.getAttributeAsBoolean(orAttrib, false);
                if ( attribVal ) // already true, so no reason to continue:
                    break;
            }
            mergedAttribs.put(orAttrib, attribVal);
        }

        return mergedAttribs;
    }

    /**
     * Check if variants can be merged into the multi-nucleotide polymorphism (MNP)
     *
     * @param genomeLocParser parse the genome locations
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @return true if variants are an SNP or MNP, on the same contig, variant location 1 is not before variant location 2, unfiltered, from the same sample set and are called,
     * false otherwise
     */
    static boolean mergeIntoMNPvalidationCheck(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2) {
	// Can only merge "simple" base strings (i.e., SNPs or MNPs, but not indels):
	final boolean vc1CanBeMerged = (vc1.isSNP() || vc1.isMNP()) && !vc1.hasAllele(Allele.SPAN_DEL);
	final boolean vc2CanBeMerged = (vc2.isSNP() || vc2.isMNP()) && !vc2.hasAllele(Allele.SPAN_DEL);
	if (!vc1CanBeMerged || !vc2CanBeMerged)
            return false;

        final GenomeLoc loc1 = GATKVariantContextUtils.getLocation(genomeLocParser, vc1);
        final GenomeLoc loc2 = GATKVariantContextUtils.getLocation(genomeLocParser, vc2);

        // Must be on same contig
        if (!loc1.onSameContig(loc2))
            return false;

        // Variant 1 location must not be before variant context 2
        if (!loc1.isBefore(loc2))
            return false;

        // Variants can not be filtered
        if (vc1.isFiltered() || vc2.isFiltered())
            return false;

        // Variants must come from the same sample set
        if (!vc1.getSampleNames().equals(vc2.getSampleNames()))
            return false;

        // All of the variant genotypes must be unfiltered and called
        if (!allGenotypesAreUnfilteredAndCalled(vc1) || !allGenotypesAreUnfilteredAndCalled(vc2))
            return false;

        return true;
    }

    static boolean allGenotypesAreUnfilteredAndCalled(VariantContext vc) {
        for (final Genotype gt : vc.getGenotypes()) {
            if (gt.isNoCall() || gt.isFiltered())
                return false;
        }

        return true;
    }

    /**
     * Check if can merge genotypes from the same sample
     *
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @return true if variants are phased or either is a homozygous, false otherwise
     */
    static boolean allSamplesAreMergeable(VariantContext vc1, VariantContext vc2) {
        // Check that each sample's genotype in vc2 is uniquely appendable onto its genotype in vc1
        for (final Genotype gt1 : vc1.getGenotypes()) {
            final Genotype gt2 = vc2.getGenotype(gt1.getSampleName());
            if ( gt2 == null ) // gt2 does not have sample name
                return false;

            if (!alleleSegregationIsKnown(gt1, gt2)) // can merge if: phased, or if either is a hom
                return false;
        }

        return true;
    }

    /**
     * Check if the allele segregation is known
     *
     * @param gt1 genotype 1
     * @param gt2 genotype 2
     * @return true if genotypes have the same number of chromosomes, haplotype, number of attributes
     * as chromosomes, and either genotype is homozygous, false otherwise
     */
    static boolean alleleSegregationIsKnown(Genotype gt1, Genotype gt2) {
        // If gt1 or gt2 do not have the same number of chromosomes, then can not be merged.
        if (gt1.getPloidy() != gt2.getPloidy())
            return false;

        // If gt1 or gt2 are homozygous, then could be merged.
        if (gt1.isHom() || gt2.isHom())
            return true;

        // If gt1 or gt2 do not have a read backed phasing haplotype, then can not be merged
        if (!gt1.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY) || !gt2.hasAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY))
            return false;

        // If gt1 or gt2 do not same number of HP attributes as chromosomes, then can not be merged.
        final String[] hp1 = (String[]) gt1.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY);
        final String[] hp2 = (String[]) gt2.getAnyAttribute(GATKVCFConstants.RBP_HAPLOTYPE_KEY);
        if (hp1.length != gt1.getPloidy() || hp2.length != gt2.getPloidy())
            return false;

        // gt1 and gt2 must have the same read backed phasing haplotype identifier attributes to be merged
	    final String[] hp1Copy = Arrays.copyOf(hp1, hp1.length);
	    final String[] hp2Copy = Arrays.copyOf(hp2, hp2.length);
        Arrays.sort(hp1Copy);
        Arrays.sort(hp2Copy);
        return (Arrays.equals(hp1Copy, hp2Copy)); // The haplotype names match (though possibly in a different order)
    }

    /**
     * Check if some samples have double alternate alleles
     *
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @return true if there is a sample with double alternate alleles, false otherwise
     */
    static boolean someSampleHasDoubleNonReferenceAllele(VariantContext vc1, VariantContext vc2) {
        for (final Genotype gt1 : vc1.getGenotypes()) {
            // gt2 from the same sample as gt1
            final Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            if ( gt2 != null ) {
                // Find the alleles with the same haplotype
                final SameHaplotypeAlleles hapAlleles = matchHaplotypeAlleles(gt1, gt2);

                // Find corresponding alternate alleles
                for (AlleleOneAndTwo all1all2 : hapAlleles.hapAlleles) {
                    if (all1all2.all1.isNonReference() && all1all2.all2.isNonReference())
                        return true;
                }
            }
        }

        return false;
    }

    /**
     * Check that alleles at vc1 and at vc2 always segregate together in all samples (including reference)
     *
     * @param vc1 variant context 1
     * @param vc2 variant context 2
     * @return true if alleles segregate together, false otherwise
     */
    static boolean doubleAllelesSegregatePerfectlyAmongSamples(VariantContext vc1, VariantContext vc2) {
        final Map<Allele, Allele> allele1ToAllele2 = new HashMap<Allele, Allele>();
        final Map<Allele, Allele> allele2ToAllele1 = new HashMap<Allele, Allele>();

        // Note the segregation of the alleles for the reference genome:
        allele1ToAllele2.put(vc1.getReference(), vc2.getReference());
        allele2ToAllele1.put(vc2.getReference(), vc1.getReference());

        // Note the segregation of the alleles for each sample (and check that it is consistent with the reference and all previous samples).
        for (final Genotype gt1 : vc1.getGenotypes()) {
            final Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            SameHaplotypeAlleles hapAlleles = matchHaplotypeAlleles(gt1, gt2);
            for (AlleleOneAndTwo all1all2 : hapAlleles.hapAlleles) {
                final Allele all1 = all1all2.all1;
                final Allele all2 = all1all2.all2;

                final Allele all1To2 = allele1ToAllele2.get(all1);
                if (all1To2 == null)
                    allele1ToAllele2.put(all1, all2);
                else if (!all1To2.equals(all2)) // all1 segregates with two different alleles at site 2
                    return false;

                final Allele all2To1 = allele2ToAllele1.get(all2);
                if (all2To1 == null)
                    allele2ToAllele1.put(all2, all1);
                else if (!all2To1.equals(all1)) // all2 segregates with two different alleles at site 1
                    return false;
            }
        }

        return true;
    }

    /**
     * Class for variants merging rules
     */
    abstract static class AlleleMergeRule {
        // vc1, vc2 are ONLY passed to allelesShouldBeMerged() if mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2) AND allSamplesAreMergeable(vc1, vc2):
        abstract public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2);

        public String toString() {
            return "all samples are mergeable";
        }
    }

    /**
     * Class for storing the alleles with the same haplotype
     */
    static class SameHaplotypeAlleles {

        /// Alleles are not in the same order
        public boolean requiresSwap;

        /// Lisgt of gthe 2 alleles with the same haplotype
        public List<AlleleOneAndTwo> hapAlleles;

        public SameHaplotypeAlleles() {
            requiresSwap = false;
            hapAlleles = new LinkedList<AlleleOneAndTwo>();
        }
    }

    /**
     * Class for holding 2 alleles
     */
    static class AlleleOneAndTwo {
        /// allele 1
        private Allele all1;
        /// allele2
        private Allele all2;

        /**
         * Constructor
         *
         * @param all1 allele 1
         * @param all2 allele 2
         */
        public AlleleOneAndTwo(Allele all1, Allele all2) {
            this.all1 = all1;
            this.all2 = all2;
        }

        /**
         * Get the hah code for alleles 1 and 2
         *
         * @return hash code for alleles 1 and 2
         */
        public int hashCode() {
            return all1.hashCode() + all2.hashCode();
        }

        /**
         * Check if equal to another 2 alleles
         *
         * @param other allele to compare to
         * @return true if equal, false otherwise
         */
        public boolean equals(Object other) {
            if (!(other instanceof AlleleOneAndTwo))
                return false;

            final AlleleOneAndTwo otherAot = (AlleleOneAndTwo) other;
            return (this.all1.equals(otherAot.all1) && this.all2.equals(otherAot.all2));
        }
    }

    /**
     * Class for merging alleles
     */
    static class MergedAllelesData {
        /// merged alleles
        private Map<AlleleOneAndTwo, Allele> mergedAlleles;

        /// bases between the alleles
        private byte[] intermediateBases;

        /// number of bases between the alleles
        private int intermediateLength;

        /**
         * Constructor
         *
         * @param intermediateBases array of bases
         * @param vc1 variant context 1
         * @param vc2 variant context 2
         */
        public MergedAllelesData(byte[] intermediateBases, VariantContext vc1, VariantContext vc2) {
            this.mergedAlleles = new HashMap<AlleleOneAndTwo, Allele>(); // implemented equals() and hashCode() for AlleleOneAndTwo
            this.intermediateBases = intermediateBases;
            this.intermediateLength = this.intermediateBases != null ? this.intermediateBases.length : 0;

            // merge the reference bases from vc1 before and vc2 after the reference (intermediate) bases
            this.ensureMergedAllele(vc1.getReference(), vc2.getReference(), true);
        }

        /**
         * Ensure that the alleles are merged. The merged allele is alternate.
         *
         * @param all1 allele 1
         * @param all2 allele 2
         * @return merged allele
         */
        public Allele ensureMergedAllele(Allele all1, Allele all2) {
            return ensureMergedAllele(all1, all2, false); // false <-> since even if all1+all2 = reference, it was already created in the constructor
        }

        /**
         * Ensure that the alleles are merged.
         * all1 is before all2, if there is a gap between them, join with the intermediate bases
         *
         * @param all1 allele 1
         * @param all2 allele 2
         * @param isRef if true, merged allele is reference, if false, merged allele is alternate
         * @return merged allele
         */
        private Allele ensureMergedAllele(Allele all1, Allele all2, boolean isRef) {
            AlleleOneAndTwo all12 = new AlleleOneAndTwo(all1, all2);
            Allele mergedAllele = mergedAlleles.get(all12);

             if (mergedAllele == null) {
                final byte[] bases1 = all1.getBases();
                final byte[] bases2 = all2.getBases();

                final byte[] mergedBases = new byte[bases1.length + intermediateLength + bases2.length];
                System.arraycopy(bases1, 0, mergedBases, 0, bases1.length);
                if (intermediateBases != null)
                    System.arraycopy(intermediateBases, 0, mergedBases, bases1.length, intermediateLength);
                System.arraycopy(bases2, 0, mergedBases, bases1.length + intermediateLength, bases2.length);

                mergedAllele = Allele.create(mergedBases, isRef);
                mergedAlleles.put(all12, mergedAllele);
            }

            return mergedAllele;
        }

        /**
         * Get all merged alleles
         *
         * @return set of merged alleles values
         */
        public Set<Allele> getAllMergedAlleles() {
            return new HashSet<Allele>(mergedAlleles.values());
        }
    }
}
