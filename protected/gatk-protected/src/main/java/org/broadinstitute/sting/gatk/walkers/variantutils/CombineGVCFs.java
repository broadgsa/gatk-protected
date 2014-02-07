/*
*  By downloading the PROGRAM you agree to the following terms of use:
*
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
*
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.*;

import java.util.*;

/**
 * Combines any number of gVCF files that were produced by the Haplotype Caller into a single joint gVCF file.
 *
 * <p>
 * CombineGVCFs is meant to be used for hierarchical merging of gVCFs that will eventually be input into GenotypeGVCFs.
 * One would use this tool when needing to genotype too large a number of individual gVCFs; instead of passing them
 * all in to GenotypeGVCFs, one would first use CombineGVCFs on smaller batches of samples and then pass these combined
 * gVCFs to GenotypeGVCFs.
 *
 * Note that this tool cannot work with just any gVCF files - they must have been produced with the Haplotype Caller
 * as part of the "single sample discovery" pipeline using the '-ERC GVCF' mode, which uses a sophisticated reference
 * model to produce accurate genotype likelihoods for every position in the target.
 *
 * <h3>Input</h3>
 * <p>
 * One or more Haplotype Caller gVCFs to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CombineGVCFs \
 *   --variant gvcf1.vcf \
 *   --variant gvcf2.vcf \
 *   -o mergeGvcf.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=1))
public class CombineGVCFs extends RodWalker<CombineGVCFs.PositionalState, LinkedList<VariantContext>> {

    protected final class PositionalState {
        final List<VariantContext> VCs;
        final byte[] refBases;
        final GenomeLoc loc;
        public PositionalState(final List<VariantContext> VCs, final byte[] refBases, final GenomeLoc loc) {
            this.VCs = VCs;
            this.refBases = refBases;
            this.loc = loc;
        }
    }

    /**
     * The gVCF files to merge together
     */
    @Input(fullName="variant", shortName = "V", doc="One or more input gVCF files", required=true)
    public List<RodBindingCollection<VariantContext>> variantCollections;
    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    @Output(doc="File to which the combined gVCF should be written")
    protected VariantContextWriter vcfWriter = null;

    private GenomeLocParser genomeLocParser;

    public void initialize() {
        // take care of the VCF headers
        final Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit());
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);

        final Set<String> samples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        final VCFHeader vcfHeader = new VCFHeader(headerLines, samples);
        vcfWriter.writeHeader(vcfHeader);

        // collect the actual rod bindings into a list for use later
        for ( final RodBindingCollection<VariantContext> variantCollection : variantCollections )
            variants.addAll(variantCollection.getRodBindings());

        genomeLocParser = getToolkit().getGenomeLocParser();
    }

    public PositionalState map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        final GenomeLoc loc = ref.getLocus();
        return new PositionalState(tracker.getValues(variants, loc), ref.getBases(), loc);
    }

    public LinkedList<VariantContext> reduceInit() {
        return new LinkedList<>();
    }

    public LinkedList<VariantContext> reduce(final PositionalState startingStates, final LinkedList<VariantContext> previousState) {
        if ( startingStates == null )
            return previousState;

        final int currentPos = startingStates.loc.getStart();

        if ( !startingStates.VCs.isEmpty() ) {
            endPreviousStates(previousState, currentPos - 1, startingStates.refBases[0]);
            previousState.addAll(startingStates.VCs);
        }

        if ( containsEndingContext(previousState, currentPos) ) {
            endPreviousStates(previousState, currentPos, startingStates.refBases.length > 1 ? startingStates.refBases[1] : (byte)'N');
        }

        return previousState;
    }

    /**
     * Does the given list of VariantContexts contain any whose context ends at the given position?
     *
     * @param VCs  list of VariantContexts
     * @param pos  the position to check against
     * @return true if there are one or more VCs that end at pos, false otherwise
     */
    private boolean containsEndingContext(final List<VariantContext> VCs, final int pos) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getEnd() == pos )
                return true;
        }
        return false;
    }

    /**
     * Disrupt the VariantContexts so that they all stop at the given pos, write them out, and put the remainder back in the list.
     *
     * @param VCs  list of VariantContexts
     * @param pos  the target ending position
     * @param refBase the reference base to use at the position AFTER pos
     */
    private void endPreviousStates(final LinkedList<VariantContext> VCs, final int pos, final byte refBase) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        final List<VariantContext> stoppedVCs = new ArrayList<>(VCs.size());

        for ( int i = VCs.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = VCs.get(i);
            if ( vc.getStart() > pos )
                continue;

            // if it was ending anyways, then just remove it as is;
            // note that for the purposes of this method, deletions are considered to be single base events (as opposed
            // to ref blocks), hence the check for the number of alleles (because we know there will always be a <NON_REF> allele)
            if ( vc.getNAlleles() > 2 || vc.getEnd() == pos ) {
                stoppedVCs.add(vc);
                VCs.remove(i);
            }
            // otherwise we need to split it into two pieces
            else {
                // the first half
                final Map<String, Object> attrs = new HashMap<>(vc.getAttributes());
                if ( attrs.containsKey(VCFConstants.END_KEY) )
                    attrs.put(VCFConstants.END_KEY, Integer.toString(pos));
                stoppedVCs.add(new VariantContextBuilder(vc).stop(pos).attributes(attrs).make());

                // the second half
                final Allele refAllele = Allele.create(refBase, true);
                final List<Allele> alleles = new ArrayList<>();
                alleles.add(refAllele);
                alleles.addAll(vc.getAlternateAlleles());
                final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
                for ( final Genotype g : vc.getGenotypes() )
                    genotypes.add(new GenotypeBuilder(g).alleles(Arrays.asList(refAllele, refAllele)).make());
                VCs.set(i, new VariantContextBuilder(vc).start(pos + 1).alleles(alleles).genotypes(genotypes).make());
            }
        }

        if ( !stoppedVCs.isEmpty() ) {
            final VariantContext mergedVC = mergeVCs(stoppedVCs);
            vcfWriter.add(mergedVC);
        }
    }

    /**
     * Combine (and re-annotate) a list of VariantContexts
     *
     * @param VCs  the VariantContexts to merge
     * @return a new VariantContext
     */
    private VariantContext mergeVCs(final List<VariantContext> VCs) {
        // we need the specialized merge if the site contains anything other than ref blocks
        if ( containsTrueAltAllele(VCs) )
            return GATKVariantContextUtils.referenceConfidenceMerge(VCs, genomeLocParser.createGenomeLoc(VCs.get(0)), null, false);

        // otherwise we can drop down to the generic simple merge
        return GATKVariantContextUtils.simpleMerge(VCs, null, VCs.size(),
                    GATKVariantContextUtils.FilteredRecordMergeType.KEEP_UNCONDITIONAL,
                    GATKVariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);
    }

    /**
     * Does the given list of VariantContexts contain any with an alternate allele other than <NON_REF>?
     *
     * @param VCs  list of VariantContexts
     * @return true if there are one or more VCs that contain a true alternate allele, false otherwise
     */
    private boolean containsTrueAltAllele(final List<VariantContext> VCs) {
        if ( VCs == null ) throw new IllegalArgumentException("The list of VariantContexts cannot be null");

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 )
                return true;
        }
        return false;
    }

    @Override
    public void onTraversalDone(final LinkedList<VariantContext> state) {
        // there shouldn't be any state left unless the user cut in the middle of a gVCF block
        if ( !state.isEmpty() )
            logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
    }
}
