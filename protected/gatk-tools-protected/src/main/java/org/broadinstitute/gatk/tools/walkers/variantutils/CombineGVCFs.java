/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
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

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;

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
public class CombineGVCFs extends RodWalker<CombineGVCFs.PositionalState, CombineGVCFs.OverallState> {

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

    protected final class OverallState {
        final LinkedList<VariantContext> VCs = new LinkedList<>();
        GenomeLoc prevPos = null;
        byte refAfterPrevPos;

        public OverallState() {}
    }

    /**
     * The gVCF files to merge together
     */
    @Input(fullName="variant", shortName = "V", doc="One or more input gVCF files", required=true)
    public List<RodBindingCollection<VariantContext>> variantCollections;
    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    @Output(doc="File to which the combined gVCF should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(fullName="convertToBasePairResolution", shortName="bpResolution", doc = "If specified, convert banded gVCFs to all-sites gVCFs", required=false)
    protected boolean USE_BP_RESOLUTION = false;

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

    public OverallState reduceInit() {
        return new OverallState();
    }

    public OverallState reduce(final PositionalState startingStates, final OverallState previousState) {
        if ( startingStates == null )
            return previousState;

        final int currentPos = startingStates.loc.getStart();

        if ( !startingStates.VCs.isEmpty() ) {
            if ( ! okayToSkipThisSite(currentPos, previousState.prevPos) )
                endPreviousStates(previousState, currentPos - 1, startingStates.refBases[0]);
            previousState.VCs.addAll(startingStates.VCs);
        }

        if ( USE_BP_RESOLUTION || containsEndingContext(previousState.VCs, currentPos) ) {
            endPreviousStates(previousState, currentPos, startingStates.refBases.length > 1 ? startingStates.refBases[1] : (byte)'N');
        }

        return previousState;
    }

    /**
     * Is it okay to skip the given position?
     *
     * @param thisPos      this position
     * @param lastPosRun   the last position for which we created a VariantContext
     * @return true if it is okay to skip this position, false otherwise
     */
    private boolean okayToSkipThisSite(final int thisPos, final GenomeLoc lastPosRun) {
        return lastPosRun != null && thisPos == lastPosRun.getStart() + 1;
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
            if ( isEndingContext(vc, pos) )
                return true;
        }
        return false;
    }

    /**
     * Does the given variant context end (in terms of reference blocks, not necessarily formally) at the given position.
     * Note that for the purposes of this method/tool, deletions are considered to be single base events (as opposed to
     * reference blocks), hence the check for the number of alleles (because we know there will always be a <NON_REF> allele).
     *
     * @param vc   the variant context
     * @param pos  the position to query against
     * @return true if this variant context "ends" at this position, false otherwise
     */
    private boolean isEndingContext(final VariantContext vc, final int pos) {
        return vc.getNAlleles() > 2 || vc.getEnd() == pos;
    }

    /**
     * Disrupt the VariantContexts so that they all stop at the given pos, write them out, and put the remainder back in the list.
     *
     * @param state   the state with list of VariantContexts
     * @param pos     the target ending position
     * @param refBase the reference base to use at the position AFTER pos
     */
    private void endPreviousStates(final OverallState state, final int pos, final byte refBase) {

        final List<VariantContext> stoppedVCs = new ArrayList<>(state.VCs.size());

        for ( int i = state.VCs.size() - 1; i >= 0; i-- ) {
            final VariantContext vc = state.VCs.get(i);
            if ( vc.getStart() <= pos ) {

                stoppedVCs.add(vc);

                // if it was ending anyways, then remove it from the future state
                if ( isEndingContext(vc, pos) )
                    state.VCs.remove(i);
            }
        }

        if ( !stoppedVCs.isEmpty() ) {
            final GenomeLoc gLoc = genomeLocParser.createGenomeLoc(stoppedVCs.get(0).getChr(), pos);

            // we need the specialized merge if the site contains anything other than ref blocks
            final VariantContext mergedVC;
            if ( containsTrueAltAllele(stoppedVCs) )
                mergedVC = ReferenceConfidenceVariantContextMerger.merge(stoppedVCs, gLoc, refBase, false);
            else
                mergedVC = referenceBlockMerge(stoppedVCs, state, pos);

            vcfWriter.add(mergedVC);
            state.prevPos = gLoc;
            state.refAfterPrevPos = refBase;
        }
    }

    /**
     * Combine a list of reference block VariantContexts.
     * We can't use GATKVariantContextUtils.simpleMerge() because it is just too slow for this sort of thing.
     *
     * @param VCs   the variant contexts to merge
     * @param state the state object
     * @param end   the end of this block (inclusive)
     * @return a new merged VariantContext
     */
    private VariantContext referenceBlockMerge(final List<VariantContext> VCs, final OverallState state, final int end) {

        final VariantContext first = VCs.get(0);

        // ref allele and start
        final Allele refAllele;
        final int start;
        if ( state.prevPos == null || !state.prevPos.getContig().equals(first.getChr()) || first.getStart() >= state.prevPos.getStart() + 1) {
            start = first.getStart();
            refAllele = first.getReference();
        } else {
            start = state.prevPos.getStart() + 1;
            refAllele = Allele.create(state.refAfterPrevPos, true);
        }

        // attributes
        final Map<String, Object> attrs = new HashMap<>(1);
        if ( !USE_BP_RESOLUTION && end != start )
            attrs.put(VCFConstants.END_KEY, Integer.toString(end));

        // genotypes
        final GenotypesContext genotypes = GenotypesContext.create();
        for ( final VariantContext vc : VCs ) {
            for ( final Genotype g : vc.getGenotypes() )
                genotypes.add(new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(g.getPloidy())).make());
        }

        return new VariantContextBuilder("", first.getChr(), start, end, Arrays.asList(refAllele, GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE)).attributes(attrs).genotypes(genotypes).make();
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
    public void onTraversalDone(final OverallState state) {
        // there shouldn't be any state left unless the user cut in the middle of a gVCF block
        if ( !state.VCs.isEmpty() )
            logger.warn("You have asked for an interval that cuts in the middle of one or more gVCF blocks. Please note that this will cause you to lose records that don't end within your interval.");
    }
}
