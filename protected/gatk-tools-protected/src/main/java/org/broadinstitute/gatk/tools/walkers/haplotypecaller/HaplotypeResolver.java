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

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;

import java.util.*;

/**
 * Haplotype-based resolution of variants in separate callsets.
 *
 * <p>
 * HaplotypeResolver is a tool that takes two VCF files and constructs haplotypes based on the variants inside them.
 * From that, it can resolve potential differences in variant calls that are inherently the same (or similar) variants.
 * Records are annotated with the set and status attributes.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * Two variant files to resolve.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A single consensus VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T HaplotypeResolver \
 *   -R reference.fasta \
 *   -V:v1 input1.vcf \
 *   -V:v2 input2.vcf \
 *   -o output.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-HaplotypeResolver.ACTIVE_WINDOW,stop= HaplotypeResolver.ACTIVE_WINDOW))
public class HaplotypeResolver extends RodWalker<Integer, Integer> {

    protected static final String INTERSECTION_SET = "intersection";
    protected static final String SAME_STATUS = "same";
    protected static final String SOME_ALLELES_MATCH_STATUS = "someAllelesMatch";
    protected static final String SAME_START_DIFFERENT_ALLELES_STATUS = "sameStartDifferentAlleles";
    protected static final String SAME_BY_HAPLOTYPE_STATUS = "sameByHaplotype";
    protected static final String ONE_ALLELE_SUBSET_OF_OTHER_STATUS = "OneAlleleSubsetOfOther";
    protected static final String OVERLAPPING_EVENTS_STATUS = "overlappingEvents";

    protected final static int MAX_DISTANCE_BETWEEN_MERGED_RECORDS = 50;
    protected final static int MAX_HAPLOTYPE_TO_CONSIDER = 1000;
    protected final static int MAX_VARIANT_SIZE_TO_CONSIDER = 100;
    protected final static int ACTIVE_WINDOW = MAX_HAPLOTYPE_TO_CONSIDER + MAX_VARIANT_SIZE_TO_CONSIDER;

    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBinding<VariantContext>> variants;

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter baseWriter = null;
    private VariantContextWriter writer;

    /**
     * Set to 'null' if you don't want the set field emitted.
     */
    @Argument(fullName="setKey", shortName="setKey", doc="Key used in the INFO key=value tag emitted describing which set the combined VCF record came from", required=false)
    protected String SET_KEY = "set";

    /**
     * Set to 'null' if you don't want the status field emitted.
     */
    @Argument(fullName="statusKey", shortName="statusKey", doc="Key used in the INFO key=value tag emitted describing the extent to which records match", required=false)
    protected String STATUS_KEY = "status";

    private final LinkedList<VCcontext> queue = new LinkedList<VCcontext>();
    private String source1, source2;
    private final List<VariantContext> sourceVCs1 = new ArrayList<VariantContext>();
    private final List<VariantContext> sourceVCs2 = new ArrayList<VariantContext>();


    private class VCcontext {
        public final Collection<VariantContext> vcs;
        public final GenomeLoc loc;
        public final ReferenceContext ref;

        public VCcontext(final Collection<VariantContext> vcs, final ReferenceContext ref) {
            this.vcs = vcs;
            this.loc = getToolkit().getGenomeLocParser().createGenomeLoc(vcs.iterator().next());
            this.ref = ref;
        }
    }

    public void initialize() {

        if ( variants.size() != 2 ) {
            throw new UserException.BadArgumentValue("variant", "this tool requires exactly 2 input variant files");
        }
        source1 = variants.get(0).getName();
        source2 = variants.get(1).getName();

        if ( SET_KEY.toLowerCase().equals("null") )
            SET_KEY = null;
        if ( STATUS_KEY.toLowerCase().equals("null") )
            STATUS_KEY = null;

        // for now, INFO and FORMAT fields are not propagated to the output VCF (so they aren't put into the header)
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        if ( SET_KEY != null )
            headerLines.add(new VCFInfoHeaderLine(SET_KEY, 1, VCFHeaderLineType.String, "Source VCF for the merged record"));
        if ( STATUS_KEY != null )
            headerLines.add(new VCFInfoHeaderLine(STATUS_KEY, 1, VCFHeaderLineType.String, "Extent to which records match"));
        final VCFHeader vcfHeader = new VCFHeader(headerLines, Collections.<String>emptySet());
        baseWriter.writeHeader(vcfHeader);
        writer = VariantContextWriterFactory.sortOnTheFly(baseWriter, ACTIVE_WINDOW);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        final Collection<VariantContext> VCs = tracker.getValues(variants, context.getLocation());
        if ( VCs.size() == 0 )
            return 0;

        final VCcontext vc = new VCcontext(VariantContextUtils.sitesOnlyVariantContexts(VCs), ref);

        // TODO -- what should we do about filtered records?

        if ( !queue.isEmpty() ) {

            final VCcontext previous = queue.getLast();
            if ( !previous.loc.onSameContig(vc.loc) ||
                    previous.loc.distance(vc.loc) > MAX_DISTANCE_BETWEEN_MERGED_RECORDS ||
                    queue.getFirst().loc.distance(vc.loc) > MAX_HAPLOTYPE_TO_CONSIDER ) {
                purgeQueue();
            }
        }

        queue.addLast(vc);
        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        if ( !queue.isEmpty() )
            purgeQueue();
        writer.close();
    }

    private void purgeQueue() {

        final ReferenceContext refContext = queue.getFirst().ref;

        // divide them up by source
        while ( !queue.isEmpty() ) {
            VCcontext context = queue.removeFirst();
            for ( final VariantContext vc: context.vcs ) {
                if ( vc.getSource().equals(source1) )
                    sourceVCs1.add(vc);
                else
                    sourceVCs2.add(vc);
            }
        }

        writeAndPurgeAllEqualVariants(sourceVCs1, sourceVCs2, SAME_STATUS);

        if ( sourceVCs1.isEmpty() ) {
            writeAll(sourceVCs2, source2, null);
        } else if ( sourceVCs2.isEmpty() ) {
            writeAll(sourceVCs1, source1, null);
        } else {
            resolveByHaplotype(refContext);
        }

        // allow for GC of the data
        sourceVCs1.clear();
        sourceVCs2.clear();
    }

    private void writeAll(final List<VariantContext> sourceVCs, final String set, final String status) {
        for ( final VariantContext vc : sourceVCs ) {
            writeOne(vc, set, status);
        }
    }

    private void writeOne(final VariantContext vc, final String set, final String status) {
        final Map<String, Object> attrs = new HashMap<>();
        if ( SET_KEY != null && set != null )
            attrs.put(SET_KEY, set);
        if ( STATUS_KEY != null && status != null )
            attrs.put(STATUS_KEY, status);
        writer.add(new VariantContextBuilder(vc).attributes(attrs).make());
    }

    private void writeAndPurgeAllEqualVariants(final List<VariantContext> sourceVCs1, final List<VariantContext> sourceVCs2, final String status) {

        int currentIndex1 = 0, currentIndex2 = 0;
        int size1 = sourceVCs1.size(), size2 = sourceVCs2.size();
        VariantContext current1 = (currentIndex1 < size1 ? sourceVCs1.get(currentIndex1): null);
        VariantContext current2 = (currentIndex2 < size2 ? sourceVCs2.get(currentIndex2): null);

        while ( current1 != null && current2 != null ) {

            final GenomeLoc loc1 = getToolkit().getGenomeLocParser().createGenomeLoc(current1);
            final GenomeLoc loc2 = getToolkit().getGenomeLocParser().createGenomeLoc(current2);

            if ( loc1.equals(loc2) ||
                    (loc1.getStart() == loc2.getStart() && (current1.getAlternateAlleles().size() > 1 || current2.getAlternateAlleles().size() > 1)) ) {
                // test the alleles
                if ( determineAndWriteOverlap(current1, current2, status) ) {
                    sourceVCs1.remove(currentIndex1);
                    sourceVCs2.remove(currentIndex2);
                    size1--;
                    size2--;
                } else {
                    currentIndex1++;
                    currentIndex2++;
                }
                current1 = (currentIndex1 < size1 ? sourceVCs1.get(currentIndex1): null);
                current2 = (currentIndex2 < size2 ? sourceVCs2.get(currentIndex2): null);
            } else if ( loc1.isBefore(loc2) ) {
                currentIndex1++;
                current1 = (currentIndex1 < size1 ? sourceVCs1.get(currentIndex1): null);
            } else {
                currentIndex2++;
                current2 = (currentIndex2 < size2 ? sourceVCs2.get(currentIndex2): null);
            }
        }
    }

    private boolean determineAndWriteOverlap(final VariantContext vc1, final VariantContext vc2, final String status) {
        final int allelesFrom1In2 = findOverlap(vc1, vc2);
        final int allelesFrom2In1 = findOverlap(vc2, vc1);
        final int totalAllelesIn1 = vc1.getAlternateAlleles().size();
        final int totalAllelesIn2 = vc2.getAlternateAlleles().size();

        final boolean allAllelesFrom1Overlap = allelesFrom1In2 == totalAllelesIn1;
        final boolean allAllelesFrom2Overlap = allelesFrom2In1 == totalAllelesIn2;

        boolean thereIsOverlap = true;

        if ( allAllelesFrom1Overlap && allAllelesFrom2Overlap ) {
            writeOne(vc1, INTERSECTION_SET, status);
        } else if ( allAllelesFrom1Overlap ) {
            writeOne(vc2, INTERSECTION_SET, source1 + "IsSubsetOf" + source2);
        } else if ( allAllelesFrom2Overlap ) {
            writeOne(vc1, INTERSECTION_SET, source2 + "IsSubsetOf" + source1);
        } else if ( allelesFrom1In2 > 0 ) {
            writeOne(vc1, INTERSECTION_SET, SOME_ALLELES_MATCH_STATUS);
        } else if ( totalAllelesIn1 > 1 || totalAllelesIn2 > 1 ) { // we don't handle multi-allelics in the haplotype-based reconstruction
            writeOne(vc1, INTERSECTION_SET, SAME_START_DIFFERENT_ALLELES_STATUS);
        } else {
            thereIsOverlap = false;
        }

        return thereIsOverlap;
    }

    private static int findOverlap(final VariantContext target, final VariantContext comparison) {
        int overlap = 0;
        for ( final Allele allele : target.getAlternateAlleles() ) {
            if ( comparison.hasAlternateAllele(allele) )
                overlap++;
        }
        return overlap;
    }

    private static final int SW_MATCH = 40;
    private static final int SW_MISMATCH = -100;
    private static final int SW_GAP = -250;
    private static final int SW_GAP_EXTEND = -13;

    private void resolveByHaplotype(final ReferenceContext refContext) {

        final byte[] source1Haplotype = generateHaplotype(sourceVCs1, refContext);
        final byte[] source2Haplotype = generateHaplotype(sourceVCs2, refContext);

        final SWPairwiseAlignment swConsensus1 = new SWPairwiseAlignment( refContext.getBases(), source1Haplotype, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );
        final SWPairwiseAlignment swConsensus2 = new SWPairwiseAlignment( refContext.getBases(), source2Haplotype, SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND );

        // protect against SW failures
        if( swConsensus1.getCigar().toString().contains("S") || swConsensus1.getCigar().getReferenceLength() < 20 ||
                swConsensus2.getCigar().toString().contains("S") || swConsensus2.getCigar().getReferenceLength() < 20 ) {
            // TODO -- handle errors appropriately
            logger.debug("Bad SW alignment; aborting at " + refContext.getLocus());
            return;
        }

        // order results by start position
        final TreeMap<Integer, VariantContext> source1Map = new TreeMap<Integer, VariantContext>(HaplotypeCallerGenotypingEngine.generateVCsFromAlignment(new Haplotype(source1Haplotype, false, 0, swConsensus1.getCigar()), refContext.getBases(), refContext.getWindow(), source1));
        final TreeMap<Integer, VariantContext> source2Map = new TreeMap<Integer, VariantContext>(HaplotypeCallerGenotypingEngine.generateVCsFromAlignment(new Haplotype(source2Haplotype, false, 0, swConsensus2.getCigar()), refContext.getBases(), refContext.getWindow(), source2));
        if ( source1Map.size() == 0 || source2Map.size() == 0 ) {
            // TODO -- handle errors appropriately
            logger.debug("No source alleles; aborting at " + refContext.getLocus());
            return;
        }

        // create lists and test for equality
        final List<VariantContext> source1Alleles = new ArrayList<VariantContext>(source1Map.values());
        final List<VariantContext> source2Alleles = new ArrayList<VariantContext>(source2Map.values());

        writeAndPurgeAllEqualVariants(source1Alleles, source2Alleles, SAME_BY_HAPLOTYPE_STATUS);
        if ( source1Alleles.isEmpty() ) {
            writeAll(source2Alleles, source2, null);
        } else if ( source2Alleles.isEmpty() ) {
            writeAll(source1Alleles, source1, null);
        } else {
            writeDifferences(source1Alleles, source2Alleles);
        }
    }

    private byte[] generateHaplotype(final List<VariantContext> sourceVCs, final ReferenceContext refContext) {

        final StringBuilder sb = new StringBuilder();

        final int startPos = refContext.getWindow().getStart();
        int currentPos = startPos;
        final byte[] reference = refContext.getBases();

        for ( final VariantContext vc : sourceVCs ) {
            // add any missing reference context
            int vcStart = vc.getStart();
            final int refAlleleLength = vc.getReference().length();
            if ( refAlleleLength == vc.getEnd() - vc.getStart() ) // this is a deletion (whereas for other events the padding base isn't part of the position)
                vcStart++;

            while ( currentPos < vcStart )
                sb.append((char)reference[currentPos++ - startPos]);

            // add the alt allele
            sb.append(vc.getAlternateAllele(0).getBaseString());

            // skip the reference allele
            currentPos += refAlleleLength;
        }
        // add any missing reference context
        final int stopPos = refContext.getWindow().getStop();
        while ( currentPos < stopPos )
            sb.append((char)reference[currentPos++ - startPos]);

        return sb.toString().getBytes();
    }

    private void writeDifferences(final List<VariantContext> source1Alleles, final List<VariantContext> source2Alleles) {
        int currentIndex1 = 0, currentIndex2 = 0;
        final int size1 = source1Alleles.size(), size2 = source2Alleles.size();
        VariantContext current1 = source1Alleles.get(0);
        VariantContext current2 = source2Alleles.get(0);

        while ( currentIndex1 < size1 || currentIndex2 < size2 ) {
            if ( current1 == null ) {
                writeOne(current2, source2, null);
                currentIndex2++;
                current2 = (currentIndex2 < size2 ? source2Alleles.get(currentIndex2): null);
            } else if ( current2 == null ) {
                writeOne(current1, source1, null);
                currentIndex1++;
                current1 = (currentIndex1 < size1 ? source1Alleles.get(currentIndex1): null);
            } else {

                final GenomeLoc loc1 = getToolkit().getGenomeLocParser().createGenomeLoc(current1);
                final GenomeLoc loc2 = getToolkit().getGenomeLocParser().createGenomeLoc(current2);

                if ( loc1.getStart() == loc2.getStart() || loc1.overlapsP(loc2) ) {
                    String status;
                    if ( loc1.getStart() == loc2.getStart() ) {
                        final String allele1 = current1.getAlternateAllele(0).getBaseString();
                        final String allele2 = current2.getAlternateAllele(0).getBaseString();
                        if ( allele1.indexOf(allele2) != -1 || allele2.indexOf(allele1) != -1 )
                            status = ONE_ALLELE_SUBSET_OF_OTHER_STATUS;
                        else
                            status = SAME_START_DIFFERENT_ALLELES_STATUS;
                    } else {
                        status = OVERLAPPING_EVENTS_STATUS;
                    }

                    writeOne(current1, INTERSECTION_SET, status);
                    currentIndex1++;
                    currentIndex2++;
                    current1 = (currentIndex1 < size1 ? source1Alleles.get(currentIndex1): null);
                    current2 = (currentIndex2 < size2 ? source2Alleles.get(currentIndex2): null);
                } else if ( loc1.isBefore(loc2) ) {
                    writeOne(current1, source1, null);
                    currentIndex1++;
                    current1 = (currentIndex1 < size1 ? source1Alleles.get(currentIndex1): null);
                } else {
                    writeOne(current2, source2, null);
                    currentIndex2++;
                    current2 = (currentIndex2 < size2 ? source2Alleles.get(currentIndex2): null);
                }
            }
        }
    }
}
