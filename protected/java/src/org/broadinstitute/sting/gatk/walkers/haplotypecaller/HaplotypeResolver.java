/*
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriterFactory;

import java.util.*;

/**
 * Haplotype-based resolution of variants in 2 different eval files.
 *
 * <p>
 * HaplotypeResolver is a tool that takes 2 VCF files and constructs haplotypes based on the variants inside them.
 * From that, it can resolve potential differences in variant calls that are inherently the same (or similar) variants.
 * Records are annotated with the set and status attributes.
 *
 * <h2>Input</h2>
 * <p>
 * 2 variant files to resolve.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A single consensus VCF.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx1g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T HaplotypeResolver \
 *   -V:v1 input1.vcf \
 *   -V:v2 input2.vcf \
 *   -o output.vcf
 * </pre>
 *
 */
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

    @Output(doc="File to which variants should be written", required=true)
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
        final Map<String, Object> attrs = new HashMap<String, Object>(vc.getAttributes());
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

    private static final double SW_MATCH = 4.0;
    private static final double SW_MISMATCH = -10.0;
    private static final double SW_GAP = -25.0;
    private static final double SW_GAP_EXTEND = -1.3;
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
        final TreeMap<Integer, VariantContext> source1Map = new TreeMap<Integer, VariantContext>(GenotypingEngine.generateVCsFromAlignment(0, swConsensus1.getCigar(), refContext.getBases(), source1Haplotype, refContext.getWindow(), source1, 0));
        final TreeMap<Integer, VariantContext> source2Map = new TreeMap<Integer, VariantContext>(GenotypingEngine.generateVCsFromAlignment(0, swConsensus2.getCigar(), refContext.getBases(), source2Haplotype, refContext.getWindow(), source2, 0));
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
