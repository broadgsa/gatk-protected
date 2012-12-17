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

import com.google.java.contract.Ensures;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardCallerArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.pairhmm.PairHMM;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Call SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. Haplotypes are evaluated using an affine gap penalty Pair HMM.
 *
 * <h2>Input</h2>
 * <p>
 * Input bam file(s) from which to make calls
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * VCF file with raw, unrecalibrated SNP and indel calls.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *   java
 *     -jar GenomeAnalysisTK.jar
 *     -T HaplotypeCaller
 *     -R reference/human_g1k_v37.fasta
 *     -I sample1.bam [-I sample2.bam ...] \
 *     --dbsnp dbSNP.vcf \
 *     -stand_call_conf [50.0] \
 *     -stand_emit_conf 10.0 \
 *     [-L targets.interval_list]
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * <h2>Caveats</h2>
 * <ul>
 * <li>The system is under active and continuous development. All outputs, the underlying likelihood model, and command line arguments are likely to change often.</li>
 * </ul>
 *
 * @author rpoplin
 * @since 8/22/11
 */

@DocumentedGATKFeature( groupName = "Variant Discovery Tools", extraDocs = {CommandLineGATK.class} )
@PartitionBy(PartitionType.LOCUS)
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ActiveRegionExtension(extension=65, maxRegion=300)
public class HaplotypeCaller extends ActiveRegionWalker<Integer, Integer> implements AnnotatorCompatible {

    /**
     * A raw, unfiltered, highly sensitive callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VariantContextWriter vcfWriter = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false)
    protected PrintStream graphWriter = null;

    /**
     * The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    @Argument(fullName = "pair_hmm_implementation", shortName = "pairHMM", doc = "The PairHMM implementation to use for genotype likelihood calculations", required = false)
    public PairHMM.HMM_IMPLEMENTATION pairHMM = PairHMM.HMM_IMPLEMENTATION.LOGLESS_CACHING;

    @Hidden
    @Argument(fullName="keepRG", shortName="keepRG", doc="Only use read from this read group when making calls (but use all reads to build the assembly)", required = false)
    protected String keepRG = null;

    @Argument(fullName="minPruning", shortName="minPruning", doc = "The minimum allowed pruning factor in assembly graph. Paths with <= X supporting kmers are pruned from the graph", required = false)
    protected int MIN_PRUNE_FACTOR = 1;

    @Advanced
    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="Flat gap continuation penalty for use in the Pair HMM", required = false)
    protected int gcpHMM = 10;

    @Argument(fullName="downsampleRegion", shortName="dr", doc="coverage, per-sample, to downsample each active region to", required = false)
    protected int DOWNSAMPLE_PER_SAMPLE_PER_REGION = 1000;

    @Argument(fullName="useAllelesTrigger", shortName="allelesTrigger", doc = "If specified, use additional trigger on variants found in an external alleles file", required=false)
    protected boolean USE_ALLELES_TRIGGER = false;

    @Advanced
    @Argument(fullName="useFilteredReadsForAnnotations", shortName="useFilteredReadsForAnnotations", doc = "If specified, use the contamination-filtered read maps for the purposes of annotating variants", required=false)
    protected boolean USE_FILTERED_READ_MAP_FOR_ANNOTATIONS = false;

    /**
     * rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();
    public RodBinding<VariantContext> getDbsnpRodBinding() { return dbsnp.dbsnp; }

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     *  as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field).
     *  Records that are filtered in the comp track will be ignored.
     *  Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Input(fullName="comp", shortName = "comp", doc="comparison VCF file", required=false)
    public List<RodBinding<VariantContext>> comps = Collections.emptyList();
    public List<RodBinding<VariantContext>> getCompRodBindings() { return comps; }

    // The following are not used by the Unified Genotyper
    public RodBinding<VariantContext> getSnpEffRodBinding() { return null; }
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return Collections.emptyList(); }
    public boolean alwaysAppendDbsnpId() { return false; }

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>(Arrays.asList(new String[]{"ClippingRankSumTest"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the -A or -G arguments,
     * so annotations will be excluded even if they are explicitly included with the other options.
     */
    @Argument(fullName="excludeAnnotation", shortName="XA", doc="One or more specific annotations to exclude", required=false)
    protected List<String> annotationsToExclude = new ArrayList<String>(Arrays.asList(new String[]{"SpanningDeletions", "TandemRepeatAnnotator"}));

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    @ArgumentCollection
    private StandardCallerArgumentCollection SCAC = new StandardCallerArgumentCollection();

    @Argument(fullName="debug", shortName="debug", doc="If specified, print out very verbose debug information about each triggering active region", required = false)
    protected boolean DEBUG;

    // the UG engines
    private UnifiedGenotyperEngine UG_engine = null;
    private UnifiedGenotyperEngine UG_engine_simple_genotyper = null;

    // the assembly engine
    private LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    private LikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    private GenotypingEngine genotypingEngine = null;

    // fasta reference reader to supplement the edges of the reference sequence
    private CachingIndexedFastaSequenceFile referenceReader;

    // reference base padding size
    private static final int REFERENCE_PADDING = 900;

    // bases with quality less than or equal to this value are trimmed off the tails of the reads
    private static final byte MIN_TAIL_QUALITY = 20;

    private ArrayList<String> samplesList = new ArrayList<String>();
    private final static double LOG_ONE_HALF = -Math.log10(2.0);
    private final static double LOG_ONE_THIRD = -Math.log10(3.0);
    private final ArrayList<VariantContext> allelesToGenotype = new ArrayList<VariantContext>();

    private final static Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private final static Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        super.initialize();

        // get all of the unique sample names
        Set<String> samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());
        samplesList.addAll( samples );
        // initialize the UnifiedGenotyper Engine which is used to call into the exact model
        final UnifiedArgumentCollection UAC = new UnifiedArgumentCollection( SCAC ); // this adapter is used so that the full set of unused UG arguments aren't exposed to the HC user
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples, VariantContextUtils.DEFAULT_PLOIDY);

        // create a UAC but with the exactCallsLog = null, so we only output the log for the HC caller itself, if requested
        UnifiedArgumentCollection simpleUAC = new UnifiedArgumentCollection(UAC);
        simpleUAC.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;
        simpleUAC.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY;
        simpleUAC.STANDARD_CONFIDENCE_FOR_CALLING = Math.min( 4.0, UAC.STANDARD_CONFIDENCE_FOR_CALLING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.STANDARD_CONFIDENCE_FOR_EMITTING = Math.min( 4.0, UAC.STANDARD_CONFIDENCE_FOR_EMITTING ); // low values used for isActive determination only, default/user-specified values used for actual calling
        simpleUAC.CONTAMINATION_FRACTION = 0.0;
        simpleUAC.exactCallsLog = null;
        UG_engine_simple_genotyper = new UnifiedGenotyperEngine(getToolkit(), simpleUAC, logger, null, null, samples, VariantContextUtils.DEFAULT_PLOIDY);

        // initialize the output VCF header
        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());
        // all callers need to add these standard annotation header lines
        VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true,
                VCFConstants.DOWNSAMPLED_KEY,
                VCFConstants.MLE_ALLELE_COUNT_KEY,
                VCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        // FILTER fields are added unconditionally as it's not always 100% certain the circumstances
        // where the filters are used.  For example, in emitting all sites the lowQual field is used
        headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));

        vcfWriter.writeHeader(new VCFHeader(headerInfo, samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        assemblyEngine = new SimpleDeBruijnAssembler( DEBUG, graphWriter );
        likelihoodCalculationEngine = new LikelihoodCalculationEngine( (byte)gcpHMM, DEBUG, pairHMM );
        genotypingEngine = new GenotypingEngine( DEBUG, annotationEngine, USE_FILTERED_READ_MAP_FOR_ANNOTATIONS );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // isActive
    //
    //---------------------------------------------------------------------------------------------------------------

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary and extended reads in the active region
    @Override
    public EnumSet<ActiveRegionReadState> desiredReadStates() {
        return EnumSet.of(
                ActiveRegionReadState.PRIMARY,
                ActiveRegionReadState.NONPRIMARY,
                ActiveRegionReadState.EXTENDED
        );
    }

    @Override
    @Ensures({"result.isActiveProb >= 0.0", "result.isActiveProb <= 1.0"})
    public ActivityProfileResult isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {

        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for( final VariantContext vc : tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()) ) {
                if( !allelesToGenotype.contains(vc) ) {
                    allelesToGenotype.add(vc); // save for later for processing during the ActiveRegion's map call. Should be folded into a RefMetaDataTracker object
                }
            }
            if( tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()).size() > 0 ) {
                return new ActivityProfileResult(ref.getLocus(), 1.0);
            }
        }

        if( USE_ALLELES_TRIGGER ) {
            return new ActivityProfileResult( ref.getLocus(), tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null ) { return new ActivityProfileResult(ref.getLocus(), 0.0); }

        final List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
        noCall.add(Allele.NO_CALL);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();
        for( final Map.Entry<String, AlignmentContext> sample : splitContexts.entrySet() ) {
            final double[] genotypeLikelihoods = new double[3]; // ref versus non-ref (any event)
            Arrays.fill(genotypeLikelihoods, 0.0);

            for( final PileupElement p : sample.getValue().getBasePileup() ) {
                final byte qual = p.getQual();
                if( p.isDeletion() || qual > (byte) 18) {
                    int AA = 0; final int AB = 1; int BB = 2;
                    if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeDeletedBase() || p.isAfterDeletedBase() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ) {
                        AA = 2;
                        BB = 0;
                        if( p.isNextToSoftClip() ) {
                            averageHQSoftClips.add(AlignmentUtils.calcNumHighQualitySoftClips(p.getRead(), (byte) 28));
                        }
                    }
                    genotypeLikelihoods[AA] += p.getRepresentativeCount() * QualityUtils.qualToProbLog10(qual);
                    genotypeLikelihoods[AB] += p.getRepresentativeCount() * MathUtils.approximateLog10SumLog10( QualityUtils.qualToProbLog10(qual) + LOG_ONE_HALF, QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD + LOG_ONE_HALF );
                    genotypeLikelihoods[BB] += p.getRepresentativeCount() * QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD;
                }
            }
            genotypes.add( new GenotypeBuilder(sample.getKey()).alleles(noCall).PL(genotypeLikelihoods).make() );
        }

        final ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add( FAKE_REF_ALLELE );
        alleles.add( FAKE_ALT_ALLELE );
        final VariantCallContext vcOut = UG_engine_simple_genotyper.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.INDEL);
        final double isActiveProb = vcOut == null ? 0.0 : QualityUtils.qualToProb( vcOut.getPhredScaledQual() );

        return new ActivityProfileResult( ref.getLocus(), isActiveProb, averageHQSoftClips.mean() > 6.0 ? ActivityProfileResult.ActivityProfileResultState.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileResult.ActivityProfileResultState.NONE, averageHQSoftClips.mean() );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer map( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion, final RefMetaDataTracker metaDataTracker ) {

        final ArrayList<VariantContext> activeAllelesToGenotype = new ArrayList<VariantContext>();

        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for( final VariantContext vc : allelesToGenotype ) {
                if( activeRegion.getLocation().overlapsP( getToolkit().getGenomeLocParser().createGenomeLoc(vc) ) ) {
                    activeAllelesToGenotype.add(vc); // do something with these VCs during GGA mode
                }
            }
            allelesToGenotype.removeAll( activeAllelesToGenotype );
        }

        if( !activeRegion.isActive ) { return 0; } // Not active so nothing to do!
        if( activeRegion.size() == 0 && UG_engine.getUAC().GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) { return 0; } // No reads here so nothing to do!
        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES && activeAllelesToGenotype.isEmpty() ) { return 0; } // No alleles found in this region so nothing to do!

        finalizeActiveRegion( activeRegion ); // merge overlapping fragments, clip adapter and low qual tails
        final Haplotype referenceHaplotype = new Haplotype(activeRegion.getActiveRegionReference(referenceReader)); // Create the reference haplotype which is the bases from the reference that make up the active region
        referenceHaplotype.setIsReference(true);
        final byte[] fullReferenceWithPadding = activeRegion.getFullReference(referenceReader, REFERENCE_PADDING);
        //int PRUNE_FACTOR = Math.max(MIN_PRUNE_FACTOR, determinePruneFactorFromCoverage( activeRegion ));
        final ArrayList<Haplotype> haplotypes = assemblyEngine.runLocalAssembly( activeRegion, referenceHaplotype, fullReferenceWithPadding, getPaddedLoc(activeRegion), MIN_PRUNE_FACTOR, activeAllelesToGenotype );
        if( haplotypes.size() == 1 ) { return 1; } // only the reference haplotype remains so nothing else to do!

        activeRegion.hardClipToActiveRegion(); // only evaluate the parts of reads that are overlapping the active region
        final List<GATKSAMRecord> filteredReads = filterNonPassingReads( activeRegion ); // filter out reads from genotyping which fail mapping quality based criteria
        if( activeRegion.size() == 0 ) { return 1; } // no reads remain after filtering so nothing else to do!

        // sort haplotypes to take full advantage of haplotype start offset optimizations in PairHMM
        Collections.sort( haplotypes, new Haplotype.HaplotypeBaseComparator() );

        // evaluate each sample's reads against all haplotypes
        final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = likelihoodCalculationEngine.computeReadLikelihoods( haplotypes, splitReadsBySample( activeRegion.getReads() ) );
        final Map<String, ArrayList<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );

        // subset down to only the best haplotypes to be genotyped in all samples ( in GGA mode use all discovered haplotypes )
        final ArrayList<Haplotype> bestHaplotypes = ( UG_engine.getUAC().GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? likelihoodCalculationEngine.selectBestHaplotypes( haplotypes, stratifiedReadMap ) : haplotypes );

        for( final VariantContext call : genotypingEngine.assignGenotypeLikelihoodsAndCallIndependentEvents( UG_engine,
                                                                                                             bestHaplotypes,
                                                                                                             samplesList,
                                                                                                             stratifiedReadMap,
                                                                                                             perSampleFilteredReadList,
                                                                                                             fullReferenceWithPadding,
                                                                                                             getPaddedLoc(activeRegion),
                                                                                                             activeRegion.getLocation(),
                                                                                                             getToolkit().getGenomeLocParser(),
                                                                                                             activeAllelesToGenotype ) ) {
            vcfWriter.add( call );
        }

        if( DEBUG ) { System.out.println("----------------------------------------------------------------------------------"); }

        return 1; // One active region was processed during this map call
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer cur, Integer sum) {
        return cur + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
        logger.info("Ran local assembly on " + result + " active regions");
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // private helper functions
    //
    //---------------------------------------------------------------------------------------------------------------

    private void finalizeActiveRegion( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion ) {
        if( DEBUG ) { System.out.println("\nAssembling " + activeRegion.getLocation() + " with " + activeRegion.size() + " reads:    (with overlap region = " + activeRegion.getExtendedLoc() + ")"); }
        final ArrayList<GATKSAMRecord> finalizedReadList = new ArrayList<GATKSAMRecord>();
        final FragmentCollection<GATKSAMRecord> fragmentCollection = FragmentUtils.create( ReadUtils.sortReadsByCoordinate(activeRegion.getReads()) );
        activeRegion.clearReads();

        // Join overlapping paired reads to create a single longer read
        finalizedReadList.addAll( fragmentCollection.getSingletonReads() );
        for( final List<GATKSAMRecord> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
            finalizedReadList.addAll( FragmentUtils.mergeOverlappingPairedFragments(overlappingPair) );
        }

        Collections.shuffle(finalizedReadList, GenomeAnalysisEngine.getRandomGenerator());

        // Loop through the reads hard clipping the adaptor and low quality tails
        for( final GATKSAMRecord myRead : finalizedReadList ) {
            final GATKSAMRecord postAdapterRead = ( myRead.getReadUnmappedFlag() ? myRead : ReadClipper.hardClipAdaptorSequence( myRead ) );
            if( postAdapterRead != null && !postAdapterRead.isEmpty() && postAdapterRead.getCigar().getReadLength() > 0 ) {
                final GATKSAMRecord clippedRead = ReadClipper.hardClipLowQualEnds( postAdapterRead, MIN_TAIL_QUALITY );
                // protect against INTERVALS with abnormally high coverage
                // BUGBUG: remove when positional downsampler is hooked up to ART/HC
                if( clippedRead.getReadLength() > 0 && activeRegion.size() < samplesList.size() * DOWNSAMPLE_PER_SAMPLE_PER_REGION ) {
                    activeRegion.add(clippedRead);
                }
            }
        }
    }

    private List<GATKSAMRecord> filterNonPassingReads( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion ) {
        final ArrayList<GATKSAMRecord> readsToRemove = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getReadLength() < 24 || rec.getMappingQuality() < 20 || BadMateFilter.hasBadMate(rec) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion ) {
        final int padLeft = Math.max(activeRegion.getReferenceLoc().getStart()-REFERENCE_PADDING, 1);
        final int padRight = Math.min(activeRegion.getReferenceLoc().getStop()+REFERENCE_PADDING, referenceReader.getSequenceDictionary().getSequence(activeRegion.getReferenceLoc().getContig()).getSequenceLength());
        return getToolkit().getGenomeLocParser().createGenomeLoc(activeRegion.getReferenceLoc().getContig(), padLeft, padRight);
    }

    private HashMap<String, ArrayList<GATKSAMRecord>> splitReadsBySample( final List<GATKSAMRecord> reads ) {
        final HashMap<String, ArrayList<GATKSAMRecord>> returnMap = new HashMap<String, ArrayList<GATKSAMRecord>>();
        for( final String sample : samplesList) {
            ArrayList<GATKSAMRecord> readList = returnMap.get( sample );
            if( readList == null ) {
                readList = new ArrayList<GATKSAMRecord>();
                returnMap.put(sample, readList);
            }
        }
        for( final GATKSAMRecord read : reads ) {
            returnMap.get(read.getReadGroup().getSample()).add(read);
        }

        return returnMap;
    }

    /*
    private int determinePruneFactorFromCoverage( final ActiveRegion activeRegion ) {
        final ArrayList<Integer> readLengthDistribution = new ArrayList<Integer>();
        for( final GATKSAMRecord read : activeRegion.getReads() ) {
            readLengthDistribution.add(read.getReadLength());
        }
        final double meanReadLength = MathUtils.average(readLengthDistribution);
        final double meanCoveragePerSample = (double) activeRegion.getReads().size() / ((double) activeRegion.getExtendedLoc().size() / meanReadLength) / (double) samplesList.size();
        int PRUNE_FACTOR = 0;
        if( meanCoveragePerSample > 8.5 ) {
            PRUNE_FACTOR = (int) Math.floor( Math.sqrt( meanCoveragePerSample - 5.0 ) );
        } else if( meanCoveragePerSample > 3.0 ) {
            PRUNE_FACTOR = 1;
        }

        if( DEBUG ) { System.out.println(String.format("Mean coverage per sample = %.1f --> prune factor = %d", meanCoveragePerSample, PRUNE_FACTOR)); }
        return PRUNE_FACTOR;
    }
    */
}