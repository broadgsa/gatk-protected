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
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

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
 *     -I input.bam
 *     -o output.raw.snps.indels.vcf
 * </pre>
 *
 * @author rpoplin
 * @since 8/22/11
 */

@PartitionBy(PartitionType.LOCUS)
@ActiveRegionExtension(extension=65, maxRegion=275)
public class HaplotypeCaller extends ActiveRegionWalker<Integer, Integer> implements AnnotatorCompatibleWalker {

    /**
     * A raw, unfiltered, highly specific callset in VCF format.
     */
    @Output(doc="File to which variants should be written", required = true)
    protected VariantContextWriter vcfWriter = null;

    @Output(fullName="graphOutput", shortName="graph", doc="File to which debug assembly graph information should be written", required = false)
    protected PrintStream graphWriter = null;

    @Argument(fullName = "assembler", shortName = "assembler", doc = "Assembler to use; currently only SIMPLE_DE_BRUIJN is available.", required = false)
    protected LocalAssemblyEngine.ASSEMBLER ASSEMBLER_TO_USE = LocalAssemblyEngine.ASSEMBLER.SIMPLE_DE_BRUIJN;

    @Argument(fullName="keepRG", shortName="keepRG", doc="keepRG", required = false)
    protected String keepRG = null;
    
    @Argument(fullName="mnpLookAhead", shortName="mnpLookAhead", doc = "The number of bases to combine together to form MNPs out of nearby consecutive SNPs on the same haplotype", required = false)
    protected int MNP_LOOK_AHEAD = 0;

    @Argument(fullName="minPruning", shortName="minPruning", doc = "The minimum allowed pruning factor in assembly graph. Paths with <= X supporting kmers are pruned from the graph", required = true)
    protected int MIN_PRUNE_FACTOR = 0;

    @Argument(fullName="genotypeFullActiveRegion", shortName="genotypeFullActiveRegion", doc = "If specified, alternate alleles are considered to be the full active region for the purposes of genotyping", required = false)
    protected boolean GENOTYPE_FULL_ACTIVE_REGION = false;

    @Argument(fullName="fullHaplotype", shortName="fullHaplotype", doc = "If specified, output the full haplotype sequence instead of converting to individual variants w.r.t. the reference", required = false)
    protected boolean OUTPUT_FULL_HAPLOTYPE_SEQUENCE = false;

    @Argument(fullName="gcpHMM", shortName="gcpHMM", doc="gcpHMM", required = false)
    protected int gcpHMM = 10;

    @Argument(fullName="downsampleRegion", shortName="dr", doc="coverage per sample to downsample each region to", required = false)
    protected int DOWNSAMPLE_PER_SAMPLE_PER_REGION = 1000;

    @Argument(fullName="useExpandedTriggerSet", shortName="expandedTriggers", doc = "If specified, use additional, experimental triggers designed to capture larger indels but which may lead to an increase in the false positive rate", required=false)
    protected boolean USE_EXPANDED_TRIGGER_SET = false;

    @Argument(fullName="useAllelesTrigger", shortName="allelesTrigger", doc = "If specified, use additional trigger on variants found in an external alleles file", required=false)
    protected boolean USE_ALLELES_TRIGGER = false;

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
    protected List<String> annotationsToExclude = new ArrayList<String>(Arrays.asList(new String[]{"HaplotypeScore", "MappingQualityZero", "SpanningDeletions", "TandemRepeatAnnotator"}));

    /**
     * Which groups of annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available groups.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    @ArgumentCollection
    private UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // the calculation arguments
    private UnifiedGenotyperEngine UG_engine = null;
    private UnifiedGenotyperEngine UG_engine_simple_genotyper = null;
    
    @Argument(fullName="debug", shortName="debug", doc="If specified, print out very verbose debug information about each triggering active region", required = false)
    protected boolean DEBUG;

    @Argument(fullName="doBanded", shortName="doBanded", doc="If specified, use the banded option", required = false)
    protected boolean doBanded = false;

    // the assembly engine
    LocalAssemblyEngine assemblyEngine = null;

    // the likelihoods engine
    LikelihoodCalculationEngine likelihoodCalculationEngine = null;

    // the genotyping engine
    GenotypingEngine genotypingEngine = null;

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    // fasta reference reader to supplement the edges of the reference sequence
    private IndexedFastaSequenceFile referenceReader;

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
        UG_engine = new UnifiedGenotyperEngine(getToolkit(), UAC.clone(), logger, null, null, samples, VariantContextUtils.DEFAULT_PLOIDY);
        UAC.OutputMode = UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY; // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.GenotypingMode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.DISCOVERY; // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.STANDARD_CONFIDENCE_FOR_CALLING = (USE_EXPANDED_TRIGGER_SET ? 0.3 : Math.max( 4.0, UAC.STANDARD_CONFIDENCE_FOR_CALLING) ); // low values used for isActive determination only, default/user-specified values used for actual calling
        UAC.STANDARD_CONFIDENCE_FOR_EMITTING = (USE_EXPANDED_TRIGGER_SET ? 0.3 : Math.max( 4.0, UAC.STANDARD_CONFIDENCE_FOR_EMITTING) ); // low values used for isActive determination only, default/user-specified values used for actual calling
        UG_engine_simple_genotyper = new UnifiedGenotyperEngine(getToolkit(), UAC, logger, null, null, samples, VariantContextUtils.DEFAULT_PLOIDY);

        // initialize the output VCF header
        annotationEngine = new VariantAnnotatorEngine(Arrays.asList(annotationClassesToUse), annotationsToUse, annotationsToExclude, this, getToolkit());

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
        // header lines for the experimental HaplotypeCaller-specific annotations
        headerInfo.add(new VCFInfoHeaderLine("NVH", 1, VCFHeaderLineType.Integer, "Number of variants found on the haplotype that contained this variant"));
        headerInfo.add(new VCFInfoHeaderLine("NumHapEval", 1, VCFHeaderLineType.Integer, "Number of haplotypes that were chosen for evaluation in this active region"));
        headerInfo.add(new VCFInfoHeaderLine("NumHapAssembly", 1, VCFHeaderLineType.Integer, "Number of haplotypes created during the assembly of this active region"));
        headerInfo.add(new VCFInfoHeaderLine("ActiveRegionSize", 1, VCFHeaderLineType.Integer, "Number of base pairs that comprise this active region"));
        headerInfo.add(new VCFInfoHeaderLine("EVENTLENGTH", 1, VCFHeaderLineType.Integer, "Max length of all the alternate alleles"));
        headerInfo.add(new VCFInfoHeaderLine("TYPE", 1, VCFHeaderLineType.String, "Type of event: SNP or INDEL"));
        headerInfo.add(new VCFInfoHeaderLine("extType", 1, VCFHeaderLineType.String, "Extended type of event: SNP, MNP, INDEL, or COMPLEX"));
        headerInfo.add(new VCFInfoHeaderLine("QDE", 1, VCFHeaderLineType.Float, "QD value divided by the number of variants found on the haplotype that contained this variant"));

        vcfWriter.writeHeader(new VCFHeader(headerInfo, samples));

        try {
            // fasta reference reader to supplement the edges of the reference sequence
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

        assemblyEngine = new SimpleDeBruijnAssembler( DEBUG, graphWriter );
        likelihoodCalculationEngine = new LikelihoodCalculationEngine( (byte)gcpHMM, DEBUG, doBanded );
        genotypingEngine = new GenotypingEngine( DEBUG, MNP_LOOK_AHEAD, OUTPUT_FULL_HAPLOTYPE_SEQUENCE );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // isActive
    //
    //---------------------------------------------------------------------------------------------------------------

    // enable deletions in the pileup
    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable non primary reads in the active region
    @Override
    public boolean wantsNonPrimaryReads() { return true; }

    @Override
    @Ensures({"result >= 0.0", "result <= 1.0"})
    public double isActive( final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context ) {

        if( UG_engine.getUAC().GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ) {
            for( final VariantContext vc : tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()) ) {
                if( !allelesToGenotype.contains(vc) ) {
                    allelesToGenotype.add(vc); // save for later for processing during the ActiveRegion's map call. Should be folded into a ReadMetaDataTracker object
                }
            }
            if( tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()).size() > 0 ) {
                return 1.0;
            }
        }

        if( USE_ALLELES_TRIGGER ) {
            return ( tracker.getValues(UG_engine.getUAC().alleles, ref.getLocus()).size() > 0 ? 1.0 : 0.0 );
        }

        if( context == null ) { return 0.0; }

        final List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
        noCall.add(Allele.NO_CALL);

        final Map<String, AlignmentContext> splitContexts = AlignmentContextUtils.splitContextBySampleName(context);
        final GenotypesContext genotypes = GenotypesContext.create(splitContexts.keySet().size());
        for( final String sample : splitContexts.keySet() ) {
            final double[] genotypeLikelihoods = new double[3]; // ref versus non-ref (any event)
            Arrays.fill(genotypeLikelihoods, 0.0);

            for( final PileupElement p : splitContexts.get(sample).getBasePileup() ) {
                final byte qual = ( USE_EXPANDED_TRIGGER_SET ?
                                        ( p.isNextToSoftClip() || p.isBeforeInsertion() || p.isAfterInsertion() ? ( p.getQual() > QualityUtils.MIN_USABLE_Q_SCORE ? p.getQual() : (byte) 20 ) : p.getQual() )
                                        : p.getQual() );
                if( p.isDeletion() || qual > (USE_EXPANDED_TRIGGER_SET ? QualityUtils.MIN_USABLE_Q_SCORE : (byte) 18) ) {
                    int AA = 0; final int AB = 1; int BB = 2;
                    if( USE_EXPANDED_TRIGGER_SET ) {
                        if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeDeletedBase() || p.isAfterDeletedBase() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ||
                                (!p.getRead().getNGSPlatform().equals(NGSPlatform.SOLID) && ((p.getRead().getReadPairedFlag() && p.getRead().getMateUnmappedFlag()) || BadMateFilter.hasBadMate(p.getRead()))) ) {
                            AA = 2;
                            BB = 0;
                        }
                    } else {
                        if( p.getBase() != ref.getBase() || p.isDeletion() || p.isBeforeDeletedBase() || p.isAfterDeletedBase() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ) {
                            AA = 2;
                            BB = 0;
                        }
                    }
                    genotypeLikelihoods[AA] += QualityUtils.qualToProbLog10(qual);
                    genotypeLikelihoods[AB] += MathUtils.approximateLog10SumLog10( QualityUtils.qualToProbLog10(qual) + LOG_ONE_HALF, QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD + LOG_ONE_HALF );
                    genotypeLikelihoods[BB] += QualityUtils.qualToErrorProbLog10(qual) + LOG_ONE_THIRD;
                }
            }
            genotypes.add( new GenotypeBuilder(sample).alleles(noCall).PL(genotypeLikelihoods).make() );
        }

        final ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add( FAKE_REF_ALLELE );
        alleles.add( FAKE_ALT_ALLELE );
        final VariantCallContext vcOut = UG_engine_simple_genotyper.calculateGenotypes(new VariantContextBuilder("HCisActive!", context.getContig(), context.getLocation().getStart(), context.getLocation().getStop(), alleles).genotypes(genotypes).make(), GenotypeLikelihoodsCalculationModel.Model.INDEL);
        return ( vcOut == null ? 0.0 : QualityUtils.qualToProb( vcOut.getPhredScaledQual() ) );
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public Integer map( final ActiveRegion activeRegion, final RefMetaDataTracker metaDataTracker ) {

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

        // evaluate each sample's reads against all haplotypes
        final HashMap<String, ArrayList<GATKSAMRecord>> perSampleReadList = splitReadsBySample( activeRegion.getReads() );
        final HashMap<String, ArrayList<GATKSAMRecord>> perSampleFilteredReadList = splitReadsBySample( filteredReads );
        likelihoodCalculationEngine.computeReadLikelihoods( haplotypes, perSampleReadList );

        // subset down to only the best haplotypes to be genotyped in all samples ( in GGA mode use all discovered haplotypes )
        final ArrayList<Haplotype> bestHaplotypes = ( UG_engine.getUAC().GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? likelihoodCalculationEngine.selectBestHaplotypes( haplotypes ) : haplotypes );

        for( final Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>> callResult :
                ( GENOTYPE_FULL_ACTIVE_REGION && UG_engine.getUAC().GenotypingMode != GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
                  ? genotypingEngine.assignGenotypeLikelihoodsAndCallHaplotypeEvents( UG_engine, bestHaplotypes, fullReferenceWithPadding, getPaddedLoc(activeRegion), activeRegion.getLocation(), getToolkit().getGenomeLocParser() )
                  : genotypingEngine.assignGenotypeLikelihoodsAndCallIndependentEvents( UG_engine, bestHaplotypes, fullReferenceWithPadding, getPaddedLoc(activeRegion), activeRegion.getLocation(), getToolkit().getGenomeLocParser(), activeAllelesToGenotype ) ) ) {
            if( DEBUG ) { System.out.println(callResult.getFirst().toStringWithoutGenotypes()); }

            final Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedReadMap = LikelihoodCalculationEngine.partitionReadsBasedOnLikelihoods( getToolkit().getGenomeLocParser(), perSampleReadList, perSampleFilteredReadList, callResult );
            final VariantContext annotatedCall = annotationEngine.annotateContext(stratifiedReadMap, callResult.getFirst());

            // add some custom annotations to the calls
            final Map<String, Object> myAttributes = new LinkedHashMap<String, Object>(annotatedCall.getAttributes());
            // Calculate the number of variants on the haplotype
            int maxNumVar = 0;
            for( final Allele allele : callResult.getFirst().getAlleles() ) {
                if( !allele.isReference() ) {
                    for( final Haplotype haplotype : callResult.getSecond().get(allele) ) {
                        final int numVar = haplotype.getEventMap().size();
                        if( numVar > maxNumVar ) { maxNumVar = numVar; }
                    }
                }
            }
            // Calculate the event length
            int maxLength = 0;
            for ( final Allele a : annotatedCall.getAlternateAlleles() ) {
                final int length = a.length() - annotatedCall.getReference().length();
                if( Math.abs(length) > Math.abs(maxLength) ) { maxLength = length; }
            }

            myAttributes.put("NVH", maxNumVar);
            myAttributes.put("NumHapEval", bestHaplotypes.size());
            myAttributes.put("NumHapAssembly", haplotypes.size());
            myAttributes.put("ActiveRegionSize", activeRegion.getLocation().size());
            myAttributes.put("EVENTLENGTH", maxLength);
            myAttributes.put("TYPE", (annotatedCall.isSNP() || annotatedCall.isMNP() ? "SNP" : "INDEL") );
            myAttributes.put("extType", annotatedCall.getType().toString() );

            //if( likelihoodCalculationEngine.haplotypeScore != null ) {
            //    myAttributes.put("HaplotypeScore", String.format("%.4f", likelihoodCalculationEngine.haplotypeScore));
            //}
            if( annotatedCall.hasAttribute("QD") ) {
                myAttributes.put("QDE", String.format("%.2f", Double.parseDouble((String)annotatedCall.getAttribute("QD")) / ((double)maxNumVar)) );
            }

            vcfWriter.add( new VariantContextBuilder(annotatedCall).attributes(myAttributes).make() );
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

    private void finalizeActiveRegion( final ActiveRegion activeRegion ) {
        if( DEBUG ) { System.out.println("\nAssembling " + activeRegion.getExtendedLoc() + " with " + activeRegion.size() + " reads:"); }
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
                if( clippedRead.getReadLength() > 0 && activeRegion.size() < samplesList.size() * DOWNSAMPLE_PER_SAMPLE_PER_REGION ) {
                    activeRegion.add(clippedRead);
                }
            }
        }
    }

    private List<GATKSAMRecord> filterNonPassingReads( final ActiveRegion activeRegion ) {
        final ArrayList<GATKSAMRecord> readsToRemove = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord rec : activeRegion.getReads() ) {
            if( rec.getReadLength() < 24 || rec.getMappingQuality() <= 20 || BadMateFilter.hasBadMate(rec) || (keepRG != null && !rec.getReadGroup().getId().equals(keepRG)) ) {
                readsToRemove.add(rec);
            }
        }
        activeRegion.removeAll( readsToRemove );
        return readsToRemove;
    }

    private GenomeLoc getPaddedLoc( final ActiveRegion activeRegion ) {
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