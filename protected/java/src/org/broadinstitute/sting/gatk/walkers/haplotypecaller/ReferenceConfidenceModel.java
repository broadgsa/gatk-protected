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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.haplotype.Haplotype;
import org.broadinstitute.sting.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.sting.utils.haplotypeBAMWriter.ReadDestination;
import org.broadinstitute.sting.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFSimpleHeaderLine;

import java.io.File;
import java.util.*;

/**
 * Code for estimating the reference confidence
 *
 * This code can estimate the probability that the data for a single sample is consistent with a
 * well-determined REF/REF diploid genotype.
 *
 * User: depristo
 * Date: 6/21/13
 * Time: 12:52 PM
 */
public class ReferenceConfidenceModel {

    //public final static String INDEL_INFORMATIVE_DEPTH = "CD"; // temporarily taking this extra genotype level information out for now
    public final static String ALTERNATE_ALLELE_STRING = "ALT"; // arbitrary alternate allele

    private final GenomeLocParser genomeLocParser;
    private final Set<String> samples;
    private final SAMFileHeader header; // TODO -- really shouldn't depend on this
    private final int indelInformativeDepthIndelSize;

    private final static boolean WRITE_DEBUGGING_BAM = false;
    private final SAMFileWriter debuggingWriter;

    private final static byte REF_MODEL_DELETION_QUAL = (byte) 30;

    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param genomeLocParser how we create genome locs
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     */
    public ReferenceConfidenceModel(final GenomeLocParser genomeLocParser,
                                    final Set<String> samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize) {
        if ( genomeLocParser == null ) throw new IllegalArgumentException("genomeLocParser cannot be null");
        if ( samples == null ) throw new IllegalArgumentException("samples cannot be null");
        if ( samples.isEmpty() ) throw new IllegalArgumentException("samples cannot be empty");
        if ( header == null ) throw new IllegalArgumentException("header cannot be empty");
        if ( indelInformativeDepthIndelSize < 0) throw new IllegalArgumentException("indelInformativeDepthIndelSize must be >= 1 but got " + indelInformativeDepthIndelSize);

        this.genomeLocParser = genomeLocParser;
        this.samples = samples;
        this.header = header;
        this.indelInformativeDepthIndelSize = indelInformativeDepthIndelSize;

        if ( WRITE_DEBUGGING_BAM ) {
            final SAMFileWriterFactory factory = new SAMFileWriterFactory();
            factory.setCreateIndex(true);
            debuggingWriter = factory.makeBAMWriter(header, false, new File("refCalc.bam"));
        } else {
            debuggingWriter = null;
        }

        initializeIndelPLCache();
    }

    /**
     * Get the VCF header lines to include when emitting reference confidence values via calculateRefConfidence
     * @return a non-null set of VCFHeaderLines
     */
    public Set<VCFHeaderLine> getVCFHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();
        // TODO - do we need a new kind of VCF Header subclass for specifying arbitrary alternate alleles?
        headerLines.add(new VCFSimpleHeaderLine(ALTERNATE_ALLELE_STRING, GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE_NAME, "Represents any possible alternative allele at this location"));
        //headerLines.add(new VCFFormatHeaderLine(INDEL_INFORMATIVE_DEPTH, 1, VCFHeaderLineType.Integer, "Number of reads at locus that are informative about an indel of size <= " + indelInformativeDepthIndelSize));
        return headerLines;
    }

    /**
     * Close down this reference model, closing down any debugging information opened during execution
     */
    public void close() {
        if ( debuggingWriter != null ) debuggingWriter.close();
    }


    /**
     * Calculate the reference confidence for a single sample given the its read data
     *
     * Returns a list of variant contexts, one for each position in the activeregion.getLoc(), each containing
     * detailed information about the certainty that the sample is hom-ref for each base in the region.
     *
     *
     *
     * @param refHaplotype the reference haplotype, used to get the reference bases across activeRegion.getLoc()
     * @param calledHaplotypes a list of haplotypes that segregate in this region, for realignment of the reads in the
     *                         stratifiedReadMap, corresponding to each reads best haplotype.  Must contain the refHaplotype.
     * @param paddedReferenceLoc the location of refHaplotype (which might be larger than activeRegion.getLoc())
     * @param activeRegion the active region we want to get the reference confidence over
     * @param stratifiedReadMap a map from a single sample to its PerReadAlleleLikelihoodMap for each haplotype in calledHaplotypes
     * @param variantCalls calls made in this region.  The return result will contain any variant call in this list in the
     *                     correct order by genomic position, and any variant in this list will stop us emitting a ref confidence
     *                     under any position it covers (for snps and insertions that is 1 bp, but for deletions its the entire ref span)
     * @return an ordered list of variant contexts that spans activeRegion.getLoc() and includes both reference confidence
     *         contexts as well as calls from variantCalls if any were provided
     */
    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final GenomeLoc paddedReferenceLoc,
                                                       final ActiveRegion activeRegion,
                                                       final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap,
                                                       final List<VariantContext> variantCalls) {
        if ( refHaplotype == null ) throw new IllegalArgumentException("refHaplotype cannot be null");
        if ( calledHaplotypes == null ) throw new IllegalArgumentException("calledHaplotypes cannot be null");
        if ( !calledHaplotypes.contains(refHaplotype)) throw new IllegalArgumentException("calledHaplotypes must contain the refHaplotype");
        if ( paddedReferenceLoc == null ) throw new IllegalArgumentException("paddedReferenceLoc cannot be null");
        if ( activeRegion == null ) throw new IllegalArgumentException("activeRegion cannot be null");
        if ( stratifiedReadMap == null ) throw new IllegalArgumentException("stratifiedReadMap cannot be null");
        if ( stratifiedReadMap.size() != 1 ) throw new IllegalArgumentException("stratifiedReadMap must contain exactly one sample but it contained " + stratifiedReadMap.size());
        if ( refHaplotype.length() != activeRegion.getExtendedLoc().size() ) throw new IllegalArgumentException("refHaplotype " + refHaplotype.length() + " and activeRegion location size " + activeRegion.getLocation().size() + " are different");

        final GenomeLoc refSpan = activeRegion.getLocation();
        final List<ReadBackedPileup> refPileups = getPileupsOverReference(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, refSpan, stratifiedReadMap);
        final byte[] ref = refHaplotype.getBases();
        final List<VariantContext> results = new ArrayList<>(refSpan.size());
        final String sampleName = stratifiedReadMap.keySet().iterator().next();

        final int globalRefOffset = refSpan.getStart() - activeRegion.getExtendedLoc().getStart();
        for ( final ReadBackedPileup pileup : refPileups ) {
            final GenomeLoc curPos = pileup.getLocation();
            final int offset = curPos.getStart() - refSpan.getStart();

            final VariantContext overlappingSite = getOverlappingVariantContext(curPos, variantCalls);
            if ( overlappingSite != null ) {
                // we have some overlapping site, add it to the list of positions
                if ( overlappingSite.getStart() == curPos.getStart() )
                    results.add(overlappingSite);
            } else {
                // otherwise emit a reference confidence variant context
                final int refOffset = offset + globalRefOffset;
                final byte refBase = ref[refOffset];
                final RefVsAnyResult homRefCalc = calcGenotypeLikelihoodsOfRefVsAny(pileup, refBase, (byte)6, null);
                homRefCalc.capByHomRefLikelihood();

                final Allele refAllele = Allele.create(refBase, true);
                final List<Allele> refSiteAlleles = Arrays.asList(refAllele, GATKVariantContextUtils.NON_REF_SYMBOLIC_ALLELE);
                final VariantContextBuilder vcb = new VariantContextBuilder("HC", curPos.getContig(), curPos.getStart(), curPos.getStart(), refSiteAlleles);
                final GenotypeBuilder gb = new GenotypeBuilder(sampleName, Arrays.asList(refAllele, refAllele));
                gb.AD(homRefCalc.AD_Ref_Any);
                gb.DP(homRefCalc.getDP());

                // genotype likelihood calculation
                final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(homRefCalc.genotypeLikelihoods);
                final int nIndelInformativeReads = calcNIndelInformativeReads(pileup, refOffset, ref, indelInformativeDepthIndelSize);
                final GenotypeLikelihoods indelGLs = getIndelPLs(nIndelInformativeReads);

                // now that we have the SNP and indel GLs, we take the one with the least confidence,
                // as this is the most conservative estimate of our certainty that we are hom-ref.
                // For example, if the SNP PLs are 0,10,100 and the indel PLs are 0,100,1000
                // we are very certain that there's no indel here, but the SNP confidence imply that we are
                // far less confident that the ref base is actually the only thing here.  So we take 0,10,100
                // as our GLs for the site.
                final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);

                gb.GQ((int) (-10 * leastConfidenceGLs.getLog10GQ(GenotypeType.HOM_REF)));
                gb.PL(leastConfidenceGLs.getAsPLs());
                //gb.attribute(INDEL_INFORMATIVE_DEPTH, nIndelInformativeReads);

                vcb.genotypes(gb.make());
                results.add(vcb.make());
//                logger.info("  => VariantContext " + vcb.make());
            }
        }

        return results;
    }

    /**
     * Get the GenotypeLikelihoods with the least strong corresponding GQ value
     * @param gl1 first to consider (cannot be null)
     * @param gl2 second to consider (cannot be null)
     * @return gl1 or gl2, whichever has the worst GQ
     */
    protected final GenotypeLikelihoods getGLwithWorstGQ(final GenotypeLikelihoods gl1, final GenotypeLikelihoods gl2) {
        return gl1.getLog10GQ(GenotypeType.HOM_REF) > gl2.getLog10GQ(GenotypeType.HOM_REF) ? gl1 : gl2;
    }

    /**
     * Get indel PLs corresponding to seeing N nIndelInformativeReads at this site
     *
     * @param nInformativeReads the number of reads that inform us about being ref without an indel at this site
     * @return non-null GenotypeLikelihoods given N
     */
    protected final GenotypeLikelihoods getIndelPLs(final int nInformativeReads) {
        return indelPLCache[nInformativeReads > MAX_N_INDEL_INFORMATIVE_READS ? MAX_N_INDEL_INFORMATIVE_READS : nInformativeReads];
    }

    protected static final int MAX_N_INDEL_INFORMATIVE_READS = 40; // more than this is overkill because GQs are capped at 99 anyway
    private static final GenotypeLikelihoods[] indelPLCache = new GenotypeLikelihoods[MAX_N_INDEL_INFORMATIVE_READS + 1];
    private static final double INDEL_ERROR_RATE = -4.5; // 10^-4.5 indel errors per bp

    private void initializeIndelPLCache() {
        for( int nInformativeReads = 0; nInformativeReads <= MAX_N_INDEL_INFORMATIVE_READS; nInformativeReads++ ) {
            final double homRef = 0.0;
            final double het    = MathUtils.LOG_ONE_HALF * nInformativeReads;
            final double homVar = INDEL_ERROR_RATE * nInformativeReads;
            indelPLCache[nInformativeReads] = GenotypeLikelihoods.fromLog10Likelihoods(new double[]{homRef, het, homVar});
        }
    }

    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
     * @return a RefVsAnyResult genotype call
     */
    public RefVsAnyResult calcGenotypeLikelihoodsOfRefVsAny(final ReadBackedPileup pileup, final byte refBase, final byte minBaseQual, final MathUtils.RunningAverage hqSoftClips) {
        final RefVsAnyResult result = new RefVsAnyResult();

        for( final PileupElement p : pileup ) {
            final byte qual = (p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual());
            if( p.isDeletion() || qual > minBaseQual ) {
                int AA = 0; final int AB = 1; int BB = 2;
                if( p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart() || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip() ) {
                    AA = 2;
                    BB = 0;
                    if( hqSoftClips != null && p.isNextToSoftClip() ) {
                        hqSoftClips.add(AlignmentUtils.calcNumHighQualitySoftClips(p.getRead(), (byte) 28));
                    }
                    result.AD_Ref_Any[1]++;
                } else {
                    result.AD_Ref_Any[0]++;
                }
                result.genotypeLikelihoods[AA] += p.getRepresentativeCount() * QualityUtils.qualToProbLog10(qual);
                result.genotypeLikelihoods[AB] += p.getRepresentativeCount() * MathUtils.approximateLog10SumLog10( QualityUtils.qualToProbLog10(qual) + MathUtils.LOG_ONE_HALF, QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG_ONE_THIRD + MathUtils.LOG_ONE_HALF );
                result.genotypeLikelihoods[BB] += p.getRepresentativeCount() * QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG_ONE_THIRD;
            }
        }

        return result;
    }

    /**
     * Get a list of pileups that span the entire active region span, in order, one for each position
     */
    private List<ReadBackedPileup> getPileupsOverReference(final Haplotype refHaplotype,
                                                           final Collection<Haplotype> calledHaplotypes,
                                                           final GenomeLoc paddedReferenceLoc,
                                                           final ActiveRegion activeRegion,
                                                           final GenomeLoc activeRegionSpan,
                                                           final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap) {

        if ( refHaplotype == null ) throw new IllegalArgumentException("refHaplotype cannot be null");
        if ( calledHaplotypes == null ) throw new IllegalArgumentException("calledHaplotypes cannot be null");
        if ( !calledHaplotypes.contains(refHaplotype)) throw new IllegalArgumentException("calledHaplotypes must contain the refHaplotype");
        if ( paddedReferenceLoc == null ) throw new IllegalArgumentException("paddedReferenceLoc cannot be null");
        if ( activeRegion == null ) throw new IllegalArgumentException("activeRegion cannot be null");
        if ( stratifiedReadMap == null ) throw new IllegalArgumentException("stratifiedReadMap cannot be null");
        if ( stratifiedReadMap.size() != 1 ) throw new IllegalArgumentException("stratifiedReadMap must contain exactly one sample but it contained " + stratifiedReadMap.size());

        List<GATKSAMRecord> realignedReads;

        if( calledHaplotypes.size() == 1 ) { // only contains ref haplotype so an optimization is to just trust the alignments to the reference haplotype as provided by the aligner
            realignedReads = activeRegion.getReads();
        } else {
            final ReadDestination.ToList realignedReadsDest = new ReadDestination.ToList(header, "FOO");
            final HaplotypeBAMWriter writer = HaplotypeBAMWriter.create(HaplotypeBAMWriter.Type.CALLED_HAPLOTYPES, realignedReadsDest);
            writer.setWriteHaplotypesAsWell(false); // don't write out reads for the haplotypes, as we only want the realigned reads themselves
            writer.setOnlyRealignInformativeReads(true);
            writer.writeReadsAlignedToHaplotypes(calledHaplotypes, paddedReferenceLoc, stratifiedReadMap);
            realignedReads = ReadUtils.sortReadsByCoordinate(realignedReadsDest.getReads());
        }

        if ( debuggingWriter != null )
            for ( final GATKSAMRecord read : realignedReads )
                debuggingWriter.addAlignment(read);

        final LocusIteratorByState libs = new LocusIteratorByState(realignedReads.iterator(), LocusIteratorByState.NO_DOWNSAMPLING,
                true, genomeLocParser, samples, false);

        final List<ReadBackedPileup> pileups = new LinkedList<>();
        final int startPos = activeRegionSpan.getStart();
        AlignmentContext next = libs.advanceToLocus(startPos, true);
        for ( int curPos = startPos; curPos <= activeRegionSpan.getStop(); curPos++ ) {
            if ( next != null && next.getLocation().getStart() == curPos ) {
                pileups.add(next.getBasePileup());
                next = libs.hasNext() ? libs.next() : null;
            } else {
                // no data, so we create empty pileups
                pileups.add(new ReadBackedPileupImpl(genomeLocParser.createGenomeLoc(activeRegionSpan.getContig(), curPos)));
            }
        }

        return pileups;
    }

    /**
     * Return the rightmost variant context in maybeOverlapping that overlaps curPos
     *
     * @param curPos non-null genome loc
     * @param maybeOverlapping a collection of variant contexts that might overlap curPos
     * @return a VariantContext, or null if none overlaps
     */
    protected final VariantContext getOverlappingVariantContext(final GenomeLoc curPos, final Collection<VariantContext> maybeOverlapping) {
        VariantContext overlaps = null;
        for ( final VariantContext vc : maybeOverlapping ) {
            if ( genomeLocParser.createGenomeLoc(vc).overlapsP(curPos) ) {
                if ( overlaps == null || vc.getStart() > overlaps.getStart() ) {
                    overlaps = vc;
                }
            }
        }
        return overlaps;
    }

    /**
     * Compute the sum of mismatching base qualities for readBases aligned to refBases at readStart / refStart
     * assuming no insertions or deletions in the read w.r.t. the reference
     *
     * @param readBases non-null bases of the read
     * @param readQuals non-null quals of the read
     * @param readStart the starting position of the read (i.e., that aligns it to a position in the reference)
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @param maxSum if the sum goes over this value, return immediately
     * @return the sum of quality scores for readBases that mismatch their corresponding ref bases
     */
    protected final int sumMismatchingQualities(final byte[] readBases,
                                                final byte[] readQuals,
                                                final int readStart,
                                                final byte[] refBases,
                                                final int refStart,
                                                final int maxSum) {
        final int n = Math.min(readBases.length - readStart, refBases.length - refStart);
        int sum = 0;

        for ( int i = 0; i < n; i++ ) {
            final byte readBase = readBases[readStart + i];
            final byte refBase  = refBases[refStart + i];
            if ( readBase != refBase ) {
                sum += readQuals[readStart + i];
                if ( sum > maxSum ) // abort early
                    return sum;
            }
        }

        return sum;
    }

    /**
     * Compute whether a read is informative to eliminate an indel of size <= maxIndelSize segregating at readStart/refStart
     *
     * @param readBases non-null bases of the read
     * @param readQuals non-null quals of the read
     * @param readStart the starting position of the read (i.e., that aligns it to a position in the reference)
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @param maxIndelSize the max indel size to consider for the read to be informative
     * @return true if read can eliminate the possibility that there's an indel of size <= maxIndelSize segregating at refStart
     */
    protected boolean isReadInformativeAboutIndelsOfSize(final byte[] readBases,
                                                         final byte[] readQuals,
                                                         final int readStart,
                                                         final byte[] refBases,
                                                         final int refStart,
                                                         final int maxIndelSize) {
        // fast exit when n bases left < maxIndelSize
        if( readBases.length - readStart < maxIndelSize || refBases.length - refStart < maxIndelSize ) {
            return false;
        }

        final int baselineMMSum = sumMismatchingQualities(readBases, readQuals, readStart, refBases, refStart, Integer.MAX_VALUE);

        // consider each indel size up to max in term, checking if an indel that deletes either the ref bases (deletion
        // or read bases (insertion) would fit as well as the origin baseline sum of mismatching quality scores
        for ( int indelSize = 1; indelSize <= maxIndelSize; indelSize++ ) {
            for ( final boolean checkInsertion : Arrays.asList(true, false) ) {
                final int readI, refI;
                if ( checkInsertion ) {
                    readI = readStart + indelSize;
                    refI = refStart;
                } else {
                    readI = readStart;
                    refI = refStart + indelSize;
                }

                final int score = sumMismatchingQualities(readBases, readQuals, readI, refBases, refI, baselineMMSum);
                if ( score <= baselineMMSum )
                    return false;
            }
        }

        return true;
    }

    /**
     * Calculate the number of indel informative reads at pileup
     *
     * @param pileup a pileup
     * @param pileupOffsetIntoRef the position of the pileup in the reference
     * @param ref the ref bases
     * @param maxIndelSize maximum indel size to consider in the informativeness calculation
     * @return an integer >= 0
     */
    protected final int calcNIndelInformativeReads(final ReadBackedPileup pileup, final int pileupOffsetIntoRef, final byte[] ref, final int maxIndelSize) {
        int nInformative = 0;
        for ( final PileupElement p : pileup ) {
            final GATKSAMRecord read = p.getRead();
            final int offset = p.getOffset();

            // doesn't count as evidence
            if ( p.isBeforeDeletionStart() || p.isBeforeInsertion() || p.isDeletion() )
                continue;

            // todo -- this code really should handle CIGARs directly instead of relying on the above tests
            if ( isReadInformativeAboutIndelsOfSize(read.getReadBases(), read.getBaseQualities(), offset, ref, pileupOffsetIntoRef, maxIndelSize) ) {
                nInformative++;
                if( nInformative > MAX_N_INDEL_INFORMATIVE_READS ) {
                    return MAX_N_INDEL_INFORMATIVE_READS;
                }
            }
        }
        return nInformative;
    }

    /**
     * Create a reference haplotype for an active region
     *
     * @param activeRegion the active region
     * @param refBases the ref bases
     * @param paddedReferenceLoc the location spanning of the refBases -- can be longer than activeRegion.getLocation()
     * @return a reference haplotype
     */
    public static Haplotype createReferenceHaplotype(final ActiveRegion activeRegion, final byte[] refBases, final GenomeLoc paddedReferenceLoc) {
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final int alignmentStart = activeRegion.getExtendedLoc().getStart() - paddedReferenceLoc.getStart();
        if ( alignmentStart < 0 ) throw new IllegalStateException("Bad alignment start in createReferenceHaplotype " + alignmentStart);
        refHaplotype.setAlignmentStartHapwrtRef(alignmentStart);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);
        return refHaplotype;
    }
}
