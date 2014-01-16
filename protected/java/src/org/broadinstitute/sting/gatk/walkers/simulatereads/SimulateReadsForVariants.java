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

package org.broadinstitute.sting.gatk.walkers.simulatereads;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;

import java.util.*;

/**
 * Generates simulated reads for variants
 *
 * <p>Given a set of variants, this tool will generate simulated reads that support the input variants.</p>
 *
 * <h3>Caveats</h3>
 * <p>For practical reasons, only bi-allelic variants that are not too close to the ends of contigs (< 1/2 read length) are supported; all others will simply be ignored.</p>
 *
 * <h3>Input</h3>
 * <p>A VCF file containing variants.</p>
 *
 * <h3>Output</h3>
 * <p>A BAM file containing simulated sequence reads that support the input variants, with the requested error rate and coverage depth.</p>
 *
 * <h3>Example</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -T SimulateReadsForVariants \
 *   -R reference.fasta \
 *   -V input_variants.vcf \
 *   -o simulated_reads.bam \
 *   --readDepth 50 \
 *   --errorRate 25
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class}, gotoDev = HelpConstants.EB)

@Reference(window=@Window(start=-200,stop=200))
public class SimulateReadsForVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();
    /**
     * The simulated reads will be written to a BAM file.
     */
    @Output(doc="Reads corresponding to variants", required=true)
    protected StingSAMFileWriter readWriter;
    /**
     * Use this argument to set the desired target read depth. See the readSamplingMode argument for options that
     * determine whether coverage distribution will be exactly this value or an approximation.
     */
    @Argument(fullName="readDepth", shortName="DP", doc="Read depth to generate", required=false, minValue = 0, minRecommendedValue = 1, maxRecommendedValue = 1000, maxValue = Integer.MAX_VALUE)
    public int readDepth = 20;
    /**
     * Errors will be generated at this rate in the simulated reads. Base qualities are therefore also assigned this value.
     */
    @Argument(fullName="errorRate", shortName="ER", doc="Base error rate (Phred-scaled)", required=false, minValue = 0, maxValue = Integer.MAX_VALUE)
    public int phredErrorRate = 20;
    /**
     * All simulated reads will be exactly this length.
     */
    @Argument(fullName="readLength", shortName="RL", doc="Read lengths (bp)", required=false, minValue = 1, maxValue = Integer.MAX_VALUE)
    public int readLength = 101;
    /**
     * The corresponding platform identifier will be specified in the simulated read group PL tag. This setting does not
     * affect the properties of the simulated reads.
     */
    @Advanced
    @Argument(fullName="rgPlatform", shortName="RGPL", doc="Sequencing platform", required=false)
    public NGSPlatform rgPlatform = NGSPlatform.ILLUMINA;
    /**
     * This determines how read sampling is achieved, and affects the coverage distribution of simulated reads.
     * CONSTANT sampling will produce uniform depth at all positions, while POISSON sampling will produce a
     * distribution of coverages around the requested value.
     */
    @Advanced
    @Argument(fullName="readSamplingMode", shortName="RSM", doc="Sampling mode", required=false)
    public ReadSamplingMode samplingMode = ReadSamplingMode.CONSTANT;
    public enum ReadSamplingMode { CONSTANT, POISSON };

    @Hidden
    @Argument(fullName = "no_pg_tag", shortName = "npt", doc ="Discard program tags, for integration tests", required=false)
    public boolean NO_PG_TAG = false;

    @Hidden
    @Argument(fullName="verbose", shortName="verbose", doc="Verbose", required=false)
    public boolean verbose = false;

    public static final String PROGRAM_RECORD_NAME = "GATK SimulateReadsForVariants";

    // variables used to store state
    private long readNameCounter = 1;
    private int halfReadLength;
    private double errorRate;
    private byte[] readQuals;
    private SAMFileHeader header = null;

    // randomness related variables
    private static final long RANDOM_SEED = 1252863495;
    private static final Random ran = GenomeAnalysisEngine.getRandomGenerator();
    private Poisson poissonRandom = null;

    // samples and read groups
    private final Map<String, SAMReadGroupRecord> sample2RG = new HashMap<String, SAMReadGroupRecord>();

    private SAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }

    private SAMReadGroupRecord createRG(String name) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(name);
        rg.setPlatform(rgPlatform.getDefaultPlatform());
        rg.setSample(name);
        return rg;
    }

    // class to store the bases, offset, and representative CIGAR of a haplotype
    private static class ArtificialHaplotype {
        public final byte[] bases;
        public final int offset;
        public final String cigar;

        public ArtificialHaplotype(final byte[] bases, final int offset, final String cigar) {
            this.bases = bases;
            this.offset = offset;
            this.cigar = cigar;
        }
    }

    @Override
    public void initialize() {

        // initialize sample -> read group map
        final List<SAMReadGroupRecord> sampleRGs = new ArrayList<SAMReadGroupRecord>();
        for ( final String sample : SampleUtils.getUniqueSamplesFromRods(getToolkit(), Arrays.asList(variantCollection.variants.getName())) ) {
            final SAMReadGroupRecord rg = createRG(sample);
            sampleRGs.add(rg);
            sample2RG.put(sample, rg);
        }

        // initialize BAM headers
        header = new SAMFileHeader();
        header.setSequenceDictionary(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setReadGroups(sampleRGs);

        final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
        if ( !NO_PG_TAG ) {
            final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
            programRecord.setProgramVersion(headerInfo.getString("org.broadinstitute.sting.gatk.version"));
            programRecord.setCommandLine(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
        }
        header.setProgramRecords(Arrays.asList(programRecord));

        readWriter.setPresorted(false);
        readWriter.writeHeader(header);

        halfReadLength = readLength / 2;
        errorRate = QualityUtils.qualToErrorProb((byte)phredErrorRate);
        readQuals = new byte[readLength];
        Arrays.fill(readQuals, (byte)phredErrorRate);
        if ( samplingMode == ReadSamplingMode.POISSON )
           poissonRandom = new Poisson(readDepth, new MersenneTwister((int)RANDOM_SEED));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        if ( ref.getLocus().getStart() < readLength || ! BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        final VariantContext vc = tracker.getFirstValue(variantCollection.variants, context.getLocation());
        if ( vc == null || !vc.isBiallelic() )
            return 0;

        if ( verbose ) logger.info(String.format("Generating reads for %s", vc));

        generateReadsForVariant(vc, ref);

        return 1;
    }

    /**
     * Contstructs an artifical haplotype given an allele and original reference context
     *
     * @param allele     the allele to model (can be reference)
     * @param refLength  the length of the reference allele
     * @param ref        the original reference context
     * @return a non-null ArtificialHaplotype
     */
    private ArtificialHaplotype constructHaplotype(final Allele allele, final int refLength, final ReferenceContext ref) {

        final byte[] haplotype = new byte[readLength];

        final int alleleLength = allele.getBases().length;
        final int halfAlleleLength = (alleleLength + 1) / 2;

        // this is how far back to move from the event to start copying bases
        final int offset = halfReadLength - halfAlleleLength;

        // copy bases before the event
        final int locusPosOnRefContext = (int)(ref.getLocus().getStart() - ref.getWindow().getStart());
        int posOnRefContext = locusPosOnRefContext - offset;
        System.arraycopy(ref.getBases(), posOnRefContext, haplotype, 0, offset);
        int copiedCount = offset;

        // copy the event bases
        System.arraycopy(allele.getBases(), 0, haplotype, copiedCount, alleleLength);
        copiedCount += alleleLength;

        // copy bases after the event
        posOnRefContext = locusPosOnRefContext + refLength;
        final int remainder = readLength - copiedCount;
        System.arraycopy(ref.getBases(), posOnRefContext, haplotype, copiedCount, remainder);

        final String cigar;
        if ( refLength == alleleLength )
            cigar = readLength + "M";
        else
            cigar = (offset + 1) + "M" + Math.abs(refLength - alleleLength) + (refLength > alleleLength ? "D" : "I") + remainder + "M";

        return new ArtificialHaplotype(haplotype, offset, cigar);
    }

    /**
     * Generates the artificial reads for a given variant
     *
     * @param vc         the (bi-allelic) variant context for which to generate artificial reads
     * @param ref        the original reference context
     */
    private void generateReadsForVariant(final VariantContext vc, final ReferenceContext ref) {

        final int refLength = vc.getReference().getBases().length;
        final ArtificialHaplotype refHap = constructHaplotype(vc.getReference(), refLength, ref);
        final ArtificialHaplotype altHap = constructHaplotype(vc.getAlternateAllele(0), refLength, ref);

        int gi = 0;
        for ( final Genotype g : vc.getGenotypes() ) {
            final int myDepth = sampleDepth();
            for ( int d = 0; d < myDepth; d++ ) {

                final ArtificialHaplotype haplotype = chooseRefHaplotype(g) ? refHap : altHap;
                final byte[] readBases = Arrays.copyOf(haplotype.bases, readLength);

                addMachineErrors(readBases, errorRate);
                writeRead(readBases, vc.getChr(), vc.getStart() - haplotype.offset, haplotype.cigar, g.getSampleName(), gi++ % 2 == 0);
            }
        }
    }

    /**
     * Decides whether or not to choose the reference haplotype, depending on the given genotype
     *
     * @param g  the genotype of the given sample
     * @return true if one should use the reference haplotype, false otherwise
     */
    private boolean chooseRefHaplotype(final Genotype g) {
        final double refP;
        if ( g.isHomRef() )     refP = 1;
        else if ( g.isHet() )   refP = 0.5;
        else                    refP = 0.0;

        return ran.nextDouble() < refP;
    }

    /**
     * Generates the artificial read depth
     *
     * @return a non-negative int
     */
    private int sampleDepth() {
        switch ( samplingMode ) {
            case CONSTANT: return readDepth;
            case POISSON: return poissonRandom.nextInt();
            default:
                throw new IllegalStateException("Unexpected DepthSamplingType " + samplingMode);
        }
    }

    /**
     * Creates and writes an artificial read given the appropriate data
     *
     * @param readBases   the bases
     * @param contig      the contig
     * @param start       the read start
     * @param cigar       the cigar string
     * @param sample      the sample name (used to get the right read group)
     * @param isNegStrand should this read be on the negative strand?
     */
    private void writeRead(final byte[] readBases, final String contig, final int start,
			   final String cigar, final String sample, final boolean isNegStrand) {
        final GATKSAMRecord read = new GATKSAMRecord(header);
        read.setBaseQualities(readQuals);
        read.setReadBases(readBases);
        read.setReadName("" + readNameCounter++);
        read.setCigarString(cigar);
        read.setReadPairedFlag(false);
        read.setAlignmentStart(start);
        read.setMappingQuality(60);
        read.setReferenceName(contig);
        read.setReadNegativeStrandFlag(isNegStrand);
        read.setAttribute("RG", sampleRG(sample).getReadGroupId());

        readWriter.addAlignment(read);
    }

    /**
     * Adds machine errors at the appropriate rate to the provided read bases
     *
     * @param readBases   the read bases
     * @param errorRate   the rate at which to produce errors
     */
    private void addMachineErrors(final byte[] readBases, final double errorRate) {
        for ( int i = 0; i < readBases.length; i++ ) {
            final double r = ran.nextDouble();
            if ( r < errorRate ) {
                byte errorBase = BaseUtils.baseIndexToSimpleBase(BaseUtils.getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex(readBases[i])));
                if ( errorBase == readBases[i] ) throw new IllegalStateException("Read and error bases are the same");
                readBases[i] = errorBase;
            }
        }
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }
}
