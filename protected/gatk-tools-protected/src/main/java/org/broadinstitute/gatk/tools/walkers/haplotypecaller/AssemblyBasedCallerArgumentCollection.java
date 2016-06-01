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

import org.broadinstitute.gatk.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.haplotypeBAMWriter.HaplotypeBAMWriter;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;

/**
 * Set of arguments for Assembly Based Callers
 *
 * @author Kristian Cibulskis &lt;kcibul@broadinstitute.org&gt;
 */
public class AssemblyBasedCallerArgumentCollection extends StandardCallerArgumentCollection {

    @Advanced
    @Argument(fullName="debug", shortName="debug", doc="Print out very verbose debug information about each triggering active region", required = false)
    public boolean DEBUG;

    @Advanced
    @Argument(fullName="useFilteredReadsForAnnotations", shortName="useFilteredReadsForAnnotations", doc = "Use the contamination-filtered read maps for the purposes of annotating variants", required=false)
    public boolean USE_FILTERED_READ_MAP_FOR_ANNOTATIONS = false;

    /**
     * The reference confidence mode makes it possible to emit variant calls in GVCF format, which includes either a per-base
     * pair (BP_RESOLUTION) or a summarized (GVCF) confidence estimate for each position being strictly homozygous-reference.
     * See http://www.broadinstitute.org/gatk/guide/article?id=2940 for more details of how this works.
     * Note that if you use <code>-ERC</code> to emit a <code>GVCF</code> or <code>BP_RESOLUTION</code> output, you either
     * need to give the output file the extension <code>.g.vcf</code> or set the parameters <code>-variant_index_type LINEAR</code>
     * and <code>-variant_index_parameter 128000</code> (with those exact values!). This has to do with index compression.
     */
    @Advanced
    @Argument(fullName="emitRefConfidence", shortName="ERC", doc="Mode for emitting reference confidence scores", required = false)
    protected ReferenceConfidenceMode emitReferenceConfidence = ReferenceConfidenceMode.NONE;

    @Override
    public AssemblyBasedCallerArgumentCollection clone() {
        return (AssemblyBasedCallerArgumentCollection) super.clone();
    }

    /**
     * The assembled haplotypes and locally realigned reads will be written as BAM to this file if requested.  This is
     * intended to be used only for troubleshooting purposes, in specific areas where you want to better understand
     * why the caller is making specific calls. Turning on this mode may result in serious performance cost for the
     * caller, so we do NOT recommend using this argument systematically as it will significantly increase runtime.
     *
     * The candidate haplotypes (called or all, depending on mode) are emitted as single reads covering the entire
     * active region, coming from sample "HC" and a special read group called "ArtificialHaplotype". This will increase
     * the pileup depth compared to what would be expected from the reads only, especially in complex regions.
     *
     * The reads are written out containing an "HC" tag (integer) that encodes which haplotype each read best matches
     * according to the haplotype caller's likelihood calculation.  The use of this tag is primarily intended
     * to allow good coloring of reads in IGV.  Simply go to "Color Alignments By > Tag" and enter "HC" to more
     * easily see which reads go with these haplotype. You can also tell IGV to group reads by sample, which will
     * separate the potential haplotypes from the reads. These features are illustrated in
     * <a href="https://www.dropbox.com/s/xvy7sbxpf13x5bp/haplotypecaller%20bamout%20for%20docs.png">this screenshot</a>.
     *
     * Note that only reads that are actually informative about the haplotypes are emitted with the HC tag.
     * By informative we mean that there's a meaningful difference in the likelihood of the read coming from one
     * haplotype compared to the next best haplotype. When coloring reads by HC tag in IGV, uninformative reads will
     * remain grey.
     *
     * Note also that not every input read is emitted to the bam in this mode. To include all trimmed, downsampled,
     * filtered and uninformative reads, add the <code>--emitDroppedReads</code> argument.
     *
     * If multiple BAMs are passed as input to the tool (as is common for MuTect2), then they will be combined in the
     * <code>-bamout</code> output and tagged with the appropriate sample names.
     */
    @Advanced
    @Output(fullName="bamOutput", shortName="bamout", doc="File to which assembled haplotypes should be written", required = false, defaultToStdout = false)
    public GATKSAMFileWriter bamWriter = null;

    /**
     * The type of <code>-bamout</code> output we want to see. This determines whether HC will write out all of the haplotypes it
     * considered (top 128 max) or just the ones that were selected as alleles and assigned to samples.
     */
    @Advanced
    @Argument(fullName="bamWriterType", shortName="bamWriterType", doc="Which haplotypes should be written to the BAM", required = false)
    public HaplotypeBAMWriter.Type bamWriterType = HaplotypeBAMWriter.Type.CALLED_HAPLOTYPES;

    /**
     * Determines whether dropped reads will be tracked and emitted when <code>-bamout</code> is specified. Use this in combination
     * with a specific interval of interest to avoid accumulating a large number of reads in the <code>-bamout</code> file.
     */
    @Advanced
    @Argument(fullName="emitDroppedReads", shortName="edr", doc="Emit reads that are dropped for filtering, trimming, realignment failure", required = false)
    public boolean emitDroppedReads = false;

    /**
     * If set, certain "early exit" optimizations in HaplotypeCaller, which aim to save compute and time by skipping
     * calculations if an ActiveRegion is determined to contain no variants, will be disabled. This is most likely to be
     * useful if you're using the <code>-bamout</code> argument to examine the placement of reads following reassembly
     * and are interested in seeing the mapping of reads in regions with no variations. Setting the <code>-forceActive</code>
     * and <code>-dontTrimActiveRegions</code> flags may also be helpful.
     */
    @Advanced
    @Argument(fullName = "disableOptimizations", shortName="disableOptimizations", doc="Don't skip calculations in ActiveRegions with no variants",
            required = false)
    public boolean disableOptimizations = false;

}
