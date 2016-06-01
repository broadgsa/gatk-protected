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

import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.commandline.Output;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Set of arguments related to the {@link ReadThreadingAssembler}
 *
 * @author Kristian Cibulskis &lt;kcibul@broadinstitute.org&gt;
 */
public class ReadThreadingAssemblerArgumentCollection {

    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the read threading assembler
    // -----------------------------------------------------------------------------------------------

    /**
     * Multiple kmer sizes can be specified, using e.g. `-kmerSize 10 -kmerSize 25`.
     */
    @Advanced
    @Argument(fullName="kmerSize", shortName="kmerSize", doc="Kmer size to use in the read threading assembler", required = false)
    public List<Integer> kmerSizes = Arrays.asList(10, 25);

    /**
     * When graph cycles are detected, the normal behavior is to increase kmer sizes iteratively until the cycles are
     * resolved. Disabling this behavior may cause the program to give up on assembling the ActiveRegion.
     */
    @Advanced
    @Argument(fullName="dontIncreaseKmerSizesForCycles", shortName="dontIncreaseKmerSizesForCycles", doc="Disable iterating over kmer sizes when graph cycles are detected", required = false)
    public boolean dontIncreaseKmerSizesForCycles = false;

    /**
     * By default, the program does not allow processing of reference sections that contain non-unique kmers. Disabling
     * this check may cause problems in the assembly graph.
     */
    @Advanced
    @Argument(fullName="allowNonUniqueKmersInRef", shortName="allowNonUniqueKmersInRef", doc="Allow graphs that have non-unique kmers in the reference", required = false)
    public boolean allowNonUniqueKmersInRef = false;

    /**
     * If fewer samples than the specified number pass the minPruning threshold for a given path, that path will be eliminated from the graph.
     */
    @Advanced
    @Argument(fullName="numPruningSamples", shortName="numPruningSamples", doc="Number of samples that must pass the minPruning threshold", required = false)
    public int numPruningSamples = 1;

    /**
     * As of version 3.3, this argument is no longer needed because dangling end recovery is now the default behavior. See GATK 3.3 release notes for more details.
     */
    @Deprecated
    @Argument(fullName="recoverDanglingHeads", shortName="recoverDanglingHeads", doc="This argument is deprecated since version 3.3", required = false)
    public boolean DEPRECATED_RecoverDanglingHeads = false;

    /**
     * By default, the read threading assembler will attempt to recover dangling heads and tails. See the `minDanglingBranchLength` argument documentation for more details.
     */
    @Hidden
    @Argument(fullName="doNotRecoverDanglingBranches", shortName="doNotRecoverDanglingBranches", doc="Disable dangling head and tail recovery", required = false)
    public boolean doNotRecoverDanglingBranches = false;

    /**
     * When constructing the assembly graph we are often left with "dangling" branches.  The assembly engine attempts to rescue these branches
     * by merging them back into the main graph.  This argument describes the minimum length of a dangling branch needed for the engine to
     * try to rescue it.  A smaller number here will lead to higher sensitivity to real variation but also to a higher number of false positives.
     */
    @Advanced
    @Argument(fullName="minDanglingBranchLength", shortName="minDanglingBranchLength", doc="Minimum length of a dangling branch to attempt recovery", required = false)
    public int minDanglingBranchLength = 4;

    /**
     * This argument is specifically intended for 1000G consensus analysis mode. Setting this flag will inject all
     * provided alleles to the assembly graph but will not forcibly genotype all of them.
     */
    @Advanced
    @Argument(fullName="consensus", shortName="consensus", doc="1000G consensus mode", required = false)
    public boolean consensusMode = false;

    /**
     * The assembly graph can be quite complex, and could imply a very large number of possible haplotypes.  Each haplotype
     * considered requires N PairHMM evaluations if there are N reads across all samples.  In order to control the
     * run of the haplotype caller we only take maxNumHaplotypesInPopulation paths from the graph, in order of their
     * weights, no matter how many paths are possible to generate from the graph.  Putting this number too low
     * will result in dropping true variation because paths that include the real variant are not even considered.
     * You can consider increasing this number when calling organisms with high heterozygosity.
     */
    @Advanced
    @Argument(fullName="maxNumHaplotypesInPopulation", shortName="maxNumHaplotypesInPopulation", doc="Maximum number of haplotypes to consider for your population", required = false)
    public int maxNumHaplotypesInPopulation = 128;

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName="errorCorrectKmers", shortName="errorCorrectKmers", doc = "Use an exploratory algorithm to error correct the kmers used during assembly", required=false)
    public boolean errorCorrectKmers = false;

    /**
     * Paths with fewer supporting kmers than the specified threshold will be pruned from the graph.
     *
     * Be aware that this argument can dramatically affect the results of variant calling and should only be used with great caution.
     * Using a prune factor of 1 (or below) will prevent any pruning from the graph, which is generally not ideal; it can make the
     * calling much slower and even less accurate (because it can prevent effective merging of "tails" in the graph).  Higher values
     * tend to make the calling much faster, but also lowers the sensitivity of the results (because it ultimately requires higher
     * depth to produce calls).
     */
    @Advanced
    @Argument(fullName="minPruning", shortName="minPruning", doc = "Minimum support to not prune paths in the graph", required = false)
    public int MIN_PRUNE_FACTOR = 2;

    @Hidden
    @Argument(fullName="debugGraphTransformations", shortName="debugGraphTransformations", doc="Write DOT formatted graph files out of the assembler for only this graph size", required = false)
    public boolean debugGraphTransformations = false;

    @Hidden
    @Argument(fullName="allowCyclesInKmerGraphToGeneratePaths", shortName="allowCyclesInKmerGraphToGeneratePaths", doc="Allow cycles in the kmer graphs to generate paths with multiple copies of the path sequenece rather than just the shortest paths", required = false)
    public boolean allowCyclesInKmerGraphToGeneratePaths = false;

    /**
     * This argument is meant for debugging and is not immediately useful for normal analysis use.
     */
    @Output(fullName="graphOutput", shortName="graph", doc="Write debug assembly graph information to this file", required = false, defaultToStdout = false)
    public PrintStream graphWriter = null;

    //---------------------------------------------------------------------------------------------------------------
    //
    // Read Error Corrector Related Parameters
    //
    // ---------------------------------------------------------------------------------------------------------------

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName="kmerLengthForReadErrorCorrection", shortName="kmerLengthForReadErrorCorrection", doc = "Use an exploratory algorithm to error correct the kmers used during assembly", required=false)
    public int kmerLengthForReadErrorCorrection = 25;

    @Hidden
    @Argument(fullName="minObservationsForKmerToBeSolid", shortName="minObservationsForKmerToBeSolid", doc = "A k-mer must be seen at least these times for it considered to be solid", required=false)
    public int minObservationsForKmerToBeSolid = 20;







}
