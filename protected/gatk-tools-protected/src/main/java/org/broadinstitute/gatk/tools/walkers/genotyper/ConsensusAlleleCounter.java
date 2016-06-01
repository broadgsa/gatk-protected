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

package org.broadinstitute.gatk.tools.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import htsjdk.variant.variantcontext.*;

import java.util.*;

/**
 * Code for determining which indels are segregating among the samples.
 *
 * This code is just a refactor of the original code from Guillermo in the UG.
 *
 * @author Mark DePristo
 * @since 3/26/12
 */
public class ConsensusAlleleCounter {
    final protected static Logger logger = Logger.getLogger(ConsensusAlleleCounter.class);
    private final int minIndelCountForGenotyping;
    private final boolean doMultiAllelicCalls;
    private final double minFractionInOneSample;

    public ConsensusAlleleCounter(final boolean doMultiAllelicCalls,
                                  final int minIndelCountForGenotyping,
                                  final double minFractionInOneSample) {
        this.minIndelCountForGenotyping = minIndelCountForGenotyping;
        this.doMultiAllelicCalls = doMultiAllelicCalls;
        this.minFractionInOneSample = minFractionInOneSample;
    }

    /**
     * Returns a list of Alleles at this locus that may be segregating
     * 
     * @param ref
     * @param contexts
     * @param contextType
     * @return
     */
    public List<Allele> computeConsensusAlleles(ReferenceContext ref,
                                                Map<String, AlignmentContext> contexts,
                                                AlignmentContextUtils.ReadOrientation contextType) {
        final Map<String, Integer> consensusIndelStrings = countConsensusAlleles(ref, contexts, contextType);
        return consensusCountsToAlleles(ref, consensusIndelStrings);
    }

    //
    // TODO -- WARNING DOESN'T WORK WITH REDUCED READS
    //
    private Map<String, Integer> countConsensusAlleles(ReferenceContext ref,
                                                       Map<String, AlignmentContext> contexts,
                                                       AlignmentContextUtils.ReadOrientation contextType) {
        final GenomeLoc loc = ref.getLocus();
        HashMap<String, Integer> consensusIndelStrings = new HashMap<String, Integer>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for ( Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            final AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedPileup indelPileup = context.getBasePileup();
            insCount += indelPileup.getNumberOfInsertionsAfterThisElement();
            delCount += indelPileup.getNumberOfDeletionsAfterThisElement();
        }

        if ( insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping )
            return Collections.emptyMap();

        for (Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            // todo -- warning, can be duplicating expensive partition here
            AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadBackedPileup indelPileup = context.getBasePileup();

            final int nIndelReads = indelPileup.getNumberOfInsertionsAfterThisElement() + indelPileup.getNumberOfDeletionsAfterThisElement();
            final int nReadsOverall = indelPileup.getNumberOfElements();

            if ( nIndelReads == 0 || (nIndelReads / (1.0 * nReadsOverall)) < minFractionInOneSample ) {
                continue;
            }

            for (PileupElement p : indelPileup) {
                final GATKSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());
                if (read == null)
                    continue;
                if (ReadUtils.is454Read(read)) {
                    continue;
                }

                if ( p.isBeforeInsertion() ) {
                    final String insertionBases = p.getBasesOfImmediatelyFollowingInsertion();
                    // edge case: ignore a deletion immediately preceding an insertion as p.getBasesOfImmediatelyFollowingInsertion() returns null [EB]
                    if ( insertionBases == null )
                        continue;

                    boolean foundKey = false;
                    // copy of hashmap into temp arrayList
                    ArrayList<Pair<String,Integer>> cList = new ArrayList<Pair<String,Integer>>();
                    for (Map.Entry<String, Integer> s : consensusIndelStrings.entrySet()) {
                        cList.add(new Pair<String, Integer>(s.getKey(), s.getValue()));
                    }

                    if (read.getAlignmentEnd() == loc.getStart()) {
                        // first corner condition: a read has an insertion at the end, and we're right at the insertion.
                        // In this case, the read could have any of the inserted bases and we need to build a consensus

                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            // case 1: current insertion is prefix of indel in hash map
                            if (s.startsWith(insertionBases)) {
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (insertionBases.startsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(insertionBases,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(insertionBases,1));

                    }
                    else if (read.getAlignmentStart() == loc.getStart()+1) {
                        // opposite corner condition: read will start at current locus with an insertion
                        for (int k=0; k < cList.size(); k++) {
                            String s = cList.get(k).getFirst();
                            int cnt = cList.get(k).getSecond();
                            if (s.endsWith(insertionBases)) {
                                // case 1: current insertion (indelString) is suffix of indel in hash map (s)
                                cList.set(k,new Pair<String, Integer>(s,cnt+1));
                                foundKey = true;
                            }
                            else if (insertionBases.endsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,new Pair<String, Integer>(insertionBases,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                            cList.add(new Pair<String, Integer>(insertionBases,1));


                    }
                    else {
                        // normal case: insertion somewhere in the middle of a read: add count to arrayList
                        int cnt = consensusIndelStrings.containsKey(insertionBases)? consensusIndelStrings.get(insertionBases):0;
                        cList.add(new Pair<String, Integer>(insertionBases,cnt+1));
                    }

                    // copy back arrayList into hashMap
                    consensusIndelStrings.clear();
                    for (Pair<String,Integer> pair : cList) {
                        consensusIndelStrings.put(pair.getFirst(),pair.getSecond());
                    }

                }
                else if ( p.isBeforeDeletionStart() ) {
                    final String deletionString = String.format("D%d",p.getLengthOfImmediatelyFollowingIndel());
                    int cnt = consensusIndelStrings.containsKey(deletionString)? consensusIndelStrings.get(deletionString):0;
                    consensusIndelStrings.put(deletionString,cnt+1);
                }
            }
        }

        return consensusIndelStrings;
    }

    private List<Allele> consensusCountsToAlleles(final ReferenceContext ref,
                                                  final Map<String, Integer> consensusIndelStrings) {
        final GenomeLoc loc = ref.getLocus();
        final Collection<VariantContext> vcs = new ArrayList<VariantContext>();
        int maxAlleleCnt = 0;
        Allele refAllele, altAllele;

        for (final Map.Entry<String, Integer> elt : consensusIndelStrings.entrySet()) {
            final String s = elt.getKey();
            final int curCnt = elt.getValue();
            int stop = 0;

            // if observed count if above minimum threshold, we will genotype this allele
            if (curCnt < minIndelCountForGenotyping)
                continue;

            if (s.startsWith("D")) {
                // get deletion length
                final int dLen = Integer.valueOf(s.substring(1));
                // get ref bases of accurate deletion
                final int startIdxInReference = 1 + loc.getStart() - ref.getWindow().getStart();
                stop = loc.getStart() + dLen;
                final byte[] refBases = Arrays.copyOfRange(ref.getBases(), startIdxInReference - 1, startIdxInReference + dLen);   // add reference padding

                if (Allele.acceptableAlleleBases(refBases, false)) {
                    refAllele = Allele.create(refBases, true);
                    altAllele = Allele.create(ref.getBase(), false);
                }
                else continue; // don't go on with this allele if refBases are non-standard
            } else {
                // insertion case
                final String insertionBases = (char)ref.getBase() + s;  // add reference padding
                if (Allele.acceptableAlleleBases(insertionBases, false)) { // don't allow N's in insertions
                    refAllele = Allele.create(ref.getBase(), true);
                    altAllele = Allele.create(insertionBases, false);
                    stop = loc.getStart();
                }
                else continue; // go on to next allele if consensus insertion has any non-standard base.
            }


            final VariantContextBuilder builder = new VariantContextBuilder().source("");
            builder.loc(loc.getContig(), loc.getStart(), stop);
            builder.alleles(Arrays.asList(refAllele, altAllele));
            builder.noGenotypes();
            if (doMultiAllelicCalls) {
                vcs.add(builder.make());
                if (vcs.size() >= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
                    break;
            } else if (curCnt > maxAlleleCnt) {
                maxAlleleCnt = curCnt;
                vcs.clear();
                vcs.add(builder.make());
            }
        }

        if (vcs.isEmpty())
            return Collections.emptyList(); // nothing else to do, no alleles passed minimum count criterion

        final VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(vcs, null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, GATKVariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);
        return mergedVC.getAlleles();
    }
}
