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

package org.broadinstitute.gatk.tools.walkers.phasing;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

// Streams in VariantContext objects and streams out VariantContexts produced by merging phased segregating polymorphisms into MNP VariantContexts

class MergeSegregatingAlternateAllelesVCFWriter implements VariantContextWriter {
    private VariantContextWriter innerWriter;

    private GenomeLocParser genomeLocParser;

    private ReferenceSequenceFile referenceFileForMNPmerging;

    private VariantContextMergeRule vcMergeRule;
    private PhasingUtils.AlleleMergeRule alleleMergeRule;

    private String useSingleSample = null;

    private boolean emitOnlyMergedRecords;

    private VCFRecord vcfrWaitingToMerge;
    private List<VCFRecord> filteredVcfrList;

    private int numRecordsAttemptToMerge;
    private int numRecordsSatisfyingMergeRule;
    private int numMergedRecords;
    private AltAlleleStatsForSamples altAlleleStats = null;

    private Logger logger;

    // Should we call innerWriter.close() in close()
    private boolean takeOwnershipOfInner;

    public MergeSegregatingAlternateAllelesVCFWriter(VariantContextWriter innerWriter, GenomeLocParser genomeLocParser, File referenceFile, VariantContextMergeRule vcMergeRule, PhasingUtils.AlleleMergeRule alleleMergeRule, String singleSample, boolean emitOnlyMergedRecords, Logger logger, boolean takeOwnershipOfInner, boolean trackAltAlleleStats) {
        this.innerWriter = innerWriter;
        this.genomeLocParser = genomeLocParser;
        try {
            this.referenceFileForMNPmerging = new CachingIndexedFastaSequenceFile(referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }        
        this.vcMergeRule = vcMergeRule;
        this.alleleMergeRule = alleleMergeRule;
        this.useSingleSample = singleSample;
        this.emitOnlyMergedRecords = emitOnlyMergedRecords;

        this.vcfrWaitingToMerge = null;
        this.filteredVcfrList = new LinkedList<VCFRecord>();
        this.numRecordsSatisfyingMergeRule = 0;
        this.numMergedRecords = 0;

        if (trackAltAlleleStats)
            this.altAlleleStats = new AltAlleleStatsForSamples();

        this.logger = logger;
        this.takeOwnershipOfInner = takeOwnershipOfInner;
    }

    public MergeSegregatingAlternateAllelesVCFWriter(VariantContextWriter innerWriter, GenomeLocParser genomeLocParser, File referenceFile, int maxGenomicDistanceForMNP, Logger logger, boolean takeOwnershipOfInner) {
        this(innerWriter, genomeLocParser, referenceFile, new DistanceMergeRule(maxGenomicDistanceForMNP, genomeLocParser), new SegregatingMNPmergeAllelesRule(), null, false, logger, takeOwnershipOfInner, true); // by default: consider all samples, emit all records, keep track of alt allele statistics
    }

    public void writeHeader(VCFHeader header) {
        if (useSingleSample != null) { // only want to output context for one sample
            Set<String> singSampSet = new TreeSet<String>();
            singSampSet.add(useSingleSample);
            header = new VCFHeader(header.getMetaDataInSortedOrder(), singSampSet);
        }

        innerWriter.writeHeader(header);
    }

    public void close() {
        stopWaitingToMerge();

        if (takeOwnershipOfInner)
            innerWriter.close();
    }

    public void add(VariantContext vc) {
        if (useSingleSample != null) { // only want to output context for one sample
            Genotype sampGt = vc.getGenotype(useSingleSample);
            if (sampGt != null) // TODO: subContextFromGenotypes() does not handle any INFO fields [AB, HaplotypeScore, MQ, etc.].  Note that even SelectVariants.subsetRecord() only handles AC,AN,AF, and DP!
                vc = vc.subContextFromSample(sampGt.getSampleName());
            else // asked for a sample that this vc does not contain, so ignore this vc:
                return;
        }

        logger.debug("Next VC input = " + GATKVariantContextUtils.getLocation(genomeLocParser, vc));
        boolean curVcIsNotFiltered = vc.isNotFiltered();

        if (vcfrWaitingToMerge == null) {
            logger.debug("NOT Waiting to merge...");

            if (!filteredVcfrList.isEmpty())
                throw new ReviewedGATKException("filteredVcfrList should be empty if not waiting to merge a vc!");

            if (curVcIsNotFiltered) { // still need to wait before can release vc
                logger.debug("Waiting for new variant " + GATKVariantContextUtils.getLocation(genomeLocParser, vc));
                vcfrWaitingToMerge = new VCFRecord(vc, false);
            }
            else if (!emitOnlyMergedRecords) { // filtered records are never merged
                logger.debug("DIRECTLY output " + GATKVariantContextUtils.getLocation(genomeLocParser, vc));
                innerWriter.add(vc);
            }
        }
        else { // waiting to merge vcfrWaitingToMerge
            logger.debug("Waiting to merge " + GATKVariantContextUtils.getLocation(genomeLocParser, vcfrWaitingToMerge.vc));

            if (!curVcIsNotFiltered) {
                if (!emitOnlyMergedRecords) { // filtered records are never merged
                    logger.debug("Caching unprocessed output " + GATKVariantContextUtils.getLocation(genomeLocParser, vc));
                    filteredVcfrList.add(new VCFRecord(vc, false));
                }
            }
            else { // waiting to merge vcfrWaitingToMerge, and curVcIsNotFiltered. So, attempt to merge them:
                numRecordsAttemptToMerge++;
                boolean shouldAttemptToMerge = vcMergeRule.shouldAttemptToMerge(vcfrWaitingToMerge.vc, vc);
                logger.debug("shouldAttemptToMerge? = " + shouldAttemptToMerge);

                /*
                TODO: -- CONSIDER THE FOLLOWING EXAMPLE: WHAT DO WE WANT HERE??? --
                If the following 3 genotypes originally exist for a sample [at sites 1, 2, and 3]:
                1/1
                0|1
                0|1

                Then, after merging the first two, we have [at sites 1 and 3]:
                1/2
                0|1

                Then, not having merged would consider sites 2 and 3 as a MNP (since it's a diploid het site with haplotypes: REF-REF and ALT-ALT)
                But, since we merged sites 1 and 2, we get that sites 1-2 and 3 are counted as two haplotypes of: ALT-REF and ALT-ALT
                 */
                if (altAlleleStats != null)
                    altAlleleStats.updateSampleStats(vcfrWaitingToMerge.vc, vc, shouldAttemptToMerge);

                boolean mergedRecords = false;
                if (shouldAttemptToMerge) {
                    numRecordsSatisfyingMergeRule++;
                    VariantContext mergedVc = PhasingUtils.mergeIntoMNP(genomeLocParser, vcfrWaitingToMerge.vc, vc, referenceFileForMNPmerging, alleleMergeRule);

                    if (mergedVc != null) {
                        mergedRecords = true;

                        Map<String, Object> addedAttribs = vcMergeRule.addToMergedAttributes(vcfrWaitingToMerge.vc, vc);
                        addedAttribs.putAll(mergedVc.getAttributes());
                        mergedVc = new VariantContextBuilder(mergedVc).attributes(addedAttribs).make();

                        vcfrWaitingToMerge = new VCFRecord(mergedVc, true);
                        numMergedRecords++;
                    }
                }

                if (!mergedRecords) {
                    stopWaitingToMerge();
                    vcfrWaitingToMerge = new VCFRecord(vc, false);
                }
                logger.debug("Merged? = " + mergedRecords);
            }
        }
    }

    private void stopWaitingToMerge() {
        if (vcfrWaitingToMerge == null) {
            if (!filteredVcfrList.isEmpty())
                throw new ReviewedGATKException("filteredVcfrList should be empty if not waiting to merge a vc!");
            return;
        }

        if (!emitOnlyMergedRecords || vcfrWaitingToMerge.resultedFromMerge)
            innerWriter.add(vcfrWaitingToMerge.vc);
        vcfrWaitingToMerge = null;

        for (VCFRecord vcfr : filteredVcfrList)
            innerWriter.add(vcfr.vc);
        filteredVcfrList.clear();
    }

    /**
     * Gets a string representation of this object.
     *
     * @return
     */
    @Override
    public String toString() {
        return getClass().getName();
    }

    private static class VCFRecord {
        public VariantContext vc;
        public boolean resultedFromMerge;

        public VCFRecord(VariantContext vc, boolean resultedFromMerge) {
            this.vc = vc;
            this.resultedFromMerge = resultedFromMerge;
        }
    }

    private class AltAlleleStats {
        public int numSuccessiveGenotypes;
        public int numSuccessiveGenotypesAttemptedToBeMerged;

        public int oneSampleMissing;
        public int atLeastOneSampleNotCalledOrFiltered;
        public int segregationUnknown;
        public int eitherNotVariant;

        public int bothInPairHaveVariant;

        public int ref_ref_pair;
        public int ref_alt_pair;
        public int alt_ref_pair;
        public int alt_alt_pair;

        public int MNPsites;
        public int CHetSites;

        public AltAlleleStats() {
            this.numSuccessiveGenotypes = 0;
            this.numSuccessiveGenotypesAttemptedToBeMerged = 0;

            this.oneSampleMissing = 0;
            this.atLeastOneSampleNotCalledOrFiltered = 0;
            this.segregationUnknown = 0;
            this.eitherNotVariant = 0;

            this.bothInPairHaveVariant = 0;

            this.ref_ref_pair = 0;
            this.ref_alt_pair = 0;
            this.alt_ref_pair = 0;
            this.alt_alt_pair = 0;

            this.MNPsites = 0;
            this.CHetSites = 0;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append("Sample missing:\t" + oneSampleMissing + "\n");
            sb.append("Not called or filtered:\t" + atLeastOneSampleNotCalledOrFiltered + "\n");

            sb.append("* Number of successive pairs of genotypes:\t" + numSuccessiveGenotypes + "\n");
            sb.append("Number of successive pairs of genotypes with " + vcMergeRule + ":\t" + numSuccessiveGenotypesAttemptedToBeMerged + "\n");

            sb.append("Unknown segregation, " + vcMergeRule + ":\t" + segregationUnknown + "\n");
            sb.append("Not variant at least one of pair, segregation known, " + vcMergeRule + ":\t" + eitherNotVariant + "\n");
            sb.append("* Variant at both, segregation known, " + vcMergeRule + ":\t" + percentageString(bothInPairHaveVariant, numSuccessiveGenotypes) + "\n");

            sb.append("[Total haplotypes at pairs:\t" + (ref_ref_pair + ref_alt_pair + alt_ref_pair + alt_alt_pair) + "\n");
            sb.append("REF-REF:\t" + ref_ref_pair + "\n");
            sb.append("REF-ALT:\t" + ref_alt_pair + "\n");
            sb.append("ALT-REF:\t" + alt_ref_pair + "\n");
            sb.append("ALT-ALT:\t" + alt_alt_pair + "]\n");

            int hetAfterHetSites = MNPsites + CHetSites;
            sb.append("* Het-Het sites (with REF allele present at each):\t" + percentageString(hetAfterHetSites, bothInPairHaveVariant) + "\n");
            sb.append("* MNPs:\t" + percentageString(MNPsites, hetAfterHetSites) + "\n");
            sb.append("Compound Hets:\t" + CHetSites + "\n");

            return sb.toString();
        }

        private String percentageString(int count, int baseCount) {
            int NUM_DECIMAL_PLACES = 1;
            String percent = new Formatter().format("%." + NUM_DECIMAL_PLACES + "f", MathUtils.percentage(count, baseCount)).toString();
            return count + " (" + percent + "%)";
        }
    }

    private class AltAlleleStatsForSamples {
        private Map<String, AltAlleleStats> sampleStats;

        public AltAlleleStatsForSamples() {
            this.sampleStats = new HashMap<String, AltAlleleStats>();
        }

        public void updateSampleStats(VariantContext vc1, VariantContext vc2, boolean shouldAttemptToMerge) {
            if (vc1.isFiltered() || vc2.isFiltered())
                return;

            Set<String> allSamples = new TreeSet<String>(vc1.getSampleNames());
            allSamples.addAll(vc2.getSampleNames());

            for (String samp : allSamples) {
                AltAlleleStats aas = sampleStats.get(samp);
                if (aas == null) {
                    aas = new AltAlleleStats();
                    sampleStats.put(samp, aas);
                }

                Genotype gt1 = vc1.getGenotype(samp);
                Genotype gt2 = vc2.getGenotype(samp);
                if (gt1 == null || gt2 == null) {
                    aas.oneSampleMissing++;
                }
                else if (gt1.isNoCall() || gt1.isFiltered() || gt2.isNoCall() || gt2.isFiltered()) {
                    aas.atLeastOneSampleNotCalledOrFiltered++;
                }
                else {
                    aas.numSuccessiveGenotypes++;

                    if (shouldAttemptToMerge) {
                        aas.numSuccessiveGenotypesAttemptedToBeMerged++;

                        if (!PhasingUtils.alleleSegregationIsKnown(gt1, gt2)) {
                            aas.segregationUnknown++;
                            logger.debug("Unknown segregation of alleles [not phased] for " + samp + " at " + GATKVariantContextUtils.getLocation(genomeLocParser, vc1) + ", " + GATKVariantContextUtils.getLocation(genomeLocParser, vc2));
                        }
                        else if (gt1.isHomRef() || gt2.isHomRef()) {
                            logger.debug("gt1.isHomRef() || gt2.isHomRef() for " + samp + " at " + GATKVariantContextUtils.getLocation(genomeLocParser, vc1) + ", " + GATKVariantContextUtils.getLocation(genomeLocParser, vc2));
                            aas.eitherNotVariant++;
                        }
                        else { // BOTH gt1 and gt2 have at least one variant allele (so either hets, or homozygous variant):
                            aas.bothInPairHaveVariant++;

                            List<Allele> site1Alleles = gt1.getAlleles();
                            List<Allele> site2Alleles = gt2.getAlleles();

                            Iterator<Allele> all2It = site2Alleles.iterator();
                            for (Allele all1 : site1Alleles) {
                                Allele all2 = all2It.next(); // this is OK, since alleleSegregationIsKnown(gt1, gt2)

                                if (all1.isReference()) {
                                    if (all2.isReference())
                                        aas.ref_ref_pair++;
                                    else
                                        aas.ref_alt_pair++;
                                }
                                else { // all1.isNonReference()
                                    if (all2.isReference())
                                        aas.alt_ref_pair++;
                                    else
                                        aas.alt_alt_pair++;
                                }
                            }

                            // Check MNPs vs. CHets:
                            if (containsRefAllele(site1Alleles) && containsRefAllele(site2Alleles)) {
                                logger.debug("HET-HET for " + samp + " at " + GATKVariantContextUtils.getLocation(genomeLocParser, vc1) + ", " + GATKVariantContextUtils.getLocation(genomeLocParser, vc2));
                                if (logger.isDebugEnabled() && !(gt1.isHet() && gt2.isHet()))
                                    throw new ReviewedGATKException("Since !gt1.isHomRef() && !gt2.isHomRef(), yet both have ref alleles, they BOTH must be hets!");

                                // There's the potential to only have REF-ALT, ALT-REF (CHet), or possibly ALT-ALT together (MNP)
                                boolean hasMNP = false;

                                all2It = site2Alleles.iterator();
                                for (Allele all1 : site1Alleles) {
                                    Allele all2 = all2It.next(); // this is OK, since alleleSegregationIsKnown(gt1, gt2)

                                    if (all1.isNonReference() && all2.isNonReference()) {
                                        hasMNP = true; // has at least one haplotype of ALT-ALT that segregates!
                                        break;
                                    }
                                }

                                if (hasMNP)
                                    aas.MNPsites++;
                                else
                                    aas.CHetSites++;
                            }
                        }
                    }
                }
            }
        }

        private boolean containsRefAllele(List<Allele> siteAlleles) {
            for (Allele all : siteAlleles) {
                if (all.isReference())
                    return true;
            }

            return false;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("-------------------------------------------------------------------------\n");
            sb.append("Per-sample alternate allele statistics [" + vcMergeRule + "]\n");
            sb.append("-------------------------------------------------------------------------");

            for (Map.Entry<String, AltAlleleStats> sampAltAllStatsEntry : sampleStats.entrySet()) {
                String samp = sampAltAllStatsEntry.getKey();
                AltAlleleStats stats = sampAltAllStatsEntry.getValue();
                sb.append("\n* Sample:\t" + samp + "\n" + stats);
            }

            return sb.toString();
        }
    }

    /**
     * Check the return from PrintStream.checkError() if underlying stream for a java.io.PrintStream
     * @return false, no error since the underlying stream is not a java.io.PrintStream
     */
    public boolean checkError(){
        return false;
    }
}



/*
 External classes:
 */

abstract class VariantContextMergeRule {
    abstract public boolean shouldAttemptToMerge(VariantContext vc1, VariantContext vc2);

    public Map<String, Object> addToMergedAttributes(VariantContext vc1, VariantContext vc2) {
        return new HashMap<String, Object>();
    }
}

class DistanceMergeRule extends VariantContextMergeRule {
    private int maxGenomicDistanceForMNP;
    private GenomeLocParser genomeLocParser;

    public DistanceMergeRule(int maxGenomicDistanceForMNP, GenomeLocParser genomeLocParser) {
        this.maxGenomicDistanceForMNP = maxGenomicDistanceForMNP;
        this.genomeLocParser = genomeLocParser;
    }

    public boolean shouldAttemptToMerge(VariantContext vc1, VariantContext vc2) {
        return minDistance(vc1, vc2) <= maxGenomicDistanceForMNP;
    }

    public String toString() {
        return "Merge distance <= " + maxGenomicDistanceForMNP;
    }

    public int minDistance(VariantContext vc1, VariantContext vc2) {
        return GATKVariantContextUtils.getLocation(genomeLocParser, vc1).minDistance(GATKVariantContextUtils.getLocation(genomeLocParser, vc2));
    }
}


class ExistsDoubleAltAlleleMergeRule extends PhasingUtils.AlleleMergeRule {
    public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2) {
        return PhasingUtils.someSampleHasDoubleNonReferenceAllele(vc1, vc2);
    }

    public String toString() {
        return super.toString() + ", some sample has a MNP of ALT alleles";
    }
}

class SegregatingMNPmergeAllelesRule extends ExistsDoubleAltAlleleMergeRule {
    public SegregatingMNPmergeAllelesRule() {
        super();
    }

    public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2) {
        // Must be interesting AND consistent:
        return super.allelesShouldBeMerged(vc1, vc2) && PhasingUtils.doubleAllelesSegregatePerfectlyAmongSamples(vc1, vc2);
    }

    public String toString() {
        return super.toString() + ", all alleles segregate consistently";
    }
}