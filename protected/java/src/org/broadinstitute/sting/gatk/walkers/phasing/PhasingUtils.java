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

package org.broadinstitute.sting.gatk.walkers.phasing;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.variant.variantcontext.*;

import java.util.*;

/**
 * [Short one sentence description of this walker]
 * <p/>
 * <p>
 * [Functionality of this walker]
 * </p>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *  </pre>
 *
 * @author Your Name
 * @since Date created
 */
class PhasingUtils {
    static VariantContext mergeIntoMNP(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile, AlleleMergeRule alleleMergeRule) {
        if (!mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2))
            return null;

        // Check that it's logically possible to merge the VCs:
        if (!allSamplesAreMergeable(vc1, vc2))
            return null;

        // Check if there's a "point" in merging the VCs (e.g., annotations could be changed)
        if (!alleleMergeRule.allelesShouldBeMerged(vc1, vc2))
            return null;

        return reallyMergeIntoMNP(vc1, vc2, referenceFile);
    }

    static VariantContext reallyMergeIntoMNP(VariantContext vc1, VariantContext vc2, ReferenceSequenceFile referenceFile) {
        int startInter = vc1.getEnd() + 1;
        int endInter = vc2.getStart() - 1;
        byte[] intermediateBases = null;
        if (startInter <= endInter) {
            intermediateBases = referenceFile.getSubsequenceAt(vc1.getChr(), startInter, endInter).getBases();
            StringUtil.toUpperCase(intermediateBases);
        }
        MergedAllelesData mergeData = new MergedAllelesData(intermediateBases, vc1, vc2); // ensures that the reference allele is added

        GenotypesContext mergedGenotypes = GenotypesContext.create();
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            List<Allele> mergedAllelesForSample = new LinkedList<Allele>();

            /* NOTE: Since merged alleles are added to mergedAllelesForSample in the SAME order as in the input VC records,
               we preserve phase information (if any) relative to whatever precedes vc1:
             */
            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next(); // this is OK, since allSamplesAreMergeable()

                Allele mergedAllele = mergeData.ensureMergedAllele(all1, all2);
                mergedAllelesForSample.add(mergedAllele);
            }

            double mergedGQ = Math.max(gt1.getLog10PError(), gt2.getLog10PError());

            Map<String, Object> mergedGtAttribs = new HashMap<String, Object>();
            PhaseAndQuality phaseQual = calcPhaseForMergedGenotypes(gt1, gt2);
            if (phaseQual.PQ != null)
                mergedGtAttribs.put(ReadBackedPhasing.PQ_KEY, phaseQual.PQ);

            Genotype mergedGt = new GenotypeBuilder(gt1.getSampleName(), mergedAllelesForSample).log10PError(mergedGQ).attributes(mergedGtAttribs).phased(phaseQual.isPhased).make();
            mergedGenotypes.add(mergedGt);
        }

        String mergedName = mergeVariantContextNames(vc1.getSource(), vc2.getSource());
        double mergedLog10PError = Math.min(vc1.getLog10PError(), vc2.getLog10PError());
        Set<String> mergedFilters = new HashSet<String>(); // Since vc1 and vc2 were unfiltered, the merged record remains unfiltered
        Map<String, Object> mergedAttribs = mergeVariantContextAttributes(vc1, vc2);

        // ids
        List<String> mergedIDs = new ArrayList<String>();
        if ( vc1.hasID() ) mergedIDs.add(vc1.getID());
        if ( vc2.hasID() ) mergedIDs.add(vc2.getID());
        String mergedID = mergedIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(VCFConstants.ID_FIELD_SEPARATOR, mergedIDs);

        VariantContextBuilder mergedBuilder = new VariantContextBuilder(mergedName, vc1.getChr(), vc1.getStart(), vc2.getEnd(), mergeData.getAllMergedAlleles()).id(mergedID).genotypes(mergedGenotypes).log10PError(mergedLog10PError).filters(mergedFilters).attributes(mergedAttribs);
        VariantContextUtils.calculateChromosomeCounts(mergedBuilder, true);
        return mergedBuilder.make();
    }

    static String mergeVariantContextNames(String name1, String name2) {
        return name1 + "_" + name2;
    }

    static Map<String, Object> mergeVariantContextAttributes(VariantContext vc1, VariantContext vc2) {
        Map<String, Object> mergedAttribs = new HashMap<String, Object>();

        List<VariantContext> vcList = new LinkedList<VariantContext>();
        vcList.add(vc1);
        vcList.add(vc2);

        String[] MERGE_OR_ATTRIBS = {VCFConstants.DBSNP_KEY};
        for (String orAttrib : MERGE_OR_ATTRIBS) {
            boolean attribVal = false;
            for (VariantContext vc : vcList) {
                attribVal = vc.getAttributeAsBoolean(orAttrib, false);
                if (attribVal) // already true, so no reason to continue:
                    break;
            }
            mergedAttribs.put(orAttrib, attribVal);
        }

        return mergedAttribs;
    }

    static boolean mergeIntoMNPvalidationCheck(GenomeLocParser genomeLocParser, VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = GATKVariantContextUtils.getLocation(genomeLocParser, vc1);
        GenomeLoc loc2 = GATKVariantContextUtils.getLocation(genomeLocParser, vc2);

        if (!loc1.onSameContig(loc2))
            throw new ReviewedStingException("Can only merge vc1, vc2 if on the same chromosome");

        if (!loc1.isBefore(loc2))
            throw new ReviewedStingException("Can only merge if vc1 is BEFORE vc2");

        if (vc1.isFiltered() || vc2.isFiltered())
            return false;

        if (!vc1.getSampleNames().equals(vc2.getSampleNames())) // vc1, vc2 refer to different sample sets
            return false;

        if (!allGenotypesAreUnfilteredAndCalled(vc1) || !allGenotypesAreUnfilteredAndCalled(vc2))
            return false;

        return true;
    }

    static boolean allGenotypesAreUnfilteredAndCalled(VariantContext vc) {
        for (final Genotype gt : vc.getGenotypes()) {
            if (gt.isNoCall() || gt.isFiltered())
                return false;
        }

        return true;
    }

    static boolean allSamplesAreMergeable(VariantContext vc1, VariantContext vc2) {
        // Check that each sample's genotype in vc2 is uniquely appendable onto its genotype in vc1:
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            if (!alleleSegregationIsKnown(gt1, gt2)) // can merge if: phased, or if either is a hom
                return false;
        }

        return true;
    }

    static boolean alleleSegregationIsKnown(Genotype gt1, Genotype gt2) {
        if (gt1.getPloidy() != gt2.getPloidy())
            return false;

        /* If gt2 is phased or hom, then could even be MERGED with gt1 [This is standard].

           HOWEVER, EVEN if this is not the case, but gt1.isHom(),
           it is trivially known that each of gt2's alleles segregate with the single allele type present in gt1.
         */
        return (gt2.isPhased() || gt2.isHom() || gt1.isHom());
    }

    static PhaseAndQuality calcPhaseForMergedGenotypes(Genotype gt1, Genotype gt2) {
        if (gt2.isPhased() || gt2.isHom())
            return new PhaseAndQuality(gt1); // maintain the phase of gt1

        if (!gt1.isHom())
            throw new ReviewedStingException("alleleSegregationIsKnown(gt1, gt2) implies: gt2.genotypesArePhased() || gt2.isHom() || gt1.isHom()");

        /* We're dealing with: gt1.isHom(), gt2.isHet(), !gt2.genotypesArePhased(); so, the merged (het) Genotype is not phased relative to the previous Genotype

           For example, if we're merging the third Genotype with the second one:
           0/1
           1|1
           0/1

           Then, we want to output:
           0/1
           1/2
         */
        return new PhaseAndQuality(gt2); // maintain the phase of gt2 [since !gt2.genotypesArePhased()]
    }

    static boolean someSampleHasDoubleNonReferenceAllele(VariantContext vc1, VariantContext vc2) {
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next(); // this is OK, since allSamplesAreMergeable()

                if (all1.isNonReference() && all2.isNonReference()) // corresponding alleles are alternate
                    return true;
            }
        }

        return false;
    }

    static boolean doubleAllelesSegregatePerfectlyAmongSamples(VariantContext vc1, VariantContext vc2) {
        // Check that Alleles at vc1 and at vc2 always segregate together in all samples (including reference):
        Map<Allele, Allele> allele1ToAllele2 = new HashMap<Allele, Allele>();
        Map<Allele, Allele> allele2ToAllele1 = new HashMap<Allele, Allele>();

        // Note the segregation of the alleles for the reference genome:
        allele1ToAllele2.put(vc1.getReference(), vc2.getReference());
        allele2ToAllele1.put(vc2.getReference(), vc1.getReference());

        // Note the segregation of the alleles for each sample (and check that it is consistent with the reference and all previous samples).
        for (final Genotype gt1 : vc1.getGenotypes()) {
            Genotype gt2 = vc2.getGenotype(gt1.getSampleName());

            List<Allele> site1Alleles = gt1.getAlleles();
            List<Allele> site2Alleles = gt2.getAlleles();

            Iterator<Allele> all2It = site2Alleles.iterator();
            for (Allele all1 : site1Alleles) {
                Allele all2 = all2It.next();

                Allele all1To2 = allele1ToAllele2.get(all1);
                if (all1To2 == null)
                    allele1ToAllele2.put(all1, all2);
                else if (!all1To2.equals(all2)) // all1 segregates with two different alleles at site 2
                    return false;

                Allele all2To1 = allele2ToAllele1.get(all2);
                if (all2To1 == null)
                    allele2ToAllele1.put(all2, all1);
                else if (!all2To1.equals(all1)) // all2 segregates with two different alleles at site 1
                    return false;
            }
        }

        return true;
    }

    abstract static class AlleleMergeRule {
        // vc1, vc2 are ONLY passed to allelesShouldBeMerged() if mergeIntoMNPvalidationCheck(genomeLocParser, vc1, vc2) AND allSamplesAreMergeable(vc1, vc2):
        abstract public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2);

        public String toString() {
            return "all samples are mergeable";
        }
    }

    static class AlleleOneAndTwo {
        private Allele all1;
        private Allele all2;

        public AlleleOneAndTwo(Allele all1, Allele all2) {
            this.all1 = all1;
            this.all2 = all2;
        }

        public int hashCode() {
            return all1.hashCode() + all2.hashCode();
        }

        public boolean equals(Object other) {
            if (!(other instanceof AlleleOneAndTwo))
                return false;

            AlleleOneAndTwo otherAot = (AlleleOneAndTwo) other;
            return (this.all1.equals(otherAot.all1) && this.all2.equals(otherAot.all2));
        }
    }

    static class MergedAllelesData {
        private Map<AlleleOneAndTwo, Allele> mergedAlleles;
        private byte[] intermediateBases;
        private int intermediateLength;

        public MergedAllelesData(byte[] intermediateBases, VariantContext vc1, VariantContext vc2) {
            this.mergedAlleles = new HashMap<AlleleOneAndTwo, Allele>(); // implemented equals() and hashCode() for AlleleOneAndTwo
            this.intermediateBases = intermediateBases;
            this.intermediateLength = this.intermediateBases != null ? this.intermediateBases.length : 0;

            this.ensureMergedAllele(vc1.getReference(), vc2.getReference(), true);
        }

        public Allele ensureMergedAllele(Allele all1, Allele all2) {
            return ensureMergedAllele(all1, all2, false); // false <-> since even if all1+all2 = reference, it was already created in the constructor
        }

        private Allele ensureMergedAllele(Allele all1, Allele all2, boolean creatingReferenceForFirstTime) {
            AlleleOneAndTwo all12 = new AlleleOneAndTwo(all1, all2);
            Allele mergedAllele = mergedAlleles.get(all12);

            if (mergedAllele == null) {
                byte[] bases1 = all1.getBases();
                byte[] bases2 = all2.getBases();

                byte[] mergedBases = new byte[bases1.length + intermediateLength + bases2.length];
                System.arraycopy(bases1, 0, mergedBases, 0, bases1.length);
                if (intermediateBases != null)
                    System.arraycopy(intermediateBases, 0, mergedBases, bases1.length, intermediateLength);
                System.arraycopy(bases2, 0, mergedBases, bases1.length + intermediateLength, bases2.length);

                mergedAllele = Allele.create(mergedBases, creatingReferenceForFirstTime);
                mergedAlleles.put(all12, mergedAllele);
            }

            return mergedAllele;
        }

        public Set<Allele> getAllMergedAlleles() {
            return new HashSet<Allele>(mergedAlleles.values());
        }
    }

    static class PhaseAndQuality {
        public boolean isPhased;
        public Double PQ = null;

        public PhaseAndQuality(Genotype gt) {
            this.isPhased = gt.isPhased();
            if (this.isPhased) {
                this.PQ = gt.getAttributeAsDouble(ReadBackedPhasing.PQ_KEY, -1);
                if ( this.PQ == -1 ) this.PQ = null;
            }
        }
    }
}
