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

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

/**
 * Compute the most likely genotype combination and phasing for trios and parent/child pairs
 *
 * <p>
 * This tool performs two functions:
 * </p>
 * <ol>
 *     <li>Compute the most likely genotype combination of trios and parent/child pairs given their genotype likelihoods and a mutation prior;</li>
 *     <li>Phase all sites were parent/child transmission can be inferred unambiguously. </li>
 * </ol>
 *
 * <p>The tool ultimately reports the genotype combination (and hence phasing) probability.</p>
 *
 * <h4>Ambiguous sites are:</h4>
 * <ul>
 *     <li>Sites where all individuals are heterozygous</li>
 *     <li>Sites where there is a Mendelian violation</li>
 * </ul>
 *
 * <h4>Missing genotypes are handled as follows:</h4>
 * <ul>
 *     <li>In parent/child pairs: If an individual genotype is missing at one site, the other one is phased if it is homozygous. No phasing probability is emitted.</li>
 *     <li>In trios: If the child is missing, parents are treated as separate individuals and phased if homozygous. No phasing probability is emitted.</li>
 *     <li>In trios: If one of the parents is missing, it is handled like a parent/child pair. Phasing is done unless both the parent and child are heterozygous and a phasing probability is emitted.</li>
 *     <li>In trios: If two individuals are missing, the remaining individual is phased if it is homozygous. No phasing probability is emitted.</li>
 * </ul>
 *
 * <h3>Input</h3>
 * <p>
 * <ul>
 *     <li>A VCF variant set containing trio(s) and/or parent/child pair(s).</li>
 *     <li>A PED pedigree file containing the description of the individuals relationships.</li>
 * </ul>
 * </p>
 *
 * <h3>Important options</h3>
 *     <ul>
 *         <li>MendelianViolationsFile: An optional argument for reporting. If a file is specified, all sites that
 *         remain in mendelian violation after being assigned the most likely genotype combination will be reported
 *         there. Information reported: chromosome, position, filter, allele count in VCF, family, transmission
 *         probability, and each individual genotype, depth, allelic depth and likelihoods.</li>
 *         <li>DeNovoPrior: Prior probability of de novo mutations. The default value of 1e-8 is fairly stringent, so if
 *         you are interested in maximizing sensitivity at the expense of specificity (i.e. are ok with seeing some false
 *         positives as long as all true positives are detected) you will need to relax this value.</li>
 *     </ul>
 *
 * <h3>Output</h3>
 * <p>
 * An VCF with genotypes recalibrated as most likely under the familial constraint and phased by descent (where non
 * ambiguous).
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T PhaseByTransmission \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -ped input.ped \
 *   -o output.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
public class PhaseByTransmission extends RodWalker<HashMap<Byte,Integer>, HashMap<Byte,Integer>> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName = "mvf",required = false,fullName = "MendelianViolationsFile", doc="File to output the mendelian violation details.")
    private PrintStream mvFile = null;

    @Argument(shortName = "prior",required = false,fullName = "DeNovoPrior", doc="Prior for de novo mutations. Default: 1e-8")
    private double deNovoPrior=1e-8;

    @Argument(shortName = "fatherAlleleFirst",required = false,fullName = "FatherAlleleFirst", doc="Ouputs the father allele as the first allele in phased child genotype. i.e. father|mother rather than mother|father.")
    private boolean fatherFAlleleFirst=false;

    @Output
    protected VariantContextWriter vcfWriter = null;

    private final String SOURCE_NAME = "PhaseByTransmission";

    public final double NO_TRANSMISSION_PROB = -1.0;

    private ArrayList<Sample> trios = new ArrayList<Sample>();

    //Matrix of priors for all genotype combinations
    private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> mvCountMatrix;

    //Matrix of allele transmission
    private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,TrioPhase>>> transmissionMatrix;

    //Metrics counters hash keys
    private final Byte NUM_TRIO_GENOTYPES_CALLED = 0;
    private final Byte NUM_TRIO_GENOTYPES_NOCALL = 1;
    private final Byte NUM_TRIO_GENOTYPES_PHASED = 2;
    private final Byte NUM_TRIO_HET_HET_HET = 3;
    private final Byte NUM_TRIO_VIOLATIONS = 4;
    private final Byte NUM_TRIO_DOUBLE_VIOLATIONS = 10;
    private final Byte NUM_PAIR_GENOTYPES_CALLED = 5;
    private final Byte NUM_PAIR_GENOTYPES_NOCALL = 6;
    private final Byte NUM_PAIR_GENOTYPES_PHASED = 7;
    private final Byte NUM_PAIR_HET_HET = 8;
    private final Byte NUM_PAIR_VIOLATIONS = 9;
    private final Byte NUM_GENOTYPES_MODIFIED = 11;

    //Random number generator
    private Random rand = new Random();

    private enum FamilyMember {
        MOTHER,
        FATHER,
        CHILD
    }

    //Stores a conceptual trio or parent/child pair genotype combination along with its phasing.
    //This combination can then be "applied" to a given trio or pair using the getPhasedGenotypes method.
    private class TrioPhase {

        //Create 2 fake alleles
        //The actual bases will never be used but the Genotypes created using the alleles will be.
        private final Allele REF = Allele.create("A",true);
        private final Allele VAR = Allele.create("A",false);
        private final Allele NO_CALL = Allele.create(".",false);
        private final String DUMMY_NAME = "DummySample";

        private EnumMap<FamilyMember,Genotype> trioPhasedGenotypes = new EnumMap<FamilyMember, Genotype>(FamilyMember.class);

        private ArrayList<Allele> getAlleles(GenotypeType genotype){
            ArrayList<Allele> alleles = new ArrayList<Allele>(2);
            if(genotype == GenotypeType.HOM_REF){
                alleles.add(REF);
                alleles.add(REF);
            }
            else if(genotype == GenotypeType.HET){
                alleles.add(REF);
                alleles.add(VAR);
            }
            else if(genotype == GenotypeType.HOM_VAR){
                alleles.add(VAR);
                alleles.add(VAR);
            }
            else{
                return null;
            }
            return alleles;
        }

        private boolean isPhasable(GenotypeType genotype){
            return genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HET || genotype == GenotypeType.HOM_VAR;
        }

        //Create a new Genotype based on information from a single individual
        //Homozygous genotypes will be set as phased, heterozygous won't be
        private void phaseSingleIndividualAlleles(GenotypeType genotype, FamilyMember familyMember){
            boolean phase = genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HOM_VAR;
            trioPhasedGenotypes.put(familyMember, makeGenotype(genotype, phase));
        }

        private Genotype makeGenotype(final GenotypeType type, boolean phase) {
            return makeGenotype(getAlleles(type), phase);
        }

        private Genotype makeGenotype(final List<Allele> alleles, boolean phase) {
            final GenotypeBuilder gb = new GenotypeBuilder(DUMMY_NAME, alleles);
            gb.phased(phase);
            return gb.make();
        }

        //Find the phase for a parent/child pair
        private void phasePairAlleles(GenotypeType parentGenotype, GenotypeType childGenotype, FamilyMember parent){

            //Special case for Het/Het as it is ambiguous
            if(parentGenotype == GenotypeType.HET && childGenotype == GenotypeType.HET){
                trioPhasedGenotypes.put(parent, makeGenotype(parentGenotype, false));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childGenotype, false));
                return;
            }

            ArrayList<Allele> parentAlleles = getAlleles(parentGenotype);
            ArrayList<Allele> childAlleles = getAlleles(childGenotype);
            ArrayList<Allele> parentPhasedAlleles = new ArrayList<Allele>(2);
            ArrayList<Allele> childPhasedAlleles = new ArrayList<Allele>(2);

            //If there is a possible phasing between the parent and child => phase
            int childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(0));
            if(childTransmittedAlleleIndex > -1){
                trioPhasedGenotypes.put(parent, makeGenotype(parentAlleles, true));
                childPhasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
                if(parent.equals(FamilyMember.MOTHER))
                        childPhasedAlleles.add(childAlleles.get(0));
                else
                        childPhasedAlleles.add(0,childAlleles.get(0));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childPhasedAlleles, true));
            }
            else if((childTransmittedAlleleIndex = childAlleles.indexOf(parentAlleles.get(1))) > -1){
                parentPhasedAlleles.add(parentAlleles.get(1));
                parentPhasedAlleles.add(parentAlleles.get(0));
                trioPhasedGenotypes.put(parent, makeGenotype(parentPhasedAlleles, true));
                childPhasedAlleles.add(childAlleles.remove(childTransmittedAlleleIndex));
                if(parent.equals(FamilyMember.MOTHER))
                    childPhasedAlleles.add(childAlleles.get(0));
                else
                    childPhasedAlleles.add(0,childAlleles.get(0));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childPhasedAlleles, true));
            }
            //This is a Mendelian Violation => Do not phase
            else{
                trioPhasedGenotypes.put(parent, makeGenotype(parentGenotype, false));
                trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childGenotype, false));
            }
        }

        //Phases a family by transmission
        private void phaseFamilyAlleles(GenotypeType mother, GenotypeType father, GenotypeType child){

            Set<ArrayList<Allele>> possiblePhasedChildGenotypes = new HashSet<ArrayList<Allele>>();
            ArrayList<Allele> motherAlleles = getAlleles(mother);
            ArrayList<Allele> fatherAlleles = getAlleles(father);
            ArrayList<Allele> childAlleles = getAlleles(child);

            //Build all possible child genotypes for the given parent's genotypes
            for (Allele momAllele : motherAlleles) {
                for (Allele fatherAllele : fatherAlleles) {
                    ArrayList<Allele> possiblePhasedChildAlleles = new ArrayList<Allele>(2);
                    possiblePhasedChildAlleles.add(momAllele);
                    possiblePhasedChildAlleles.add(fatherAllele);
                    possiblePhasedChildGenotypes.add(possiblePhasedChildAlleles);
                }
            }

            for (ArrayList<Allele> childPhasedAllelesAlleles : possiblePhasedChildGenotypes) {
                int firstAlleleIndex = childPhasedAllelesAlleles.indexOf(childAlleles.get(0));
                int secondAlleleIndex = childPhasedAllelesAlleles.lastIndexOf(childAlleles.get(1));
                //If a possible combination has been found, create the genotypes
                if (firstAlleleIndex != secondAlleleIndex && firstAlleleIndex > -1 && secondAlleleIndex > -1) {
                    //Create mother's genotype
                    ArrayList<Allele> motherPhasedAlleles = new ArrayList<Allele>(2);
                    motherPhasedAlleles.add(childPhasedAllelesAlleles.get(0));
                    if(motherAlleles.get(0) != motherPhasedAlleles.get(0))
                        motherPhasedAlleles.add(motherAlleles.get(0));
                    else
                        motherPhasedAlleles.add(motherAlleles.get(1));
                    trioPhasedGenotypes.put(FamilyMember.MOTHER, makeGenotype(motherPhasedAlleles, true));

                    //Create father's genotype
                    ArrayList<Allele> fatherPhasedAlleles = new ArrayList<Allele>(2);
                    fatherPhasedAlleles.add(childPhasedAllelesAlleles.get(1));
                    if(fatherAlleles.get(0) != fatherPhasedAlleles.get(0))
                        fatherPhasedAlleles.add(fatherAlleles.get(0));
                    else
                        fatherPhasedAlleles.add(fatherAlleles.get(1));
                    trioPhasedGenotypes.put(FamilyMember.FATHER, makeGenotype(fatherPhasedAlleles,true));

                    //Create child's genotype
                    trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(childPhasedAllelesAlleles,true));

                    //Once a phased combination is found; exit
                    return;
                }
            }

            //If this is reached then no phasing could be found
            trioPhasedGenotypes.put(FamilyMember.MOTHER, makeGenotype(mother,false));
            trioPhasedGenotypes.put(FamilyMember.FATHER, makeGenotype(father,false));
            trioPhasedGenotypes.put(FamilyMember.CHILD, makeGenotype(child,false));
        }

        /*  Constructor: Creates a conceptual trio genotype combination from the given genotypes.
            If one or more genotypes are set as NO_CALL or UNAVAILABLE, it will phase them like a pair
            or single individual.
        */
        public TrioPhase(GenotypeType mother, GenotypeType father, GenotypeType child){

            //Take care of cases where one or more family members are no call
            if(!isPhasable(child)){
                phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
            }
            else if(!isPhasable(mother)){
                phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
                if(!isPhasable(father)){
                    phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                    phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
                }
                else
                    phasePairAlleles(father, child, FamilyMember.FATHER);
            }
            else if(!isPhasable(father)){
                phasePairAlleles(mother, child, FamilyMember.MOTHER);
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
            }
            //Special case for Het/Het/Het as it is ambiguous
            else if(mother == GenotypeType.HET && father  == GenotypeType.HET && child == GenotypeType.HET){
                phaseSingleIndividualAlleles(mother, FamilyMember.MOTHER);
                phaseSingleIndividualAlleles(father, FamilyMember.FATHER);
                phaseSingleIndividualAlleles(child, FamilyMember.CHILD);
            }
            //All family members have genotypes and at least one of them is not Het
            else{
                phaseFamilyAlleles(mother, father, child);
            }

            //If child should phased genotype should be father first, then swap the alleles
            if(fatherFAlleleFirst && trioPhasedGenotypes.get(FamilyMember.CHILD).isPhased()){
                ArrayList<Allele> childAlleles = new ArrayList<Allele>(trioPhasedGenotypes.get(FamilyMember.CHILD).getAlleles());
                childAlleles.add(childAlleles.remove(0));
                trioPhasedGenotypes.put(FamilyMember.CHILD,makeGenotype(childAlleles,true));
            }

        }

        /**
         * Applies the trio genotype combination to the given trio.
         * @param ref: Reference allele
         * @param alt: Alternate allele
         * @param motherGenotype: Genotype of the mother to phase using this trio genotype combination
         * @param fatherGenotype: Genotype of the father to phase using this trio genotype combination
         * @param childGenotype: Genotype of the child to phase using this trio genotype combination
         * @param transmissionProb: Probability for this trio genotype combination to be correct (pass NO_TRANSMISSION_PROB if unavailable)
         * @param phasedGenotypes: An ArrayList<Genotype> to which the newly phased genotypes are added in the following order: Mother, Father, Child
         */
        public void getPhasedGenotypes(Allele ref, Allele alt, Genotype motherGenotype, Genotype fatherGenotype, Genotype childGenotype, double transmissionProb,ArrayList<Genotype> phasedGenotypes){
            phasedGenotypes.add(getPhasedGenotype(ref,alt,motherGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.MOTHER)));
            phasedGenotypes.add(getPhasedGenotype(ref,alt,fatherGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.FATHER)));
            phasedGenotypes.add(getPhasedGenotype(ref,alt,childGenotype,transmissionProb,this.trioPhasedGenotypes.get(FamilyMember.CHILD)));
        }

        private Genotype getPhasedGenotype(Allele refAllele, Allele altAllele, Genotype genotype, double transmissionProb, Genotype phasedGenotype){

            int phredScoreTransmission = -1;
            if(transmissionProb != NO_TRANSMISSION_PROB){
                double dphredScoreTransmission = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1 - (transmissionProb)));
                phredScoreTransmission = dphredScoreTransmission < Byte.MAX_VALUE ? (byte)dphredScoreTransmission : Byte.MAX_VALUE;
            }
           //Handle null, missing and unavailable genotypes
           //Note that only cases where a null/missing/unavailable genotype was passed in the first place can lead to a null/missing/unavailable
           //genotype so it is safe to return the original genotype in this case.
           //In addition, if the phasing confidence is 0, then return the unphased, original genotypes.
           if(phredScoreTransmission ==0 || genotype == null || !isPhasable(genotype.getType()))
               return genotype;

           //Add the transmission probability
           Map<String, Object> genotypeAttributes = new HashMap<String, Object>();
           genotypeAttributes.putAll(genotype.getExtendedAttributes());
           if(transmissionProb>NO_TRANSMISSION_PROB)
                genotypeAttributes.put(GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY, phredScoreTransmission);

            ArrayList<Allele> phasedAlleles = new ArrayList<Allele>(2);
            for(Allele allele : phasedGenotype.getAlleles()){
                if(allele.isReference())
                    phasedAlleles.add(refAllele);
                else if(allele.isNonReference())
                    phasedAlleles.add(altAllele);
                    //At this point there should not be any other alleles left
                else
                    throw new UserException(String.format("BUG: Unexpected allele: %s. Please report.",allele.toString()));

            }

            //Compute the new Log10Error if the genotype is different from the original genotype
            double log10Error;
            if(genotype.getType() == phasedGenotype.getType())
                log10Error = genotype.getLog10PError();
            else
                log10Error = genotype.getLikelihoods().getLog10GQ(phasedGenotype.getType());

            return new GenotypeBuilder(genotype).alleles(phasedAlleles)
                    .log10PError(log10Error)
                    .attributes(genotypeAttributes)
                    .phased(phasedGenotype.isPhased()).make();
        }


    }

    /**
     * Parse the familial relationship specification, build the transmission matrices and initialize VCF writer
     */
    public void initialize() {
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //Get the trios from the families passed as ped
        setTrios(vcfSamples);
        if(trios.size()<1)
            throw new UserException.BadInput("No PED file passed or no *non-skipped* trios found in PED file. Aborted.");


        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(GATKVCFUtils.getHeaderFields(this.getToolkit()));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        buildMatrices();

        if(mvFile != null)
            mvFile.println("CHROM\tPOS\tAC\tFAMILY\tTP\tMOTHER_GT\tMOTHER_DP\tMOTHER_AD\tMOTHER_PL\tFATHER_GT\tFATHER_DP\tFATHER_AD\tFATHER_PL\tCHILD_GT\tCHILD_DP\tCHILD_AD\tCHILD_PL");

    }

    /**
     * Select trios and parent/child pairs only
     */
    private void setTrios(Set<String> vcfSamples){

        Map<String,Set<Sample>> families = this.getSampleDB().getFamilies(vcfSamples);
        Set<Sample> family;
        ArrayList<Sample> parents;
        for(Map.Entry<String,Set<Sample>> familyEntry : families.entrySet()){
            family = familyEntry.getValue();

            // Since getFamilies(vcfSamples) above still returns parents of samples in the VCF even if those parents are not in the VCF, need to subset down here:
            Set<Sample> familyMembersInVCF = new TreeSet<Sample>();
            for(Sample familyMember : family){
                if (vcfSamples.contains(familyMember.getID())) {
                    familyMembersInVCF.add(familyMember);
                }
            }
            family = familyMembersInVCF;

            if(family.size()<2 || family.size()>3){
                logger.info(String.format("Caution: Family %s has %d members; At the moment Phase By Transmission only supports trios and parent/child pairs. Family skipped.",familyEntry.getKey(),family.size()));
            }
            else{
                for(Sample familyMember : family){
                    parents = familyMember.getParents();
                    if(parents.size()>0){
                        if(family.containsAll(parents))
                            this.trios.add(familyMember);
                        else
                            logger.info(String.format("Caution: Child %s of family %s skipped as info is not provided as a complete trio nor a parent/child pair; At the moment Phase By Transmission only supports trios and parent/child pairs. Child skipped.", familyMember.getID(), familyEntry.getKey()));
                    }
                }
            }

        }



    }

    //Create the transmission matrices
    private void buildMatrices(){
        mvCountMatrix = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>>(GenotypeType.class);
        transmissionMatrix = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,TrioPhase>>>(GenotypeType.class);
        for(GenotypeType mother : GenotypeType.values()){
            mvCountMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>(GenotypeType.class));
            transmissionMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,TrioPhase>>(GenotypeType.class));
            for(GenotypeType father : GenotypeType.values()){
                mvCountMatrix.get(mother).put(father,new EnumMap<GenotypeType, Integer>(GenotypeType.class));
                transmissionMatrix.get(mother).put(father,new EnumMap<GenotypeType,TrioPhase>(GenotypeType.class));
                for(GenotypeType child : GenotypeType.values()){
                    mvCountMatrix.get(mother).get(father).put(child, getCombinationMVCount(mother, father, child));
                    transmissionMatrix.get(mother).get(father).put(child,new TrioPhase(mother,father,child));
                }
            }
        }
    }

    //Returns the number of Mendelian Violations for a given genotype combination.
    //If one of the parents genotype is missing, it will consider it as a parent/child pair
    //If the child genotype or both parents genotypes are missing, 0 is returned.
    private int getCombinationMVCount(GenotypeType mother, GenotypeType father, GenotypeType child){

        //Child is no call => No MV
        if(child == GenotypeType.NO_CALL || child == GenotypeType.UNAVAILABLE)
            return 0;
        //Add parents with genotypes for the evaluation
        ArrayList<GenotypeType> parents = new ArrayList<GenotypeType>();
        if (!(mother == GenotypeType.NO_CALL || mother == GenotypeType.UNAVAILABLE))
            parents.add(mother);
        if (!(father == GenotypeType.NO_CALL || father == GenotypeType.UNAVAILABLE))
            parents.add(father);

        //Both parents no calls => No MV
        if (parents.isEmpty())
            return 0;

        //If at least one parent had a genotype, then count the number of ref and alt alleles that can be passed
        int parentsNumRefAlleles = 0;
        int parentsNumAltAlleles = 0;

        for(GenotypeType parent : parents){
            if(parent == GenotypeType.HOM_REF){
                parentsNumRefAlleles++;
            }
            else if(parent == GenotypeType.HET){
                parentsNumRefAlleles++;
                parentsNumAltAlleles++;
            }
            else if(parent == GenotypeType.HOM_VAR){
                parentsNumAltAlleles++;
            }
        }

        //Case Child is HomRef
        if(child == GenotypeType.HOM_REF){
            if(parentsNumRefAlleles == parents.size())
                return 0;
            else return (parents.size()-parentsNumRefAlleles);
        }

        //Case child is HomVar
        if(child == GenotypeType.HOM_VAR){
            if(parentsNumAltAlleles == parents.size())
                return 0;
            else return parents.size()-parentsNumAltAlleles;
        }

        //Case child is Het
        if(child == GenotypeType.HET && ((parentsNumRefAlleles > 0 && parentsNumAltAlleles > 0) || parents.size()<2))
            return 0;

        //MV
        return 1;
    }

    //Given two trio genotypes combinations, returns the number of different genotypes between the two combinations.
    private int countFamilyGenotypeDiff(GenotypeType motherOriginal,GenotypeType fatherOriginal,GenotypeType childOriginal,GenotypeType motherNew,GenotypeType fatherNew,GenotypeType childNew){
        int count = 0;
        if(motherOriginal!=motherNew)
            count++;
        if(fatherOriginal!=fatherNew)
            count++;
        if(childOriginal!=childNew)
            count++;
        return count;
    }

    //Get a Map of genotype likelihoods.
    //In case of null, unavailable or no call, all likelihoods are 1/3.
    private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(Genotype genotype){
        if(genotype == null || !genotype.isCalled() || genotype.getLikelihoods() == null){
            EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
            likelihoods.put(GenotypeType.HOM_REF,1.0/3.0);
            likelihoods.put(GenotypeType.HET,1.0/3.0);
            likelihoods.put(GenotypeType.HOM_VAR,1.0/3.0);
            return likelihoods;
        }
        return genotype.getLikelihoods().getAsMap(true);
    }

    //Returns the GenotypeType; returns UNVAILABLE if given null
    private GenotypeType getTypeSafeNull(Genotype genotype){
        if(genotype == null)
            return GenotypeType.UNAVAILABLE;
        return genotype.getType();
    }


    /**
     * Phases the genotypes of the given trio. If one of the parents is null, it is considered a parent/child pair.
     * @param ref: Reference allele
     * @param alt: Alternative allele
     * @param mother: Mother's genotype
     * @param father: Father's genotype
     * @param child: Child's genotype
     * @param finalGenotypes: An ArrayList<Genotype> that will be added the genotypes phased by transmission in the following order: Mother, Father, Child
     * @return
     */
    private int phaseTrioGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child,ArrayList<Genotype> finalGenotypes) {

        //Check whether it is  a pair or trio
        //Always assign the first parent as the parent having genotype information in pairs
        //Always assign the mother as the first parent in trios
        int parentsCalled = 0;
        Map<GenotypeType,Double> firstParentLikelihoods;
        Map<GenotypeType,Double> secondParentLikelihoods;
        ArrayList<GenotypeType> bestFirstParentGenotype = new ArrayList<GenotypeType>();
        ArrayList<GenotypeType> bestSecondParentGenotype = new ArrayList<GenotypeType>();
        ArrayList<GenotypeType> bestChildGenotype = new ArrayList<GenotypeType>();
        GenotypeType pairSecondParentGenotype = null;
        if(mother == null || !mother.isCalled()){
            firstParentLikelihoods = getLikelihoodsAsMapSafeNull(father);
            secondParentLikelihoods = getLikelihoodsAsMapSafeNull(mother);
            bestFirstParentGenotype.add(getTypeSafeNull(father));
            bestSecondParentGenotype.add(getTypeSafeNull(mother));
            pairSecondParentGenotype = mother == null ? GenotypeType.UNAVAILABLE : mother.getType();
            if(father != null && father.isCalled())
                parentsCalled = 1;
        }
        else{
            firstParentLikelihoods = getLikelihoodsAsMapSafeNull(mother);
            secondParentLikelihoods = getLikelihoodsAsMapSafeNull(father);
            bestFirstParentGenotype.add(getTypeSafeNull(mother));
            bestSecondParentGenotype.add(getTypeSafeNull(father));
            if(father == null || !father.isCalled()){
                parentsCalled = 1;
                pairSecondParentGenotype = father == null ? GenotypeType.UNAVAILABLE : father.getType();
            }else{
                parentsCalled = 2;
            }
        }
        Map<GenotypeType,Double> childLikelihoods = getLikelihoodsAsMapSafeNull(child);
        bestChildGenotype.add(getTypeSafeNull(child));

        //Prior vars
        double bestConfigurationLikelihood = 0.0;
        double norm = 0.0;
        int configuration_index =0;
        ArrayList<Integer> bestMVCount = new ArrayList<Integer>();
        bestMVCount.add(0);

        //Get the most likely combination
        //Only check for most likely combination if at least a parent and the child have genotypes
        if(child.isCalled() && parentsCalled > 0){
            int mvCount;
            int cumulativeMVCount = 0;
            double configurationLikelihood = 0;
            for(Map.Entry<GenotypeType,Double> childGenotype : childLikelihoods.entrySet()){
                for(Map.Entry<GenotypeType,Double> firstParentGenotype : firstParentLikelihoods.entrySet()){
                    for(Map.Entry<GenotypeType,Double> secondParentGenotype : secondParentLikelihoods.entrySet()){
                        mvCount = mvCountMatrix.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(childGenotype.getKey());
                        //For parent/child pairs, sum over the possible genotype configurations of the missing parent
                        if(parentsCalled<2){
                            cumulativeMVCount += mvCount;
                            configurationLikelihood += mvCount>0 ? Math.pow(deNovoPrior,mvCount)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue() : (1.0-11*deNovoPrior)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue();
                        }
                        //Evaluate configurations of trios
                        else{
                            configurationLikelihood =  mvCount>0 ? Math.pow(deNovoPrior,mvCount)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue() : (1.0-11*deNovoPrior)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue();
                            norm += configurationLikelihood;
                            //Keep this combination if
                            //It has a better likelihood
                            //Or it has the same likelihood but requires less changes from original genotypes
                            if (configurationLikelihood > bestConfigurationLikelihood){
                                bestConfigurationLikelihood = configurationLikelihood;
                                bestMVCount.clear();
                                bestMVCount.add(mvCount);
                                bestFirstParentGenotype.clear();
                                bestFirstParentGenotype.add(firstParentGenotype.getKey());
                                bestSecondParentGenotype.clear();
                                bestSecondParentGenotype.add(secondParentGenotype.getKey());
                                bestChildGenotype.clear();
                                bestChildGenotype.add(childGenotype.getKey());
                            }
                            else if(configurationLikelihood == bestConfigurationLikelihood) {
                                bestFirstParentGenotype.add(firstParentGenotype.getKey());
                                bestSecondParentGenotype.add(secondParentGenotype.getKey());
                                bestChildGenotype.add(childGenotype.getKey());
                                bestMVCount.add(mvCount);
                            }
                        }
                    }
                    //Evaluate configurations of parent/child pairs
                    if(parentsCalled<2){
                        norm += configurationLikelihood;
                        //Keep this combination if
                        //It has a better likelihood
                        //Or it has the same likelihood but requires less changes from original genotypes
                        if (configurationLikelihood > bestConfigurationLikelihood){
                            bestConfigurationLikelihood = configurationLikelihood;
                            bestMVCount.clear();
                            bestMVCount.add(cumulativeMVCount/3);
                            bestChildGenotype.clear();
                            bestFirstParentGenotype.clear();
                            bestSecondParentGenotype.clear();
                            bestChildGenotype.add(childGenotype.getKey());
                            bestFirstParentGenotype.add(firstParentGenotype.getKey());
                            bestSecondParentGenotype.add(pairSecondParentGenotype);
                        }
                        else if(configurationLikelihood == bestConfigurationLikelihood) {
                            bestFirstParentGenotype.add(firstParentGenotype.getKey());
                            bestSecondParentGenotype.add(pairSecondParentGenotype);
                            bestChildGenotype.add(childGenotype.getKey());
                            bestMVCount.add(cumulativeMVCount/3);
                        }
                        configurationLikelihood = 0;
                    }
                }
            }

            //normalize the best configuration probability
            bestConfigurationLikelihood = bestConfigurationLikelihood / norm;

            //In case of multiple equally likely combinations, take a random one
            if(bestFirstParentGenotype.size()>1){
                configuration_index = rand.nextInt(bestFirstParentGenotype.size()-1);
            }

        }
        else{
            bestConfigurationLikelihood = NO_TRANSMISSION_PROB;
        }

        TrioPhase phasedTrioGenotypes;
        if(parentsCalled < 2 && mother == null || !mother.isCalled())
            phasedTrioGenotypes = transmissionMatrix.get(bestSecondParentGenotype.get(configuration_index)).get(bestFirstParentGenotype.get(configuration_index)).get(bestChildGenotype.get(configuration_index));
        else
            phasedTrioGenotypes = transmissionMatrix.get(bestFirstParentGenotype.get(configuration_index)).get(bestSecondParentGenotype.get(configuration_index)).get(bestChildGenotype.get(configuration_index));

        //Return the phased genotypes
        phasedTrioGenotypes.getPhasedGenotypes(ref,alt,mother,father,child,bestConfigurationLikelihood,finalGenotypes);
        return bestMVCount.get(configuration_index);

    }


    private void updatePairMetricsCounters(Genotype parent, Genotype child, int mvCount, HashMap<Byte,Integer> counters){

        //Increment metrics counters
        if(parent.isCalled() && child.isCalled()){
            counters.put(NUM_PAIR_GENOTYPES_CALLED,counters.get(NUM_PAIR_GENOTYPES_CALLED)+1);
            if(parent.isPhased())
                counters.put(NUM_PAIR_GENOTYPES_PHASED,counters.get(NUM_PAIR_GENOTYPES_PHASED)+1);
            else{
                counters.put(NUM_PAIR_VIOLATIONS,counters.get(NUM_PAIR_VIOLATIONS)+mvCount);
                if(parent.isHet() && child.isHet())
                    counters.put(NUM_PAIR_HET_HET,counters.get(NUM_PAIR_HET_HET)+1);
            }
        }else{
            counters.put(NUM_PAIR_GENOTYPES_NOCALL,counters.get(NUM_PAIR_GENOTYPES_NOCALL)+1);
        }

    }

    private void updateTrioMetricsCounters(Genotype mother, Genotype father, Genotype child, int mvCount, HashMap<Byte,Integer> counters){

        //Increment metrics counters
        if(mother.isCalled() && father.isCalled() && child.isCalled()){
            counters.put(NUM_TRIO_GENOTYPES_CALLED,counters.get(NUM_TRIO_GENOTYPES_CALLED)+1);
            if(mother.isPhased())
                counters.put(NUM_TRIO_GENOTYPES_PHASED,counters.get(NUM_TRIO_GENOTYPES_PHASED)+1);

            else{
                if(mvCount > 0){
                    if(mvCount >1)
                        counters.put(NUM_TRIO_DOUBLE_VIOLATIONS,counters.get(NUM_TRIO_DOUBLE_VIOLATIONS)+1);
                    else
                        counters.put(NUM_TRIO_VIOLATIONS,counters.get(NUM_TRIO_VIOLATIONS)+1);
                }
                else if(mother.isHet() && father.isHet() && child.isHet())
                    counters.put(NUM_TRIO_HET_HET_HET,counters.get(NUM_TRIO_HET_HET_HET)+1);

            }
        }else{
            counters.put(NUM_TRIO_GENOTYPES_NOCALL,counters.get(NUM_TRIO_GENOTYPES_NOCALL)+1);
        }
    }

    /**
     * For each variant in the file, determine the phasing for the child and replace the child's genotype with the trio's genotype
     *
     * @param tracker  the reference meta-data tracker
     * @param ref      the reference context
     * @param context  the alignment context
     * @return null
     */
    @Override
    public HashMap<Byte,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        HashMap<Byte,Integer> metricsCounters = new HashMap<Byte, Integer>(10);
        metricsCounters.put(NUM_TRIO_GENOTYPES_CALLED,0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_NOCALL,0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_PHASED,0);
        metricsCounters.put(NUM_TRIO_HET_HET_HET,0);
        metricsCounters.put(NUM_TRIO_VIOLATIONS,0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_CALLED,0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_NOCALL,0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_PHASED,0);
        metricsCounters.put(NUM_PAIR_HET_HET,0);
        metricsCounters.put(NUM_PAIR_VIOLATIONS,0);
        metricsCounters.put(NUM_TRIO_DOUBLE_VIOLATIONS,0);
        metricsCounters.put(NUM_GENOTYPES_MODIFIED,0);

        String mvfLine;

        if (tracker == null)
            return metricsCounters;

        final VariantContext vc = tracker.getFirstValue(variantCollection.variants, context.getLocation());
        if ( vc == null )
            return metricsCounters;

        if ( !vc.isBiallelic() ) {
            vcfWriter.add(vc);
            return metricsCounters;
        }

        final VariantContextBuilder builder = new VariantContextBuilder(vc);

        final GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());
        for (Sample sample : trios) {
            Genotype mother = vc.getGenotype(sample.getMaternalID());
            Genotype father = vc.getGenotype(sample.getPaternalID());
            Genotype child = vc.getGenotype(sample.getID());

            //Keep only trios and parent/child pairs
            if(mother == null && father == null || child == null)
                continue;

            ArrayList<Genotype> trioGenotypes = new ArrayList<Genotype>(3);
            final int mvCount = phaseTrioGenotypes(vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), mother, father, child,trioGenotypes);

            Genotype phasedMother = trioGenotypes.get(0);
            Genotype phasedFather = trioGenotypes.get(1);
            Genotype phasedChild = trioGenotypes.get(2);

            //Fill the genotype map with the new genotypes and increment metrics counters
            genotypesContext.replace(phasedChild);
            if(mother != null){
                genotypesContext.replace(phasedMother);
                if(father != null){
                    genotypesContext.replace(phasedFather);
                    updateTrioMetricsCounters(phasedMother,phasedFather,phasedChild,mvCount,metricsCounters);
                    mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                            vc.getChr(),vc.getStart(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),sample.getFamilyID(),
                            phasedMother.getExtendedAttribute(GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY),phasedMother.getGenotypeString(),phasedMother.getDP(),printAD(phasedMother.getAD()),
                            phasedMother.getLikelihoodsString(), phasedFather.getGenotypeString(),phasedFather.getDP(),printAD(phasedFather.getAD()),phasedFather.getLikelihoodsString(),
                            phasedChild.getGenotypeString(),phasedChild.getDP(),printAD(phasedChild.getAD()),phasedChild.getLikelihoodsString());
                    if(!(phasedMother.getType()==mother.getType() && phasedFather.getType()==father.getType() && phasedChild.getType()==child.getType()))
                        metricsCounters.put(NUM_GENOTYPES_MODIFIED,metricsCounters.get(NUM_GENOTYPES_MODIFIED)+1);
                }
                else{
                    updatePairMetricsCounters(phasedMother,phasedChild,mvCount,metricsCounters);
                    if(!(phasedMother.getType()==mother.getType() && phasedChild.getType()==child.getType()))
                        metricsCounters.put(NUM_GENOTYPES_MODIFIED,metricsCounters.get(NUM_GENOTYPES_MODIFIED)+1);
                    mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s:%s:%s:%s\t.\t.\t.\t.\t%s\t%s\t%s\t%s",
                            vc.getChr(),vc.getStart(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),sample.getFamilyID(),
                            phasedMother.getExtendedAttribute(GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY),phasedMother.getGenotypeString(),phasedMother.getDP(),printAD(phasedMother.getAD()),phasedMother.getLikelihoodsString(),
                            phasedChild.getGenotypeString(),phasedChild.getDP(),printAD(phasedChild.getAD()),phasedChild.getLikelihoodsString());
                }
            }
            else{
                genotypesContext.replace(phasedFather);
                updatePairMetricsCounters(phasedFather,phasedChild,mvCount,metricsCounters);
                if(!(phasedFather.getType()==father.getType() && phasedChild.getType()==child.getType()))
                    metricsCounters.put(NUM_GENOTYPES_MODIFIED,metricsCounters.get(NUM_GENOTYPES_MODIFIED)+1);
                mvfLine =   String.format("%s\t%d\t%s\t%s\t%s\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                        vc.getChr(),vc.getStart(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),sample.getFamilyID(),
                        phasedFather.getExtendedAttribute(GATKVCFConstants.TRANSMISSION_PROBABILITY_KEY),phasedFather.getGenotypeString(),phasedFather.getDP(),printAD(phasedFather.getAD()),phasedFather.getLikelihoodsString(),
                        phasedChild.getGenotypeString(),phasedChild.getDP(),printAD(phasedChild.getAD()),phasedChild.getLikelihoodsString());
            }

            //Report violation if set so
            //TODO: ADAPT FOR PAIRS TOO!!
            if(mvCount>0 && mvFile != null && !vc.isFiltered())
                mvFile.println(mvfLine);
        }

        builder.genotypes(genotypesContext);
        vcfWriter.add(builder.make());

        return metricsCounters;
    }

    private static String printAD(final int[] AD) {
        if ( AD == null || AD.length == 0 )
            return ".";
        final StringBuilder sb = new StringBuilder();
        sb.append(AD[0]);
        for ( int i = 1; i < AD.length; i++) {
            sb.append(",");
            sb.append(AD[i]);
        }
        return sb.toString();
    }

    /**
     * Initializes the reporting counters.
     *
     * @return All counters initialized to 0
     */
    @Override
    public HashMap<Byte,Integer> reduceInit() {
        HashMap<Byte,Integer> metricsCounters = new HashMap<Byte, Integer>(10);
        metricsCounters.put(NUM_TRIO_GENOTYPES_CALLED,0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_NOCALL,0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_PHASED,0);
        metricsCounters.put(NUM_TRIO_HET_HET_HET,0);
        metricsCounters.put(NUM_TRIO_VIOLATIONS,0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_CALLED,0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_NOCALL,0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_PHASED,0);
        metricsCounters.put(NUM_PAIR_HET_HET,0);
        metricsCounters.put(NUM_PAIR_VIOLATIONS,0);
        metricsCounters.put(NUM_TRIO_DOUBLE_VIOLATIONS,0);
        metricsCounters.put(NUM_GENOTYPES_MODIFIED,0);

        return metricsCounters;
    }

    /**
     * Adds the value of the site phased to the reporting counters.
     *
     * @param value Site values
     * @param sum   accumulator for the reporting counters
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public HashMap<Byte,Integer> reduce(HashMap<Byte,Integer> value, HashMap<Byte,Integer> sum) {
        sum.put(NUM_TRIO_GENOTYPES_CALLED,value.get(NUM_TRIO_GENOTYPES_CALLED)+sum.get(NUM_TRIO_GENOTYPES_CALLED));
        sum.put(NUM_TRIO_GENOTYPES_NOCALL,value.get(NUM_TRIO_GENOTYPES_NOCALL)+sum.get(NUM_TRIO_GENOTYPES_NOCALL));
        sum.put(NUM_TRIO_GENOTYPES_PHASED,value.get(NUM_TRIO_GENOTYPES_PHASED)+sum.get(NUM_TRIO_GENOTYPES_PHASED));
        sum.put(NUM_TRIO_HET_HET_HET,value.get(NUM_TRIO_HET_HET_HET)+sum.get(NUM_TRIO_HET_HET_HET));
        sum.put(NUM_TRIO_VIOLATIONS,value.get(NUM_TRIO_VIOLATIONS)+sum.get(NUM_TRIO_VIOLATIONS));
        sum.put(NUM_PAIR_GENOTYPES_CALLED,value.get(NUM_PAIR_GENOTYPES_CALLED)+sum.get(NUM_PAIR_GENOTYPES_CALLED));
        sum.put(NUM_PAIR_GENOTYPES_NOCALL,value.get(NUM_PAIR_GENOTYPES_NOCALL)+sum.get(NUM_PAIR_GENOTYPES_NOCALL));
        sum.put(NUM_PAIR_GENOTYPES_PHASED,value.get(NUM_PAIR_GENOTYPES_PHASED)+sum.get(NUM_PAIR_GENOTYPES_PHASED));
        sum.put(NUM_PAIR_HET_HET,value.get(NUM_PAIR_HET_HET)+sum.get(NUM_PAIR_HET_HET));
        sum.put(NUM_PAIR_VIOLATIONS,value.get(NUM_PAIR_VIOLATIONS)+sum.get(NUM_PAIR_VIOLATIONS));
        sum.put(NUM_TRIO_DOUBLE_VIOLATIONS,value.get(NUM_TRIO_DOUBLE_VIOLATIONS)+sum.get(NUM_TRIO_DOUBLE_VIOLATIONS));
        sum.put(NUM_GENOTYPES_MODIFIED,value.get(NUM_GENOTYPES_MODIFIED)+sum.get(NUM_GENOTYPES_MODIFIED));

        return sum;
    }


    /**
     * Reports statistics on the phasing by transmission process.
     * @param result Accumulator with all counters.
     */
    @Override
    public void onTraversalDone(HashMap<Byte,Integer> result) {
        logger.info("Number of complete trio-genotypes: " + result.get(NUM_TRIO_GENOTYPES_CALLED));
        logger.info("Number of trio-genotypes containing no call(s): " + result.get(NUM_TRIO_GENOTYPES_NOCALL));
        logger.info("Number of trio-genotypes phased: " + result.get(NUM_TRIO_GENOTYPES_PHASED));
        logger.info("Number of resulting Het/Het/Het trios: " + result.get(NUM_TRIO_HET_HET_HET));
        logger.info("Number of remaining single mendelian violations in trios: " + result.get(NUM_TRIO_VIOLATIONS));
        logger.info("Number of remaining double mendelian violations in trios: " + result.get(NUM_TRIO_DOUBLE_VIOLATIONS));
        logger.info("Number of complete pair-genotypes: " + result.get(NUM_PAIR_GENOTYPES_CALLED));
        logger.info("Number of pair-genotypes containing no call(s): " + result.get(NUM_PAIR_GENOTYPES_NOCALL));
        logger.info("Number of pair-genotypes phased: " + result.get(NUM_PAIR_GENOTYPES_PHASED));
        logger.info("Number of resulting Het/Het pairs: " + result.get(NUM_PAIR_HET_HET));
        logger.info("Number of remaining mendelian violations in pairs: " + result.get(NUM_PAIR_VIOLATIONS));
        logger.info("Number of genotypes updated: " + result.get(NUM_GENOTYPES_MODIFIED));

    }
}
