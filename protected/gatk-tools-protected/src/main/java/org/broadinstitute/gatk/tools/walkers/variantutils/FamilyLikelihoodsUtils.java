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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.*;

import java.util.*;

/**
 * FamilyLikelihoodsUtils code is based on PhaseByTransmission with added posterior probability calculations
 */

public class FamilyLikelihoodsUtils {
    private static Logger logger = Logger.getLogger(FamilyLikelihoodsUtils.class);

    //Matrix of priors for all genotype combinations
    final private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> mvCountMatrix =
            new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>>(GenotypeType.class);

    //Matrix of allele transmission
    final private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,TrioGenotypes>>> transmissionMatrix =
             new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,TrioGenotypes>>>(GenotypeType.class);

    final double[] configurationLikelihoodsMatrix = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //27 is # of trio genotype combos, initialize to zero

    ArrayList<Sample> trios = new ArrayList<Sample>();

    //Random number generator
    final private Random rand = new GenomeAnalysisEngine().getRandomGenerator();

    private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
    private final String PHRED_SCALED_POSTERIORS_KEY = "PP";

    public final double NO_TRANSMISSION_PROB = -1.0;

    private double deNovoPrior = 1e-8;

    private enum FamilyMember {
        MOTHER,
        FATHER,
        CHILD
    }

    //Stores a conceptual trio or parent/child pair genotype combination
    //This combination can then be "applied" to a given trio or pair using the getUpdatedGenotypes method.
    private class TrioGenotypes {

        //Create 2 fake alleles
        //The actual bases will never be used but the Genotypes created using the alleles will be.
        private final Allele REF = Allele.create("A",true);
        private final Allele VAR = Allele.create("A",false);
        private final Allele NO_CALL = Allele.create(".",false);
        private final String DUMMY_NAME = "DummySample";
        private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> MVcountMatrix;

        private final EnumMap<FamilyMember,Genotype> familyGenotypes = new EnumMap<FamilyMember, Genotype>(FamilyMember.class);

        /*  Constructor: Creates a conceptual trio genotype combination from the given genotypes.
        */
        public TrioGenotypes(GenotypeType mother, GenotypeType father, GenotypeType child){
            familyGenotypes.put(FamilyMember.MOTHER, makeGenotype(mother));
            familyGenotypes.put(FamilyMember.FATHER, makeGenotype(father));
            familyGenotypes.put(FamilyMember.CHILD, makeGenotype(child));
        }

        private ArrayList<Allele> getAlleles(GenotypeType genotype){
            final ArrayList<Allele> alleles = new ArrayList<Allele>(2);
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

        public void setMVcountMatrix(EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> inputMat) {
            MVcountMatrix = inputMat;
        }

        private boolean hasCalledGT(GenotypeType genotype){
            return genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HET || genotype == GenotypeType.HOM_VAR;
        }

        //TODO: this was stupid stuff for phasing -- let's take it out
        private Genotype makeGenotype(final GenotypeType type) {
            return makeGenotype(getAlleles(type));
        }

        private Genotype makeGenotype(final List<Allele> alleles) {
            final GenotypeBuilder gb = new GenotypeBuilder(DUMMY_NAME, alleles);
            return gb.make();
        }

        /**
         * Applies the trio genotype combination to the given trio.
         * @param ref: Reference allele
         * @param alt: Alternate allele
         * @param motherGenotype: Genotype of the mother in this trio genotype combination
         * @param fatherGenotype: Genotype of the father in this trio genotype combination
         * @param childGenotype: Genotype of the child in this trio genotype combination
         * @param transmissionProb: Probability for this trio genotype combination to be correct (pass NO_TRANSMISSION_PROB if unavailable)
         * @param configurationLikelihoodsMatrix: (Non-normalized) likelihoods for each trio genotype combination, for use in generating new PLs
         * @param updatedGenotypes: An ArrayList<Genotype> to which the newly updated genotypes are added in the following order: Mother, Father, Child
         */
        public void getUpdatedGenotypes(final Allele ref, final Allele alt, final Genotype motherGenotype, final Genotype fatherGenotype, final Genotype childGenotype, final double transmissionProb, double[] configurationLikelihoodsMatrix, final ArrayList<Genotype> updatedGenotypes){
            //default to flat priors in case input genotypes are not called
            double[] motherPosteriors = {1,1,1};
            double[] fatherPosteriors = {1,1,1};
            double[] childPosteriors = {1,1,1};

            //genotypes here can be no call
            boolean fatherIsCalled = fatherGenotype != null && hasCalledGT(fatherGenotype.getType());
            boolean motherIsCalled = motherGenotype != null && hasCalledGT(motherGenotype.getType());
            boolean childIsCalled = childGenotype != null && hasCalledGT(childGenotype.getType());

            if (fatherIsCalled && childIsCalled) {
                motherPosteriors = getPosteriors(FamilyMember.MOTHER);
            }
            updatedGenotypes.add(getUpdatedGenotype(ref, alt, motherGenotype, transmissionProb, motherPosteriors));

            if (motherIsCalled && childIsCalled) {
                fatherPosteriors = getPosteriors(FamilyMember.FATHER);
            }
            updatedGenotypes.add(getUpdatedGenotype(ref, alt, fatherGenotype, transmissionProb, fatherPosteriors));

            if (motherIsCalled && fatherIsCalled) {
                childPosteriors = getPosteriors(FamilyMember.CHILD);
            }
            updatedGenotypes.add(getUpdatedGenotype(ref, alt, childGenotype, transmissionProb, childPosteriors));
        }

        private Genotype getUpdatedGenotype(Allele refAllele, Allele altAllele, Genotype genotype, double transmissionProb, double[] normalizedPosteriors){

            int phredScoreTransmission = -1;
            if(transmissionProb != NO_TRANSMISSION_PROB){
                double dphredScoreTransmission = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1 - (transmissionProb)));
                phredScoreTransmission = dphredScoreTransmission < Byte.MAX_VALUE ? (byte)dphredScoreTransmission : Byte.MAX_VALUE;
            }
            //Handle null, missing and unavailable genotypes
            //Note that only cases where a null/missing/unavailable genotype was passed in the first place can lead to a null/missing/unavailable
            //genotype so it is safe to return the original genotype in this case.
            //In addition, if the genotype configuration confidence is 0, then return the original genotypes.
            if(phredScoreTransmission ==0 || genotype == null || !hasCalledGT(genotype.getType()))
                return genotype;

            //Add the transmission probability
            final Map<String, Object> genotypeAttributes = new HashMap<String, Object>();
            genotypeAttributes.putAll(genotype.getExtendedAttributes());
            if(transmissionProb>NO_TRANSMISSION_PROB)
                genotypeAttributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, phredScoreTransmission);

            final ArrayList<Allele> usedAlleles = new ArrayList<Allele>(2);
            usedAlleles.add(refAllele);
            usedAlleles.add(altAllele);

            final GenotypeBuilder builder = new GenotypeBuilder(genotype);

            final double[] log10Posteriors = MathUtils.toLog10(normalizedPosteriors);

            //note that there will there be times when posteriors don't agree with genotype predicted by configuration likelihoods
            GATKVariantContextUtils.updateGenotypeAfterSubsetting(usedAlleles, builder,
                    GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, log10Posteriors, usedAlleles);



            builder.attribute(PHRED_SCALED_POSTERIORS_KEY,
                    Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(log10Posteriors).getAsPLs()));
            builder.attributes(genotypeAttributes);
            return builder.make();
        }

        //marginalize over the configurationLikelihoodsMatrix and normalize to get the posteriors
        private double[] getPosteriors(FamilyMember recalcInd) {
            double marginalOverChangedHR, marginalOverChangedHET, marginalOverChangedHV;
            marginalOverChangedHR = marginalOverChangedHET = marginalOverChangedHV = 0;
            final double[] recalcPosteriors = new double[3];

            GenotypeType[] calledTypes = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};

            switch (recalcInd) {
                case MOTHER:
                    for(final GenotypeType father : calledTypes) {
                        for(final GenotypeType child : calledTypes) {
                            GenotypeType mother;
                            mother = GenotypeType.HOM_REF;
                            marginalOverChangedHR += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                            mother = GenotypeType.HET;
                            marginalOverChangedHET += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                            mother = GenotypeType.HOM_VAR;
                            marginalOverChangedHV += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                        }
                    }
                    break;
                case FATHER:
                    for(final GenotypeType mother : calledTypes){
                        for (final GenotypeType child : calledTypes){
                            GenotypeType father;
                            father = GenotypeType.HOM_REF;
                            marginalOverChangedHR += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                            father = GenotypeType.HET;
                            marginalOverChangedHET += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                            father = GenotypeType.HOM_VAR;
                            marginalOverChangedHV += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                        }
                    }
                    break;
                case CHILD:
                    for(final GenotypeType mother : calledTypes){
                        for (final GenotypeType father: calledTypes){
                            GenotypeType child;
                            child = GenotypeType.HOM_REF;
                            marginalOverChangedHR += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                            child = GenotypeType.HET;
                            marginalOverChangedHET += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                            child = GenotypeType.HOM_VAR;
                            marginalOverChangedHV += configurationLikelihoodsMatrix[getLikelihoodIndex(mother,father,child,false)];
                        }
                    }
                    break;
                default:
                    throw new UserException(String.format("%d does not indicate a valid trio individual -- use 0 for mother, 1 for father, 2 for child",recalcInd));
            }
            recalcPosteriors[0] = marginalOverChangedHR;
            recalcPosteriors[1] = marginalOverChangedHET;
            recalcPosteriors[2] = marginalOverChangedHV;

            final double[] normalizedPosteriors = MathUtils.normalizeFromRealSpace(recalcPosteriors);

            return normalizedPosteriors;
        }

    }

    public void initialize(double DNprior, Set<String> vcfSamples, Map<String,Set<Sample>> families){
        this.deNovoPrior = DNprior;
        buildMatrices();
        trios = setTrios(vcfSamples, families);
    }

    public GenotypesContext calculatePosteriorGLs(VariantContext vc){
        final GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());

        for (final Sample sample : trios) {
            Genotype mother = vc.getGenotype(sample.getMaternalID());
            Genotype father = vc.getGenotype(sample.getPaternalID());
            Genotype child = vc.getGenotype(sample.getID());

            //Keep only trios and parent/child pairs
            if(mother == null && father == null || child == null) {
                logger.warn("null genos in var "+vc.toStringDecodeGenotypes());
                continue;
            }

            final ArrayList<Genotype> trioGenotypes = new ArrayList<Genotype>(3);
            final int mvCount = updateFamilyGenotypes(vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), mother, father, child, trioGenotypes);

            Genotype updatedMother = trioGenotypes.get(0);
            Genotype updatedFather = trioGenotypes.get(1);
            Genotype updatedChild = trioGenotypes.get(2);

            genotypesContext.replace(updatedChild);
            genotypesContext.replace(updatedFather);
            genotypesContext.replace(updatedMother);
        }

    return genotypesContext;
    }

    /**
     * Select trios and parent/child pairs only
     */
    private ArrayList<Sample> setTrios(Set<String> vcfSamples, Map<String,Set<Sample>> families){

        Set<Sample> family;
        ArrayList<Sample> parents;
        final ArrayList<Sample> trios = new ArrayList<Sample>();
        for(final Map.Entry<String,Set<Sample>> familyEntry : families.entrySet()){
            family = familyEntry.getValue();

            // Since getFamilies(vcfSamples) above still returns parents of samples in the VCF even if those parents are not in the VCF, need to subset down here:
            final Set<Sample> familyMembersInVCF = new TreeSet<Sample>();
            for(final Sample familyMember : family){
                if (vcfSamples.contains(familyMember.getID())) {
                    familyMembersInVCF.add(familyMember);
                }
            }
            family = familyMembersInVCF;

            if(family.size() == 3){
                for(final Sample familyMember : family){
                    parents = familyMember.getParents();
                    if(parents.size()>0){
                        if(family.containsAll(parents))
                            trios.add(familyMember);
                    }
                }
            }

        }
        return trios;
    }

    //Create the transmission matrices
    //TODO: pass in the real genotypes so we have that info
    private void buildMatrices(){
        for(final GenotypeType mother : GenotypeType.values()){
            mvCountMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>(GenotypeType.class));
            transmissionMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,TrioGenotypes>>(GenotypeType.class));
            for(final GenotypeType father : GenotypeType.values()){
                mvCountMatrix.get(mother).put(father,new EnumMap<GenotypeType, Integer>(GenotypeType.class));
                transmissionMatrix.get(mother).put(father,new EnumMap<GenotypeType,TrioGenotypes>(GenotypeType.class));
                for(final GenotypeType child : GenotypeType.values()){
                    mvCountMatrix.get(mother).get(father).put(child, getCombinationMVCount(mother, father, child));
                    transmissionMatrix.get(mother).get(father).put(child,new TrioGenotypes(mother,father,child));
                }
            }
        }
    }

    //Returns the number of Mendelian Violations for a given genotype combination.
    //If one of the parents' genotypes is missing, it will consider it as a parent/child pair
    //If the child genotype or both parents genotypes are missing, 0 is returned.
    private int getCombinationMVCount(GenotypeType mother, GenotypeType father, GenotypeType child){

        //Child is no call => No MV
        if(child == GenotypeType.NO_CALL || child == GenotypeType.UNAVAILABLE)
            return 0;
        //Add parents with genotypes for the evaluation
        final ArrayList<GenotypeType> parents = new ArrayList<GenotypeType>();
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

        for(final GenotypeType parent : parents){
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

    /**
     * Updates the genotypes of the given trio. If one of the parents is null, it is considered a parent/child pair.
     * @param ref: Reference allele
     * @param alt: Alternative allele
     * @param mother: Mother's genotype from vc input
     * @param father: Father's genotype from vc input
     * @param child: Child's genotype from vc input
     * @param finalGenotypes: An ArrayList<Genotype> containing the updated genotypes
     * @return
     */
    private int updateFamilyGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child, ArrayList<Genotype> finalGenotypes) {

        //Check whether it is  a pair or trio
        //Always assign the first parent as the parent having genotype information in pairs
        //Always assign the mother as the first parent in trios
        int parentsCalled = 0;
        Map<GenotypeType,Double> firstParentLikelihoods;
        Map<GenotypeType,Double> secondParentLikelihoods;
        final ArrayList<GenotypeType> bestFirstParentGenotype = new ArrayList<GenotypeType>();
        final ArrayList<GenotypeType> bestSecondParentGenotype = new ArrayList<GenotypeType>();
        final ArrayList<GenotypeType> bestChildGenotype = new ArrayList<GenotypeType>();
        GenotypeType pairSecondParentGenotype = null;
        boolean parentsAreFlipped = false; //usually mother comes first, like for indexing of transmissionMatrix
        final int INVALID_INDEX = -1;

        //if only one parent is called, make uncalled parent the secondParent
        if(mother == null || !mother.isCalled()){
            firstParentLikelihoods = getLikelihoodsAsMapSafeNull(father);
            secondParentLikelihoods = getLikelihoodsAsMapSafeNull(mother);
            bestFirstParentGenotype.add(getTypeSafeNull(father));
            bestSecondParentGenotype.add(getTypeSafeNull(mother));
            pairSecondParentGenotype = mother == null ? GenotypeType.UNAVAILABLE : mother.getType();
            parentsAreFlipped = true;
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
        final ArrayList<Integer> bestMVCount = new ArrayList<Integer>();
        bestMVCount.add(0);

        //Get the most likely combination
        //Only check for most likely combination if at least a parent and the child have genotypes
        int matInd;
        if(child.isCalled() && parentsCalled > 0){
            int mvCount;
            int cumulativeMVCount = 0;
            double configurationLikelihood = 0;
            for(final Map.Entry<GenotypeType,Double> childGenotype :
                    childLikelihoods.entrySet()){
                for(final Map.Entry<GenotypeType,Double> firstParentGenotype :
                        firstParentLikelihoods.entrySet()){
                    for(final Map.Entry<GenotypeType,Double> secondParentGenotype :
                            secondParentLikelihoods.entrySet()){
                        mvCount = mvCountMatrix.get(firstParentGenotype.getKey()).get(secondParentGenotype.getKey()).get(childGenotype.getKey());
                        //For parent/child pairs, sum over the possible genotype configurations of the missing parent
                        if(parentsCalled<2){
                            cumulativeMVCount += mvCount;
                            configurationLikelihood += mvCount>0 ? Math.pow(deNovoPrior,mvCount)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue() : (1.0-10*deNovoPrior-deNovoPrior*deNovoPrior)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue();
                        }
                        //Evaluate configurations of trios
                        else{
                            configurationLikelihood =  mvCount>0 ? Math.pow(deNovoPrior,mvCount)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue() : (1.0-10*deNovoPrior-deNovoPrior*deNovoPrior)*firstParentGenotype.getValue()*secondParentGenotype.getValue()*childGenotype.getValue();
                            norm += configurationLikelihood;
                            matInd = getLikelihoodIndex(firstParentGenotype.getKey(), secondParentGenotype.getKey(),childGenotype.getKey(),parentsAreFlipped);
                            if (matInd > INVALID_INDEX)  //still a slim chance of a MIXED GT
                                configurationLikelihoodsMatrix[matInd] = configurationLikelihood;
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
                        matInd = getLikelihoodIndex(firstParentGenotype.getKey(), GenotypeType.HOM_REF,childGenotype.getKey(), parentsAreFlipped);
                        if (matInd > INVALID_INDEX)
                            configurationLikelihoodsMatrix[matInd] = configurationLikelihood;
                        matInd = getLikelihoodIndex(firstParentGenotype.getKey(), GenotypeType.HET,childGenotype.getKey(),parentsAreFlipped);
                        if (matInd > INVALID_INDEX)
                            configurationLikelihoodsMatrix[matInd] = configurationLikelihood;
                        matInd = getLikelihoodIndex(firstParentGenotype.getKey(), GenotypeType.HOM_VAR,childGenotype.getKey(),parentsAreFlipped);
                        if (matInd > INVALID_INDEX)
                            configurationLikelihoodsMatrix[matInd] = configurationLikelihood;

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

        TrioGenotypes updatedTrioGenotypes;
        if(parentsCalled < 2 && mother == null || !mother.isCalled())
            updatedTrioGenotypes = transmissionMatrix.get(bestSecondParentGenotype.get(configuration_index)).get(bestFirstParentGenotype.get(configuration_index)).get(bestChildGenotype.get(configuration_index));
        else
            updatedTrioGenotypes = transmissionMatrix.get(bestFirstParentGenotype.get(configuration_index)).get(bestSecondParentGenotype.get(configuration_index)).get(bestChildGenotype.get(configuration_index));

        //Return the updated genotypes
        updatedTrioGenotypes.setMVcountMatrix(mvCountMatrix);
        updatedTrioGenotypes.getUpdatedGenotypes(ref, alt, mother, father, child, bestConfigurationLikelihood, configurationLikelihoodsMatrix, finalGenotypes);
        return bestMVCount.get(configuration_index);

    }

    //Get a Map of genotype likelihoods, normalized from log10-space.
    //In case of null, unavailable or no call, all likelihoods are 1/3.
    private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(Genotype genotype){
        if (genotype != null && genotype.isCalled() && genotype.hasExtendedAttribute(PHRED_SCALED_POSTERIORS_KEY)) {
            final EnumMap<GenotypeType,Double> likelihoodsMap = new EnumMap<GenotypeType, Double>(GenotypeType.class);
            Object GPfromVCF = genotype.getExtendedAttribute(PHRED_SCALED_POSTERIORS_KEY);
            //parse the GPs into a vector of probabilities
            final String[] likelihoodsAsStringVector = ((String)GPfromVCF).split(",");
            final double[] likelihoodsAsVector = new double[likelihoodsAsStringVector.length];
            for ( int i = 0; i < likelihoodsAsStringVector.length; i++ ) {
                likelihoodsAsVector[i] = Double.parseDouble(likelihoodsAsStringVector[i]) / -10.0;
            }
            double[] likelihoods = GeneralUtils.normalizeFromLog10(likelihoodsAsVector);
            likelihoodsMap.put(GenotypeType.HOM_REF,likelihoods[GenotypeType.HOM_REF.ordinal()-1]);
            likelihoodsMap.put(GenotypeType.HET,likelihoods[GenotypeType.HET.ordinal()-1]);
            likelihoodsMap.put(GenotypeType.HOM_VAR, likelihoods[GenotypeType.HOM_VAR.ordinal() - 1]);
            return likelihoodsMap;
        }

        if(genotype == null || !genotype.isCalled() || genotype.getLikelihoods() == null){
            final EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
            likelihoods.put(GenotypeType.HOM_REF,1.0/3.0);
            likelihoods.put(GenotypeType.HET,1.0/3.0);
            likelihoods.put(GenotypeType.HOM_VAR,1.0/3.0);
            return likelihoods;
        }
        return genotype.getLikelihoods().getAsMap(true);
    }

    //Returns the GenotypeType; returns UNAVAILABLE if given null
    private GenotypeType getTypeSafeNull(Genotype genotype){
        if(genotype == null)
            return GenotypeType.UNAVAILABLE;
        return genotype.getType();
    }

    private int getLikelihoodIndex(GenotypeType firstParent, GenotypeType secondParent, GenotypeType child, boolean parentsAreFlipped){
        int childInd = genotypeTypeValue(child);
        int motherInd;
        int fatherInd;
        final int NUM_CALLED_GENOTYPETYPES = 3;
        final int INVALID = -1;
        if (parentsAreFlipped)
        {
            motherInd = genotypeTypeValue(secondParent);
            fatherInd = genotypeTypeValue(firstParent);
        }
        else {
            motherInd = genotypeTypeValue(firstParent);
            fatherInd = genotypeTypeValue(secondParent);
        }


        if (childInd == INVALID || motherInd == INVALID || fatherInd == INVALID) //any of the genotypes are NO_CALL, UNAVAILABLE or MIXED
            return INVALID;

        //index into array playing the part of a 3x3x3 matrix (where 3=NUM_CALLED_GENOTYPETYPES)
        return motherInd*NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES + fatherInd*NUM_CALLED_GENOTYPETYPES + childInd;
    }

    private int genotypeTypeValue(GenotypeType input){
        if (input == GenotypeType.HOM_REF) return 0;
        if (input == GenotypeType.HET) return 1;
        if (input == GenotypeType.HOM_VAR) return 2;
        return -1;
    }

}