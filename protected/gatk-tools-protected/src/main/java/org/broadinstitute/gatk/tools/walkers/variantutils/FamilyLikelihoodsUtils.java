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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.*;

import java.util.Arrays;
import java.util.*;

/**
 * FamilyLikelihoodsUtils code is based on PhaseByTransmission with added posterior probability calculations
 */

public class FamilyLikelihoodsUtils {
    private static Logger logger = Logger.getLogger(FamilyLikelihoodsUtils.class);

    //Matrix of priors for all genotype combinations
    final private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>> mvCountMatrix =
            new EnumMap<>(GenotypeType.class);

    final int NUM_CALLED_GENOTYPETYPES = 3; //HOM_REF, HET, and HOM_VAR

    double[] configurationLikelihoodsMatrix = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];

    ArrayList<Sample> trios = new ArrayList<>();

    public final double NO_JOINT_VALUE = -1.0;

    private double deNovoPrior = 1e-8;

    private final double ONE_THIRD = 0.333333333333333333;
    private final double LOG10_OF_ONE_THIRD = -0.4771213;

    private enum FamilyMember {
        MOTHER,
        FATHER,
        CHILD
    }

    /**
     * Applies the trio genotype combination to the given trio.
     * @param motherGenotype: Original genotype of the mother
     * @param fatherGenotype: Original genotype of the father
     * @param childGenotype: Original genotype of the child
     * @param updatedGenotypes: An ArrayList<Genotype> to which the newly updated genotypes are added in the following order: Mother, Father, Child
     */
    public void getUpdatedGenotypes(final VariantContext vc, final Genotype motherGenotype, final Genotype fatherGenotype, final Genotype childGenotype, final ArrayList<Genotype> updatedGenotypes){
        //genotypes here can be no call
        boolean fatherIsCalled = fatherGenotype != null && hasCalledGT(fatherGenotype.getType()) && fatherGenotype.hasLikelihoods();
        boolean motherIsCalled = motherGenotype != null && hasCalledGT(motherGenotype.getType()) && motherGenotype.hasLikelihoods();
        boolean childIsCalled = childGenotype != null && hasCalledGT(childGenotype.getType()) && childGenotype.hasLikelihoods();

        //default to posteriors equal to likelihoods (flat priors) in case input genotypes are not called
        double[] uninformativeLikelihoods = {ONE_THIRD, ONE_THIRD, ONE_THIRD};

        double[] motherLikelihoods = motherIsCalled? GeneralUtils.normalizeFromLog10(motherGenotype.getLikelihoods().getAsVector()) : uninformativeLikelihoods;
        double[] fatherLikelihoods = fatherIsCalled? GeneralUtils.normalizeFromLog10(fatherGenotype.getLikelihoods().getAsVector()) : uninformativeLikelihoods;
        double[] childLikelihoods = childIsCalled? GeneralUtils.normalizeFromLog10(childGenotype.getLikelihoods().getAsVector()) : uninformativeLikelihoods;

        //these are also in log10 space
        double[] motherLog10Posteriors = getPosteriors(FamilyMember.MOTHER);
        double[] fatherLog10Posteriors = getPosteriors(FamilyMember.FATHER);
        double[] childLog10Posteriors = getPosteriors(FamilyMember.CHILD);

        double[] motherPosteriors = GeneralUtils.normalizeFromLog10(motherLog10Posteriors);
        double[] fatherPosteriors = GeneralUtils.normalizeFromLog10(fatherLog10Posteriors);
        double[] childPosteriors = GeneralUtils.normalizeFromLog10(childLog10Posteriors);


        double jointPosteriorProbability =  -1;
        //jointTrioLikelihood is combined likelihoods (before prior) of best configuration after applying prior
        double jointTrioLikelihood = -1;
        if(childIsCalled && motherIsCalled && fatherIsCalled) {
            jointTrioLikelihood = motherLikelihoods[MathUtils.maxElementIndex(motherPosteriors)]*fatherLikelihoods[MathUtils.maxElementIndex(fatherPosteriors)]*childLikelihoods[MathUtils.maxElementIndex(childPosteriors)];
            jointPosteriorProbability = MathUtils.arrayMax(motherPosteriors)*MathUtils.arrayMax(fatherPosteriors)*MathUtils.arrayMax(childPosteriors);
        }

        updatedGenotypes.add(getUpdatedGenotype(vc, motherGenotype, jointTrioLikelihood, jointPosteriorProbability, motherLog10Posteriors));
        updatedGenotypes.add(getUpdatedGenotype(vc, fatherGenotype, jointTrioLikelihood, jointPosteriorProbability, fatherLog10Posteriors));
        updatedGenotypes.add(getUpdatedGenotype(vc, childGenotype, jointTrioLikelihood, jointPosteriorProbability, childLog10Posteriors));
    }

    private Genotype getUpdatedGenotype(final VariantContext vc, final Genotype genotype, final double jointLikelihood, final double jointPosteriorProb, final double[] log10Posteriors){
        //Don't update null, missing or unavailable genotypes
        if(genotype == null || !hasCalledGT(genotype.getType()))
            return genotype;

        int phredScaledJL = -1;
        int phredScaledJP = -1;
        if(jointLikelihood != NO_JOINT_VALUE){
            double dphredScaledJL = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1-jointLikelihood));
            phredScaledJL = dphredScaledJL < Byte.MAX_VALUE ? (byte)dphredScaledJL : Byte.MAX_VALUE;
        }
        if(jointPosteriorProb != NO_JOINT_VALUE){
            double dphredScaledJP = QualityUtils.phredScaleLog10ErrorRate(Math.log10(1-jointPosteriorProb));
            phredScaledJP = dphredScaledJP < Byte.MAX_VALUE ? (byte)dphredScaledJP : Byte.MAX_VALUE;
        }

        //Add the joint trio calculations
        final Map<String, Object> genotypeAttributes = new HashMap<>();
        genotypeAttributes.putAll(genotype.getExtendedAttributes());
        genotypeAttributes.put(GATKVCFConstants.JOINT_LIKELIHOOD_TAG_NAME, phredScaledJL);
        genotypeAttributes.put(GATKVCFConstants.JOINT_POSTERIOR_TAG_NAME, phredScaledJP);

        final GenotypeBuilder builder = new GenotypeBuilder(genotype);

        //final double[] log10Posteriors = MathUtils.toLog10(normalizedPosteriors);

        //update genotype types based on posteriors
        GATKVariantContextUtils.updateGenotypeAfterSubsetting(vc.getAlleles(), builder,
                GATKVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, log10Posteriors, vc.getAlleles());

        builder.attribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY,
                Utils.listFromPrimitives(GenotypeLikelihoods.fromLog10Likelihoods(log10Posteriors).getAsPLs()));
        builder.attributes(genotypeAttributes);
        return builder.make();
    }

    //marginalize over the configurationLikelihoodsMatrix and normalize to get the posteriors
    private double[] getPosteriors(final FamilyMember recalcInd) {
        double[] marginalOverChangedHR = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];
        double[] marginalOverChangedHET = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];
        double[] marginalOverChangedHV = new double[NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES];
        final double[] recalcPosteriors = new double[NUM_CALLED_GENOTYPETYPES];

        final GenotypeType[] calledTypes = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};
        int counter = 0;

        switch (recalcInd) {
            case MOTHER:
                for(final GenotypeType father : calledTypes) {
                    for(final GenotypeType child : calledTypes) {
                        GenotypeType mother;
                        mother = GenotypeType.HOM_REF;
                        marginalOverChangedHR[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        mother = GenotypeType.HET;
                        marginalOverChangedHET[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        mother = GenotypeType.HOM_VAR;
                        marginalOverChangedHV[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        counter++;
                    }
                }
                break;
            case FATHER:
                for(final GenotypeType mother : calledTypes){
                    for (final GenotypeType child : calledTypes){
                        GenotypeType father;
                        father = GenotypeType.HOM_REF;
                        marginalOverChangedHR[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        father = GenotypeType.HET;
                        marginalOverChangedHET[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        father = GenotypeType.HOM_VAR;
                        marginalOverChangedHV[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        counter++;
                    }
                }
                break;
            case CHILD:
                for(final GenotypeType mother : calledTypes){
                    for (final GenotypeType father: calledTypes){
                        GenotypeType child;
                        child = GenotypeType.HOM_REF;
                        marginalOverChangedHR[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        child = GenotypeType.HET;
                        marginalOverChangedHET[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        child = GenotypeType.HOM_VAR;
                        marginalOverChangedHV[counter] = configurationLikelihoodsMatrix[getLikelihoodMatrixIndex(mother, father, child)];
                        counter++;
                    }
                }
                break;
            default:
                throw new UserException(String.format("%d does not indicate a valid trio FamilyMember -- use 0 for mother, 1 for father, 2 for child",recalcInd.ordinal()));
        }

        recalcPosteriors[0] = MathUtils.log10sumLog10(marginalOverChangedHR,0);
        recalcPosteriors[1] = MathUtils.log10sumLog10(marginalOverChangedHET,0);
        recalcPosteriors[2] = MathUtils.log10sumLog10(marginalOverChangedHV,0);

        return MathUtils.normalizeFromLog10(recalcPosteriors,true,true);
    }

    public void initialize(final double DNprior, final Set<String> vcfSamples, final Map<String,Set<Sample>> families){
        this.deNovoPrior = DNprior;
        Arrays.fill(configurationLikelihoodsMatrix,0);
        buildMatrices();
        trios = setTrios(vcfSamples, families);
    }

    public GenotypesContext calculatePosteriorGLs(final VariantContext vc){
        final GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());

        for (final Sample sample : trios) {
            Genotype mother = vc.getGenotype(sample.getMaternalID());
            Genotype father = vc.getGenotype(sample.getPaternalID());
            Genotype child = vc.getGenotype(sample.getID());

            //Keep only trios and parent/child pairs
            if(mother == null && father == null || child == null) {
                logger.warn("Null genotypes in variant: "+vc.toStringDecodeGenotypes());
                continue;
            }

            final ArrayList<Genotype> trioGenotypes = new ArrayList<>(3);
            updateFamilyGenotypes(vc, mother, father, child, trioGenotypes);

            //replace uses sample names to match genotypes, so order doesn't matter
            if (trioGenotypes.size() > 0) {
                genotypesContext.replace(trioGenotypes.get(0));
                genotypesContext.replace(trioGenotypes.get(1));
                genotypesContext.replace(trioGenotypes.get(2));
            }
        }

        return genotypesContext;
    }

    /**
     * Select trios and parent/child pairs only
     */
    private ArrayList<Sample> setTrios(Set<String> vcfSamples, Map<String,Set<Sample>> families){
        Set<Sample> family;
        ArrayList<Sample> parents;
        final ArrayList<Sample> trios = new ArrayList<>();
        for(final Map.Entry<String,Set<Sample>> familyEntry : families.entrySet()){
            family = familyEntry.getValue();

            // Since getFamilies(vcfSamples) above still returns parents of samples in the VCF even if those parents are not in the VCF, need to subset down here:
            final Set<Sample> familyMembersInVCF = new TreeSet<>();
            for(final Sample familyMember : family){
                if (vcfSamples.contains(familyMember.getID())) {
                    familyMembersInVCF.add(familyMember);
                }
            }
            family = familyMembersInVCF;

            if(family.size() == 3){
                for(final Sample familyMember : family){
                    parents = familyMember.getParents();
                    if(parents.size()==2){
                        if(family.containsAll(parents))
                            trios.add(familyMember);
                    }
                }
            }

        }
        return trios;
    }

    //Create a lookup matrix to find the number of MVs for each family genotype combination
    private void buildMatrices(){
        for(final GenotypeType mother : GenotypeType.values()){
            mvCountMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,Integer>>(GenotypeType.class));
            for(final GenotypeType father : GenotypeType.values()){
                mvCountMatrix.get(mother).put(father,new EnumMap<GenotypeType, Integer>(GenotypeType.class));
                for(final GenotypeType child : GenotypeType.values()){
                    mvCountMatrix.get(mother).get(father).put(child, getCombinationMVCount(mother, father, child));
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
        final ArrayList<GenotypeType> parents = new ArrayList<>();
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
     * @param vc: Input variant context
     * @param mother: Mother's genotype from vc input
     * @param father: Father's genotype from vc input
     * @param child: Child's genotype from vc input
     * @param finalGenotypes: An ArrayList<Genotype> containing the updated genotypes
     */
    private void updateFamilyGenotypes(VariantContext vc, Genotype mother, Genotype father, Genotype child, ArrayList<Genotype> finalGenotypes) {

        //If one of the parents is not called, fill in with uninformative likelihoods
        Map<GenotypeType,Double> motherLikelihoods = getLikelihoodsAsMapSafeNull(mother);
        Map<GenotypeType,Double> fatherLikelihoods = getLikelihoodsAsMapSafeNull(father);
        Map<GenotypeType,Double> childLikelihoods = getLikelihoodsAsMapSafeNull(child);

        //if the child isn't called or neither parent is called, there's no extra inheritance information in that trio so return
        if (!hasCalledGT(child.getType()) || (!hasCalledGT(mother.getType()) && !hasCalledGT(father.getType())))
            return;

        //Fill the configurationLikelihoodsMatrix for each genotype combination
        int matInd;
        int mvCount;
        double jointLikelihood;
        double mvCoeff;
        double configurationLikelihood;
        for(final Map.Entry<GenotypeType,Double> childGenotype :
                childLikelihoods.entrySet()){
            for(final Map.Entry<GenotypeType,Double> motherGenotype :
                    motherLikelihoods.entrySet()){
                for(final Map.Entry<GenotypeType,Double> fatherGenotype :
                        fatherLikelihoods.entrySet()){
                    mvCount = mvCountMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey());
                    jointLikelihood = motherGenotype.getValue()+fatherGenotype.getValue()+childGenotype.getValue();
                    mvCoeff = mvCount>0 ? Math.pow(deNovoPrior,mvCount) : (1.0-10*deNovoPrior-deNovoPrior*deNovoPrior);
                    configurationLikelihood =  Math.log10(mvCoeff) + jointLikelihood;
                    matInd = getLikelihoodMatrixIndex(motherGenotype.getKey(), fatherGenotype.getKey(), childGenotype.getKey());
                    configurationLikelihoodsMatrix[matInd] = configurationLikelihood;
                }
            }
        }

        getUpdatedGenotypes(vc, mother, father, child, finalGenotypes);
    }

    //Get a Map of genotype (log10)likelihoods
    private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(Genotype genotype){
        final EnumMap<GenotypeType,Double> likelihoodsMap = new EnumMap<>(GenotypeType.class);
        double[] likelihoods;

        if (genotype != null && hasCalledGT(genotype.getType()) && genotype.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY)) {
            Object GPfromVCF = genotype.getExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
            //parse the GPs into a vector of probabilities
            final String[] likelihoodsAsStringVector = ((String)GPfromVCF).split(",");
            final double[] likelihoodsAsVector = new double[likelihoodsAsStringVector.length];
            for ( int i = 0; i < likelihoodsAsStringVector.length; i++ ) {
                likelihoodsAsVector[i] = Double.parseDouble(likelihoodsAsStringVector[i]) / -10.0;
            }
            //keep in log10 space for large GQs
            likelihoods = GeneralUtils.normalizeFromLog10(likelihoodsAsVector, true, true);
        }

        //In case of null, unavailable or no call, all likelihoods are log10(1/3)
        else if(genotype == null || !hasCalledGT(genotype.getType()) || genotype.getLikelihoods() == null){
            likelihoods = new double[3];
            likelihoods[0] = LOG10_OF_ONE_THIRD;
            likelihoods[1] = LOG10_OF_ONE_THIRD;
            likelihoods[2] = LOG10_OF_ONE_THIRD;
        }

        //No posteriors in VC, use PLs
        else
            likelihoods = GeneralUtils.normalizeFromLog10(genotype.getLikelihoods().getAsVector(),true,true);

        likelihoodsMap.put(GenotypeType.HOM_REF,likelihoods[genotypeTypeToValue(GenotypeType.HOM_REF)]);
        likelihoodsMap.put(GenotypeType.HET,likelihoods[genotypeTypeToValue(GenotypeType.HET)]);
        likelihoodsMap.put(GenotypeType.HOM_VAR, likelihoods[genotypeTypeToValue(GenotypeType.HOM_VAR)]);
        return likelihoodsMap;
    }

    private int getLikelihoodMatrixIndex(GenotypeType mother, GenotypeType father, GenotypeType child){
        int childInd = genotypeTypeToValue(child);
        int motherInd;
        int fatherInd;
        final int INVALID = -1;
        motherInd = genotypeTypeToValue(mother);
        fatherInd = genotypeTypeToValue(father);

        if (childInd == INVALID || motherInd == INVALID || fatherInd == INVALID) //any of the genotypes are NO_CALL, UNAVAILABLE or MIXED
            return INVALID;

        //index into array playing the part of a 3x3x3 matrix (where 3=NUM_CALLED_GENOTYPETYPES)
        return motherInd*NUM_CALLED_GENOTYPETYPES*NUM_CALLED_GENOTYPETYPES + fatherInd*NUM_CALLED_GENOTYPETYPES + childInd;
    }

    private int genotypeTypeToValue(GenotypeType input){
        if (input == GenotypeType.HOM_REF) return 0;
        if (input == GenotypeType.HET) return 1;
        if (input == GenotypeType.HOM_VAR) return 2;
        return -1;
    }

    //this excludes mixed genotypes, whereas the htsjdk Genotype.isCalled() will return true if the GenotypeType is mixed
    private boolean hasCalledGT(GenotypeType genotype){
        return genotype == GenotypeType.HOM_REF || genotype == GenotypeType.HET || genotype == GenotypeType.HOM_VAR;
    }

}