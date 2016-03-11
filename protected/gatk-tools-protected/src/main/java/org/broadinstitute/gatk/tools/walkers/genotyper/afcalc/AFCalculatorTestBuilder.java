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

package org.broadinstitute.gatk.tools.walkers.genotyper.afcalc;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.gatk.tools.walkers.genotyper.AFPriorProvider;
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.UnifiedGenotypingEngine;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.variant.HomoSapiensConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AFCalculatorTestBuilder {
    final static Allele A = Allele.create("A", true);
    final static Allele C = Allele.create("C");
    final static Allele G = Allele.create("G");
    final static Allele T = Allele.create("T");
    final static Allele AA = Allele.create("AA");
    final static Allele AT = Allele.create("AT");
    final static Allele AG = Allele.create("AG");

    static int sampleNameCounter = 0;

    final int nSamples;
    final int numAltAlleles;
    final AFCalculatorImplementation modelType;
    final PriorType priorType;

    public AFCalculatorTestBuilder(final int nSamples, final int numAltAlleles,
                                   final AFCalculatorImplementation modelType, final PriorType priorType) {
        this.nSamples = nSamples;
        this.numAltAlleles = numAltAlleles;
        this.modelType = modelType;
        this.priorType = priorType;
    }

    @Override
    public String toString() {
        return String.format("AFCalcTestBuilder nSamples=%d nAlts=%d model=%s prior=%s", nSamples, numAltAlleles, modelType, priorType);
    }

    public enum PriorType {
        flat,
        human
    }

    public int getNumAltAlleles() {
        return numAltAlleles;
    }

    public int getnSamples() {
        return nSamples;
    }

    public AFCalculator makeModel() {
        return AFCalculatorFactory.createCalculator(modelType, nSamples, getNumAltAlleles(), HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    public double[] makePriors() {
        final int nPriorValues = 2*nSamples+1;
        final double human_theta = 0.001;

        switch ( priorType ) {
            case flat:
                return MathUtils.normalizeFromLog10(new double[nPriorValues], true);  // flat priors

            //TODO break dependency with human... avoid special reference to this species.
            case human:

                final AFPriorProvider log10priorProvider = GenotypingEngine.composeAlleleFrequencyPriorProvider(2*nSamples, human_theta, new ArrayList<Double>());
                final double[] humanPriors = log10priorProvider.forTotalPloidy(2*nSamples);
                return humanPriors;
            default:
                throw new RuntimeException("Unexpected type " + priorType);
        }
    }

    public VariantContext makeACTest(final List<Integer> ACs, final int nNonInformative, final int nonTypePL) {
        return makeACTest(ArrayUtils.toPrimitive(ACs.toArray(new Integer[]{})), nNonInformative, nonTypePL);
    }

    public VariantContext makeACTest(final int[] ACs, final int nNonInformative, final int nonTypePL) {
        final int nChrom = nSamples * 2;

        final int[] nhet = new int[numAltAlleles];
        final int[] nhomvar = new int[numAltAlleles];

        for ( int i = 0; i < ACs.length; i++ ) {
            final double p = ACs[i] / (1.0 * nChrom);
            nhomvar[i] = (int)Math.floor((nSamples - nNonInformative) * p * p);
            nhet[i] = ACs[i] - 2 * nhomvar[i];

            if ( nhet[i] < 0 )
                throw new IllegalStateException("Bug! nhet[i] < 0");
        }

        final long calcAC = MathUtils.sum(nhet) + 2 * MathUtils.sum(nhomvar);
        if ( calcAC != MathUtils.sum(ACs) )
            throw new IllegalStateException("calculated AC " + calcAC + " not equal to desired AC " + Utils.join(",", ACs));

        return makeACTest(nhet, nhomvar, nNonInformative, nonTypePL);
    }

    public VariantContext makeACTest(final int[] nhet, final int[] nhomvar, final int nNonInformative, final int nonTypePL) {
        List<Genotype> samples = new ArrayList<Genotype>(nSamples);

        for ( int altI = 0; altI < nhet.length; altI++ ) {
            for ( int i = 0; i < nhet[altI]; i++ )
                samples.add(makePL(GenotypeType.HET, nonTypePL, altI+1));
            for ( int i = 0; i < nhomvar[altI]; i++ )
                samples.add(makePL(GenotypeType.HOM_VAR, nonTypePL, altI+1));
        }

        final Genotype nonInformative = makeNonInformative();
        samples.addAll(Collections.nCopies(nNonInformative, nonInformative));

        final int nRef = Math.max((int) (nSamples - nNonInformative - MathUtils.sum(nhet) - MathUtils.sum(nhomvar)), 0);
        samples.addAll(Collections.nCopies(nRef, makePL(GenotypeType.HOM_REF, nonTypePL, 0)));

        samples = samples.subList(0, nSamples);

        if ( samples.size() > nSamples )
            throw new IllegalStateException("too many samples");

        VariantContextBuilder vcb = new VariantContextBuilder("x", "1", 1, 1, getAlleles());
        vcb.genotypes(samples);
        return vcb.make();
    }

    public List<Allele> getAlleles() {
        return Arrays.asList(A, C, G, T, AA, AT, AG).subList(0, numAltAlleles+1);
    }

    public List<Allele> getAlleles(final GenotypeType type, final int altI) {
        switch (type) {
            case HOM_REF: return Arrays.asList(getAlleles().get(0), getAlleles().get(0));
            case HET:     return Arrays.asList(getAlleles().get(0), getAlleles().get(altI));
            case HOM_VAR: return Arrays.asList(getAlleles().get(altI), getAlleles().get(altI));
            default: throw new IllegalArgumentException("Unexpected type " + type);
        }
    }

    public Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }

    private int numPLs() {
        return GenotypeLikelihoods.numLikelihoods(numAltAlleles+1, 2);
    }

    public Genotype makeNonInformative() {
        final int[] nonInformativePLs = new int[GenotypeLikelihoods.numLikelihoods(numAltAlleles, 2)];
        return makePL(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), nonInformativePLs);
    }

    public Genotype makePL(final GenotypeType type, final int nonTypePL, final int altI) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(getAlleles(type, altI));

        final int[] pls = new int[numPLs()];
        Arrays.fill(pls, nonTypePL);

        int index = 0;
        switch ( type ) {
            case HOM_REF: index = GenotypeLikelihoods.calculatePLindex(0, 0); break;
            case HET:     index = GenotypeLikelihoods.calculatePLindex(0, altI); break;
            case HOM_VAR: index = GenotypeLikelihoods.calculatePLindex(altI, altI); break;
        }
        pls[index] = 0;
        gb.PL(pls);

        return gb.make();
    }
}