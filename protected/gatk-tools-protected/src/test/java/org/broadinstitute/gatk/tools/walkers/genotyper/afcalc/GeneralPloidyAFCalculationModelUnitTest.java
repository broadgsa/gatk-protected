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

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.tools.walkers.genotyper.GeneralPloidyGenotypeLikelihoods;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 3/28/12
 * Time: 7:44 AM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralPloidyAFCalculationModelUnitTest extends BaseTest {

    static double[] AA1, AB1, BB1;
    static double[] AA2, AB2, AC2, BB2, BC2, CC2;
    static double[] A4_1, B4_1, C4_1, D4_1, E4_1,F4_1;
    static double[] A4_400, B4_310, C4_220, D4_130, E4_121, F4_013;
    static final int numSamples = 4;
    static final int samplePloidy = 4;   // = 2*samplesPerPool

    @BeforeSuite
    public void before() {
        // legacy diploid cases
        AA1 = new double[]{-5.0, -20.0, -20.0};
        AB1 = new double[]{-20.0, 0.0, -20.0};
        BB1 = new double[]{-20.0, -20.0, 0.0};

        // diploid, nAlleles = 3. Ordering is [2 0 0] [1 1 0] [0 2 0] [1 0 1] [0 1 1] [0 0 2], ie AA AB BB AC BC CC
        AA2 = new double[]{0.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        AB2 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0, -20.0};
        AC2 = new double[]{-20.0, -20.0, -20.0, 0.0, -20.0, -20.0};
        BB2 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0, -20.0};
        BC2 = new double[]{-20.0, -20.0, -20.0, -20.0, 0.0, -20.0};
        CC2 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, 0.0};
        
        // pool (i.e. polyploid cases)
        // NAlleles = 2, ploidy=4
        // ordering is [4 0] [3 1] [2 2 ] [1 3] [0 4]

        A4_1 = new double[]{-3.0, -20.0, -20.0, -20.0, -20.0};
        B4_1 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0};
        C4_1 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0};
        D4_1 = new double[]{-20.0, -20.0, 0.0,   0.0, -20.0};
        E4_1 = new double[]{-20.0, -20.0, 0.0,   0.0, -20.0};
        F4_1 = new double[]{-20.0, -20.0, -20.0,   -20.0, 0.0};

        // NAlleles = 3, ploidy = 4
        // ordering is [4 0 0] [3 1 0] [2 2 0] [1 3 0] [0 4 0] [3 0 1] [2 1 1] [1 2 1] [0 3 1] [2 0 2] [1 1 2] [0 2 2] [1 0 3] [0 1 3] [0 0 4]
        A4_400 = new double[]{0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        B4_310 = new double[]{-20.0, 0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        C4_220 = new double[]{-20.0, -20.0, 0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        D4_130 = new double[]{-20.0, -20.0, -20.0,   0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        E4_121 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0,   0.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
        F4_013 = new double[]{-20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, 0.0, -20.0};

    }

    private class GetGLsTest extends TestDataProvider {
        GenotypesContext GLs;
        int numAltAlleles;
        String name;
        int ploidy;
        private GetGLsTest(String name, int numAltAlleles, int ploidy, Genotype... arg) {
            super(GetGLsTest.class, name);
            GLs = GenotypesContext.create(arg);
            this.name = name;
            this.numAltAlleles = numAltAlleles;
            this.ploidy = ploidy;
        }

        public String toString() {
            return String.format("%s input=%s", super.toString(), GLs);
        }
    }

    private static Genotype createGenotype(String name, double[] gls, int ploidy) {
        Allele[] alleles = new Allele[ploidy];
        
        for (int i=0; i < ploidy; i++)
            alleles[i] = Allele.NO_CALL;
        
        return new GenotypeBuilder(name, Arrays.asList(alleles)).PL(gls).make();
    }                              

    @DataProvider(name = "getGLs")
    public Object[][] createGLsData() {

        // bi-allelic diploid case
        new GetGLsTest("B0", 1, 2, createGenotype("AA1", AA1,2), createGenotype("AA2", AA1,2), createGenotype("AA3", AA1,2));
        new GetGLsTest("B1", 1, 2, createGenotype("AA1", AA1,2), createGenotype("AA2", AA1,2), createGenotype("AB", AB1,2));
        new GetGLsTest("B2", 1, 2, createGenotype("AA1", AA1,2), createGenotype("BB", BB1,2), createGenotype("AA2", AA1,2));
        new GetGLsTest("B3a", 1, 2, createGenotype("AB", AB1,2), createGenotype("AA", AA1,2), createGenotype("BB", BB1,2));
        new GetGLsTest("B3b", 1, 2, createGenotype("AB1", AB1,2), createGenotype("AB2", AB1,2), createGenotype("AB3", AB1,2));
        new GetGLsTest("B4", 1, 2, createGenotype("BB1", BB1,2), createGenotype("BB2", BB1,2), createGenotype("AA", AA1,2));
        new GetGLsTest("B5", 1, 2, createGenotype("BB1", BB1,2), createGenotype("AB", AB1,2), createGenotype("BB2", BB1,2));
        new GetGLsTest("B6", 1, 2, createGenotype("BB1", BB1,2), createGenotype("BB2", BB1,2), createGenotype("BB3", BB1,2));

        // tri-allelic diploid case
        new GetGLsTest("B1C0", 2, 2, createGenotype("AA1", AA2,2), createGenotype("AA2", AA2,2), createGenotype("AB", AB2,2));
        new GetGLsTest("B0C1", 2, 2, createGenotype("AA1", AA2,2), createGenotype("AA2", AA2,2), createGenotype("AC", AC2,2));
        new GetGLsTest("B1C1a", 2,2,  createGenotype("AA", AA2,2), createGenotype("AB", AB2,2), createGenotype("AC", AC2,2));
        new GetGLsTest("B1C1b", 2,2,  createGenotype("AA1", AA2,2), createGenotype("AA2", AA2,2), createGenotype("BC", BC2,2));
        new GetGLsTest("B2C1", 2, 2, createGenotype("AB1", AB2,2), createGenotype("AB2", AB2,2), createGenotype("AC", AC2,2));
        new GetGLsTest("B3C2a", 2, 2, createGenotype("AB", AB2,2), createGenotype("BC1", BC2,2), createGenotype("BC2", BC2,2));
        new GetGLsTest("B3C2b", 2, 2, createGenotype("AB", AB2,2), createGenotype("BB", BB2,2), createGenotype("CC", CC2,2));

        // bi-allelic pool case
        new GetGLsTest("P0", 1, samplePloidy, createGenotype("A4_1", A4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P1", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("B4_1", B4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P2a", 1,samplePloidy, createGenotype("A4_1", A4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P2b", 1, samplePloidy,createGenotype("B4_1", B4_1,samplePloidy), createGenotype("B4_1", B4_1,samplePloidy), createGenotype("A4_1", A4_1,samplePloidy));
        new GetGLsTest("P4", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy));
        new GetGLsTest("P6", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("F4_1", F4_1,samplePloidy), createGenotype("C4_1", C4_1,samplePloidy));
        new GetGLsTest("P8", 1, samplePloidy,createGenotype("A4_1", A4_1,samplePloidy), createGenotype("F4_1", F4_1,samplePloidy), createGenotype("F4_1", F4_1,samplePloidy));

        // multi-allelic pool case
        new GetGLsTest("B1C3", 2, samplePloidy,createGenotype("A4_400", A4_400,samplePloidy), createGenotype("A4_400", A4_400,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B3C9", 2, samplePloidy,createGenotype("F4_013", F4_013,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B6C0", 2, samplePloidy,createGenotype("B4_310", B4_310,samplePloidy), createGenotype("C4_220", C4_220,samplePloidy), createGenotype("D4_130", D4_130,samplePloidy));
        new GetGLsTest("B6C4", 2, samplePloidy,createGenotype("D4_130", D4_130,samplePloidy), createGenotype("E4_121", E4_121,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B4C7", 2, samplePloidy,createGenotype("F4_013", F4_013,samplePloidy), createGenotype("E4_121", E4_121,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy));
        new GetGLsTest("B2C3", 2, samplePloidy,createGenotype("A4_400", A4_400,samplePloidy), createGenotype("F4_013", F4_013,samplePloidy), createGenotype("B4_310", B4_310,samplePloidy));

        return GetGLsTest.getTests(GetGLsTest.class);
    }

    @Test(dataProvider = "getGLs")
    public void testGLs(GetGLsTest cfg) {
        final int len = GeneralPloidyGenotypeLikelihoods.getNumLikelihoodElements(1 + cfg.numAltAlleles, cfg.ploidy * cfg.GLs.size());
        double[] priors = new double[len];  // flat priors

        final GeneralPloidyExactAFCalculator calc = new GeneralPloidyExactAFCalculator();
        calc.combineSinglePools(cfg.GLs, cfg.ploidy,cfg.numAltAlleles + 1, priors);
        int nameIndex = 1;

        for ( int allele = 0; allele < cfg.numAltAlleles; allele++, nameIndex+=2 ) {
            int expectedAlleleCount = Integer.valueOf(cfg.name.substring(nameIndex, nameIndex + 1));
            int calculatedAlleleCount = calc.getAltAlleleCountOfMAP(allele);
            Assert.assertEquals(calculatedAlleleCount, expectedAlleleCount);
        }
    }


}
