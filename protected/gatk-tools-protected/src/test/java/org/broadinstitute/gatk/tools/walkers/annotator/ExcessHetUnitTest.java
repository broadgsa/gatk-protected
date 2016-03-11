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

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: mshand
 * Date: 9/3/15
 */
public class ExcessHetUnitTest {
    private static double DELTA_PRECISION = .001;
    private Allele Aref, T, C;
    private int[] hetPLs, homRefPLs;

    @BeforeSuite
    public void setup() {
        // alleles
        Aref = Allele.create("A", true);
        T = Allele.create("T");
        C = Allele.create("C");

        // simulating 20 reads with Q30 base qualities
        hetPLs = new int[]{240, 0, 240};
        homRefPLs = new int[]{0, 60, 600};
    }

    private Genotype makeGwithPLs(String sample, Allele a1, Allele a2, double[] pls) {
        Genotype gt = new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
        if (pls != null && pls.length > 0) {
            Assert.assertNotNull(gt.getPL());
            Assert.assertTrue(gt.getPL().length > 0);
            for (int i : gt.getPL()) {
                Assert.assertTrue(i >= 0);
            }
            Assert.assertNotEquals(Arrays.toString(gt.getPL()), "[0]");
        }
        return gt;
    }

    private Genotype makeG(String sample, Allele a1, Allele a2, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... genotypes) {
        int start = 10;
        int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles)
                .genotypes(Arrays.asList(genotypes))
                .filters((String) null)
                .make();
    }

    @Test
    public void testExcessHetForMultiallelicVC() {
        //make sure that compound gets (with no ref) don't add to het count
        VariantContext test1 = makeVC("1", Arrays.asList(Aref, T, C),
                makeG("s1", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
                makeG("s4", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s6", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s8", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s10", Aref, T, 2530, 0, 7099, 366, 3056, 14931));

        final double EHresult1 = new ExcessHet().calculateEH(test1, test1.getGenotypes());
        Assert.assertEquals(EHresult1, 5.85, DELTA_PRECISION, "Pass");

        //make sure that hets with different alternate alleles all get counted
        VariantContext test2 = makeVC("2", Arrays.asList(Aref, T, C),
                makeG("s1", Aref, C, 4878, 1623, 11297, 0, 7970, 8847),
                makeG("s2", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s3", Aref, T, 3382, 0, 6364, 1817, 5867, 12246),
                makeG("s4", Aref, T, 2488, 0, 9110, 3131, 9374, 12505),
                makeG("s5", Aref, C, 4530, 2006, 18875, 0, 6847, 23949),
                makeG("s6", Aref, T, 5325, 0, 18692, 389, 16014, 24570),
                makeG("s7", Aref, T, 2936, 0, 29743, 499, 21979, 38630),
                makeG("s8", Aref, T, 6902, 0, 8976, 45, 5844, 9061),
                makeG("s9", Aref, T, 5732, 0, 10876, 6394, 11408, 17802),
                makeG("s10", Aref, T, 2780, 0, 25045, 824, 23330, 30939));

        final double EHresult2 = new ExcessHet().calculateEH(test2, test2.getGenotypes());
        Assert.assertEquals(EHresult2, 25.573, DELTA_PRECISION, "Pass");
    }

    @Test
    public void testSingletonVsCommonAllele() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 10000;
        for (int i = 0; i < numHomRefGTs; i++)
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));

        allGTs.add(makeG("het0", Aref, T, hetPLs));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHsingleton = new ExcessHet().calculateEH(singleton, singleton.getGenotypes());

        final int targetNumHetGTs = 20;
        for (int i = numHetGTs; i < targetNumHetGTs; i++)
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext common = makeVC("common", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHcommon = new ExcessHet().calculateEH(common, common.getGenotypes());

        Assert.assertTrue(Math.abs(EHsingleton) < Math.abs(EHcommon), String.format("singleton=%f common=%f", EHsingleton, EHcommon));
    }

    @Test
    public void testLargeCohorts() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 1000000;
        for (int i = 0; i < numHomRefGTs; i++)
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));

        allGTs.add(makeG("het0", Aref, T, hetPLs));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHsingleton = new ExcessHet().calculateEH(singleton, singleton.getGenotypes());

        for (int i = numHetGTs; i < 100; i++) {
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));
            numHetGTs++;
        }

        final VariantContext hundredton = makeVC("hundredton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHhundredton = new ExcessHet().calculateEH(hundredton, hundredton.getGenotypes());

        Assert.assertTrue(Math.abs(EHsingleton) < Math.abs(EHhundredton), String.format("singleton=%f hundredton=%f", EHsingleton, EHhundredton));

        for (int i = numHetGTs; i < numHomRefGTs; i++)
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext common = makeVC("common", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHcommon = new ExcessHet().calculateEH(common, common.getGenotypes());

        Assert.assertTrue(Math.abs(EHhundredton) < Math.abs(EHcommon), String.format("hundredton=%f common=%f", EHhundredton, EHcommon));
    }

    @Test
    public void testAllHetsForLargeCohorts() {

        final int numGTs = 1000000;

        final List<Genotype> singletonGTs = new ArrayList<>();
        for (int i = 0; i < numGTs; i++)
            singletonGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));

        singletonGTs.add(makeG("het0", Aref, T, hetPLs));

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), singletonGTs.toArray(new Genotype[singletonGTs.size()]));
        final double EHsingleton = new ExcessHet().calculateEH(singleton, singleton.getGenotypes());

        final List<Genotype> allHetGTs = new ArrayList<>();
        for (int i = 0; i < numGTs; i++)
            allHetGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext allHet = makeVC("allHet", Arrays.asList(Aref, T), allHetGTs.toArray(new Genotype[allHetGTs.size()]));
        final double EHHets = new ExcessHet().calculateEH(allHet, allHet.getGenotypes());

        Assert.assertTrue(Math.abs(EHsingleton) < Math.abs(EHHets), String.format("singleton=%f allHets=%f", EHsingleton, EHHets));
    }

    @DataProvider(name = "smallSets")
    public Object[][] counts() {
        return new Object[][]{
                {1, 0, 0, .5},
                {1, 1, 0, .5},
                {1, 1, 1, .7},
                {4, 0, 0, .114},
                {2, 1, 1, .571},
                {0, 2, 2, .957},
                {1, 1, 40, .982},
                {3, 0, 39, .482},
        };
    }


    @Test(dataProvider = "smallSets")
    public void smallSets(int hetCount, int homrefCount, int homvarCount, double expected) {
        double actual = new ExcessHet().exactTest(new int[]{homrefCount, hetCount, homvarCount});
        Assert.assertEquals(actual, expected, DELTA_PRECISION, "Pass");
    }
}
