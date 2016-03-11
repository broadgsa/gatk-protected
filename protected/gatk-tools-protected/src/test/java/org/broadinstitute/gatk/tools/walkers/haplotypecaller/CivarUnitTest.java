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

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: valentin
 * Date: 8/7/13
 * Time: 5:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class CivarUnitTest extends BaseTest {


    @Test(dataProvider="validCivarExamples")
    public void testValidCivarInstanciation(final String civarString) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertNotNull(civar);
    }


    @Test(dataProvider="expectedElementLengths")
    public void testValidCivarElementLength(final String civarString, final int expected) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.elements().size(), expected);
    }


    @Test(dataProvider="expectedElementSizes")
    public void testValidCivarElementSizes(final String civarString, final int[] expected) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.elements().size(),expected.length);
        for (int i  = 0; i < expected.length; i++) {
            Assert.assertEquals(civar.elements().get(i).size(),expected[i]);
        }
    }

    @Test(dataProvider="expectedElementOperators")
    public void testValidCivarElementOperators(final String civarString, final String expected) {

        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.elements().size(),expected.length());
        for (int i  = 0; i < expected.length(); i++) {
            Assert.assertEquals(civar.elements().get(i).operator().charValue,expected.charAt(i));
        }
    }

    @Test(dataProvider="expectedMinimumSequenceLength")
    public void testValidCivarMinimumSequenceLength(final String civarString, final int expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.minimumTemplateSequenceSize(),expected);
    }

    @Test(dataProvider="expectedHasVariation")
    public void testValidCivarHasVariation(final String civarString, final boolean expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.hasVariation(),expected);
    }


    @Test(dataProvider="invalidCivarExamples", expectedExceptions = {IllegalArgumentException.class})
    public void testInvalidInstanciation(final String civarString) {

        final Civar civar = Civar.fromCharSequence(civarString);
    }

    @Test(dataProvider="unrolledTestDataIsUnrolledExamples")
    public void testInUnrolled(final String civarString, final boolean expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.isUnrolled(),expected);
    }

    @Test(dataProvider="unrolledTestDataUnrolledCivarExamples")
    public void testValidCivarUnrolling(final String civarString, final String[] expected) {
        Set<String> expectedSet = new HashSet<>();
        expectedSet.addAll(Arrays.asList(expected));

        final Civar civar = Civar.fromCharSequence(civarString);
        java.util.List<Civar> unrolledList = civar.unroll();
        Assert.assertEquals(unrolledList.size(),expected.length);
        for (int i  = 0; i < expected.length; i++) {
            Assert.assertTrue(expectedSet.contains(unrolledList.get(i).toString()),
                    "Unrolled civar " + unrolledList.get(i).toString() + " not present in expected Set: " +
                            Arrays.toString(expected) + ". Unrolled set is: " + Arrays.toString(unrolledList.toArray()));
        }
    }

    @Test(dataProvider="applyToDataExamples")
    public void testValidCivarUnrolling(final String civarString, final String before, final String expectedAfter) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.applyTo(before),expectedAfter);
    }

    @Test(dataProvider="optionizeDataExamples")
    public void testValidOptionizeAll(final String civarString, final String expected) {
        final Civar civar = Civar.fromCharSequence(civarString);
        Assert.assertEquals(civar.optionalizeAll().toString(),expected);
    }

    @DataProvider(name="validCivarExamples")
    public Iterator<Object[]> validCivarExamples() {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < VALID_CIVAR_EXAMPLES.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { VALID_CIVAR_EXAMPLES[i++][0] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="expectedHasVariation")
    public Iterator<Object[]> expectedHasVariation () {
        return validCivarExamples(5);
    }

    @DataProvider(name="expectedMinimumSequenceLength")
    public Iterator<Object[]> expectedMinimumSequenceLength () {
        return validCivarExamples(4);
    }

    @DataProvider(name="expectedElementOperators")
    public Iterator<Object[]> expectedElementOperators() {
        return validCivarExamples(3);
    }

    @DataProvider(name="expectedElementSizes")
    public Iterator<Object[]> expectedElementSizes() {
        return validCivarExamples(2);
    }

    @DataProvider(name="expectedElementLengths")
    public Iterator<Object[]> expectedElementLengths() {
        return validCivarExamples(1);
    }

    public Iterator<Object[]> validCivarExamples(final int field) {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < VALID_CIVAR_EXAMPLES.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { VALID_CIVAR_EXAMPLES[i][0], VALID_CIVAR_EXAMPLES[i++][field] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="unrolledTestDataIsUnrolledExamples")
    public Iterator<Object[]> unrolledTestDataIsUnrolledExamples() {
        return unrolledTestDataExamples(1);
    }

    @DataProvider(name="unrolledTestDataUnrolledCivarExamples")
    public Iterator<Object[]> unrolledTestDataUnrolledCivarExamples() {
        return unrolledTestDataExamples(2);
    }

    public Iterator<Object[]> unrolledTestDataExamples(final int field) {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < UNROLLED_TEST_DATA.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { UNROLLED_TEST_DATA[i][0], UNROLLED_TEST_DATA[i++][field] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="optionizeDataExamples")
    public Iterator<Object[]> optionizeDataExamples() {
        return optionizeDataExamples(1);
    }

    public Iterator<Object[]> optionizeDataExamples(final int field) {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < OPTIONIZED_TEST_DATA.length;
            }

            @Override
            public Object[] next() {
                return new Object[] { OPTIONIZED_TEST_DATA[i][0], OPTIONIZED_TEST_DATA[i++][field] };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="applyToDataExamples")
    public Iterator<Object[]> applyToDataExamples() {
        return new Iterator<Object[]>() {

            int i = 0;

            @Override
            public boolean hasNext() {
                return i < APPLY_TO_TEST_DATA.length;
            }

            @Override
            public Object[] next() {
                return APPLY_TO_TEST_DATA[i++];
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name="invalidCivarExamples")
    public Object[][] invalidCivarExamples() {
            return INVALID_CIVAR_EXAMPLES;
    }

    // columns : Civar string, number of elements.
    private static final Object[][] INVALID_CIVAR_EXAMPLES = new Object[][] {
            {"(100="},
            {"*=)"},
            {"10(=2T30="},
            {"2*=2T/3*="},
            {"3I(acc)"},
            {"a"},
            {")"},
            {"100&=1"},
            {"?100="},

    };


    private static final Object[][] VALID_CIVAR_EXAMPLES = new Object[][] {
            {"100=", 1, ints(100), "=", 100, false },
            {"*=", 1 , ints(1), "=", 0, false },
            {"10=2T30=", 3, ints(10,2,30), "=T=",42 , true},
            {"*=2T3*=", 3, ints(1,2,3), "=T=",2 , true},
            {"3Iacc",1 , ints(3), "I", 0, true},
            {"Ia",1, ints(1), "I", 0, true},
            {"10D",1, ints(10), "D", 10, true},
            {"*", 1, ints(1), "=", 0, false},
            {"*D", 1, ints(1), "D", 0, true},
            {"10(1D)10=",3, ints(10,1,10), "=(=", 21, true},
            {"1*",1, ints(1), "=", 0, false},
            {"1*2*",2, ints(1,2), "==", 0, false},
            {"*11",2, ints(1,11), "==", 11, false},
            {"100=1T100=", 3, ints(100,1,100), "=T=", 201, true},
            {"100=3Iacg101=", 3, ints(100,3,101), "=I=", 201, true},
            {"100=30Igctcggatgccttgcggggctccagagtcc101=", 3 , ints(100,30,101), "=I=", 201, true},
            {"99=3D99=", 3, ints(99,3,99), "=D=", 201, true},
            {"84=30D84=", 3, ints(84,30,84), "=D=", 198, true},
            {"91=1T9=3Iacg100=", 5, ints(91,1,9,3,100), "=T=I=", 201, true},
            {"71=1T29=3Iacg100=",5, ints(71,1,29,3,100), "=T=I=",201, true},
            {"75=1T8=1T8=1T8=1T8=1T75=", 11, ints(75,1,8,1,8,1,8,1,8,1,75), "=T=T=T=T=T=",187, true},
            {"75=1T?8=", 3, ints(75,1,8), "=T=", 84, true}
    };

    private static final Object[][] UNROLLED_TEST_DATA = new Object[][] {
            { "10=1D10=", true, strs( "10=1D10=") },
            { "10=(1D)10=", true, strs( "10=(1D)10=") },
            { "10=1D?10=", false, strs("10=1=10=", "10=1D10=") },
            { "10=1D?10=3Iacg?10=", false , strs("10=1=10=0=10=","10=1=10=3Iacg10=", "10=1D10=0=10=", "10=1D10=3Iacg10=") },
            {  "10=1D?10=" , false, strs("10=1D10=","10=1=10=") },
            {  "100=1T?100=" , false, strs("100=1T100=","100=1=100=") },
            {  "100=3Iacg?101=" , false, strs("100=3Iacg101=","100=0=101=") },
            {  "100=30Igctcggatgccttgcggggctccagagtcc?101=", false ,strs("100=30Igctcggatgccttgcggggctccagagtcc101=", "100=0=101=") },
            {  "99=3D?99=", false , strs("99=3D99=","99=3=99=") },
            {  "84=30D?84=", false, strs("84=30D84=", "84=30=84=")},
            {  "91=1T?9=3Iacg?100=", false, strs("91=1T9=3Iacg100=", "91=1=9=3Iacg100=", "91=1=9=0=100=", "91=1T9=0=100=") },
            {  "71=1T?29=3Iacg?100=", false , strs("71=1T29=3Iacg100=","71=1=29=3Iacg100=","71=1=29=0=100=", "71=1T29=0=100=") },
           // {  "75=1T?8=1T?8=1T?8=1T?8=1T?75=", false, },
            {  "75=1T?8=", false, strs("75=1T8=","75=1=8=") }
    };

    private static final Object[][] OPTIONIZED_TEST_DATA = new Object[][] {
            { "10=1D10=", "10=1D?10=" },
            {"100=1T100=","100=1T?100=" },
            {"100=3Iacg101=", "100=3Iacg?101=" },
            {"100=30Igctcggatgccttgcggggctccagagtcc101=","100=30Igctcggatgccttgcggggctccagagtcc?101="},
            {"99=3D99=", "99=3D?99="},
            {"84=30D84=", "84=30D?84="},
            {"91=1T9=3Iacg100=", "91=1T?9=3Iacg?100="},
            {"71=1T29=3Iacg100=","71=1T?29=3Iacg?100="},
            {"75=1T8=1T8=1T8=1T8=1T75=", "75=1T?8=1T?8=1T?8=1T?8=1T?75="},
            {"75=1T?8=", "75=1T?8="}
    };

    private static final Object[][] APPLY_TO_TEST_DATA = new Object[][] {
            {"3=1D3=", "ACTAACT", "ACTACT" },
            {"*=1C*=","ACTTACT", "ACTAACT" },
            {"4=3Iacg3=","ACTGACT","ACTGACGACT" },
            {"*=30Igctcggatgccttgcggggctccagagtcc*=","AA","AGCTCGGATGCCTTGCGGGGCTCCAGAGTCCA"},
            {"*=3D*=", "ACTTTAC","ACAC"},
            {"1=30D1=", "AGCTCGGATGCCTTGCGGGGCTCCAGAGTCCA","AA"},
    };


    private static int[] ints(final int ... iii) {
        return iii;
    }

    private static String[] strs(final String ... sss) {
        return sss;
    }

}
