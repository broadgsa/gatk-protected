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

package org.broadinstitute.gatk.utils.collections;

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

/**
 * Unit tests for {@link CountSet}
 */
public class CountSetUnitTest extends BaseTest {

    @Test(dataProvider="capacities")
    public void testSize(final int capacity) {
        final CountSet empty = new CountSet(capacity);
        Assert.assertEquals(empty.size(), 0);
        CountSet nonEmpty = new CountSet(capacity);
        for (int i = 0; i < capacity * 3; i++) {
            nonEmpty.add(i);
            Assert.assertEquals(nonEmpty.size(),i + 1);
        }
    }

    @Test
    public void testSingleValueAdd() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final HashSet<Integer> reasuranceSet = new HashSet<>(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(500);
            boolean expectedResult = reasuranceSet.add(newInt);
            boolean result = subject.add(newInt);
            Assert.assertEquals(result,expectedResult);
            Assert.assertEquals(subject.size(),reasuranceSet.size());
        }
        for (final int j : reasuranceSet)
            Assert.assertTrue(subject.contains(j));
        for (int j = 0; j < 501; j++)
            Assert.assertEquals(subject.contains(j),reasuranceSet.contains(j));
    }

    @Test
    public void testToIntArray() {
        final CountSet subject = new CountSet(10);
        subject.addAll(1,4,7);
        final int[] intArray = subject.toIntArray();
        Assert.assertEquals(intArray.length,3);
        Assert.assertEquals(intArray[0],1);
        Assert.assertEquals(intArray[1],4);
        Assert.assertEquals(intArray[2],7);
    }

    @Test
    public void testCopyTo() {
        final CountSet subject = new CountSet(10);
        subject.addAll(1,4,7);
        final int[] intArray = new int[3];
        subject.copyTo(intArray);
        Assert.assertEquals(intArray[0],1);
        Assert.assertEquals(intArray[1],4);
        Assert.assertEquals(intArray[2],7);
    }

    @Test
    public void testSetToSingleValue() {
        final CountSet subject = new CountSet(10);
        subject.setTo(-31);
        Assert.assertEquals(subject.size(),1);
        Assert.assertEquals(subject.min(),-31);
        Assert.assertEquals(subject.max(),-31);
        Assert.assertTrue(subject.contains(-31));
        Assert.assertFalse(subject.contains(-21));
    }

    @Test void testSetToArrayOfValues() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(Integer.MAX_VALUE) * (rnd.nextBoolean() ? -1 : 1);
            values[i] = newInt;
        }
        subject.setTo(values);
        Arrays.sort(values);
        Assert.assertEquals(subject.size(),REPEATS);
        Assert.assertEquals(subject.min(),values[0]);
        Assert.assertEquals(subject.max(),values[REPEATS - 1]);
    }

    @Test
    public void testMinMax() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(Integer.MAX_VALUE) * (rnd.nextBoolean() ? -1 : 1);
            values[i] = newInt;
        }
        subject.addAll(values);
        Arrays.sort(values);
        Assert.assertEquals(subject.min(),values[0]);
        Assert.assertEquals(subject.max(),values[REPEATS - 1]);
    }

    @Test
    public void testIncrease() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final HashSet<Integer> reasuranceSet = new HashSet<>(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        final Integer[] valueWrappers = new Integer[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(500);
            values[i] = newInt;
            valueWrappers[i] = newInt;
        }

        subject.incAll(3);

        for (final int j : reasuranceSet)
            Assert.assertTrue(subject.contains(j+3));
        for (int j = 0; j < 501; j++)
            Assert.assertEquals(subject.contains(j+3),reasuranceSet.contains(j));

    }

    @Test
    public void testArrayValueAdd() {
        final int CAPACITY = 10;
        final CountSet subject = new CountSet(CAPACITY);
        final HashSet<Integer> reasuranceSet = new HashSet<>(CAPACITY);
        final int REPEATS = 1000;
        final Random rnd = new Random(13);
        final int[] values = new int[REPEATS];
        final Integer[] valueWrappers = new Integer[REPEATS];
        for (int i = 0; i < REPEATS; i++) {
            int newInt = rnd.nextInt(500);
            values[i] = newInt;
            valueWrappers[i] = newInt;
        }

        boolean expectedResult = reasuranceSet.addAll(Arrays.asList(valueWrappers));
        boolean result = subject.addAll(values);
        Assert.assertEquals(result,expectedResult);
        Assert.assertEquals(subject.size(),reasuranceSet.size());

        for (final int j : reasuranceSet)
            Assert.assertTrue(subject.contains(j));
        for (int j = 0; j < 501; j++)
            Assert.assertEquals(subject.contains(j),reasuranceSet.contains(j));

    }

    @Test
    public void testAddRange() {
        final CountSet subject = new CountSet(10);
        subject.addRange(10,21);
        Assert.assertEquals(subject.size(),12);
        for (int i = 10; i < 22; i++)
            Assert.assertTrue(subject.contains(i));
        for (int i = -1; i < 10; i++)
            Assert.assertFalse(subject.contains(i));
        for (int i = 22; i < 31; i++)
            Assert.assertFalse(subject.contains(i));
    }

    @DataProvider(name="capacities")
    public Iterator<Object[]> capacities() {
        final int MIN = 0;
        final int MAX = 255;
        return new Iterator<Object[]>() {
                private int current = MIN;


            @Override
            public boolean hasNext() {
                return current < MAX;
            }

            @Override
            public Object[] next() {
                return new Object[] { Integer.valueOf(current++) };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }
}
