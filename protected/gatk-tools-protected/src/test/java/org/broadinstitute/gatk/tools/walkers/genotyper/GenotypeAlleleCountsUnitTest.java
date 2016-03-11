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

package org.broadinstitute.gatk.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Test {@link org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeAlleleCounts}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeAlleleCountsUnitTest {

    @Test(dataProvider="ploidyData")
    public void testFirst(final int ploidy) {
        final GenotypeAlleleCounts subject = GenotypeAlleleCounts.first(ploidy);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.ploidy(), ploidy);
        Assert.assertEquals(subject.distinctAlleleCount(),1);
        Assert.assertEquals(subject.alleleCountAt(0),ploidy);
        Assert.assertEquals(subject.alleleCountFor(0),ploidy);
        Assert.assertEquals(subject.alleleRankFor(0),0);
        Assert.assertEquals(subject.alleleRankFor(1),-2);
        Assert.assertTrue(subject.containsAllele(0));
        Assert.assertFalse(subject.containsAllele(1));
        Assert.assertEquals(subject.alleleIndexAt(0),0);
        Assert.assertEquals(subject.maximumAlleleIndex(),0);
        Assert.assertEquals(subject.minimumAlleleIndex(),0);
        Assert.assertTrue(subject.compareTo(subject) == 0);
        Assert.assertTrue(subject.equals(subject));
        Assert.assertEquals(subject.index(),0);
        Assert.assertEquals(subject.asAlleleList(testAlleles),Collections.nCopies(ploidy,testAlleles.get(0)));
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            expected[0] = ploidy;
            Assert.assertEquals(subject.alleleCountsByIndex(maximumAlleleIndex),expected);
        }
    }

    @Test(dataProvider = "ploidyData",dependsOnMethods = "testFirst")
    public void testNext(final int ploidy) {
        if (ploidy == 0)
            testNextZeroPloidy();
        else if (ploidy == 1)
            testNextOnePloidy();
        else
            testPloidyTwoOrMore(ploidy);
    }

    @Test(dataProvider = "ploidyData",dependsOnMethods = "testNext")
    public void testIncrease(final int ploidy) {
        if (ploidy == 0)
            testNextZeroPloidyIncrease();
        else if (ploidy == 1)
            testNextOnePloidyIncrease();
        else
            testPloidyTwoOrMoreIncrease(ploidy);
    }

    private void testNextZeroPloidy() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(0);
        final GenotypeAlleleCounts next = first.next();
        Assert.assertEquals(first,next);
        Assert.assertEquals(first.compareTo(next),0);
        Assert.assertEquals(next.compareTo(first), 0);
        Assert.assertEquals(next.distinctAlleleCount(),0);
        Assert.assertEquals(next.ploidy(),0);
        Assert.assertEquals(next.index(),0);
        Assert.assertEquals(next.asAlleleList(testAlleles),Collections.EMPTY_LIST);
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= 10; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex),expected);
        }
    }

    private void testNextOnePloidy() {
       final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(1);
       GenotypeAlleleCounts current = first;

       while (!current.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts next = current.next();
            Assert.assertEquals(next.minimumAlleleIndex(),next.maximumAlleleIndex());
            Assert.assertEquals(next.minimumAlleleIndex(),current.minimumAlleleIndex() + 1);
            Assert.assertEquals(next.alleleCountAt(0),1);
            Assert.assertEquals(next.alleleIndexAt(0),next.minimumAlleleIndex());
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex()),0);
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex() + 1),-2);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex()),1);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex()+1),0);
            Assert.assertEquals(next.ploidy(),1);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));

            Assert.assertEquals(next.index(), current.index() + 1);
            Assert.assertEquals(next.ploidy(),current.ploidy());

            Assert.assertEquals(next.asAlleleList(testAlleles),Collections.singletonList(testAlleles.get(next.maximumAlleleIndex())));

           for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
               final int[] expected = new int[maximumAlleleIndex + 1];
               if (maximumAlleleIndex >= current.minimumAlleleIndex() + 1) expected[current.minimumAlleleIndex() + 1] = 1;
               Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex),expected);
           }
           current = next;
       }
    }

    private void testPloidyTwoOrMore(final int ploidy) {
        if (ploidy < 2)
            throw new IllegalArgumentException();

        GenotypeAlleleCounts current = GenotypeAlleleCounts.first(ploidy);

        while (!current.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts next = current.next();
            if (current.distinctAlleleCount() == 1) {
                Assert.assertEquals(next.maximumAlleleIndex(),current.maximumAlleleIndex() + 1);
                Assert.assertEquals(next.distinctAlleleCount(), 2 );
                Assert.assertEquals(next.minimumAlleleIndex(), 0 );
            } else {
                Assert.assertEquals(next.maximumAlleleIndex(),current.maximumAlleleIndex());
                Assert.assertEquals(next.minimumAlleleIndex(),current.alleleCountAt(0) > 1 ? 0
                    : current.alleleCountAt(0) == 1 ? current.minimumAlleleIndex() + 1 : current.minimumAlleleIndex());
            }

            // Checking on 0's new count and current.minAllele + 1 alleles.
            Assert.assertEquals(next.alleleCountFor(0),current.alleleCountFor(current.minimumAlleleIndex()) - 1);
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex() + 1),
                    current.alleleCountFor(current.minimumAlleleIndex() + 1) + 1);

            // Checks current.minAllele count
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex()),
                    current.minimumAlleleIndex() == 0 ? current.alleleCountAt(0) - 1 : 0);

            int totalCountSum = 0;
            final int[] expectedAlleleCountsByIndex = new int[Math.max(MAXIMUM_ALLELE_INDEX,next.maximumAlleleIndex()) + 1];
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                final int count = next.alleleCountAt(i);
                final int index = next.alleleIndexAt(i);
                expectedAlleleCountsByIndex[index] = count;
                // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
                Assert.assertEquals(next.alleleCountFor(index),count);
                totalCountSum += count;
                // Check on counts of, in theory, unaffected allele counts.
                if (index > current.minimumAlleleIndex() + 1)
                    Assert.assertEquals(next.alleleCountFor(index),current.alleleCountFor(index));
            }
            Assert.assertTrue(Arrays.equals(next.alleleCountsByIndex(Math.max(MAXIMUM_ALLELE_INDEX,next.maximumAlleleIndex())),expectedAlleleCountsByIndex));
            Assert.assertEquals(totalCountSum,ploidy);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));
            Assert.assertEquals(next.index(),current.index() + 1);
            Assert.assertEquals(next.ploidy(),ploidy);

            //Check asAlleleList.
            final List<Allele> expectedList = new ArrayList<>(ploidy);
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                for (int j = 0; j < next.alleleCountAt(i); j++) {
                    expectedList.add(testAlleles.get(next.alleleIndexAt(i)));
                }
            }
            Assert.assertEquals(next.asAlleleList(testAlleles),expectedList);

            current = next;
        }
    }

    private void testNextZeroPloidyIncrease() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(0);
        final GenotypeAlleleCounts next = first.clone();
        next.increase();
        Assert.assertEquals(first,next);
        Assert.assertEquals(first.compareTo(next),0);
        Assert.assertEquals(next.compareTo(first), 0);
        Assert.assertEquals(next.distinctAlleleCount(),0);
        Assert.assertEquals(next.ploidy(),0);
        Assert.assertEquals(next.index(),0);
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= 10; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex),expected);
        }
    }

    private void testNextOnePloidyIncrease() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(1);
        GenotypeAlleleCounts next = first;

        while (!next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts current = next.clone();
            next.increase();
            Assert.assertEquals(next.minimumAlleleIndex(),next.maximumAlleleIndex());
            Assert.assertEquals(next.minimumAlleleIndex(),current.minimumAlleleIndex() + 1);
            Assert.assertEquals(next.alleleCountAt(0),1);
            Assert.assertEquals(next.alleleIndexAt(0),next.minimumAlleleIndex());
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex()),0);
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex() + 1),-2);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex()),1);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex()+1),0);
            Assert.assertEquals(next.ploidy(),1);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));

            Assert.assertEquals(next.index(), current.index() + 1);
            Assert.assertEquals(next.ploidy(),current.ploidy());

            for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
                final int[] expected = new int[maximumAlleleIndex + 1];
                if (maximumAlleleIndex >= current.minimumAlleleIndex() + 1) expected[current.minimumAlleleIndex() + 1] = 1;
                Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex),expected);
            }
        }
    }

    private void testPloidyTwoOrMoreIncrease(final int ploidy) {
        if (ploidy < 2)
            throw new IllegalArgumentException();

        GenotypeAlleleCounts next = GenotypeAlleleCounts.first(ploidy);

        while (!next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts current = next.clone();
            next.increase();
            if (current.distinctAlleleCount() == 1) {
                Assert.assertEquals(next.maximumAlleleIndex(),current.maximumAlleleIndex() + 1);
                Assert.assertEquals(next.distinctAlleleCount(), 2 );
                Assert.assertEquals(next.minimumAlleleIndex(), 0 );
            } else {
                Assert.assertEquals(next.maximumAlleleIndex(),current.maximumAlleleIndex());
                Assert.assertEquals(next.minimumAlleleIndex(),current.alleleCountAt(0) > 1 ? 0
                        : current.alleleCountAt(0) == 1 ? current.minimumAlleleIndex() + 1 : current.minimumAlleleIndex());
            }

            // Checking on 0's new count and current.minAllele + 1 alleles.
            Assert.assertEquals(next.alleleCountFor(0),current.alleleCountFor(current.minimumAlleleIndex()) - 1);
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex() + 1),
                    current.alleleCountFor(current.minimumAlleleIndex() + 1) + 1);

            // Checks current.minAllele count
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex()),
                    current.minimumAlleleIndex() == 0 ? current.alleleCountAt(0) - 1 : 0);

            int totalCountSum = 0;
            final int[] expectedAlleleCountsByIndex = new int[Math.max(MAXIMUM_ALLELE_INDEX,next.maximumAlleleIndex()) + 1];
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                final int count = next.alleleCountAt(i);
                final int index = next.alleleIndexAt(i);
                expectedAlleleCountsByIndex[index] = count;
                // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
                Assert.assertEquals(next.alleleCountFor(index),count);
                totalCountSum += count;
                // Check on counts of, in theory, unaffected allele counts.
                if (index > current.minimumAlleleIndex() + 1)
                    Assert.assertEquals(next.alleleCountFor(index),current.alleleCountFor(index));
            }
            Assert.assertTrue(Arrays.equals(next.alleleCountsByIndex(Math.max(MAXIMUM_ALLELE_INDEX,next.maximumAlleleIndex())),expectedAlleleCountsByIndex));
            Assert.assertEquals(totalCountSum,ploidy);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));
            Assert.assertEquals(next.index(),current.index() + 1);
            Assert.assertEquals(next.ploidy(),ploidy);
        }
    }

    private static final int MAXIMUM_ALLELE_INDEX = 10;

    private static final List<Allele> testAlleles;

    static {
        final StringBuilder sb = new StringBuilder(51);
        testAlleles = new ArrayList<>(51);
        sb.append('A');
        for (int i = 0; i <= 50; i++) {
            testAlleles.add(Allele.create(sb.toString().getBytes(), i == 0));
            sb.append('A');
        }
    }

    private static final int[] PLOIDY = new int[] { 1, 2, 3, 7, 10};

    @DataProvider(name="ploidyData")
    public Object[][] ploidyData() {
        final Object[][] result = new Object[PLOIDY.length][];
        for (int i = 0; i < PLOIDY.length; i++)
            result[i] = new Object[] { PLOIDY[i ]};
        return result;
    }
}
