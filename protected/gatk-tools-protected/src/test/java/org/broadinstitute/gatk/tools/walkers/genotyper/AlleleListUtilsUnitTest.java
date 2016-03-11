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
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.genotyper.AlleleList;
import org.broadinstitute.gatk.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.gatk.utils.genotyper.AlleleListUtils;
import org.broadinstitute.gatk.utils.genotyper.IndexedAlleleList;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test {@link org.broadinstitute.gatk.utils.genotyper.AlleleListUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class AlleleListUtilsUnitTest {

    @Test(dataProvider = "singleAlleleListData")
    public void testAsList(final List<Allele> alleles1) {
         final Allele[] uniqueAlleles = new LinkedHashSet<>(alleles1).toArray(new Allele[0]);
         final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
         final List<Allele> asList = AlleleListUtils.asList(alleleList);
         final Allele[] asListArray = asList.toArray(new Allele[asList.size()]);
         Assert.assertTrue(Arrays.equals(uniqueAlleles,asListArray));
    }

    @Test(dataProvider = "singleAlleleListData")
    public void testIndexOfReference(final List<Allele> alleles1) {
        final Allele[] uniqueAlleles = new LinkedHashSet<>(alleles1).toArray(new Allele[0]);
        for (int i = 0; i < uniqueAlleles.length; i++) {
            final Allele[] actualAlleles = uniqueAlleles.clone();
            actualAlleles[i] = Allele.create(actualAlleles[i].getBases(),true);
            final AlleleList<Allele> alleleList = new IndexedAlleleList<>(actualAlleles);
            Assert.assertEquals(AlleleListUtils.indexOfReference(alleleList),i);
        }
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(uniqueAlleles);
        Assert.assertEquals(AlleleListUtils.indexOfReference(alleleList),-1);
    }

    @Test(dataProvider = "twoAlleleListData", dependsOnMethods={"testAsList"})
    public void testEquals(final List<Allele> alleles1, final List<Allele> alleles2) {
        final AlleleList<Allele> alleleList1 = new IndexedAlleleList<Allele>(alleles1);
        final AlleleList<Allele> alleleList2 = new IndexedAlleleList<Allele>(alleles2);
        Assert.assertTrue(AlleleListUtils.equals(alleleList1,alleleList1));
        Assert.assertTrue(AlleleListUtils.equals(alleleList2,alleleList2));
        Assert.assertEquals(AlleleListUtils.equals(alleleList1, alleleList2),
                Arrays.equals(AlleleListUtils.asList(alleleList1).toArray(new Allele[alleleList1.alleleCount()]),
                        AlleleListUtils.asList(alleleList2).toArray(new Allele[alleleList2.alleleCount()]))
        );
        Assert.assertEquals(AlleleListUtils.equals(alleleList1,alleleList2),
                            AlleleListUtils.equals(alleleList2,alleleList1));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods= "testEquals" )
    public void testSelfPermutation(final List<Allele> alleles1) {
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        final AlleleListPermutation<Allele> selfPermutation = AlleleListUtils.permutation(originalAlleleList,originalAlleleList);
        Assert.assertEquals(selfPermutation.fromSize(),originalAlleleList.alleleCount());
        Assert.assertEquals(selfPermutation.toSize(),originalAlleleList.alleleCount());
        Assert.assertTrue(selfPermutation.isNonPermuted());
        Assert.assertFalse(selfPermutation.isPartial());
        for (int i = 0; i < originalAlleleList.alleleCount(); i++) {
            Assert.assertEquals(selfPermutation.fromIndex(i), i);
            Assert.assertEquals(selfPermutation.toIndex(i),i);
            Assert.assertEquals(selfPermutation.fromList(),selfPermutation.toList());
            AlleleListUnitTester.assertAlleleList(originalAlleleList, selfPermutation.fromList());
        }
        Assert.assertTrue(AlleleListUtils.equals(selfPermutation,originalAlleleList));
    }

    private final Random rnd = Utils.getRandomGenerator();

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods = "testEquals")
    public void testSubsetPermutation(final List<Allele> alleles1) {
        final List<Allele> subsetAlleles = new ArrayList<>(alleles1.size());
        for (final Allele allele : alleles1)
            if (rnd.nextBoolean()) subsetAlleles.add(allele);
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        final AlleleList<Allele> targetAlleleList = new IndexedAlleleList<>(subsetAlleles);
        final AlleleListPermutation<Allele> subset = AlleleListUtils.permutation(originalAlleleList,targetAlleleList);
        if (originalAlleleList.alleleCount() == targetAlleleList.alleleCount())
            throw new SkipException("no real subset");
        Assert.assertTrue(subset.isPartial());
        Assert.assertFalse(subset.isNonPermuted());
        Assert.assertEquals(subset.fromSize(),originalAlleleList.alleleCount());
        Assert.assertEquals(subset.toSize(),targetAlleleList.alleleCount());
        AlleleListUnitTester.assertAlleleList(originalAlleleList,subset.fromList());
        AlleleListUnitTester.assertAlleleList(targetAlleleList,subset.toList());

        for (int i = 0; i < targetAlleleList.alleleCount(); i++)
            Assert.assertEquals(subset.fromIndex(i), originalAlleleList.alleleIndex(targetAlleleList.alleleAt(i)));

        for (int j = 0; j < originalAlleleList.alleleCount(); j++) {
            final Allele allele = originalAlleleList.alleleAt(j);
            Assert.assertEquals(subset.toIndex(j),targetAlleleList.alleleIndex(allele));
        }

        Assert.assertTrue(AlleleListUtils.equals(subset,targetAlleleList));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods = {"testAsList","testEquals"})
    public void testShufflePermutation(final List<Allele> alleles1) {
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        if (originalAlleleList.alleleCount() <= 1)
            throw new SkipException("non-shuffle allele-list");

        final Allele[] targetAlleleArray = AlleleListUtils.asList(originalAlleleList).toArray(new Allele[originalAlleleList.alleleCount()]);
        final int[] fromIndex = new int[targetAlleleArray.length];
        for (int i = 0; i < fromIndex.length; i++)
            fromIndex[i] = i;

        for (int i = 0; i < targetAlleleArray.length - 1; i++) {
            final int swapIndex = rnd.nextInt(targetAlleleArray.length - i - 1);
            final int otherIndex = fromIndex[swapIndex + i + 1];
            final Allele other = targetAlleleArray[swapIndex + i + 1];
            fromIndex[swapIndex + i + 1] = fromIndex[i];
            fromIndex[i] = otherIndex;
            targetAlleleArray[swapIndex + i + 1] = targetAlleleArray[i];
            targetAlleleArray[i] = other;
        }
        final AlleleList<Allele> targetAlleleList = new IndexedAlleleList<>(targetAlleleArray);

        final AlleleListPermutation<Allele> permutation = AlleleListUtils.permutation(originalAlleleList,targetAlleleList);
        Assert.assertFalse(permutation.isNonPermuted());
        AlleleListUnitTester.assertAlleleList(originalAlleleList,permutation.fromList());
        AlleleListUnitTester.assertAlleleList(targetAlleleList,permutation.toList());
        Assert.assertFalse(permutation.isPartial());
        Assert.assertEquals(permutation.fromSize(),originalAlleleList.alleleCount());
        Assert.assertEquals(permutation.toSize(),targetAlleleList.alleleCount());
        for (int i = 0; i < permutation.fromSize(); i++) {
            Assert.assertEquals(permutation.toIndex(i),targetAlleleList.alleleIndex(originalAlleleList.alleleAt(i)));
            Assert.assertEquals(permutation.fromIndex(i),originalAlleleList.alleleIndex(targetAlleleList.alleleAt(i)));
            Assert.assertEquals(permutation.fromIndex(i),fromIndex[i]);
        }
        Assert.assertTrue(AlleleListUtils.equals(permutation,targetAlleleList));

    }


    private List<Allele>[] alleleLists;

    @BeforeClass
    public void setUp() {
        alleleLists = new List[ALLELE_COUNT.length * MAX_ALLELE_LENGTH.length];
        int nextIndex = 0;
        for (int i = 0; i < ALLELE_COUNT.length; i++)
            for (int j = 0; j < MAX_ALLELE_LENGTH.length; j++)
                alleleLists[nextIndex++] = Arrays.asList(AlleleListUnitTester.generateRandomAlleles(ALLELE_COUNT[i], MAX_ALLELE_LENGTH[j]));
    }

    private static final int[] ALLELE_COUNT = { 0, 1, 5, 10, 20};

    private static final int[] MAX_ALLELE_LENGTH = { 1, 2, 3, 10 };

    @DataProvider(name="singleAlleleListData")
    public Object[][] singleAlleleListData() {
        final Object[][] result = new Object[alleleLists.length][];
        for (int i = 0; i < alleleLists.length; i++)
            result[i] = new Object[] { alleleLists[i]};
        return result;
    }

    @DataProvider(name="twoAlleleListData")
    public Object[][] twoAlleleListData() {
        final Object[][] result = new Object[alleleLists.length * alleleLists.length][];
        int index = 0;
        for (int i = 0; i < alleleLists.length; i++)
            for (int j = 0; j < alleleLists.length; j++)
                result[index++] = new Object[] { alleleLists[i], alleleLists[j]};
        return result;
    }







}
