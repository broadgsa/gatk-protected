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

package org.broadinstitute.gatk.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.HaplotypeGraph;
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests for {@link HaplotypeGraph}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class HaplotypeGraphUnitTest extends BaseTest {

   @Test(dataProvider="buildByStringDataProvider")
   public void testBuildByString(final String string, final int kmerSize, final int vertexCount, final int edgeCount) {
       final HaplotypeGraph g = new HaplotypeGraph(string);
       Assert.assertEquals(g.getKmerSize(),kmerSize,g.toString());
       Assert.assertEquals(g.vertexSet().size(),vertexCount,g.toString());
       Assert.assertEquals(g.edgeSet().size(),edgeCount,g.toString());
   }

   @Test(dataProvider="equalTestDataProvider")
   public void testEquals(final HaplotypeGraph one, final HaplotypeGraph two, final boolean expected) {
      Assert.assertEquals(one.equals(two),expected);
   }

   @Test(dataProvider="equalTestDataProvider")
   public void testEqualReciprocal(final HaplotypeGraph one, final HaplotypeGraph two, final boolean expected) {
       Assert.assertEquals(two.equals(one),expected);
   }


   @Test(dataProvider="equalTestDataProvider")
   public void testReflexibeEquals(final HaplotypeGraph one, final HaplotypeGraph two,
                                   @SuppressWarnings("unused") final boolean expected) {
       Assert.assertTrue(one.equals(one));
       Assert.assertTrue(two.equals(two));
   }

   @Test(dataProvider="mergingCommonChainsDataProvider")
   public void testMergingCommonChains(final HaplotypeGraph actual, HaplotypeGraph expected) {


        final Map<Kmer,MultiDeBruijnVertex> beforeMap = new HashMap<>(actual.uniqueKmerMap());
        actual.mergeCommonChains();
        final Map<Kmer,MultiDeBruijnVertex> afterMap = new HashMap<>(actual.uniqueKmerMap());
        final Map<Kmer,MultiDeBruijnVertex> mergedMap = new HashMap<>(expected.uniqueKmerMap());

        Assert.assertEquals(actual, expected,""+actual.vertexSet() + " EDGES " + actual.edgeSet());
        Assert.assertEquals(beforeMap.size(),afterMap.size());
        Assert.assertEquals(afterMap.size(),mergedMap.size());
        for (final Kmer k : beforeMap.keySet()) {
            Assert.assertTrue(afterMap.containsKey(k));
            Assert.assertTrue(mergedMap.containsKey(k));
            final byte[] seq1 = beforeMap.get(k).getSequence();
            final byte[] seq2 = afterMap.get(k).getSequence();
            final byte[] seq3 = mergedMap.get(k).getSequence();
            Assert.assertEquals(seq1.length,seq2.length);
            Assert.assertEquals(seq2.length,seq3.length);
            for (int i = 0; i < seq3.length; i++) {
                final byte bk = k.base(i);
                final byte b1 = seq1[i];
                final byte b2 = seq2[i];
                final byte b3 = seq3[i];
                final byte theByte = b1 == 'N' || b2 == 'N' || b3 == 'N' ? (byte)'N' : b1;
                if (theByte == 'N') continue;
                Assert.assertEquals(b1,b2);
                Assert.assertEquals(b2,b3);
                Assert.assertEquals(bk,b1);
            }
        }
   }


   @DataProvider(name="mergingCommonChainsDataProvider")
   public Iterator<Object[]> mergingCommonChainsDataProvider() {
      final List<Object[]> list = new LinkedList<>();
      for (int i = 0; i < MERGING_COMMON_CHAINS_DATA.length; i += 2) {
          final HaplotypeGraph before = new HaplotypeGraph(MERGING_COMMON_CHAINS_DATA[i]);
          final HaplotypeGraph after = new HaplotypeGraph(MERGING_COMMON_CHAINS_DATA[i+1]);
          list.add(new Object[] { before , after});
      }
      return list.iterator();
   }

   @DataProvider(name="equalTestDataProvider")
   public Iterator<Object[]> equalsTestDataProvider() {
      final List<Object[]> result = new LinkedList<>();
      for (int i = 0; i < EQUAL_TEST_DATA.length; i += 3) {
          final HaplotypeGraph g1 = new HaplotypeGraph(EQUAL_TEST_DATA[i]);
          final HaplotypeGraph g2 = new HaplotypeGraph(EQUAL_TEST_DATA[i+1]);
          final boolean outcome = Boolean.parseBoolean(EQUAL_TEST_DATA[i+2]);
          result.add(new Object[] { g1, g2, outcome});
      }
      return result.iterator();
   }

   @DataProvider(name="buildByStringDataProvider")
   public Iterator<Object[]> buildByStringDataProvider() {
      return Arrays.asList(BUILD_BY_STRING_TEST_DATA).iterator();
   }

   private static final Object[][] BUILD_BY_STRING_TEST_DATA = new Object[][] {
           {"[ks=3]{REF: ACT}",3,1,0},
           {"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                         "{ (3) -> A -> G ->          (2) }" +
                                "{  (1) -> A -> G ->  (2) }",3,8,9},
           {"[ks=3]{REF: ACT -> C(1) -> G}{ACT -> C(1) -> G}{ACT -> C(1) -> G}",3,5,4} ,
           {"[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> C -> G -> T }" +
                              "{A(1) -> T -> A(2) }",3,8,8}  ,
           {"[ks=3]{REF: ACT -> A -> T(2) -> C -> A -> G -> T -> A -> C -> G -> T -> A(1) -> T}" +
                   "{ ACT -> A -> T(2) -> C -> A -> G -> T -> A -> C -> G -> T -> A(1) -> T}",3,15,14} ,
           {"[ks=3]{REF: ACT -> A -> T    -> C -> A -> G -> T -> A -> C -> G -> T -> A    -> T}",3,13,12},
           {"[ks=3]{REF: ACT -> A -> T(1) }" +
                   "{ ACT -> A -> T(1) }",3,5,4},
           {"[ks=3]{REF: TTT -> A(1) -> C -> T(2)}{ A(1) -> T(2) } ",3,4,4}
   };

   private static final String[] EQUAL_TEST_DATA = new String[] {
           "[ks=3]{REF: ACT}","[ks=3]{REF: ACT}", "true",
           "[ks=3]{REF: TCA}","[ks=3]{REF: ACT}", "false",
           "[ks=4]{REF: ACTG}","[ks=3]{REF: ACT}", "false",
           "[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                        "{ (3) -> A -> G ->          (2) }" +
                               "{  (1) -> A -> G ->  (2) }"
                   ,"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                                                                 "{  (1) -> A -> G ->  (2) }" +
                                                          "{ (3) -> A -> G ->          (2) }", "true",
           "[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                        "{ (3) -> A -> T ->          (2) }" +
                              "{  (1) -> A -> G ->  (2) }"
                    ,"[ks=3]{REF: ACT(3) -> T(1)      ->      G(2) -> A}" +
                     "{  (1) -> A -> G ->  (2) }" +
                     "{ (3) -> A -> G ->          (2) }", "true",
           "[ks=3]{REF: ACT -> G -> C(2) }{ ACT -> T -> C(2) }","[ks=3]{REF: ACT -> T -> C(2) }{ ACT -> G -> C(2) }","false",

   };

   private static final String[] MERGING_COMMON_CHAINS_DATA = new String[] {  // pairs before and after.
           "[ks=3]{REF: ACT -> A(1) -> G -> A -> G(2) -> T }" +
                             "{A(1) -> T -> A -> G(2) }",
           "[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> G -> T }" +
                             "{A(1) -> T -> A(2) }",

           "[ks=3]{REF: ACT -> A(1) -> G -> A -> C -> G(2) -> T }" +
                             "{A(1) -> T -> A -> C -> G(2) }",
           "[ks=3]{REF: ACT -> A(1) -> G -> A(2) -> C -> G -> T }" +
                             "{A(1) -> T -> A(2) }",

           "[ks=3]{REF: ACT -> A -> T(1) -> C -> A -> G -> T -> A -> C(2) -> G -> T -> A}" +
                                 "{ T(1) -> A -> A -> G -> T -> A -> C(2) }",
           "[ks=3]{REF: ACT -> A -> T(1) -> C -> A(2) -> G -> T -> A -> C -> G -> T -> A}" +
                                 "{ T(1) -> A -> A(2) } ",

//           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A(1)}" +
//                     "{ ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A(1)}" ,
//           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A}" ,

           "[ks=3]{REF: ACT -> A -> T(1) }" +
                     "{ AGT -> A -> T(1) }" ,
           "[ks=3]{REF: ACT -> A(1) -> T }" +
                   "{ AGT -> A(1)  }"  ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T}" ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T}" ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T(1)}" + "{ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T(1)}" ,
           "[ks=3]{REF: ACT -> A -> T -> C -> A -> G -> T -> A -> C -> G -> T -> A -> T}" ,

           "[ks=3]{REF: TTT -> T -> T -> T -> T -> T -> T -> T -> T -> T -> T(1) -> T -> T -> T(2) -> T -> T}"
                  + "{  TTT -> T -> T -> T -> T -> T -> T -> T -> T -> T -> T(1) -> G -> T -> T -> T -> T -> T(2) -> T -> T}",
           "[ks=3]{REF: TTT -> T -> T -> T -> T -> T -> T -> T -> T -> T -> T(1) -> T(2) -> T -> T -> T -> T}"
                  +                                                      "{ T(1) -> G -> T -> T -> T(2) }",

           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C -> C -> T(2)}" +
                                 "{ G(1) -> T -> C -> C -> T(2)}" +
                                 "{ G(1) -> G -> C -> C -> T(2)}" +
                                 "{ G(1) -> C -> T(2)} ",
           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C(2) -> C(3) -> T }" +
                                 "{  G(1) -> T -> C(2) }" +
                                 "{  G(1) -> G -> C(2) }" +
                                 "{  G(1) -> C(3) }",

           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C -> G}{ TTT -> T -> G(1) -> G -> C -> G}",
           "[ks=3]{REF: TTT -> T -> G(1) -> A -> C -> G}{ G(1) -> G -> C -> G}",

   };

}
