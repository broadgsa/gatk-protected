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

package org.broadinstitute.gatk.utils.diffengine;


// the imports for unit testing.

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.diffengine.DiffElement;
import org.broadinstitute.gatk.utils.diffengine.DiffEngine;
import org.broadinstitute.gatk.utils.diffengine.DiffNode;
import org.broadinstitute.gatk.utils.diffengine.Difference;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffEngineUnitTest extends BaseTest {
    DiffEngine engine;

    @BeforeClass(enabled = true)
    public void createDiffEngine() {
        engine = new DiffEngine();
    }

    // --------------------------------------------------------------------------------
    //
    // Difference testing routines
    //
    // --------------------------------------------------------------------------------

    private class DifferenceTest extends TestDataProvider {
        public DiffElement tree1, tree2;
        public List<String> differences;

        private DifferenceTest(String tree1, String tree2) {
            this(tree1, tree2, Collections.<String>emptyList());
        }

        private DifferenceTest(String tree1, String tree2, String difference) {
            this(tree1, tree2, Arrays.asList(difference));
        }

        private DifferenceTest(String tree1, String tree2, List<String> differences) {
            super(DifferenceTest.class);
            this.tree1 = DiffNode.fromString(tree1);
            this.tree2 = DiffNode.fromString(tree2);
            this.differences = differences;
        }

        public String toString() {
            return String.format("tree1=%s tree2=%s diff=%s",
                    tree1.toOneLineString(), tree2.toOneLineString(), differences);
        }
    }

    @DataProvider(name = "trees")
    public Object[][] createTrees() {
        new DifferenceTest("A=X", "A=X");
        new DifferenceTest("A=X", "A=Y", "A:X!=Y");
        new DifferenceTest("A=X", "B=X", Arrays.asList("A:X!=MISSING", "B:MISSING!=X"));
        new DifferenceTest("A=(X=1)", "B=(X=1)", Arrays.asList("A:(X=1)!=MISSING", "B:MISSING!=(X=1)"));
        new DifferenceTest("A=(X=1)", "A=(X=1)");
        new DifferenceTest("A=(X=1 Y=2)", "A=(X=1 Y=2)");
        new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2 B=(Z=3))");
        new DifferenceTest("A=(X=1)", "A=(X=2)", "A.X:1!=2");
        new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2 B=(Z=4))", "A.B.Z:3!=4");
        new DifferenceTest("A=(X=1)", "A=(X=1 Y=2)", "A.Y:MISSING!=2");
        new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2)", "A.B:(Z=3)!=MISSING");
        return DifferenceTest.getTests(DifferenceTest.class);
    }

    @Test(enabled = true, dataProvider = "trees")
    public void testDiffs(DifferenceTest test) {
        logger.warn("Test tree1: " + test.tree1.toOneLineString());
        logger.warn("Test tree2: " + test.tree2.toOneLineString());

        List<Difference> diffs = engine.diff(test.tree1, test.tree2);
        logger.warn("Test expected diff : " + test.differences);
        logger.warn("Observed diffs     : " + diffs);
    }

    // --------------------------------------------------------------------------------
    //
    // Low-level routines for summarizing differences
    //
    // --------------------------------------------------------------------------------

    @Test(enabled = true)
    public void testLongestCommonPostfix() {
        testLongestCommonPostfixHelper("A", "A", 1);
        testLongestCommonPostfixHelper("A", "B", 0);
        testLongestCommonPostfixHelper("A.B", "A.B", 2);
        testLongestCommonPostfixHelper("A.B.C", "A.B.C", 3);
        testLongestCommonPostfixHelper("A.B.C", "X.B.C", 2);
        testLongestCommonPostfixHelper("A.B.C", "X.Y.C", 1);
        testLongestCommonPostfixHelper("A.B.C", "X.Y.Z", 0);
        testLongestCommonPostfixHelper("A.B.C", "A.X.C", 1);
        testLongestCommonPostfixHelper("A.B.C", "A.X.Z", 0);
        testLongestCommonPostfixHelper("A.B.C", "A.B.Z", 0);
    }

    public void testLongestCommonPostfixHelper(String p1, String p2, int expected) {
        String[] parts1 = p1.split("\\.");
        String[] parts2 = p2.split("\\.");
        int obs = DiffEngine.longestCommonPostfix(parts1, parts2);
        Assert.assertEquals(obs, expected, "p1=" + p1 + " p2=" + p2 + " failed");
    }

    @Test(enabled = true, dependsOnMethods = "testLongestCommonPostfix")
    public void testSummarizePath() {
        testSummarizePathHelper("A", "A", "A");
        testSummarizePathHelper("A", "B", "*");
        testSummarizePathHelper("A.B", "A.B", "A.B");
        testSummarizePathHelper("A.B", "X.B", "*.B");
        testSummarizePathHelper("A.B", "X.Y", "*.*");
        testSummarizePathHelper("A.B.C", "A.B.C", "A.B.C");
        testSummarizePathHelper("A.B.C", "X.B.C", "*.B.C");
        testSummarizePathHelper("A.B.C", "X.Y.C", "*.*.C");
        testSummarizePathHelper("A.B.C", "X.Y.Z", "*.*.*");
        testSummarizePathHelper("A.B.C", "A.X.C", "*.*.C");
        testSummarizePathHelper("A.B.C", "A.X.Z", "*.*.*");
        testSummarizePathHelper("A.B.C", "A.B.Z", "*.*.*");
    }

    public void testSummarizePathHelper(String p1, String p2, String expected) {
        String[] parts1 = DiffEngine.diffNameToPath(p1);
        String[] parts2 = DiffEngine.diffNameToPath(p2);
        int obs = DiffEngine.longestCommonPostfix(parts1, parts2);
        String path = DiffEngine.summarizedPath(parts2, obs);
        Assert.assertEquals(path, expected, "p1=" + p1 + " p2=" + p2 + " failed");
    }

    // --------------------------------------------------------------------------------
    //
    // High-level difference summary
    //
    // --------------------------------------------------------------------------------

    private class SummarizeDifferenceTest extends TestDataProvider {
        List<String> diffs = new ArrayList<String>();
        List<String> expecteds = new ArrayList<String>();

        public SummarizeDifferenceTest() { super(SummarizeDifferenceTest.class); }

        public SummarizeDifferenceTest addDiff(String... diffsToAdd) {
            diffs.addAll(Arrays.asList(diffsToAdd));
            return this;
        }

        public SummarizeDifferenceTest addSummary(String... expectedSummary) {
            expecteds.addAll(Arrays.asList(expectedSummary));
            return this;
        }

        public String toString() {
            return String.format("diffs=%s => expected=%s", diffs, expecteds);
        }

        public void test() {
            List<String[]> diffPaths = new ArrayList<String[]>(diffs.size());
            for ( String diff : diffs ) { diffPaths.add(DiffEngine.diffNameToPath(diff)); }

            List<Difference> sumDiffs = engine.summarizedDifferencesOfPathsFromString(diffs);

            Assert.assertEquals(sumDiffs.size(), expecteds.size(), "Unexpected number of summarized differences: " + sumDiffs);

            for ( int i = 0; i < sumDiffs.size(); i++ ) {
                Difference sumDiff = sumDiffs.get(i);
                String expected = expecteds.get(i);
                String[] pathCount = expected.split(":");
                String path = pathCount[0];
                int count = Integer.valueOf(pathCount[1]);
                Assert.assertEquals(sumDiff.getPath(), path, "Unexpected path at: " + expected + " obs=" + sumDiff + " all=" + sumDiffs);
                Assert.assertEquals(sumDiff.getCount(), count, "Unexpected counts at: " + expected + " obs=" + sumDiff + " all=" + sumDiffs);
            }
        }
    }

    @DataProvider(name = "summaries")
    public Object[][] createSummaries() {
        new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
        new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
        new SummarizeDifferenceTest().addDiff("A", "A", "A").addSummary("A:3");
        new SummarizeDifferenceTest().addDiff("A", "A", "A", "B").addSummary("A:3", "B:1");
        new SummarizeDifferenceTest().addDiff("A", "A", "A", "B", "B").addSummary("A:3", "B:2");
        new SummarizeDifferenceTest().addDiff("A", "A", "A", "B", "B", "C").addSummary("A:3", "B:2", "C:1");
        new SummarizeDifferenceTest().addDiff("A.X", "A.X").addSummary("A.X:2");
        new SummarizeDifferenceTest().addDiff("A.X", "A.X", "B.X").addSummary("*.X:3", "A.X:2", "B.X:1");
        new SummarizeDifferenceTest().addDiff("A.X", "A.X", "B.X", "B.X").addSummary("*.X:4", "A.X:2", "B.X:2");
        new SummarizeDifferenceTest().addDiff("A.B.C", "X.B.C").addSummary("*.B.C:2", "A.B.C:1", "X.B.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "X.Y.C", "X.Y.C").addSummary("*.*.C:3", "X.Y.C:2", "A.B.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "A.X.C", "X.Y.C").addSummary("*.*.C:3", "A.B.C:1", "A.X.C:1", "X.Y.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "A.X.C", "B.X.C").addSummary("*.*.C:3", "*.X.C:2", "A.B.C:1", "A.X.C:1", "B.X.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "A.X.C", "B.X.C", "B.X.C").addSummary("*.*.C:4", "*.X.C:3", "B.X.C:2", "A.B.C:1", "A.X.C:1");

        return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
    }


    @Test(enabled = true, dependsOnMethods = "testSummarizePath", dataProvider = "summaries")
    public void testSummarizeDifferences(SummarizeDifferenceTest test) {
        test.test();
    }
}