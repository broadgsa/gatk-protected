/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 1/31/13
 */

public class KBestPathsUnitTest {
    @DataProvider(name = "BasicBubbleDataProvider")
    public Object[][] makeBasicBubbleDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final int refBubbleLength : Arrays.asList(1, 5, 10) ) {
            for ( final int altBubbleLength : Arrays.asList(1, 5, 10) ) {
                tests.add(new Object[]{refBubbleLength, altBubbleLength});
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BasicBubbleDataProvider")
    public void testBasicBubbleData(final int refBubbleLength, final int altBubbleLength) {
        // Construct the assembly graph
        DeBruijnAssemblyGraph graph = new DeBruijnAssemblyGraph();
        final int KMER_LENGTH = 3;
        final String preRef = "ATGG";
        final String postRef = new String(Utils.dupBytes((byte) 'A', KMER_LENGTH-1)) + "GGGGC";

        DeBruijnVertex v = new DeBruijnVertex(preRef.getBytes(), KMER_LENGTH);
        DeBruijnVertex v2Ref = new DeBruijnVertex(Utils.dupBytes((byte) 'A', refBubbleLength+KMER_LENGTH-1), KMER_LENGTH);
        DeBruijnVertex v2Alt = new DeBruijnVertex(ArrayUtils.addAll(Utils.dupBytes((byte) 'A', altBubbleLength + KMER_LENGTH - 1 - 1), Utils.dupBytes((byte) 'T',1)), KMER_LENGTH);
        DeBruijnVertex v3 = new DeBruijnVertex(postRef.getBytes(), KMER_LENGTH);

        graph.addVertex(v);
        graph.addVertex(v2Ref);
        graph.addVertex(v2Alt);
        graph.addVertex(v3);
        graph.addEdge(v, v2Ref, new DeBruijnEdge(true, 10));
        graph.addEdge(v2Ref, v3, new DeBruijnEdge(true, 10));
        graph.addEdge(v, v2Alt, new DeBruijnEdge(false, 5));
        graph.addEdge(v2Alt, v3, new DeBruijnEdge(false, 5));

        // Construct the test path
        KBestPaths.Path path = new KBestPaths.Path(v, graph);
        path = new KBestPaths.Path(path, graph.getEdge(v, v2Alt));
        path = new KBestPaths.Path(path, graph.getEdge(v2Alt, v3));

        // Construct the actual cigar string implied by the test path
        Cigar expectedCigar = new Cigar();
        expectedCigar.add(new CigarElement(preRef.length(), CigarOperator.M));
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength,CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength,CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength,CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(postRef.length() - (KMER_LENGTH - 1), CigarOperator.M));

        Assert.assertEquals(path.calculateCigar().toString(), AlignmentUtils.consolidateCigar(expectedCigar).toString(), "Cigar string mismatch");
    }


    @DataProvider(name = "TripleBubbleDataProvider")
    public Object[][] makeTripleBubbleDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final int refBubbleLength : Arrays.asList(1, 5, 10) ) {
            for ( final int altBubbleLength : Arrays.asList(1, 5, 10) ) {
                for ( final boolean offRefBeginning : Arrays.asList(false) ) {
                    for ( final boolean offRefEnding : Arrays.asList(true, false) ) {
                        tests.add(new Object[]{refBubbleLength, altBubbleLength, offRefBeginning, offRefEnding});
                    }
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TripleBubbleDataProvider")
    public void testTripleBubbleData(final int refBubbleLength, final int altBubbleLength, final boolean offRefBeginning, final boolean offRefEnding) {
        // Construct the assembly graph
        DeBruijnAssemblyGraph graph = new DeBruijnAssemblyGraph();
        final int KMER_LENGTH = 3;
        final String preAltOption = "ATCGATCGATCGATCGATCG";
        final String postAltOption = "CCCC";
        final String preRef = "ATGG";
        final String postRef = new String(Utils.dupBytes((byte) 'A', KMER_LENGTH-1)) + "GGCCG";
        final String midRef1 = new String(Utils.dupBytes((byte) 'A', KMER_LENGTH-1)) + "TTCCT";
        final String midRef2 = new String(Utils.dupBytes((byte) 'A', KMER_LENGTH-1)) + "CCCAAAAAAAAAAAA";

        DeBruijnVertex preV = new DeBruijnVertex(preAltOption.getBytes(), KMER_LENGTH);
        DeBruijnVertex v = new DeBruijnVertex(preRef.getBytes(), KMER_LENGTH);
        DeBruijnVertex v2Ref = new DeBruijnVertex(Utils.dupBytes((byte) 'A', refBubbleLength+KMER_LENGTH-1), KMER_LENGTH);
        DeBruijnVertex v2Alt = new DeBruijnVertex(ArrayUtils.addAll(Utils.dupBytes((byte) 'A', altBubbleLength + KMER_LENGTH - 1 - 1), Utils.dupBytes((byte) 'T',1)), KMER_LENGTH);
        DeBruijnVertex v4Ref = new DeBruijnVertex(Utils.dupBytes((byte) 'C', refBubbleLength+KMER_LENGTH-1), KMER_LENGTH);
        DeBruijnVertex v4Alt = new DeBruijnVertex(ArrayUtils.addAll(Utils.dupBytes((byte) 'C', altBubbleLength + KMER_LENGTH - 1 - 1), Utils.dupBytes((byte) 'T',1)), KMER_LENGTH);
        DeBruijnVertex v6Ref = new DeBruijnVertex(Utils.dupBytes((byte) 'G', refBubbleLength+KMER_LENGTH-1), KMER_LENGTH);
        DeBruijnVertex v6Alt = new DeBruijnVertex(ArrayUtils.addAll(Utils.dupBytes((byte) 'G', altBubbleLength + KMER_LENGTH - 1 - 1), Utils.dupBytes((byte) 'T',1)), KMER_LENGTH);
        DeBruijnVertex v3 = new DeBruijnVertex(midRef1.getBytes(), KMER_LENGTH);
        DeBruijnVertex v5 = new DeBruijnVertex(midRef2.getBytes(), KMER_LENGTH);
        DeBruijnVertex v7 = new DeBruijnVertex(postRef.getBytes(), KMER_LENGTH);
        DeBruijnVertex postV = new DeBruijnVertex(postAltOption.getBytes(), KMER_LENGTH);

        graph.addVertex(preV);
        graph.addVertex(v);
        graph.addVertex(v2Ref);
        graph.addVertex(v2Alt);
        graph.addVertex(v3);
        graph.addVertex(v4Ref);
        graph.addVertex(v4Alt);
        graph.addVertex(v5);
        graph.addVertex(v6Ref);
        graph.addVertex(v6Alt);
        graph.addVertex(v7);
        graph.addVertex(postV);
        graph.addEdge(preV, v, new DeBruijnEdge(false, 1));
        graph.addEdge(v, v2Ref, new DeBruijnEdge(true, 10));
        graph.addEdge(v2Ref, v3, new DeBruijnEdge(true, 10));
        graph.addEdge(v, v2Alt, new DeBruijnEdge(false, 5));
        graph.addEdge(v2Alt, v3, new DeBruijnEdge(false, 5));
        graph.addEdge(v3, v4Ref, new DeBruijnEdge(true, 10));
        graph.addEdge(v4Ref, v5, new DeBruijnEdge(true, 10));
        graph.addEdge(v3, v4Alt, new DeBruijnEdge(false, 5));
        graph.addEdge(v4Alt, v5, new DeBruijnEdge(false, 5));
        graph.addEdge(v5, v6Ref, new DeBruijnEdge(true, 11));
        graph.addEdge(v6Ref, v7, new DeBruijnEdge(true, 11));
        graph.addEdge(v5, v6Alt, new DeBruijnEdge(false, 55));
        graph.addEdge(v6Alt, v7, new DeBruijnEdge(false, 55));
        graph.addEdge(v7, postV, new DeBruijnEdge(false, 1));

        // Construct the test path
        KBestPaths.Path path = new KBestPaths.Path( (offRefBeginning ? preV : v), graph);
        if( offRefBeginning ) {
            path = new KBestPaths.Path(path, graph.getEdge(preV, v));
        }
        path = new KBestPaths.Path(path, graph.getEdge(v, v2Alt));
        path = new KBestPaths.Path(path, graph.getEdge(v2Alt, v3));
        path = new KBestPaths.Path(path, graph.getEdge(v3, v4Ref));
        path = new KBestPaths.Path(path, graph.getEdge(v4Ref, v5));
        path = new KBestPaths.Path(path, graph.getEdge(v5, v6Alt));
        path = new KBestPaths.Path(path, graph.getEdge(v6Alt, v7));
        if( offRefEnding ) {
            path = new KBestPaths.Path(path, graph.getEdge(v7,postV));
        }

        // Construct the actual cigar string implied by the test path
        Cigar expectedCigar = new Cigar();
        if( offRefBeginning ) {
            expectedCigar.add(new CigarElement(preAltOption.length(), CigarOperator.I));
        }
        expectedCigar.add(new CigarElement(preRef.length() - (KMER_LENGTH - 1), CigarOperator.M));
        // first bubble
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength,CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength,CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength,CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(midRef1.length() - (KMER_LENGTH - 1), CigarOperator.M));
        // second bubble is ref path
        expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        expectedCigar.add(new CigarElement(midRef2.length() - (KMER_LENGTH - 1), CigarOperator.M));
        // third bubble
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength,CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength,CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength,CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(postRef.length() - (KMER_LENGTH - 1), CigarOperator.M));
        if( offRefEnding ) {
            expectedCigar.add(new CigarElement(postAltOption.length() - (KMER_LENGTH - 1), CigarOperator.I));
        }

        Assert.assertEquals(path.calculateCigar().toString(), AlignmentUtils.consolidateCigar(expectedCigar).toString(), "Cigar string mismatch");
    }
}
