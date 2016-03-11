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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.manager;


// the imports for unit testing.


import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.util.*;


public class StratificationManagerUnitTest extends BaseTest {
    @BeforeClass
    public void init() throws FileNotFoundException {
    }

    // --------------------------------------------------------------------------------
    //
    // Basic tests Provider
    //
    // --------------------------------------------------------------------------------

    private class StratificationStatesTestProvider extends TestDataProvider {
        final List<List<Object>> allStates = new ArrayList<List<Object>>();
        final List<IntegerStratifier> asSetOfStates = new ArrayList<IntegerStratifier>();
        final int nStates;
        
        public StratificationStatesTestProvider(final List<Integer> ... allStates) {
            super(StratificationStatesTestProvider.class);
            
            for ( List<Integer> states : allStates ) {
                this.allStates.add(new ArrayList<Object>(states));
            }

            for ( List<Object> states : this.allStates ) { 
                asSetOfStates.add(new IntegerStratifier(states));
            }
            this.nStates = Utils.nCombinations(allStates);

            setName(getName());
        }

        private String getName() {
            StringBuilder b = new StringBuilder();
            int c = 1;
            for ( List<Object> state : allStates )
                b.append(String.format("%d = [%s] ", c++, Utils.join(",", state)));
            return b.toString();
        }
        
        public List<IntegerStratifier> getStateSpaceList() {
            return asSetOfStates;
        }
        
        public ArrayList<Integer> values() {
            final ArrayList<Integer> l = new ArrayList<Integer>();
            for ( int i = 0; i < nStates; i++ )
                l.add(i);
            return l;
        }
            
        public Queue<List<Object>> getAllCombinations() {
            return getAllCombinations(new LinkedList<List<Object>>(allStates));
        }

        private Queue<List<Object>> getAllCombinations(Queue<List<Object>> states) {
            if ( states.isEmpty() ) 
                return new LinkedList<List<Object>>();
            else {
                List<Object> head = states.poll();
                Queue<List<Object>> substates = getAllCombinations(states);
                Queue<List<Object>> newStates = new LinkedList<List<Object>>();
                for ( final Object e : head) {
                    if ( substates.isEmpty() ) {
                        newStates.add(new LinkedList<Object>(Collections.singleton(e)));
                    } else {
                        for ( final List<Object> state : substates ) {
                            List<Object> newState = new LinkedList<Object>();
                            newState.add(e);
                            newState.addAll(state);
                            newStates.add(newState);
                        }
                    }
                }
                return newStates;
            }
        }
    }

    private class IntegerStratifier implements Stratifier {
        final List<Object> integers;

        private IntegerStratifier(final List<Object> integers) {
            this.integers = integers;
        }
        
        @Override
        public List<Object> getAllStates() {
            return integers;
        }
    }

    @DataProvider(name = "StratificationStatesTestProvider")
    public Object[][] makeStratificationStatesTestProvider() {
        new StratificationStatesTestProvider(Arrays.asList(0));
        new StratificationStatesTestProvider(Arrays.asList(0, 1));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3), Arrays.asList(4, 5));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4), Arrays.asList(5, 6));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4, 5), Arrays.asList(6));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4, 5), Arrays.asList(6, 7));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3), Arrays.asList(4, 5), Arrays.asList(6, 7));
        return StratificationStatesTestProvider.getTests(StratificationStatesTestProvider.class);
    }
    
    private final StratificationManager<IntegerStratifier, Integer> createManager(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> manager = new StratificationManager<IntegerStratifier, Integer>(cfg.getStateSpaceList());
        List<Integer> values = cfg.values();
        for ( int i = 0; i < cfg.nStates; i++ )
            manager.set(i, values.get(i));
        
        Assert.assertEquals(manager.values(), values, "Values not equal");
        
        return manager;
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testLeafCount(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        
        Assert.assertEquals(stratificationManager.size(), cfg.nStates);
        
        int nLeafs = 0;
        for ( final StratNode node : stratificationManager.getRoot() ) {
            if ( node.isLeaf() )
                nLeafs++;
        }
        Assert.assertEquals(nLeafs, cfg.nStates, "Unexpected number of leaves");
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testKeys(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        final Set<Integer> seenKeys = new HashSet<Integer>(cfg.nStates);
        for ( final StratNode node : stratificationManager.getRoot() ) {
            if ( node.isLeaf() ) {
                Assert.assertFalse(seenKeys.contains(node.getKey()), "Already seen the key");
                seenKeys.add(node.getKey());
            }
        }
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testFindSingleKeys(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        final Set<Integer> seenKeys = new HashSet<Integer>(cfg.nStates);
        for ( List<Object> state : cfg.getAllCombinations() ) {
            final int key = stratificationManager.getKey(state);
            Assert.assertFalse(seenKeys.contains(key), "Already saw state mapping to this key");
            Assert.assertTrue(stratificationManager.containsKey(state));
            seenKeys.add(key);

            // test value
            Assert.assertEquals(stratificationManager.get(key), cfg.values().get(key));
            Assert.assertEquals(stratificationManager.get(state), cfg.values().get(key));

            state.set(0, 12345); // not present
            Assert.assertEquals(stratificationManager.getKey(state), -1);
            Assert.assertFalse(stratificationManager.containsKey(state));
        }
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testFindMultipleKeys(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        final List<List<Object>> states = new ArrayList<List<Object>>(cfg.allStates);
        final Set<Integer> keys = stratificationManager.getKeys(states);
        Assert.assertEquals(keys.size(), cfg.nStates, "Find all states didn't find all of the expected unique keys");

        final Queue<List<Object>> combinations = cfg.getAllCombinations();
        while ( ! combinations.isEmpty() ) {
            List<Object> first = combinations.poll();
            List<Object> second = combinations.peek();
            if ( second != null ) {
                List<List<Object>> combined = StratificationManager.combineStates(first, second);
                int nExpectedKeys = Utils.nCombinations(combined);

                final int key1 = stratificationManager.getKey(first);
                final int key2 = stratificationManager.getKey(second);
                final Set<Integer> keysCombined = stratificationManager.getKeys(combined);
            
                Assert.assertTrue(keysCombined.contains(key1), "couldn't find key in data set");
                Assert.assertTrue(keysCombined.contains(key2), "couldn't find key in data set");
                
                Assert.assertEquals(keysCombined.size(), nExpectedKeys);
            }
        }
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testMapSet(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> stratificationManager = createManager(cfg);
        stratificationManager.set(0, -1);
        Assert.assertEquals((int)stratificationManager.get(0), -1);
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testStratifierByKey(StratificationStatesTestProvider cfg) {
        final StratificationManager<IntegerStratifier, Integer> manager = createManager(cfg);
        for ( int key = 0; key < cfg.nStates; key++ ) {
            List<Pair<IntegerStratifier, Object>> stratsAndStates = manager.getStratsAndStatesForKey(key);
            final List<Object> strats = manager.getStatesForKey(key);
            Assert.assertEquals((int)manager.get(strats), key, "Key -> strats -> key failed to return same key");

            for ( int i = 0; i < strats.size(); i++ ) {
                Assert.assertEquals(stratsAndStates.get(i).getSecond(), strats.get(i), "Strats and StratsAndStates differ");
            }
        }
    }
}