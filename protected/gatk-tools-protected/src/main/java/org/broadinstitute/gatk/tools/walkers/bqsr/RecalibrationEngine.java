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

package org.broadinstitute.gatk.tools.walkers.bqsr;

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.engine.recalibration.ReadCovariates;
import org.broadinstitute.gatk.engine.recalibration.RecalDatum;
import org.broadinstitute.gatk.engine.recalibration.RecalUtils;
import org.broadinstitute.gatk.engine.recalibration.RecalibrationTables;
import org.broadinstitute.gatk.utils.collections.NestedIntegerArray;
import org.broadinstitute.gatk.utils.recalibration.*;
import org.broadinstitute.gatk.engine.recalibration.covariates.Covariate;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;

public class RecalibrationEngine {
    final protected Covariate[] covariates;
    final private int numReadGroups;
    final private PrintStream maybeLogStream;
    final private boolean lowMemoryMode;

    /**
     * Has finalizeData() been called?
     */
    private boolean finalized = false;

    /**
     * The final (merged, etc) recalibration tables, suitable for downstream analysis.
     */
    private RecalibrationTables finalRecalibrationTables = null;

    private final List<RecalibrationTables> recalibrationTablesList = new LinkedList<RecalibrationTables>();

    private final ThreadLocal<RecalibrationTables> threadLocalTables = new ThreadLocal<RecalibrationTables>() {
        private synchronized RecalibrationTables makeAndCaptureTable() {
            final RecalibrationTables newTable = new RecalibrationTables(covariates, numReadGroups, maybeLogStream);
            recalibrationTablesList.add(newTable);
            return newTable;
        }

        @Override
        protected synchronized RecalibrationTables initialValue() {
            if ( lowMemoryMode ) {
                return recalibrationTablesList.isEmpty() ? makeAndCaptureTable() : recalibrationTablesList.get(0);
            } else {
                return makeAndCaptureTable();
            }
        }
    };

    /**
     * Get a recalibration table suitable for updating the underlying RecalDatums
     *
     * May return a thread-local version, or a single version, depending on the initialization
     * arguments of this instance.
     *
     * @return updated tables
     */
    protected RecalibrationTables getUpdatableRecalibrationTables() {
        return threadLocalTables.get();
    }

    /**
     * Initialize the recalibration engine
     *
     * Called once before any calls to updateDataForRead are made.  The engine should prepare itself
     * to handle any number of updateDataForRead calls containing ReadRecalibrationInfo containing
     * keys for each of the covariates provided.
     *
     * The engine should collect match and mismatch data into the recalibrationTables data.
     *
     * @param covariates an array of the covariates we'll be using in this engine, order matters
     * @param numReadGroups the number of read groups we should use for the recalibration tables
     * @param maybeLogStream an optional print stream for logging calls to the nestedhashmap in the recalibration tables
     */
    public RecalibrationEngine(final Covariate[] covariates, final int numReadGroups, final PrintStream maybeLogStream, final boolean enableLowMemoryMode) {
        if ( covariates == null ) throw new IllegalArgumentException("Covariates cannot be null");
        if ( numReadGroups < 1 ) throw new IllegalArgumentException("numReadGroups must be >= 1 but got " + numReadGroups);

        this.covariates = covariates.clone();
        this.numReadGroups = numReadGroups;
        this.maybeLogStream = maybeLogStream;
        this.lowMemoryMode = enableLowMemoryMode;
    }

    /**
     * Update the recalibration statistics using the information in recalInfo
     * @param recalInfo data structure holding information about the recalibration values for a single read
     */
    @Requires("recalInfo != null")
    public void updateDataForRead( final ReadRecalibrationInfo recalInfo ) {
        final GATKSAMRecord read = recalInfo.getRead();
        final ReadCovariates readCovariates = recalInfo.getCovariatesValues();
        final RecalibrationTables tables = getUpdatableRecalibrationTables();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = tables.getQualityScoreTable();

        for( int offset = 0; offset < read.getReadBases().length; offset++ ) {
            if( ! recalInfo.skip(offset) ) {

                for (final EventType eventType : EventType.values()) {
                    final int[] keys = readCovariates.getKeySet(offset, eventType);
                    final int eventIndex = eventType.ordinal();
                    final byte qual = recalInfo.getQual(eventType, offset);
                    final double isError = recalInfo.getErrorFraction(eventType, offset);

                    RecalUtils.incrementDatumOrPutIfNecessary(qualityScoreTable, qual, isError, keys[0], keys[1], eventIndex);

                    for (int i = 2; i < covariates.length; i++) {
                        if (keys[i] < 0)
                            continue;

                        RecalUtils.incrementDatumOrPutIfNecessary(tables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
                    }
                }
            }
        }
    }


    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to updateDataForRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public void finalizeData() {
        if ( finalized ) throw new IllegalStateException("FinalizeData() has already been called");

        // merge all of the thread-local tables
        finalRecalibrationTables = mergeThreadLocalRecalibrationTables();

        final NestedIntegerArray<RecalDatum> byReadGroupTable = finalRecalibrationTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> byQualTable = finalRecalibrationTables.getQualityScoreTable();

        // iterate over all values in the qual table
        for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[2];
            final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
            final RecalDatum qualDatum = leaf.value;

            if ( rgDatum == null ) {
                // create a copy of qualDatum, and initialize byReadGroup table with it
                byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            } else {
                // combine the qual datum with the existing datum in the byReadGroup table
                rgDatum.combine(qualDatum);
            }
        }

        finalized = true;
    }

    /**
     * Merge all of the thread local recalibration tables into a single one.
     *
     * Reuses one of the recalibration tables to hold the merged table, so this function can only be
     * called once in the engine.
     *
     * @return the merged recalibration table
     */
    @Requires("! finalized")
    private RecalibrationTables mergeThreadLocalRecalibrationTables() {
        if ( recalibrationTablesList.isEmpty() ) {
            recalibrationTablesList.add( new RecalibrationTables(covariates, numReadGroups, maybeLogStream) );
        }

        RecalibrationTables merged = null;
        for ( final RecalibrationTables table : recalibrationTablesList ) {
            if ( merged == null )
                // fast path -- if there's only only one table, so just make it the merged one
                merged = table;
            else {
                merged.combine(table);
            }
        }

        return merged;
    }

    /**
     * Get the final recalibration tables, after finalizeData() has been called
     *
     * This returns the finalized recalibration table collected by this engine.
     *
     * It is an error to call this function before finalizeData has been called
     *
     * @return the finalized recalibration table collected by this engine
     */
    public RecalibrationTables getFinalRecalibrationTables() {
        if ( ! finalized ) throw new IllegalStateException("Cannot get final recalibration tables until finalizeData() has been called");
        return finalRecalibrationTables;
    }
}
