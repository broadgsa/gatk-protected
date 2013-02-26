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

package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.bqsr.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;

import java.io.*;
import java.util.*;

/**
 * This class has all the static functionality for reading a recalibration report file into memory. 
 *
 * @author carneiro
 * @since 3/26/12
 */
public class RecalibrationReport {
    private QuantizationInfo quantizationInfo; // histogram containing the counts for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables; // quick access reference to the tables
    private final Covariate[] requestedCovariates; // list of all covariates to be used in this calculation
    private final HashMap<String, Integer> optionalCovariateIndexes;

    private final GATKReportTable argumentTable; // keep the argument table untouched just for output purposes
    private final RecalibrationArgumentCollection RAC; // necessary for quantizing qualities with the same parameter

    private final int[] tempRGarray = new int[2];
    private final int[] tempQUALarray = new int[3];
    private final int[] tempCOVarray = new int[4];

    public RecalibrationReport(final File RECAL_FILE) {
        final GATKReport report = new GATKReport(RECAL_FILE);

        argumentTable = report.getTable(RecalUtils.ARGUMENT_REPORT_TABLE_TITLE);
        RAC = initializeArgumentCollectionTable(argumentTable);

        GATKReportTable quantizedTable = report.getTable(RecalUtils.QUANTIZED_REPORT_TABLE_TITLE);
        quantizationInfo = initializeQuantizationTable(quantizedTable);

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalUtils.initializeCovariates(RAC); // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getFirst();
        ArrayList<Covariate> optionalCovariates = covariates.getSecond();
        requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        optionalCovariateIndexes = new HashMap<String, Integer>(optionalCovariates.size());
        int covariateIndex = 0;
        for (final Covariate covariate : requiredCovariates)
            requestedCovariates[covariateIndex++] = covariate;
        for (final Covariate covariate : optionalCovariates) {
            requestedCovariates[covariateIndex] = covariate;
            final String covariateName = covariate.getClass().getSimpleName().split("Covariate")[0]; // get the name of the covariate (without the "covariate" part of it) so we can match with the GATKReport
            optionalCovariateIndexes.put(covariateName, covariateIndex-2);
            covariateIndex++;
        }

        for (Covariate cov : requestedCovariates)
            cov.initialize(RAC); // initialize any covariate member variables using the shared argument collection

        recalibrationTables = new RecalibrationTables(requestedCovariates, countReadGroups(report.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE)));

        parseReadGroupTable(report.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE), recalibrationTables.getReadGroupTable());

        parseQualityScoreTable(report.getTable(RecalUtils.QUALITY_SCORE_REPORT_TABLE_TITLE), recalibrationTables.getQualityScoreTable());

        parseAllCovariatesTable(report.getTable(RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE), recalibrationTables);

    }

    /**
     * Counts the number of unique read groups in the table
     *
     * @param reportTable the GATKReport table containing data for this table
     * @return the number of unique read groups
     */
    private int countReadGroups(final GATKReportTable reportTable) {
        Set<String> readGroups = new HashSet<String>();
        for ( int i = 0; i < reportTable.getNumRows(); i++ )
            readGroups.add(reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME).toString());
        return readGroups.size();
    }

    /**
    * Combines two recalibration reports by adding all observations and errors
    *
    * Note: This method DOES NOT recalculate the empirical qualities and quantized qualities. You have to recalculate
    * them after combining. The reason for not calculating it is because this function is intended for combining a
    * series of recalibration reports, and it only makes sense to calculate the empirical qualities and quantized
    * qualities after all the recalibration reports have been combined. Having the user recalculate when appropriate,
    * makes this method faster
    *
    * Note2: The empirical quality reported, however, is recalculated given its simplicity.
    *
    * @param other the recalibration report to combine with this one
    */
    public void combine(final RecalibrationReport other) {
        for ( int tableIndex = 0; tableIndex < recalibrationTables.numTables(); tableIndex++ ) {
            final NestedIntegerArray<RecalDatum> myTable = recalibrationTables.getTable(tableIndex);
            final NestedIntegerArray<RecalDatum> otherTable = other.recalibrationTables.getTable(tableIndex);
            RecalUtils.combineTables(myTable, otherTable);
        }
    }

    public QuantizationInfo getQuantizationInfo() {
        return quantizationInfo;
    }

    public RecalibrationTables getRecalibrationTables() {
        return recalibrationTables;
    }

    public Covariate[] getRequestedCovariates() {
        return requestedCovariates;
    }

    /**
     * Compiles the list of keys for the Covariates table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param recalibrationTables    the recalibration tables
\     */
    private void parseAllCovariatesTable(final GATKReportTable reportTable, final RecalibrationTables recalibrationTables) {
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempCOVarray[0] = requestedCovariates[0].keyFromValue(rg);
            final Object qual = reportTable.get(i, RecalUtils.QUALITY_SCORE_COLUMN_NAME);
            tempCOVarray[1] = requestedCovariates[1].keyFromValue(qual);

            final String covName = (String)reportTable.get(i, RecalUtils.COVARIATE_NAME_COLUMN_NAME);
            final int covIndex = optionalCovariateIndexes.get(covName);
            final Object covValue = reportTable.get(i, RecalUtils.COVARIATE_VALUE_COLUMN_NAME);
            tempCOVarray[2] = requestedCovariates[RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal() + covIndex].keyFromValue(covValue);

            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempCOVarray[3] = event.ordinal();

            recalibrationTables.getTable(RecalibrationTables.TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal() + covIndex).put(getRecalDatum(reportTable, i, false), tempCOVarray);
        }
    }

    /**
     *
     * Compiles the list of keys for the QualityScore table and uses the shared parsing utility to produce the actual table
     * @param reportTable            the GATKReport table containing data for this table
     * @param qualTable               the map representing this table
     */
    private void parseQualityScoreTable(final GATKReportTable reportTable, final NestedIntegerArray<RecalDatum> qualTable) {
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempQUALarray[0] = requestedCovariates[0].keyFromValue(rg);
            final Object qual = reportTable.get(i, RecalUtils.QUALITY_SCORE_COLUMN_NAME);
            tempQUALarray[1] = requestedCovariates[1].keyFromValue(qual);
            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempQUALarray[2] = event.ordinal();

            qualTable.put(getRecalDatum(reportTable, i, false), tempQUALarray);
        }
    }

    /**
     * Compiles the list of keys for the ReadGroup table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param rgTable                the map representing this table
     */
    private void parseReadGroupTable(final GATKReportTable reportTable, final NestedIntegerArray<RecalDatum> rgTable) {
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            final Object rg = reportTable.get(i, RecalUtils.READGROUP_COLUMN_NAME);
            tempRGarray[0] = requestedCovariates[0].keyFromValue(rg);
            final EventType event = EventType.eventFrom((String)reportTable.get(i, RecalUtils.EVENT_TYPE_COLUMN_NAME));
            tempRGarray[1] = event.ordinal();

            rgTable.put(getRecalDatum(reportTable, i, true), tempRGarray);
        }
    }

    private double asDouble(final Object o) {
        if ( o instanceof Double )
            return (Double)o;
        else if ( o instanceof Integer )
            return (Integer)o;
        else if ( o instanceof Long )
            return (Long)o;
        else
            throw new ReviewedStingException("Object " + o + " is expected to be either a double, long or integer but it's not either: " + o.getClass());
    }

    private long asLong(final Object o) {
        if ( o instanceof Long )
            return (Long)o;
        else if ( o instanceof Integer )
            return ((Integer)o).longValue();
        else if ( o instanceof Double )
            return ((Double)o).longValue();
        else
            throw new ReviewedStingException("Object " + o + " is expected to be a long but it's not: " + o.getClass());
    }

    private RecalDatum getRecalDatum(final GATKReportTable reportTable, final int row, final boolean hasEstimatedQReportedColumn) {
        final long nObservations = asLong(reportTable.get(row, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME));
        final double nErrors = asDouble(reportTable.get(row, RecalUtils.NUMBER_ERRORS_COLUMN_NAME));
        //final double empiricalQuality = asDouble(reportTable.get(row, RecalUtils.EMPIRICAL_QUALITY_COLUMN_NAME));

        // the estimatedQreported column only exists in the ReadGroup table
        final double estimatedQReported = hasEstimatedQReportedColumn ?
                (Double) reportTable.get(row, RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME) : // we get it if we are in the read group table
                Byte.parseByte((String) reportTable.get(row, RecalUtils.QUALITY_SCORE_COLUMN_NAME)); // or we use the reported quality if we are in any other table

        final RecalDatum datum = new RecalDatum(nObservations, nErrors, (byte)1);
        datum.setEstimatedQReported(estimatedQReported);
        //datum.setEmpiricalQuality(empiricalQuality); // don't set the value here because we will want to recompute with a different conditional Q score prior value
        return datum;
    }

    /**
     * Parses the quantization table from the GATK Report and turns it into a map of original => quantized quality scores
     *
     * @param table the GATKReportTable containing the quantization mappings
     * @return an ArrayList with the quantization mappings from 0 to MAX_SAM_QUAL_SCORE
     */
    private QuantizationInfo initializeQuantizationTable(GATKReportTable table) {
        final Byte[] quals  = new Byte[QualityUtils.MAX_SAM_QUAL_SCORE + 1];
        final Long[] counts = new Long[QualityUtils.MAX_SAM_QUAL_SCORE + 1];
        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final byte originalQual = (byte)i;
            final Object quantizedObject = table.get(i, RecalUtils.QUANTIZED_VALUE_COLUMN_NAME);
            final Object countObject = table.get(i, RecalUtils.QUANTIZED_COUNT_COLUMN_NAME);
            final byte quantizedQual = Byte.parseByte(quantizedObject.toString());
            final long quantizedCount = Long.parseLong(countObject.toString());
            quals[originalQual] = quantizedQual;
            counts[originalQual] = quantizedCount;
        }
        return new QuantizationInfo(Arrays.asList(quals), Arrays.asList(counts));
    }

    /**
     * Parses the arguments table from the GATK Report and creates a RAC object with the proper initialization values
     *
     * @param table the GATKReportTable containing the arguments and its corresponding values
     * @return a RAC object properly initialized with all the objects in the table
     */
    private RecalibrationArgumentCollection initializeArgumentCollectionTable(GATKReportTable table) {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        for ( int i = 0; i < table.getNumRows(); i++ ) {
            final String argument = table.get(i, "Argument").toString();
            Object value = table.get(i, RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
            if (value.equals("null"))
                value = null; // generic translation of null values that were printed out as strings | todo -- add this capability to the GATKReport

            if (argument.equals("covariate") && value != null)
                RAC.COVARIATES = value.toString().split(",");

            else if (argument.equals("standard_covs"))
                RAC.DO_NOT_USE_STANDARD_COVARIATES = Boolean.parseBoolean((String) value);

            else if (argument.equals("solid_recal_mode"))
                RAC.SOLID_RECAL_MODE = RecalUtils.SOLID_RECAL_MODE.recalModeFromString((String) value);

            else if (argument.equals("solid_nocall_strategy"))
                RAC.SOLID_NOCALL_STRATEGY = RecalUtils.SOLID_NOCALL_STRATEGY.nocallStrategyFromString((String) value);

            else if (argument.equals("mismatches_context_size"))
                RAC.MISMATCHES_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("indels_context_size"))
                RAC.INDELS_CONTEXT_SIZE = Integer.parseInt((String) value);

            else if (argument.equals("mismatches_default_quality"))
                RAC.MISMATCHES_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("insertions_default_quality"))
                RAC.INSERTIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("deletions_default_quality"))
                RAC.DELETIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);

            else if (argument.equals("maximum_cycle_value"))
                RAC.MAXIMUM_CYCLE_VALUE = Integer.parseInt((String) value);

            else if (argument.equals("low_quality_tail"))
                RAC.LOW_QUAL_TAIL = Byte.parseByte((String) value);

            else if (argument.equals("default_platform"))
                RAC.DEFAULT_PLATFORM = (String) value;

            else if (argument.equals("force_platform"))
                RAC.FORCE_PLATFORM = (String) value;

            else if (argument.equals("quantizing_levels"))
                RAC.QUANTIZING_LEVELS = Integer.parseInt((String) value);

            else if (argument.equals("recalibration_report"))
                RAC.existingRecalibrationReport = (value == null) ? null : new File((String) value);

            else if (argument.equals("plot_pdf_file"))
                RAC.RECAL_PDF_FILE = (value == null) ? null : new File((String) value);

            else if (argument.equals("binary_tag_name"))
                RAC.BINARY_TAG_NAME = (value == null) ? null : (String) value;

            else if (argument.equals("sort_by_all_columns"))
                RAC.SORT_BY_ALL_COLUMNS = Boolean.parseBoolean((String) value);
        }

        return RAC;
    }

    /**
     * this functionality avoids recalculating the empirical qualities, estimated reported quality
     * and quantization of the quality scores during every call of combine(). Very useful for the BQSRGatherer.
     */
    public void calculateQuantizedQualities() {
        quantizationInfo = new QuantizationInfo(recalibrationTables, RAC.QUANTIZING_LEVELS);
    }

    public void output(PrintStream output) {
        RecalUtils.outputRecalibrationReport(argumentTable, quantizationInfo, recalibrationTables, requestedCovariates, output, RAC.SORT_BY_ALL_COLUMNS);
    }

    public RecalibrationArgumentCollection getRAC() {
        return RAC;
    }

    public Covariate[] getCovariates() {
        return requestedCovariates;
    }
}
