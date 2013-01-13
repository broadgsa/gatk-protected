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

import net.sf.samtools.SAMTag;
import net.sf.samtools.SAMUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.NestedIntegerArray;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;

/**
 * Utility methods to facilitate on-the-fly base quality score recalibration.
 *
 * User: carneiro and rpoplin
 * Date: 2/4/12
 */

public class BaseRecalibration {
    private static Logger logger = Logger.getLogger(BaseRecalibration.class);
    private final static boolean TEST_CACHING = false;

    private final QuantizationInfo quantizationInfo; // histogram containing the map for qual quantization (calculated after recalibration is done)
    private final RecalibrationTables recalibrationTables;
    private final Covariate[] requestedCovariates; // list of all covariates to be used in this calculation

    private final boolean disableIndelQuals;
    private final int preserveQLessThan;
    private final boolean emitOriginalQuals;

    private final NestedIntegerArray<Double> globalDeltaQs;
    private final NestedIntegerArray<Double> deltaQReporteds;


    /**
     * Constructor using a GATK Report file
     * 
     * @param RECAL_FILE         a GATK Report file containing the recalibration information
     * @param quantizationLevels number of bins to quantize the quality scores
     * @param disableIndelQuals  if true, do not emit base indel qualities
     * @param preserveQLessThan  preserve quality scores less than this value
     */
    public BaseRecalibration(final File RECAL_FILE, final int quantizationLevels, final boolean disableIndelQuals, final int preserveQLessThan, final boolean emitOriginalQuals) {
        RecalibrationReport recalibrationReport = new RecalibrationReport(RECAL_FILE);

        recalibrationTables = recalibrationReport.getRecalibrationTables();
        requestedCovariates = recalibrationReport.getRequestedCovariates();
        quantizationInfo = recalibrationReport.getQuantizationInfo();
        if (quantizationLevels == 0) // quantizationLevels == 0 means no quantization, preserve the quality scores
            quantizationInfo.noQuantization();
        else if (quantizationLevels > 0 && quantizationLevels != quantizationInfo.getQuantizationLevels()) // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
            quantizationInfo.quantizeQualityScores(quantizationLevels);

        this.disableIndelQuals = disableIndelQuals;
        this.preserveQLessThan = preserveQLessThan;
        this.emitOriginalQuals = emitOriginalQuals;

        logger.info("Calculating cached tables...");

        //
        // Create a NestedIntegerArray<Double> that maps from rgKey x errorModel -> double,
        // where the double is the result of this calculation.  The entire calculation can
        // be done upfront, on initialization of this BaseRecalibration structure
        //
        final NestedIntegerArray<RecalDatum> byReadGroupTable = recalibrationTables.getReadGroupTable();
        globalDeltaQs = new NestedIntegerArray<Double>( byReadGroupTable.getDimensions() );
        logger.info("Calculating global delta Q table...");
        for ( NestedIntegerArray.Leaf<RecalDatum> leaf : byReadGroupTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[1];
            final double globalDeltaQ = calculateGlobalDeltaQ(rgKey, EventType.eventFrom(eventIndex));
            globalDeltaQs.put(globalDeltaQ, rgKey, eventIndex);
        }


        // The calculation of the deltaQ report is constant.  key[0] and key[1] are the read group and qual, respectively
        // and globalDeltaQ is a constant for the read group.  So technically the delta Q reported is simply a lookup
        // into a matrix indexed by rgGroup, qual, and event type.
        // the code below actually creates this cache with a NestedIntegerArray calling into the actual
        // calculateDeltaQReported code.
        final NestedIntegerArray<RecalDatum> byQualTable = recalibrationTables.getQualityScoreTable();
        deltaQReporteds = new NestedIntegerArray<Double>( byQualTable.getDimensions() );
        logger.info("Calculating delta Q reported table...");
        for ( NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int qual = leaf.keys[1];
            final int eventIndex = leaf.keys[2];
            final EventType event = EventType.eventFrom(eventIndex);
            final double globalDeltaQ = getGlobalDeltaQ(rgKey, event);
            final double deltaQReported = calculateDeltaQReported(rgKey, qual, event, globalDeltaQ, (byte)qual);
            deltaQReporteds.put(deltaQReported, rgKey, qual, eventIndex);
        }

        logger.info("done calculating cache");
    }

    /**
     * Recalibrates the base qualities of a read
     *
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     *
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     *
     * Given the full recalibration table, we perform the following preprocessing steps:
     *
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     *
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param read the read to recalibrate
     */
    public void recalibrateRead(final GATKSAMRecord read) {
        if (emitOriginalQuals && read.getAttribute(SAMTag.OQ.name()) == null) { // Save the old qualities if the tag isn't already taken in the read
            try {
                read.setAttribute(SAMTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
            } catch (IllegalArgumentException e) {
                throw new UserException.MalformedBAM(read, "illegal base quality encountered; " + e.getMessage());
            }
        }

        final ReadCovariates readCovariates = RecalUtils.computeCovariates(read, requestedCovariates);
        final int readLength = read.getReadLength();

        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
            if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
                read.setBaseQualities(null, errorModel);
                continue;
            }

            final byte[] quals = read.getBaseQualities(errorModel);

            // get the keyset for this base using the error model
            final int[][] fullReadKeySet = readCovariates.getKeySet(errorModel);

            // the rg key is constant over the whole read, the global deltaQ is too
            final int rgKey = fullReadKeySet[0][0];

            final double globalDeltaQ = getGlobalDeltaQ(rgKey, errorModel);

            for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read
                final byte origQual = quals[offset];

                // only recalibrate usable qualities (the original quality will come from the instrument -- reported quality)
                if ( origQual >= preserveQLessThan ) {
                    // get the keyset for this base using the error model
                    final int[] keySet = fullReadKeySet[offset];
                    final double deltaQReported = getDeltaQReported(keySet[0], keySet[1], errorModel, globalDeltaQ);
                    final double deltaQCovariates = calculateDeltaQCovariates(recalibrationTables, keySet, errorModel, globalDeltaQ, deltaQReported, origQual);

                    // calculate the recalibrated qual using the BQSR formula
                    double recalibratedQualDouble = origQual + globalDeltaQ + deltaQReported + deltaQCovariates;

                    // recalibrated quality is bound between 1 and MAX_QUAL
                    final byte recalibratedQual = QualityUtils.boundQual(MathUtils.fastRound(recalibratedQualDouble), QualityUtils.MAX_RECALIBRATED_Q_SCORE);

                    // return the quantized version of the recalibrated quality
                    final byte recalibratedQualityScore = quantizationInfo.getQuantizedQuals().get(recalibratedQual);

                    quals[offset] = recalibratedQualityScore;
                }
            }

            // finally update the base qualities in the read
            read.setBaseQualities(quals, errorModel);
        }
    }

    private double getGlobalDeltaQ(final int rgKey, final EventType errorModel) {
        final Double cached = globalDeltaQs.get(rgKey, errorModel.ordinal());

        if ( TEST_CACHING ) {
            final double calcd = calculateGlobalDeltaQ(rgKey, errorModel);
            if ( calcd != cached )
                throw new IllegalStateException("calculated " + calcd + " and cached " + cached + " global delta q not equal at " + rgKey + " / " + errorModel);
        }

        return cachedWithDefault(cached);
    }

    private double getDeltaQReported(final int rgKey, final int qualKey, final EventType errorModel, final double globalDeltaQ) {
        final Double cached = deltaQReporteds.get(rgKey, qualKey, errorModel.ordinal());

        if ( TEST_CACHING ) {
            final double calcd = calculateDeltaQReported(rgKey, qualKey, errorModel, globalDeltaQ, (byte)qualKey);
            if ( calcd != cached )
                throw new IllegalStateException("calculated " + calcd + " and cached " + cached + " global delta q not equal at " + rgKey + " / " + qualKey + " / " + errorModel);
        }

        return cachedWithDefault(cached);
    }

    /**
     * @param d a Double (that may be null) that is the result of a delta Q calculation
     * @return a double == d if d != null, or 0.0 if it is
     */
    private double cachedWithDefault(final Double d) {
        return d == null ? 0.0 : d;
    }

    /**
     * Note that this calculation is a constant for each rgKey and errorModel.  We need only
     * compute this value once for all data.
     *
     * @param rgKey        read group key
     * @param errorModel   the event type
     * @return global delta Q
     */
    private double calculateGlobalDeltaQ(final int rgKey, final EventType errorModel) {
        double result = 0.0;

        final RecalDatum empiricalQualRG = recalibrationTables.getReadGroupTable().get(rgKey, errorModel.ordinal());

        if (empiricalQualRG != null) {
            final double globalDeltaQEmpirical = empiricalQualRG.getEmpiricalQuality();
            final double aggregrateQReported = empiricalQualRG.getEstimatedQReported();
            result = globalDeltaQEmpirical - aggregrateQReported;
        }

        return result;
    }

    private double calculateDeltaQReported(final int rgKey, final int qualKey, final EventType errorModel, final double globalDeltaQ, final byte qualFromRead) {
        double result = 0.0;

        final RecalDatum empiricalQualQS = recalibrationTables.getQualityScoreTable().get(rgKey, qualKey, errorModel.ordinal());
        if (empiricalQualQS != null) {
            final double deltaQReportedEmpirical = empiricalQualQS.getEmpiricalQuality();
            result = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
        }

        return result;
    }

    private double calculateDeltaQCovariates(final RecalibrationTables recalibrationTables, final int[] key, final EventType errorModel, final double globalDeltaQ, final double deltaQReported, final byte qualFromRead) {
        double result = 0.0;

        // for all optional covariates
        for (int i = 2; i < requestedCovariates.length; i++) {
            if (key[i] < 0)
                continue;

            result += calculateDeltaQCovariate(recalibrationTables.getTable(i),
                    key[0], key[1], key[i], errorModel,
                    globalDeltaQ, deltaQReported, qualFromRead);
        }

        return result;
    }

    private double calculateDeltaQCovariate(final NestedIntegerArray<RecalDatum> table,
                                            final int rgKey,
                                            final int qualKey,
                                            final int tableKey,
                                            final EventType errorModel,
                                            final double globalDeltaQ,
                                            final double deltaQReported,
                                            final byte qualFromRead) {
        final RecalDatum empiricalQualCO = table.get(rgKey, qualKey, tableKey, errorModel.ordinal());
        if (empiricalQualCO != null) {
            final double deltaQCovariateEmpirical = empiricalQualCO.getEmpiricalQuality();
            return deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported);
        } else {
            return 0.0;
        }
    }
}
