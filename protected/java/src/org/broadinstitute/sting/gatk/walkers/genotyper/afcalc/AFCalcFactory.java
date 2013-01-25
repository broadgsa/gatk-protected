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

package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.lang.reflect.Constructor;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Factory to make AFCalculations
 */
public class AFCalcFactory {
    /**
     * Enumeration of usable AF calculation, their constraints (i.e. ploidy).
     *
     * Note that the order these occur in the enum is the order of preference, so
     * the first value is taken over the second when multiple calculations satisfy
     * the needs of the request (i.e., considering ploidy).
     */
    public enum Calculation {
        /** default implementation */
        EXACT_INDEPENDENT(IndependentAllelesDiploidExactAFCalc.class, 2, -1),

        /** reference implementation of multi-allelic EXACT model.  Extremely slow for many alternate alleles */
        EXACT_REFERENCE(ReferenceDiploidExactAFCalc.class, 2, -1),

        /** original biallelic exact model, for testing only */
        EXACT_ORIGINAL(OriginalDiploidExactAFCalc.class, 2, 2),

        /** implementation that supports any sample ploidy */
        EXACT_GENERAL_PLOIDY("GeneralPloidyExactAFCalc", -1, -1);

        /**
         * Must be a name because we look this up dynamically
         */
        public final String className;
        public final int maxAltAlleles;
        public final int requiredPloidy;

        private Calculation(final String className, final int requiredPloidy, final int maxAltAlleles) {
            this.className = className;
            this.requiredPloidy = requiredPloidy;
            this.maxAltAlleles = maxAltAlleles;
        }

        private Calculation(final Class clazz, final int requiredPloidy, final int maxAltAlleles) {
            this(clazz.getSimpleName(), requiredPloidy, maxAltAlleles);
        }

        public boolean usableForParams(final int requestedPloidy, final int requestedMaxAltAlleles) {
            return (requiredPloidy == -1 || requiredPloidy == requestedPloidy)
                    && (maxAltAlleles == -1 || maxAltAlleles >= requestedMaxAltAlleles);
        }

        public static Calculation getDefaultModel() { return EXACT_INDEPENDENT; }
    }

    private static final Map<String, Class<? extends AFCalc>> afClasses;
    static {
        afClasses = new PluginManager<AFCalc>(AFCalc.class).getPluginsByName();
    }

    private AFCalcFactory() {

    }

    private static Class<? extends AFCalc> getClassByName(final String name) {
        for ( final Class<? extends AFCalc> clazz : afClasses.values() ) {
            if ( clazz.getSimpleName().contains(name) ) {
                return clazz;
            }
        }

        return null;
    }

    /**
     * Create a new AFCalc based on the parameters in the UAC
     *
     * @param UAC the UnifiedArgumentCollection containing the command-line parameters for the caller
     * @param nSamples the number of samples we will be using
     * @param logger an optional (can be null) logger to override the default in the model
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final UnifiedArgumentCollection UAC,
                                      final int nSamples,
                                      final Logger logger) {
        final int maxAltAlleles = UAC.MAX_ALTERNATE_ALLELES;
        if ( ! UAC.AFmodel.usableForParams(UAC.samplePloidy, maxAltAlleles) ) {
            logger.info("Requested ploidy " + UAC.samplePloidy + " maxAltAlleles " + maxAltAlleles + " not supported by requested model " + UAC.AFmodel + " looking for an option");
            final List<Calculation> supportingCalculations = new LinkedList<Calculation>();
            for ( final Calculation calc : Calculation.values() ) {
                if ( calc.usableForParams(UAC.samplePloidy, maxAltAlleles) )
                    supportingCalculations.add(calc);
            }

            if ( supportingCalculations.isEmpty() )
                throw new UserException("no AFCalculation model found that supports ploidy of " + UAC.samplePloidy + " and max alt alleles " + maxAltAlleles);
            else if ( supportingCalculations.size() > 1 )
                logger.debug("Warning, multiple supporting AFCalcs found " + Utils.join(",", supportingCalculations) + " choosing first arbitrarily");
            else
                UAC.AFmodel = supportingCalculations.get(0);
            logger.info("Selecting model " + UAC.AFmodel);
        }

        final AFCalc calc = createAFCalc(UAC.AFmodel, nSamples, UAC.MAX_ALTERNATE_ALLELES, UAC.samplePloidy);

        if ( logger != null ) calc.setLogger(logger);
        if ( UAC.exactCallsLog != null ) calc.enableProcessLog(UAC.exactCallsLog);

        return calc;
    }

    /**
     * Create a new AFCalc, choosing the best implementation based on the given parameters, assuming
     * that we will only be requesting bi-allelic variants to diploid genotypes
     *
     * @param nSamples the number of samples we'll be using
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final int nSamples) {
        return createAFCalc(chooseBestCalculation(nSamples, 2, 1), nSamples, 2, 2);
    }

    /**
     * Create a new AFCalc that supports maxAltAlleles for all variants and diploid genotypes
     *
     * @param calc the calculation we'd like to use
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles for both SNPs and indels
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final Calculation calc, final int nSamples, final int maxAltAlleles) {
        return createAFCalc(calc, nSamples, maxAltAlleles, 2);
    }

    /**
     * Create a new AFCalc, choosing the best implementation based on the given parameters
     *
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles to consider for SNPs
     * @param ploidy the sample ploidy.  Must be consistent with the calc
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final int nSamples, final int maxAltAlleles, final int ploidy) {
        return createAFCalc(chooseBestCalculation(nSamples, ploidy, maxAltAlleles), nSamples, maxAltAlleles, ploidy);
    }

    /**
     * Choose the best calculation for nSamples and ploidy
     *
     * @param nSamples
     * @param ploidy
     * @param maxAltAlleles
     * @return
     */
    private static Calculation chooseBestCalculation(final int nSamples, final int ploidy, final int maxAltAlleles) {
        for ( final Calculation calc : Calculation.values() ) {
            if ( calc.usableForParams(ploidy, maxAltAlleles) ) {
                return calc;
            }
        }

        throw new IllegalStateException("no calculation found that supports nSamples " + nSamples + " ploidy " + ploidy + " and maxAltAlleles " + maxAltAlleles);
    }

    /**
     * Create a new AFCalc
     *
     * @param calc the calculation to use
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles to consider for SNPs
     * @param ploidy the sample ploidy.  Must be consistent with the calc
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final Calculation calc, final int nSamples, final int maxAltAlleles, final int ploidy) {
        if ( calc == null ) throw new IllegalArgumentException("Calculation cannot be null");
        if ( nSamples < 0 ) throw new IllegalArgumentException("nSamples must be greater than zero " + nSamples);
        if ( maxAltAlleles < 1 ) throw new IllegalArgumentException("maxAltAlleles must be greater than zero " + maxAltAlleles);
        if ( ploidy < 1 ) throw new IllegalArgumentException("sample ploidy must be greater than zero " + ploidy);

        if ( ! calc.usableForParams(ploidy, maxAltAlleles) )
            throw new IllegalArgumentException("AFCalc " + calc + " does not support requested ploidy " + ploidy);

        final Class<? extends AFCalc> afClass = getClassByName(calc.className);
        if ( afClass == null )
            throw new IllegalArgumentException("Unexpected AFCalc " + calc);

        try {
            Object args[] = new Object[]{nSamples, maxAltAlleles, ploidy};
            Constructor c = afClass.getDeclaredConstructor(int.class, int.class, int.class);
            return (AFCalc)c.newInstance(args);
        } catch (Exception e) {
            throw new ReviewedStingException("Could not instantiate AFCalc " + calc, e);
        }
    }

    protected static List<AFCalc> createAFCalcs(final List<Calculation> calcs, final int nSamples, final int maxAltAlleles, final int ploidy) {
        final List<AFCalc> AFCalcs = new LinkedList<AFCalc>();

        for ( final Calculation calc : calcs )
            AFCalcs.add(createAFCalc(calc, nSamples, maxAltAlleles, ploidy));

        return AFCalcs;
    }
}
