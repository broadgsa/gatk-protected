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

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Class of tests to detect strand bias.
 */
//TODO: will eventually implement ReducibleAnnotation -- see RMSAnnotation.java for an example of an abstract ReducibleAnnotation
public abstract class StrandBiasTest extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {
    private final static Logger logger = Logger.getLogger(StrandBiasTest.class);
    private static boolean stratifiedPerReadAlleleLikelihoodMapWarningLogged = false;
    private static boolean inputVariantContextWarningLogged = false;
    private static boolean getTableFromSamplesWarningLogged = false;
    private static boolean decodeSBBSWarningLogged = false;

    protected static final int ARRAY_DIM = 2;
    protected static final int ARRAY_SIZE = ARRAY_DIM * ARRAY_DIM;

    @Override
    public void initialize(final AnnotatorCompatible walker, final GenomeAnalysisEngine toolkit, final Set<VCFHeaderLine> headerLines) {
        // Does the VCF header contain strand bias (SB) by sample annotation?
        for ( final VCFHeaderLine line : headerLines) {
            if ( line instanceof VCFFormatHeaderLine) {
                final VCFFormatHeaderLine formatline = (VCFFormatHeaderLine)line;
                if ( formatline.getID().equals(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY) ) {
                    logger.warn("StrandBiasBySample annotation exists in input VCF header. Attempting to use StrandBiasBySample " +
                            "values to calculate strand bias annotation values. If no sample has the SB genotype annotation, annotation may still fail.");
                    return;
                }
            }
        }

        // Are there reads from a SAM/BAM file?
        if (toolkit.getReadsDataSource().getReaderIDs().isEmpty())
            logger.warn("No StrandBiasBySample annotation or read data was found. Strand bias annotations will not be output.");
        else
            logger.info("SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values.");
    }

    @Override
    //template method for calculating strand bias annotations using the three different methods
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String,AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {

        // do not process if not a variant
        if ( !vc.isVariant() )
            return null;

        // if the genotype and strand bias are provided, calculate the annotation from the Genotype (GT) field
        if ( vc.hasGenotypes() ) {
            for (final Genotype g : vc.getGenotypes()) {
                if (g.hasAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                    return calculateAnnotationFromGTfield(vc.getGenotypes());
                }
            }
        }

        // if a the variant is a snp and has stratified contexts, calculate the annotation from the stratified contexts
        //stratifiedContexts can come come from VariantAnnotator, but will be empty if no reads were provided
        if (vc.isSNP() && stratifiedContexts != null  && !stratifiedContexts.isEmpty()) {
            return calculateAnnotationFromStratifiedContexts(stratifiedContexts, vc);
        }

        // calculate the annotation from the stratified per read likelihood map
        // stratifiedPerReadAllelelikelihoodMap can come from HaplotypeCaller call to VariantAnnotatorEngine
        else if (stratifiedPerReadAlleleLikelihoodMap != null) {
            return calculateAnnotationFromLikelihoodMap(stratifiedPerReadAlleleLikelihoodMap, vc);
        }
        else {
            // for non-snp variants, we need per-read likelihoods.
            // for snps, we can get same result from simple pileup
            // for indels that do not have a computed strand bias (SB) or strand bias by sample (SBBS)
            return null;
        }
    }

    protected abstract Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes);

    protected abstract Map<String, Object> calculateAnnotationFromStratifiedContexts(final Map<String, AlignmentContext> stratifiedContexts,
                                                                                     final VariantContext vc);

    protected abstract Map<String, Object> calculateAnnotationFromLikelihoodMap(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                                                final VariantContext vc);

    /**
     * Create the contingency table by retrieving the per-sample strand bias annotation and adding them together
     * @param genotypes the genotypes from which to pull out the per-sample strand bias annotation
     * @param minCount minimum threshold for the sample strand bias counts for each ref and alt.
     *                 If both ref and alt counts are above minCount the whole sample strand bias is added to the resulting table
     * @return the table used for several strand bias tests, will be null if none of the genotypes contain the per-sample SB annotation
     */
    protected int[][] getTableFromSamples( final GenotypesContext genotypes, final int minCount ) {
        if( genotypes == null ) {
            if ( !getTableFromSamplesWarningLogged ) {
                logger.warn("Genotypes cannot be null.");
                getTableFromSamplesWarningLogged = true;
            }
            return null;
        }

        final int[] sbArray = {0,0,0,0}; // reference-forward-reverse -by- alternate-forward-reverse
        boolean foundData = false;

        for( final Genotype g : genotypes ) {
            if( g.isNoCall() || ! g.hasAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY) )
                continue;

            foundData = true;
            int[] data;
            if ( g.getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY).getClass().equals(String.class)) {
                final String sbbsString = (String) g.getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
                data = encodeSBBS(sbbsString);
            } else if (g.getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY).getClass().equals(ArrayList.class)) {
                ArrayList sbbsList = (ArrayList) g.getAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
                data = encodeSBBS(sbbsList);
            } else
                throw new IllegalArgumentException("Unexpected " + GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY + " type");

            if ( passesMinimumThreshold(data, minCount) ) {
                for( int index = 0; index < sbArray.length; index++ ) {
                    sbArray[index] += data[index];
                }
            }
        }

        return ( foundData ? decodeSBBS(sbArray) : null );
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    protected static int[][] getSNPContingencyTable(final Map<String, AlignmentContext> stratifiedContexts,
                                                  final Allele ref,
                                                  final List<Allele> allAlts,
                                                  final int minQScoreToConsider,
                                                  final int minCount ) {
        int[][] table = new int[ARRAY_DIM][ARRAY_DIM];

        for (final Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            final int[] myTable = new int[ARRAY_SIZE];
            for (final PileupElement p : sample.getValue().getBasePileup()) {

                if ( ! isUsableBase(p) ) // ignore deletions and bad MQ
                    continue;

                if ( p.getQual() < minQScoreToConsider || p.getMappingQual() < minQScoreToConsider )
                    continue;

                updateTable(myTable, Allele.create(p.getBase(), false), p.getRead(), ref, allAlts);
            }

            if ( passesMinimumThreshold( myTable, minCount ) ) {
                copyToMainTable(myTable, table);
            }

        }

        return table;
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    public static int[][] getContingencyTable( final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                               final VariantContext vc,
                                               final int minCount) {
        if( stratifiedPerReadAlleleLikelihoodMap == null ) {
            if ( !stratifiedPerReadAlleleLikelihoodMapWarningLogged ) {
                logger.warn("stratifiedPerReadAlleleLikelihoodMap cannot be null");
                stratifiedPerReadAlleleLikelihoodMapWarningLogged = true;
            }
            return null;
        }
        if( vc == null ) {
            if ( !inputVariantContextWarningLogged ) {
                logger.warn("input vc cannot be null");
                inputVariantContextWarningLogged = true;
            }
            return null;
        }

        final Allele ref = vc.getReference();
        final Allele alt = vc.getAltAlleleWithHighestAlleleCount();
        final List<Allele> allAlts = vc.getAlternateAlleles();
        final int[][] table = new int[ARRAY_DIM][ARRAY_DIM];

        for (final PerReadAlleleLikelihoodMap maps : stratifiedPerReadAlleleLikelihoodMap.values() ) {
            final int[] myTable = new int[ARRAY_SIZE];
            for (final Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : maps.getLikelihoodReadMap().entrySet()) {
                final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
                final GATKSAMRecord read = el.getKey();
                updateTable(myTable, mostLikelyAllele.getAlleleIfInformative(), read, ref, allAlts);
            }
            if ( passesMinimumThreshold(myTable, minCount) )
                copyToMainTable(myTable, table);
        }

        return table;
    }

    /**
     * Helper method to copy the per-sample table to the main table
     *
     * @param perSampleTable   per-sample table (single dimension)
     * @param mainTable        main table (two dimensions)
     */
    private static void copyToMainTable(final int[] perSampleTable, final int[][] mainTable) {
        mainTable[0][0] += perSampleTable[0];
        mainTable[0][1] += perSampleTable[1];
        mainTable[1][0] += perSampleTable[2];
        mainTable[1][1] += perSampleTable[3];
    }


    /**
     * Can the base in this pileup element be used in comparative tests?
     *
     * @param p the pileup element to consider
     *
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    private static boolean isUsableBase(final PileupElement p) {
        return !( p.isDeletion() ||
                p.getMappingQual() == 0 ||
                p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                ((int) p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE);
    }

    private static void updateTable(final int[] table, final Allele allele, final GATKSAMRecord read, final Allele ref, final List<Allele> allAlts) {

        final boolean matchesRef = allele.equals(ref, true);
        final boolean matchesAnyAlt = allAlts.contains(allele);

        if ( matchesRef || matchesAnyAlt ) {
            final int offset = matchesRef ? 0 : ARRAY_DIM;

            if ( read.isStrandless() ) {
                // a strandless read counts as observations on both strand, at 50% weight, with a minimum of 1
                // (the 1 is to ensure that a strandless read always counts as an observation on both strands, even
                // if the read is only seen once, because it's a merged read or other)
                table[offset]++;
                table[offset + 1]++;
            } else {
                // a normal read with an actual strand
                final boolean isFW = !read.getReadNegativeStrandFlag();
                table[offset + (isFW ? 0 : 1)]++;
            }
        }
    }

    /**
     * Does this strand data array pass the minimum threshold for inclusion?
     *
     * @param data  the array
     * @param minCount The minimum threshold of counts in the array
     * @return true if it passes the minimum threshold, false otherwise
     */
    protected static boolean passesMinimumThreshold(final int[] data, final int minCount) {
        // the ref and alt totals must be greater than MIN_COUNT
        return data[0] + data[1] + data[2] + data[3] > minCount;
    }

    /**
     * Helper function to parse the genotype annotation into the SB annotation array
     * @param string the string that is returned by genotype.getAnnotation("SB")
     * @return the array used by the per-sample Strand Bias annotation
     */
    private static int[] encodeSBBS( final String string ) {
        final int[] array = new int[ARRAY_SIZE];
        final StringTokenizer tokenizer = new StringTokenizer(string, ",", false);
        for( int index = 0; index < ARRAY_SIZE; index++ ) {
            array[index] = Integer.parseInt(tokenizer.nextToken());
        }
        return array;
    }

    /**
     * Helper function to parse the genotype annotation into the SB annotation array
     * @param arrayList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    private static int[] encodeSBBS( final ArrayList<Integer> arrayList ) {
        final int[] array = new int[ARRAY_SIZE];
        int index = 0;
        for ( Integer item : arrayList )
            array[index++] = item.intValue();

        return array;
    }

    /**
     * Helper function to turn the  SB annotation array into a contingency table
     * @param array the array used by the per-sample Strand Bias annotation
     * @return the table used by the StrandOddsRatio annotation
     */
    private static int[][] decodeSBBS( final int[] array ) {
        if(array.length != ARRAY_SIZE) {
            if ( !decodeSBBSWarningLogged ) {
                logger.warn("Expecting a length = " +  ARRAY_SIZE + " strand bias array.");
                decodeSBBSWarningLogged = true;
            }
            return null;
        }
        final int[][] table = new int[ARRAY_DIM][ARRAY_DIM];
        table[0][0] = array[0];
        table[0][1] = array[1];
        table[1][0] = array[2];
        table[1][1] = array[3];
        return table;
    }
}
