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
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ReducibleAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific implementation of strand bias annotations
 */
public abstract class AS_StrandBiasTest extends StrandBiasTest implements ReducibleAnnotation {
    private final static Logger logger = Logger.getLogger(StrandBiasTest.class);
    protected final String splitDelim = "\\|"; //String.split takes a regex, so we need to escape the pipe
    protected final String printDelim = "|";
    protected final String reducedDelim = ",";
    protected AnnotatorCompatible callingWalker;
    protected final int MIN_COUNT = 2;
    protected static final double MIN_PVALUE = 1E-320;
    protected final int FORWARD = 0;
    protected final int REVERSE = 1;
    protected final ArrayList<Integer> ZERO_LIST = new ArrayList<>();

    @Override
    public void initialize(final AnnotatorCompatible walker, final GenomeAnalysisEngine toolkit, final Set<VCFHeaderLine> headerLines) {
        if (!AnnotationUtils.walkerSupportsAlleleSpecificAnnotations(walker))
            logger.warn("Allele-specific annotations can only be used with HaplotypeCaller, CombineGVCFs and GenotypeGVCFs -- no data will be output");
        callingWalker = walker;
        ZERO_LIST.add(0,0);
        ZERO_LIST.add(1,0);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        if (AnnotationUtils.walkerRequiresRawData(callingWalker))
            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
        else
            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_SB_TABLE_KEY; }

    public Map<String, Object> annotateRawData(final RefMetaDataTracker tracker,
                                               final AnnotatorCompatible walker,
                                               final ReferenceContext ref,
                                               final Map<String, AlignmentContext> stratifiedContexts,
                                               final VariantContext vc,
                                               final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {

        //for allele-specific annotations we only call from HC and we only use perReadAlleleLikelihoodMap
        if ( perReadAlleleLikelihoodMap == null)
            return new HashMap<>();

        // calculate the annotation from the stratified per read likelihood map
        // stratifiedPerReadAllelelikelihoodMap can come from HaplotypeCaller call to VariantAnnotatorEngine
        else if (perReadAlleleLikelihoodMap != null) {
            final HashMap<String, Object> annotations = new HashMap<>();
            final ReducibleAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            calculateRawData(vc, perReadAlleleLikelihoodMap, myData);
            final Map<Allele, List<Integer>> perAlleleValues = myData.getAttributeMap();
            final String annotationString = makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
            annotations.put(getRawKeyName(), annotationString);
            return annotations;
        }
        else {
            // for non-snp variants, we need per-read likelihoods.
            // for snps, we can get same result from simple pileup
            // for indels that do not have a computed strand bias (SB) or strand bias by sample (SBBS)
            return null;
        }
    }

    protected void parseRawDataString(ReducibleAnnotationData<List<Integer>> myData) {
        final String rawDataString = myData.getRawData();
        String[] rawDataPerAllele;
        String[] rawListEntriesAsStringVector;
        Map<Allele, List<Integer>> perAlleleValues = new HashMap<>();
        //Initialize maps
        for (Allele current : myData.getAlleles()) {
            perAlleleValues.put(current, new LinkedList<Integer>());
        }
        //rawDataPerAllele is the list of values for each allele (each of variable length)
        rawDataPerAllele = rawDataString.split(splitDelim);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            String alleleData = rawDataPerAllele[i];
            if (alleleData.isEmpty())
                continue;
            List<Integer> alleleList = perAlleleValues.get(myData.getAlleles().get(i));
            rawListEntriesAsStringVector = alleleData.split(",");
            //Read counts will only ever be integers
            for (String s : rawListEntriesAsStringVector) {
                if (!s.isEmpty())
                    alleleList.add(Integer.parseInt(s.trim()));
            }
        }
        myData.setAttributeMap(perAlleleValues);
    }

    @Override
    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<? extends ReducibleAnnotationData> annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        ReducibleAnnotationData combinedData = new AlleleSpecificAnnotationData(vcAlleles, null);

        for (final ReducibleAnnotationData currentValue : annotationList) {
           parseRawDataString(currentValue);
           combineAttributeMap(currentValue, combinedData);
        }
        final Map<String, Object> annotations = new HashMap<>();
        final String annotationString = makeRawAnnotationString(vcAlleles, combinedData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    protected void combineAttributeMap(final ReducibleAnnotationData<List<Integer>> toAdd, final ReducibleAnnotationData<List<Integer>> combined) {
        for (final Allele a : combined.getAlleles()) {
            if (toAdd.hasAttribute(a) && toAdd.getAttribute(a) != null) {
                if (combined.getAttribute(a) != null) {
                    combined.getAttribute(a).set(0, (int) combined.getAttribute(a).get(0) + (int) toAdd.getAttribute(a).get(0));
                    combined.getAttribute(a).set(1, (int) combined.getAttribute(a).get(1) + (int) toAdd.getAttribute(a).get(1));
                }
                else {
                    List<Integer> alleleData = new ArrayList<>();
                    alleleData.add(0, toAdd.getAttribute(a).get(0));
                    alleleData.add(1, toAdd.getAttribute(a).get(1));
                    combined.putAttribute(a,alleleData);
                }
            }
        }
    }

    protected String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleValues) {
        String annotationString = "";
        for (final Allele a : vcAlleles) {
            if (!annotationString.isEmpty())
                annotationString += printDelim;
            List<Integer> alleleValues = perAlleleValues.get(a);
            if (alleleValues == null)
                alleleValues = ZERO_LIST;
            annotationString += encode(alleleValues);
        }
        return annotationString;
    }

    protected String encode(List<Integer> alleleValues) {
        String annotationString = "";
        for (int j =0; j < alleleValues.size(); j++) {
            annotationString += alleleValues.get(j);
            if (j < alleleValues.size()-1)
                annotationString += ",";
        }
        return annotationString;
    }



    protected String makeReducedAnnotationString(VariantContext vc, Map<Allele,Double> perAltsStrandCounts) {
        String annotationString = "";
        for (Allele a : vc.getAlternateAlleles()) {
            if (!annotationString.isEmpty())
                annotationString += reducedDelim;
            if (!perAltsStrandCounts.containsKey(a))
                logger.warn("ERROR: VC allele not found in annotation alleles -- maybe there was trimming?");
            else
                annotationString += String.format("%.3f", perAltsStrandCounts.get(a));
        }
        return annotationString;
    }

    /**
     *
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    @Override
    public  Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        Map<String, Object> annotations = new HashMap<>();
        if (!vc.hasAttribute(getRawKeyName()))
            return new HashMap<>();
        String rawRankSumData = vc.getAttributeAsString(getRawKeyName(),null);
        if (rawRankSumData == null)
            return new HashMap<>();

        AlleleSpecificAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(originalVC.getAlleles(), rawRankSumData);
        parseRawDataString(myData);

        Map<Allele, Double> perAltRankSumResults = calculateReducedData(myData);

        String annotationString = makeReducedAnnotationString(vc, perAltRankSumResults);
        annotations.put(getKeyNames().get(0), annotationString);
        return annotations;
    }

    @Override
    public void calculateRawData(final VariantContext vc, Map<String, PerReadAlleleLikelihoodMap> pralm, final ReducibleAnnotationData rawAnnotations) {
        if(pralm == null)
            return;

        getStrandCountsFromLikelihoodMap(vc, pralm, rawAnnotations, MIN_COUNT);
    }

    protected abstract Map<Allele,Double> calculateReducedData(final AlleleSpecificAnnotationData<List<Integer>> combinedData );

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    public void getStrandCountsFromLikelihoodMap( final VariantContext vc,
                                                            final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                            final ReducibleAnnotationData perAlleleValues,
                                                            final int minCount) {
        if( stratifiedPerReadAlleleLikelihoodMap == null )
            return;
        if( vc == null )
            return;

        final Allele ref = vc.getReference();
        final List<Allele> allAlts = vc.getAlternateAlleles();

        for (final PerReadAlleleLikelihoodMap maps : stratifiedPerReadAlleleLikelihoodMap.values() ) {
            final ReducibleAnnotationData<List<Integer>> sampleTable = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            for (final Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : maps.getLikelihoodReadMap().entrySet()) {
                final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
                final GATKSAMRecord read = el.getKey();
                updateTable(mostLikelyAllele.getAlleleIfInformative(), read, ref, allAlts, sampleTable);
            }
            //for each sample (value in stratified PRALM), only include it if there are >minCount informative reads
            if ( passesMinimumThreshold(sampleTable, minCount) )
                combineAttributeMap(sampleTable, perAlleleValues);
        }
    }

    private void updateTable(final Allele bestAllele, final GATKSAMRecord read, final Allele ref, final List<Allele> allAlts, final ReducibleAnnotationData<List<Integer>> perAlleleValues) {

        final boolean matchesRef = bestAllele.equals(ref, true);
        final boolean matchesAnyAlt = allAlts.contains(bestAllele);

        //for uninformative reads
        if(bestAllele.isNoCall())
            return;

        //can happen if a read's most likely allele has been removed when --max_alternate_alleles is exceeded
        if (!( matchesRef || matchesAnyAlt ))
            return;

        final List<Integer> alleleStrandCounts;
        if (perAlleleValues.hasAttribute(bestAllele) && perAlleleValues.getAttribute(bestAllele) != null)
            alleleStrandCounts = perAlleleValues.getAttribute(bestAllele);
        else {
            alleleStrandCounts = new ArrayList<>();
            alleleStrandCounts.add(0,0);
            alleleStrandCounts.add(1,0);
            }
        if (read.isStrandless()) {
            // a strandless read counts as observations on both strand, at 50% weight, with a minimum of 1
            // (the 1 is to ensure that a strandless read always counts as an observation on both strands, even
            // if the read is only seen once, because it's a merged read or other)
            alleleStrandCounts.set(FORWARD, alleleStrandCounts.get(FORWARD)+1);
            alleleStrandCounts.set(REVERSE, alleleStrandCounts.get(REVERSE)+1);
        } else {
            // a normal read with an actual strand
            final boolean isFW = !read.getReadNegativeStrandFlag();
            if (isFW)
                alleleStrandCounts.set(FORWARD, alleleStrandCounts.get(FORWARD)+1);
            else
                alleleStrandCounts.set(REVERSE, alleleStrandCounts.get(REVERSE)+1);
        }
        perAlleleValues.putAttribute(bestAllele, alleleStrandCounts);
    }

    /**
     * Does this strand data array pass the minimum threshold for inclusion?
     *
     * @param sampleTable  the per-allele fwd/rev read counts for a single sample
     * @param minCount The minimum threshold of counts in the array
     * @return true if it passes the minimum threshold, false otherwise
     */
    protected boolean passesMinimumThreshold(final ReducibleAnnotationData<List<Integer>> sampleTable, final int minCount) {
        // the read total must be greater than MIN_COUNT
        int readTotal = 0;
        for (final List<Integer> alleleValues : sampleTable.getAttributeMap().values()) {
            if (alleleValues != null) {
                readTotal += alleleValues.get(FORWARD);
                readTotal += alleleValues.get(REVERSE);
            }
        }
        return readTotal > minCount;
    }


    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        return new HashMap<>();
    }

    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromStratifiedContexts(final Map<String, AlignmentContext> stratifiedContexts,
                                                                            final VariantContext vc){
        return new HashMap<>();
    }

    @Override
    //This just calls the non-allele-specific code in StrandBiasTest.java
    protected abstract Map<String, Object> calculateAnnotationFromLikelihoodMap(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                                       final VariantContext vc);

}
