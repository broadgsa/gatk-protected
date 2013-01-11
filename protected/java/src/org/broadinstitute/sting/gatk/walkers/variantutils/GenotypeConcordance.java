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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.PrintStream;
import java.util.*;

/**
 * A simple walker for performing genotype concordance calculations between two callsets
 */
public class GenotypeConcordance extends RodWalker<Pair<VariantContext,VariantContext>,ConcordanceMetrics> {

    @Input(fullName="eval",shortName="eval",doc="The variants and genotypes to evaluate",required=true)
    RodBinding<VariantContext> evalBinding;

    @Input(fullName="comp",shortName="comp",doc="The variants and genotypes to compare against",required=true)
    RodBinding<VariantContext> compBinding;

    @Argument(fullName="ignoreFilters",doc="Filters will be ignored",required=false)
    boolean ignoreFilters = false;

    @Output
    PrintStream out;

    List<String> evalSamples;
    List<String> compSamples;

    // todo -- integration test coverage
    // todo -- deal with occurrences like:
    //     Eval: 20   4000     A     C
    //     Eval: 20   4000     A    AC
    //     Comp: 20   4000     A     C
    //  currently this results in a warning and skipping
    // todo -- extend to multiple eval, multiple comp
    // todo -- table with "proportion of overlapping sites" (not just eval/comp margins)


    public ConcordanceMetrics reduceInit() {
        Map<String,VCFHeader> headerMap = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(evalBinding,compBinding));
        VCFHeader evalHeader = headerMap.get(evalBinding.getName());
        evalSamples = evalHeader.getGenotypeSamples();
        VCFHeader compHeader = headerMap.get(compBinding.getName());
        compSamples = compHeader.getGenotypeSamples();
        return new ConcordanceMetrics(evalHeader,compHeader);
    }


    public Pair<VariantContext,VariantContext> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Pair<VariantContext,VariantContext> evalCompPair = null;
        if ( tracker != null && (
                tracker.getValues(evalBinding,ref.getLocus()).size() > 0 ||
                tracker.getValues(compBinding,ref.getLocus()).size() > 0 ) ) {

            List<VariantContext> eval = tracker.getValues(evalBinding,ref.getLocus());
            List<VariantContext> comp = tracker.getValues(compBinding,ref.getLocus());
            if ( eval.size() > 1 || comp.size() > 1 ) {
                logger.warn("Eval or Comp Rod at position "+ref.getLocus().toString()+" has multiple records. Site will be skipped.");
                return evalCompPair;
            }
            // if a rod is missing, explicitly create a variant context with 'missing' genotypes. Slow, but correct.
            // note that if there is no eval rod there must be a comp rod, and also the reverse
            VariantContext evalContext = eval.size() == 1 ? eval.get(0) : createEmptyContext(ref,comp.get(0),evalSamples);
            VariantContext compContext = comp.size() == 1 ? comp.get(0) : createEmptyContext(ref,eval.get(0),compSamples);
            evalContext = filterGenotypes(evalContext,ignoreFilters);
            compContext = filterGenotypes(compContext,ignoreFilters);
            evalCompPair = new Pair<VariantContext, VariantContext>(evalContext,compContext);
        }

        return evalCompPair;
    }

    public ConcordanceMetrics reduce(Pair<VariantContext,VariantContext> evalComp, ConcordanceMetrics metrics) {
        if ( evalComp != null )
            metrics.update(evalComp.getFirst(),evalComp.getSecond());
        return metrics;
    }

    public void onTraversalDone(ConcordanceMetrics metrics) {
        GATKReport report = new GATKReport();
        GATKReportTable concordanceCounts = new GATKReportTable("GenotypeConcordance_Counts","Per-sample concordance tables: comparison counts",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceEvalProportions = new GATKReportTable("GenotypeConcordance_EvalProportions", "Per-sample concordance tables: proportions of genotypes called in eval",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceCompProportions = new GATKReportTable("GenotypeConcordance_CompProportions", "Per-sample concordance tables: proportions of genotypes called in comp",2+GenotypeType.values().length*GenotypeType.values().length);
        GATKReportTable concordanceSummary = new GATKReportTable("GenotypeConcordance_Summary","Per-sample summary statistics: NRS and NRD",2);
        GATKReportTable siteConcordance = new GATKReportTable("SiteConcordance_Summary","Site-level summary statistics",ConcordanceMetrics.SiteConcordanceType.values().length);
        concordanceCompProportions.addColumn("Sample","%s");
        concordanceCounts.addColumn("Sample","%s");
        concordanceEvalProportions.addColumn("Sample","%s");
        concordanceSummary.addColumn("Sample","%s");
        for ( GenotypeType evalType : GenotypeType.values() ) {
            for ( GenotypeType compType : GenotypeType.values() ) {
                String colKey = String.format("%s_%s", evalType.toString(), compType.toString());
                concordanceCounts.addColumn(colKey,"%d");
                if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                    concordanceEvalProportions.addColumn(colKey,"%.3f");
                if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                    concordanceCompProportions.addColumn(colKey,"%.3f");
            }
        }
        concordanceEvalProportions.addColumn("Mismatching_Alleles","%.3f");
        concordanceCompProportions.addColumn("Mismatching_Alleles","%.3f");
        concordanceCounts.addColumn("Mismatching_Alleles","%d");
        concordanceSummary.addColumn("Non-Reference Sensitivity","%.3f");
        concordanceSummary.addColumn("Non-Reference Discrepancy","%.3f");
        for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
            siteConcordance.addColumn(type.toString(),"%d");
        }

        for ( Map.Entry<String,ConcordanceMetrics.GenotypeConcordanceTable> entry : metrics.getPerSampleGenotypeConcordance().entrySet() ) {
            ConcordanceMetrics.GenotypeConcordanceTable table = entry.getValue();
            concordanceEvalProportions.set(entry.getKey(),"Sample",entry.getKey());
            concordanceCompProportions.set(entry.getKey(),"Sample",entry.getKey());
            concordanceCounts.set(entry.getKey(),"Sample",entry.getKey());
            for ( GenotypeType evalType : GenotypeType.values() ) {
                for ( GenotypeType compType : GenotypeType.values() ) {
                    String colKey = String.format("%s_%s",evalType.toString(),compType.toString());
                    int count = table.get(evalType, compType);
                    concordanceCounts.set(entry.getKey(),colKey,count);
                    if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                        concordanceEvalProportions.set(entry.getKey(),colKey,( (double) count)/table.getnEvalGenotypes(evalType));
                    if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                        concordanceCompProportions.set(entry.getKey(),colKey,( (double) count)/table.getnCompGenotypes(compType));
                }
            }
            concordanceEvalProportions.set(entry.getKey(),"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes());
            concordanceCompProportions.set(entry.getKey(),"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes());
            concordanceCounts.set(entry.getKey(),"Mismatching_Alleles",table.getnMismatchingAlt());
        }

        String rowKey = "ALL";
        concordanceCompProportions.set(rowKey,"Sample",rowKey);
        concordanceEvalProportions.set(rowKey,"Sample",rowKey);
        concordanceCounts.set(rowKey,"Sample",rowKey);
        ConcordanceMetrics.GenotypeConcordanceTable table = metrics.getOverallGenotypeConcordance();
        for ( GenotypeType evalType : GenotypeType.values() ) {
            for ( GenotypeType compType : GenotypeType.values() ) {
                String colKey = String.format("%s_%s",evalType.toString(),compType.toString());
                int count = table.get(evalType,compType);
                concordanceCounts.set(rowKey,colKey,count);
                if ( evalType == GenotypeType.HET || evalType == GenotypeType.HOM_REF || evalType == GenotypeType.HOM_VAR)
                    concordanceEvalProportions.set(rowKey,colKey,( (double) count)/table.getnEvalGenotypes(evalType));
                if ( compType == GenotypeType.HET || compType == GenotypeType.HOM_VAR || compType == GenotypeType.HOM_REF )
                    concordanceCompProportions.set(rowKey,colKey,( (double) count)/table.getnCompGenotypes(compType));
            }
        }
        concordanceEvalProportions.set(rowKey,"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledEvalGenotypes());
        concordanceCompProportions.set(rowKey,"Mismatching_Alleles", ( (double) table.getnMismatchingAlt() )/table.getnCalledCompGenotypes());
        concordanceCounts.set(rowKey,"Mismatching_Alleles",table.getnMismatchingAlt());

        for ( Map.Entry<String,Double> nrsEntry : metrics.getPerSampleNRS().entrySet() ) {
            concordanceSummary.set(nrsEntry.getKey(),"Sample",nrsEntry.getKey());
            concordanceSummary.set(nrsEntry.getKey(),"Non-Reference Sensitivity",nrsEntry.getValue());
        }
        for ( Map.Entry<String,Double> nrdEntry : metrics.getPerSampleNRD().entrySet() ) {
            concordanceSummary.set(nrdEntry.getKey(),"Non-Reference Discrepancy",nrdEntry.getValue());
        }
        concordanceSummary.set("ALL","Sample","ALL");
        concordanceSummary.set("ALL","Non-Reference Sensitivity",metrics.getOverallNRS());
        concordanceSummary.set("ALL","Non-Reference Discrepancy",metrics.getOverallNRD());

        for (ConcordanceMetrics.SiteConcordanceType type : ConcordanceMetrics.SiteConcordanceType.values() ) {
            siteConcordance.set("Comparison",type.toString(),metrics.getOverallSiteConcordance().get(type));
        }

        report.addTable(concordanceCompProportions);
        report.addTable(concordanceEvalProportions);
        report.addTable(concordanceCounts);
        report.addTable(concordanceSummary);
        report.addTable(siteConcordance);

        report.print(out);
    }

    public VariantContext createEmptyContext(ReferenceContext ref, VariantContext other, List<String> samples) {
        VariantContextBuilder builder = new VariantContextBuilder();
        // set the alleles to be the same
        builder.alleles(other.getAlleles());
        builder.loc(other.getChr(),other.getStart(),other.getEnd());
        // set all genotypes to empty
        List<Genotype> genotypes = new ArrayList<Genotype>(samples.size());
        for ( String sample : samples )
            genotypes.add(GenotypeBuilder.create(sample, new ArrayList<Allele>(0)));
        builder.genotypes(genotypes);
        return builder.make();
    }

    public VariantContext filterGenotypes(VariantContext context, boolean ignoreSiteFilter) {
        // placeholder method for genotype-level filtering. However if the site itself is filtered,
        // and such filters are not ignored, the genotype-level data should be altered to reflect this
        if ( ! context.isFiltered() || ignoreSiteFilter ) {
            // todo -- add genotype-level jexl filtering here
            return context;
        }
        VariantContextBuilder builder = new VariantContextBuilder();
        builder.alleles(Arrays.asList(context.getReference()));
        builder.loc(context.getChr(),context.getStart(),context.getEnd());
        List<Genotype> newGeno = new ArrayList<Genotype>(context.getNSamples());
        for ( Genotype g : context.getGenotypes().iterateInSampleNameOrder() ) {
            newGeno.add(GenotypeBuilder.create(g.getSampleName(),new ArrayList<Allele>()));
        }
        builder.genotypes(newGeno);
        return builder.make();
    }
}