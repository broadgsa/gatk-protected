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

package org.broadinstitute.gatk.tools.walkers.cancer.contamination;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.Utils;

import java.io.*;
import java.util.*;

/**
 * our contamination results object; this object aggregates the results of the contamination run over lanes, samples,
 * or whatever other divisor we've used on the read data
 */
public class ContaminationResults {

    public static class ContaminationData implements Comparable<ContaminationData> {
        private GenomeLoc site;
        private long basesMatching = 0l;
        private long basesMismatching = 0l;
        private double mismatchFraction = -1d;
        private double[] bins;
        private double p;

        public long getBasesMatching() {
            return basesMatching;
        }

        public long getBasesMismatching() {
            return basesMismatching;
        }

        public double getMismatchFraction() {
            return mismatchFraction;
        }

        public double[] getBins() {
            return bins;
        }

        public double getP() {
            return p;
        }

        public ContaminationData(GenomeLoc site, long basesMatching, long basesMismatching, double[] bins) {
            this.site = site;
            this.basesMatching = basesMatching;
            this.basesMismatching = basesMismatching;
            this.bins = bins;
            long totalBases = this.basesMatching + this.basesMismatching;
            if (totalBases != 0) {
                this.mismatchFraction = (double)this.basesMismatching / (double) totalBases;
            }

            int a = (int) this.getBasesMismatching() + 1;
            int b = (int) this.getBasesMatching() + 1;
            BetaDistribution dist = new BetaDistributionImpl(a,b);
            try {
                this.p = 1.0d - dist.cumulativeProbability(0.5d);
            } catch (MathException me) {
                throw new RuntimeException("Error! - " + me.toString());
            }

        }

        public int compareTo(ContaminationData other) {
            return -Double.compare(this.getP(), other.getP());
        }

        @Override
        public String toString() {
            return "ContaminationData{" +
                    "site=" + site +
                    ", basesMatching=" + basesMatching +
                    ", basesMismatching=" + basesMismatching +
                    ", mismatchFraction=" + mismatchFraction +
                    '}';
        }
    }


    // what precision are we using in our calculations
    private final double precision;

    // a map of our contamination targets and their stats
    // key: aggregation entity ("META", sample name, or lane name)
    // value: ContaminationStats (whcih
    private Map<String,Map<String, ContaminationStats>> stats = new HashMap<String,Map<String, ContaminationStats>>();

    public ContaminationResults(double precision) {
        this.precision = precision;
    }


    Map<String, Map<String, List<ContaminationData>>> storedData = new HashMap<String, Map<String, List<ContaminationData>>>();

    /**
     * add to the stats
     *
     * @param newAggregationStats a mapping of the stat name to their statistics collected
     */
    public void add(Map<String, Map<String, ContaminationStats>> newAggregationStats) {

        // for each aggregation level
        for (String aggregationKey : newAggregationStats.keySet()) {
            Map<String, ContaminationStats> populationContaminationStats = newAggregationStats.get(aggregationKey);


            // a new way of doing this... store all the data until the end...
            if (!storedData.containsKey(aggregationKey)) { storedData.put(aggregationKey, new HashMap<String, List<ContaminationData>>()); }
            for (String pop : populationContaminationStats.keySet()) {
                ContaminationStats newStats = populationContaminationStats.get(pop);

                // if it exists... just merge it
                if (!storedData.get(aggregationKey).containsKey(pop)) {
                    storedData.get(aggregationKey).put(pop, new ArrayList<ContaminationData>());
                }

                double[] newData = new double[newStats.getContamination().getBins().length];
                System.arraycopy(newStats.getContamination().getBins(),0,newData,0,newStats.getContamination().getBins().length);
                storedData.get(aggregationKey).get(pop).add(new ContaminationData(newStats.getSite(), newStats.getBasesMatching(), newStats.getBasesMismatching(), newData));
            }



            // merge the sets
            if (stats.containsKey(aggregationKey)) {

                // and for each population
                for (String pop : populationContaminationStats.keySet()) {
                    ContaminationStats newStats = populationContaminationStats.get(pop);

                    // if it exists... just merge it
                    if (stats.get(aggregationKey).containsKey(pop)) {
                        stats.get(aggregationKey).get(pop).add(newStats);
                    } else {
                        stats.get(aggregationKey).put(pop, newStats);
                    }
                }
            } else {
                stats.put(aggregationKey, populationContaminationStats);
            }
        }
    }

    /**
     * output the contamination data, and return the contamination data
     * @param out the output source
     * @return the contamination value
     */
    public void outputReport(double precision, PrintStream out, double fractionToTrim, double trimInterval, double betaThreshold) {
        out.println("name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites");

        for (Map.Entry<String,Map<String, ContaminationStats>> entry : stats.entrySet()) {
            for (ContaminationStats stats : entry.getValue().values()) {
                String aggregationLevel = entry.getKey();
                String population = stats.getContamination().getPopultationName();

                List<ContaminationData> newStats = storedData.get(aggregationLevel).get(population);
                String pm = "%3." + Math.round(Math.log10(1/precision)) +"f";

                int bins = newStats.iterator().next().getBins().length;
                int maxTrim = (int) Math.floor((double)(newStats.size()) * fractionToTrim);

                // sort the collection
                Collections.sort(newStats);

                List<ContaminationData> data = new ArrayList<ContaminationData>(newStats);

                // trim sites with > 95% p of being > 0.5 f (based on beta distribution)
                int trimmed = 0;
                for(Iterator<ContaminationData> i = data.iterator(); trimmed < maxTrim && i.hasNext();) {
                        ContaminationData x = i.next();
                        if (x.getP() >= betaThreshold) {
                            System.out.println("Trimming " + x.toString() + " with p(f>=0.5) >= " + betaThreshold + " with a value of  " + x.getP());
                            i.remove();
                            trimmed++;
                        }
                }

                double[][] matrix = new double[bins][data.size()];

                for (int i = 0; i<bins; i++) {
                    for (int j=0; j<data.size(); j++) {
                        matrix[i][j] = data.get(j).getBins()[i];
                    }
                }

                // now perform the sum
                double[] output = new double[bins];
                for (int i = 0; i<bins; i++) {
                    double[] binData = matrix[i];

                    // remove the top and bottom
                    output[i] = 0;
                    for (int x = 0; x < binData.length; x++) {
                        output[i] += binData[x];
                    }
                }
                double[] newTrimmedStats = output;

                // get the confidence interval, at the set width
                ContaminationEstimate.ConfidenceInterval newInterval = new ContaminationEstimate.ConfidenceInterval(newTrimmedStats, 0.95);

                out.println(
                        String.format("%s\t%s\t%s\t"+pm+"\t"+pm+"\t"+pm+"\t"+pm+"\t"+"%d",
                                aggregationLevel,
                                population,
                                "n/a",
                                newInterval.getContamination(),
                                (newInterval.getStop() - newInterval.getStart()),
                                newInterval.getStart(),
                                newInterval.getStop(),
                                data.size())
                );

            }
        }
    }


    public void writeCurves(PrintStream out) {
        boolean outputBins = false;
        for (Map.Entry<String, Map<String, ContaminationStats>> entry : stats.entrySet()) {
            for (ContaminationStats stats : entry.getValue().values()) {
                if (!outputBins) {
                    String[] bins = new String[stats.getContamination().getBins().length];
                    for (int index = 0; index < stats.getContamination().getBins().length; index++)
                        bins[index] = String.valueOf(100.0 * (1 - (double) index / stats.getContamination().getBins().length));
                    outputBins = true;
                    out.print("name,pop,");
                    out.println(Utils.join(",",bins));
                }
                String[] bins = new String[stats.getContamination().getBins().length];
                int index = 0;
                for (double value : stats.getContamination().getBins())
                    bins[index++] = String.valueOf(value);
                out.print(entry.getKey()+",\""+stats.getContamination().getPopultationName()+"\",");
                out.println(Utils.join(",", bins));
            }
        }
    }

    public Map<String, Map<String, ContaminationStats>> getStats() {
        return Collections.unmodifiableMap(stats);
    }

    public void setStats(Map<String, Map<String,ContaminationStats>> stats) {
        this.stats = stats;
    }
}