package org.broadinstitute.gatk.tools.walkers.coverage;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.*;


/**
 * This is based on the CallableLoci walker, designed to report sites with more than 2 alleles, where
 * the third allele is above the supplied threshold.  The original purpose was to find positions of potential duplications
 * or poor genome quality (i.e. reads from more than one duplicated gene are mapping to a single site)
 *
 * <h3>Usage example</h3>
 * <pre>
 *  java -jar GenomeAnalysisTK.jar \
 *     -T MultipleAllelesAtLoci \
 *     -R reference.fasta \
 *     -I bams.list \
 *     -o output.bed
 * </pre>
 * <p/>
 * would produce a BED file that looks like:
 * <p/>
 * <pre>
 *     chr20 10000000 10000001 2 Subj1:A/0.13;Subj12:A/0.23;
 *     chr20 10000865 10000866 2 Subj4:T/0.12
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@By(DataSource.REFERENCE)
public class MultipleAllelesAtLoci extends LocusWalker<MultipleAllelesAtLoci.SiteBaseCounter, Integer> {
    @Output
    PrintStream out;

    /**
     * Reads with MAPQ > minMappingQuality are treated as usable for variation detection, contributing to the PASS
     * state.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth.", required = false)
    byte minMappingQuality = 10;

    /**
     * Bases with less than minBaseQuality are viewed as not sufficiently high quality to contribute to the PASS state
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth.", required = false)
    byte minBaseQuality = 20;

    /**
     * If the number of QC+ bases (on reads with MAPQ > minMappingQuality and with base quality > minBaseQuality) exceeds this
     * value and is less than maxDepth the site is considered PASS.
     */
    @Advanced
    @Argument(fullName = "minDepth", shortName = "minDepth", doc = "Minimum QC+ read depth before a locus is considered callable", required = false)
    int minDepth = 4;

    @Argument(fullName = "minBasePct", shortName = "minBasePct", doc = "If a given site has more than three alleles present, any alleles beyond the top two that are above this threshold are reported.", required = false)
    double minBasePct = 0.1;

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    @Override
    public boolean includeReadsWithDeletionAtLoci() {
        return true;
    }

    @Override
    public void initialize() {

    }

    protected static class SiteBaseCounter implements HasGenomeLocation {
        final public GenomeLocParser genomeLocParser;
        public GenomeLoc loc;
        final List<FlaggedSite> flaggedSites;

        public SiteBaseCounter(GenomeLocParser genomeLocParser, GenomeLoc loc) {
            this.genomeLocParser = genomeLocParser;
            this.loc = loc;
            this.flaggedSites = new ArrayList<>();
        }

        public GenomeLoc getLocation() {
            return loc;
        }

        public void merge(SiteBaseCounter counter2) {
            if (!loc.equals(counter2.loc)) {
                throw new RuntimeException(String.format("SiteBaseCounters do not have the same locus: %s: %d, %s: %d", loc.getContig(), loc.getStart(), counter2.loc.getContig(), counter2.loc.getStart()));
            }

            flaggedSites.addAll(counter2.flaggedSites);
        }

        public void doPrint(PrintStream out) {
            if (!flaggedSites.isEmpty()) {
                Set<String> totalSubjects = new HashSet<>();
                List<String> comments = new ArrayList<>();
                for (FlaggedSite f : flaggedSites) {
                    totalSubjects.add(f.sample);
                    comments.add(String.format("%s: [%s/%f/%d/%d]", f.sample, (char)(f.base.byteValue()), f.pct, f.baseCount, f.depth));
                }

                out.println(String.format("%s\t%d\t%d\t%s\t%s\t%s\t%s", loc.getContig(), loc.getStart()-1, loc.getStop(), "MultiAllelicSite", totalSubjects.size(), "+", StringUtils.join(comments, ";")));
            }
        }

        public void addSite(FlaggedSite site) {
            flaggedSites.add(site);
        }
    }

    private static class FlaggedSite {
        public String sample;
        public Byte base;
        public double pct;
        public long depth;
        public long baseCount;

        public FlaggedSite(String sample, Byte base, double pct, long depth, long baseCount) {
            this.sample = sample;
            this.base = base;
            this.pct = pct;
            this.depth = depth;
            this.baseCount = baseCount;
        }
    }

    @Override
    public SiteBaseCounter map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        SiteBaseCounter counter = new SiteBaseCounter(getToolkit().getGenomeLocParser(), context.getLocation());

        Map<DoCOutputType.Partition,Map<String,int[]>> counts = CoverageUtils.getBaseCountsByPartition(context,minMappingQuality,Integer.MAX_VALUE,minBaseQuality,Byte.MAX_VALUE,CoverageUtils.CountPileupType.COUNT_READS, EnumSet.of(DoCOutputType.Partition.sample));
        for (String sn : getSampleDB().getSampleNames()) {
            int[] baseCounts = counts.get(DoCOutputType.Partition.sample) != null ? counts.get(DoCOutputType.Partition.sample).get(sn) : null;
            if (baseCounts == null) {
                //logger.warn(String.format("no base counts found for sample: %s, position: %s %d", sn, context.getLocation().getContig(), context.getLocation().getStart()));
                continue;
            }

            int total = sumOfArray(baseCounts);
            if (total < minDepth) {
                continue;
            }

            Map<Integer, Double> baseCountPct = getBaseCountPctMap(baseCounts, total);
            List<Map.Entry<Integer, Double>> vals = sortByValue(baseCountPct);
            if (vals.size() <= 2) {
                //nothing to do
            }
            else {
                //accept the top two bases
                vals.remove(0);
                vals.remove(0);

                //iterate remaining by pct
                for (Map.Entry<Integer, Double> entry : vals) {
                    if (entry.getValue() >= minBasePct) {
                        try {
                            counter.addSite(new FlaggedSite(sn, BaseUtils.EXTENDED_BASES[entry.getKey()], entry.getValue(), total, baseCounts[entry.getKey()]));
                        }
                        catch (IndexOutOfBoundsException e) {
                            logger.warn(String.format("index of out bounds, %d, %d, %s, %d", entry.getKey(), baseCounts.length, context.getLocation().getContig(), context.getLocation().getStart()));
                            throw e;
                        }
                    }
                }
            }
        }

        // TODO: consider moving
        // print to BED
        counter.doPrint(out);

        return counter;
    }

    public static <K, V extends Comparable<? super V>> List<Map.Entry<K, V>> sortByValue( Map<K, V> map ) {
        List<Map.Entry<K, V>> list =  new LinkedList<Map.Entry<K, V>>( map.entrySet() );
        Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
            public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 )
            {
                return (o1.getValue()).compareTo( o2.getValue() );
            }
        } );

        Collections.reverse(list);

        return list;
    }

    private Map<Integer, Double> getBaseCountPctMap(int[] baseCounts, int total)
    {
        Map<Integer, Double> ret = new HashMap<>();

        for (int i = 0; i < baseCounts.length; i++)
        {
            ret.put(i,  total == 0 ? 0 : ((double)baseCounts[i]) / total);
        }

        return ret;
    }

    private int sumOfArray(int[] array)
    {
        int ret = 0;
        for (Integer i : array)
        {
            ret += i;
        }

        return ret;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(SiteBaseCounter counter1, Integer i) {
        return 0;
    }

    @Override
    public void onTraversalDone(Integer i) {

    }
}