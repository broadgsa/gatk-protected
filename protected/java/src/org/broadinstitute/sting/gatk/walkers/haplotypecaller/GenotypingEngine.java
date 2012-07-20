/*
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFAlleleClipper;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

public class GenotypingEngine {

    private final boolean DEBUG;
    private final int MNP_LOOK_AHEAD;
    private final boolean OUTPUT_FULL_HAPLOTYPE_SEQUENCE;
    private final static List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
    private final static Allele SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele.create("<UNASSEMBLED_EVENT>", false);

    public GenotypingEngine( final boolean DEBUG, final int MNP_LOOK_AHEAD, final boolean OUTPUT_FULL_HAPLOTYPE_SEQUENCE ) {
        this.DEBUG = DEBUG;
        this.MNP_LOOK_AHEAD = MNP_LOOK_AHEAD;
        this.OUTPUT_FULL_HAPLOTYPE_SEQUENCE = OUTPUT_FULL_HAPLOTYPE_SEQUENCE;
        noCall.add(Allele.NO_CALL);
    }

    // This function is the streamlined approach, currently not being used
    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    public List<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> assignGenotypeLikelihoodsAndCallHaplotypeEvents( final UnifiedGenotyperEngine UG_engine, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc,
                                                                                                                             final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser ) {
        // Prepare the list of haplotype indices to genotype
        final ArrayList<Allele> allelesToGenotype = new ArrayList<Allele>();

        for( final Haplotype h : haplotypes ) {
            allelesToGenotype.add( Allele.create(h.getBases(), h.isReference()) );
        }
        final int numHaplotypes = haplotypes.size();

        // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
        final GenotypesContext genotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
        for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
            final double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes+1) / 2];
            final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(haplotypes, sample);
            int glIndex = 0;
            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                }
            }
            genotypes.add(new GenotypeBuilder(sample, noCall).PL(genotypeLikelihoods).make());
        }
        final VariantCallContext call = UG_engine.calculateGenotypes(new VariantContextBuilder().loc(activeRegionWindow).alleles(allelesToGenotype).genotypes(genotypes).make(), UG_engine.getUAC().GLmodel);
        if( call == null ) { return Collections.emptyList(); } // exact model says that the call confidence is below the specified confidence threshold so nothing to do here

        // Prepare the list of haplotypes that need to be run through Smith-Waterman for output to VCF
        final ArrayList<Haplotype> haplotypesToRemove = new ArrayList<Haplotype>();
        for( final Haplotype h : haplotypes ) {
            if( call.getAllele(h.getBases()) == null ) { // exact model removed this allele from the list so no need to run SW and output to VCF
                haplotypesToRemove.add(h);
            }
        }
        haplotypes.removeAll(haplotypesToRemove);

        if( OUTPUT_FULL_HAPLOTYPE_SEQUENCE ) {
            final List<Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>>> returnVCs = new ArrayList<Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>>>();
            // set up the default 1-to-1 haplotype mapping object
            final HashMap<Allele,ArrayList<Haplotype>> haplotypeMapping = new HashMap<Allele,ArrayList<Haplotype>>();
            for( final Haplotype h : haplotypes ) {
                final ArrayList<Haplotype> list = new ArrayList<Haplotype>();
                list.add(h);
                haplotypeMapping.put(call.getAllele(h.getBases()), list);
            }
            returnVCs.add( new Pair<VariantContext, HashMap<Allele, ArrayList<Haplotype>>>(call,haplotypeMapping) );
            return returnVCs;
        }

        final ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> returnCalls = new ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>>();

        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        int count = 0;
        if( DEBUG ) { System.out.println("=== Best Haplotypes ==="); }
        for( final Haplotype h : haplotypes ) {
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + h.getCigar() );
            }
            // Walk along the alignment and turn any difference from the reference into an event
            h.setEventMap( generateVCsFromAlignment( h.getAlignmentStartHapwrtRef(), h.getCigar(), ref, h.getBases(), refLoc, "HC" + count++, MNP_LOOK_AHEAD ) );
            startPosKeySet.addAll(h.getEventMap().keySet());
        }
        
        // Create the VC merge priority list
        final ArrayList<String> priorityList = new ArrayList<String>();
        for( int iii = 0; iii < haplotypes.size(); iii++ ) {
            priorityList.add("HC" + iii);
        }
        
        // Walk along each position in the key set and create each event to be outputted
        for( final int loc : startPosKeySet ) {
            if( loc >= activeRegionWindow.getStart() && loc <= activeRegionWindow.getStop() ) {
                final ArrayList<VariantContext> eventsAtThisLoc = new ArrayList<VariantContext>();
                for( final Haplotype h : haplotypes ) {
                    final HashMap<Integer,VariantContext> eventMap = h.getEventMap();
                    final VariantContext vc = eventMap.get(loc);
                    if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                        eventsAtThisLoc.add(vc);
                    }
                }
                
                // Create the allele mapping object which maps the original haplotype alleles to the alleles present in just this event
                final ArrayList<ArrayList<Haplotype>> alleleMapper = createAlleleMapper( loc, eventsAtThisLoc, haplotypes );

                // Merge the event to find a common reference representation
                final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, eventsAtThisLoc, priorityList, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);

                final HashMap<Allele, ArrayList<Haplotype>> alleleHashMap = new HashMap<Allele, ArrayList<Haplotype>>();
                int aCount = 0;
                for( final Allele a : mergedVC.getAlleles() ) {
                    alleleHashMap.put(a, alleleMapper.get(aCount++)); // BUGBUG: needs to be cleaned up and merged with alleleMapper
                }

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    //System.out.println("Event/haplotype allele mapping = " + alleleMapper);
                }

                // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
                final GenotypesContext myGenotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
                for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
                    final int myNumHaplotypes = alleleMapper.size();
                    final double[] genotypeLikelihoods = new double[myNumHaplotypes * (myNumHaplotypes+1) / 2];
                    final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(sample, alleleMapper);
                    int glIndex = 0;
                    for( int iii = 0; iii < myNumHaplotypes; iii++ ) {
                        for( int jjj = 0; jjj <= iii; jjj++ ) {
                            genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                        }
                    }

                    // using the allele mapping object translate the haplotype allele into the event allele
                    final Genotype g = new GenotypeBuilder(sample)
                            .alleles(findEventAllelesInSample(mergedVC.getAlleles(), call.getAlleles(), call.getGenotype(sample).getAlleles(), alleleMapper, haplotypes))
                            .phased(loc != startPosKeySet.first())
                            .PL(genotypeLikelihoods).make();
                    myGenotypes.add(g);
                }
                returnCalls.add( new Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>(
                                 new VariantContextBuilder(mergedVC).log10PError(call.getLog10PError()).genotypes(myGenotypes).make(), alleleHashMap) );
            }
        }
        return returnCalls;
    }

    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    public List<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> assignGenotypeLikelihoodsAndCallIndependentEvents( final UnifiedGenotyperEngine UG_engine, final ArrayList<Haplotype> haplotypes, final byte[] ref, final GenomeLoc refLoc,
                                                                                                                               final GenomeLoc activeRegionWindow, final GenomeLocParser genomeLocParser, final ArrayList<VariantContext> activeAllelesToGenotype ) {

        final ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>> returnCalls = new ArrayList<Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>>();

        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        int count = 0;
        if( DEBUG ) { System.out.println("=== Best Haplotypes ==="); }
        for( final Haplotype h : haplotypes ) {
            // Walk along the alignment and turn any difference from the reference into an event
            h.setEventMap( generateVCsFromAlignment( h, h.getAlignmentStartHapwrtRef(), h.getCigar(), ref, h.getBases(), refLoc, "HC" + count++, MNP_LOOK_AHEAD ) );
            if( activeAllelesToGenotype.isEmpty() ) { startPosKeySet.addAll(h.getEventMap().keySet()); }
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + h.getCigar() );
                System.out.println( "> Left and right breaks = (" + h.leftBreakPoint + " , " + h.rightBreakPoint + ")");
                System.out.println( ">> Events = " + h.getEventMap());
            }
        }
        // Create the VC merge priority list
        final ArrayList<String> priorityList = new ArrayList<String>();
        for( int iii = 0; iii < haplotypes.size(); iii++ ) {
            priorityList.add("HC" + iii);
        }

        cleanUpSymbolicUnassembledEvents( haplotypes, priorityList );
        if( activeAllelesToGenotype.isEmpty() && haplotypes.get(0).getSampleKeySet().size() >= 3 ) { // if not in GGA mode and have at least 3 samples try to create MNP and complex events by looking at LD structure
            mergeConsecutiveEventsBasedOnLD( haplotypes, startPosKeySet, ref, refLoc );
        }
        if( !activeAllelesToGenotype.isEmpty() ) { // we are in GGA mode!
            for( final VariantContext compVC : activeAllelesToGenotype ) {
                startPosKeySet.add( compVC.getStart() );
            }
        }


        // Walk along each position in the key set and create each event to be outputted
        for( final int loc : startPosKeySet ) {
            if( loc >= activeRegionWindow.getStart() && loc <= activeRegionWindow.getStop() ) {
                final ArrayList<VariantContext> eventsAtThisLoc = new ArrayList<VariantContext>();
                if( activeAllelesToGenotype.isEmpty() ) {
                    for( final Haplotype h : haplotypes ) {
                        final HashMap<Integer,VariantContext> eventMap = h.getEventMap();
                        final VariantContext vc = eventMap.get(loc);
                        if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                            eventsAtThisLoc.add(vc);
                        }
                    }
                } else { // we are in GGA mode!
                    for( final VariantContext compVC : activeAllelesToGenotype ) {
                        if( compVC.getStart() == loc ) {
                            priorityList.clear();
                            int alleleCount = 0;
                            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                                HashSet<Allele> alleleSet = new HashSet<Allele>(2);
                                alleleSet.add(compVC.getReference());
                                alleleSet.add(compAltAllele);
                                priorityList.add("Allele" + alleleCount);
                                eventsAtThisLoc.add(new VariantContextBuilder(compVC).alleles(alleleSet).source("Allele"+alleleCount).make());
                                alleleCount++;
                            }
                        }
                    }
                }

                if( eventsAtThisLoc.isEmpty() ) { continue; }

                // Create the allele mapping object which maps the original haplotype alleles to the alleles present in just this event
                final ArrayList<ArrayList<Haplotype>> alleleMapper = createAlleleMapper( loc, eventsAtThisLoc, haplotypes );

                // Merge the event to find a common reference representation
                final VariantContext mergedVC = VariantContextUtils.simpleMerge(genomeLocParser, eventsAtThisLoc, priorityList, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);
                if( mergedVC == null ) { continue; }

                final HashMap<Allele, ArrayList<Haplotype>> alleleHashMap = new HashMap<Allele, ArrayList<Haplotype>>();
                int aCount = 0;
                for( final Allele a : mergedVC.getAlleles() ) {
                    alleleHashMap.put(a, alleleMapper.get(aCount++)); // BUGBUG: needs to be cleaned up and merged with alleleMapper
                }

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    //System.out.println("Event/haplotype allele mapping = " + alleleMapper);
                }

                // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
                final GenotypesContext genotypes = GenotypesContext.create(haplotypes.get(0).getSampleKeySet().size());
                for( final String sample : haplotypes.get(0).getSampleKeySet() ) { // BUGBUG: assume all haplotypes saw the same samples
                    final int numHaplotypes = alleleMapper.size();
                    final double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes+1) / 2];
                    final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(sample, alleleMapper);
                    int glIndex = 0;
                    for( int iii = 0; iii < numHaplotypes; iii++ ) {
                        for( int jjj = 0; jjj <= iii; jjj++ ) {
                            genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                        }
                    }
                    genotypes.add( new GenotypeBuilder(sample).alleles(noCall).PL(genotypeLikelihoods).make() );
                }
                final VariantCallContext call = UG_engine.calculateGenotypes(new VariantContextBuilder(mergedVC).genotypes(genotypes).make(), UG_engine.getUAC().GLmodel);

                if( call != null ) {
                    returnCalls.add( new Pair<VariantContext, HashMap<Allele,ArrayList<Haplotype>>>(call, alleleHashMap) );
                }
            }
        }
        return returnCalls;
    }

    protected static void cleanUpSymbolicUnassembledEvents( final ArrayList<Haplotype> haplotypes, final ArrayList<String> priorityList ) {
        final ArrayList<Haplotype> haplotypesToRemove = new ArrayList<Haplotype>();
        final ArrayList<String> stringsToRemove = new ArrayList<String>();
        for( final Haplotype h : haplotypes ) {
            for( final VariantContext vc : h.getEventMap().values() ) {
                if( vc.isSymbolic() ) {
                    for( final Haplotype h2 : haplotypes ) {
                        for( final VariantContext vc2 : h2.getEventMap().values() ) {
                            if( vc.getStart() == vc2.getStart() && vc2.isIndel() ) {
                                haplotypesToRemove.add(h);
                                stringsToRemove.add(vc.getSource());
                                break;
                            }
                        }
                    }
                }
            }
        }
        haplotypes.removeAll(haplotypesToRemove);
        priorityList.removeAll(stringsToRemove);
    }

    protected void mergeConsecutiveEventsBasedOnLD( final ArrayList<Haplotype> haplotypes, final TreeSet<Integer> startPosKeySet, final byte[] ref, final GenomeLoc refLoc ) {
        final int MAX_SIZE_TO_COMBINE = 10;
        final double MERGE_EVENTS_R2_THRESHOLD = 0.95;
        if( startPosKeySet.size() <= 1 ) { return; }

        boolean mapWasUpdated = true;
        while( mapWasUpdated ) {
            mapWasUpdated = false;

            // loop over the set of start locations and consider pairs that start near each other
            final Iterator<Integer> iter = startPosKeySet.iterator();
            int thisStart = iter.next();
            while( iter.hasNext() ) {
                final int nextStart = iter.next();
                if( nextStart - thisStart < MAX_SIZE_TO_COMBINE) {
                    boolean isBiallelic = true;
                    VariantContext thisVC = null;
                    VariantContext nextVC = null;
                    int x11 = 0;
                    int x12 = 0;
                    int x21 = 0;
                    int x22 = 0;

                    for( final Haplotype h : haplotypes ) {
                        // only make complex substitutions out of consecutive biallelic sites
                        final VariantContext thisHapVC = h.getEventMap().get(thisStart);
                        if( thisHapVC != null && !thisHapVC.isSymbolic() ) { // something was found at this location on this haplotype
                            if( thisVC == null ) {
                                thisVC = thisHapVC;
                            } else if( !thisHapVC.hasSameAllelesAs( thisVC ) ) {
                                isBiallelic = false;
                                break;
                            }
                        }
                        final VariantContext nextHapVC = h.getEventMap().get(nextStart);
                        if( nextHapVC != null && !nextHapVC.isSymbolic() ) { // something was found at the next location on this haplotype
                            if( nextVC == null ) {
                                nextVC = nextHapVC;
                            } else if( !nextHapVC.hasSameAllelesAs( nextVC ) ) {
                                isBiallelic = false;
                                break;
                            }
                        }
                        // count up the co-occurrences of the events for the R^2 calculation
                        // BUGBUG: use haplotype likelihoods per-sample to make this more accurate
                        if( thisHapVC == null ) {
                            if( nextHapVC == null ) { x11++; }
                            else { x12++; }
                        } else {
                            if( nextHapVC == null ) { x21++; }
                            else { x22++; }
                        }
                    }
                    if( thisVC == null || nextVC == null ) {
                        continue;
                        //throw new ReviewedStingException("StartPos TreeSet has an entry for an event that is found on no haplotype. start pos = " + thisStart + ", next pos = " + nextStart);
                    }
                    if( isBiallelic ) {
                        final double R2 = calculateR2LD( x11, x12, x21, x22 );
                        if( DEBUG ) {
                            System.out.println("Found consecutive biallelic events with R^2 = " + String.format("%.4f", R2));
                            System.out.println("-- " + thisVC);
                            System.out.println("-- " + nextVC);
                        }
                        if( R2 > MERGE_EVENTS_R2_THRESHOLD ) {

                            final VariantContext mergedVC = createMergedVariantContext(thisVC, nextVC, ref, refLoc);

                            // remove the old event from the eventMap on every haplotype and the start pos key set, replace with merged event
                            for( final Haplotype h : haplotypes ) {
                                final HashMap<Integer, VariantContext> eventMap = h.getEventMap();
                                if( eventMap.containsKey(thisStart) && eventMap.containsKey(nextStart) ) {
                                    eventMap.remove(thisStart);
                                    eventMap.remove(nextStart);
                                    eventMap.put(mergedVC.getStart(), mergedVC);
                                }
                            }
                            startPosKeySet.add(mergedVC.getStart());
                            boolean containsStart = false;
                            boolean containsNext = false;
                            for( final Haplotype h : haplotypes ) {
                                final HashMap<Integer, VariantContext> eventMap = h.getEventMap();
                                if( eventMap.containsKey(thisStart) ) { containsStart = true; }
                                if( eventMap.containsKey(nextStart) ) { containsNext = true; }
                            }
                            if(!containsStart) { startPosKeySet.remove(thisStart); }
                            if(!containsNext) { startPosKeySet.remove(nextStart); }

                            if( DEBUG ) { System.out.println("====> " + mergedVC); }
                            mapWasUpdated = true;
                            break; // break out of tree set iteration since it was just updated, start over from the beginning and keep merging events
                        }
                    }
                }
                thisStart = nextStart;
            }
        }
    }

    // BUGBUG: make this merge function more general
    protected static VariantContext createMergedVariantContext( final VariantContext thisVC, final VariantContext nextVC, final byte[] ref, final GenomeLoc refLoc ) {
        final int thisStart = thisVC.getStart();
        final int nextStart = nextVC.getStart();
        byte[] refBases = ( thisVC.hasReferenceBaseForIndel() ? new byte[]{ thisVC.getReferenceBaseForIndel() } : new byte[]{} );
        byte[] altBases = ( thisVC.hasReferenceBaseForIndel() ? new byte[]{ thisVC.getReferenceBaseForIndel() } : new byte[]{} );
        refBases = ArrayUtils.addAll(refBases, thisVC.getReference().getBases());
        altBases = ArrayUtils.addAll(altBases, thisVC.getAlternateAllele(0).getBases());
        for( int locus = thisStart + refBases.length; locus < nextStart; locus++ ) {
            final byte refByte = ref[locus - refLoc.getStart()];
            refBases = ArrayUtils.add(refBases, refByte);
            altBases = ArrayUtils.add(altBases, refByte);
        }
        if( nextVC.hasReferenceBaseForIndel() ) {
            refBases = ArrayUtils.add(refBases, nextVC.getReferenceBaseForIndel());
            altBases = ArrayUtils.add(altBases, nextVC.getReferenceBaseForIndel());
        }
        refBases = ArrayUtils.addAll(refBases, nextVC.getReference().getBases());
        altBases = ArrayUtils.addAll(altBases, nextVC.getAlternateAllele(0).getBases());

        int iii = 0;
        if( refBases.length == altBases.length && VCFAlleleClipper.needsPadding(thisVC) ) { // special case of insertion + deletion of same length creates an MNP --> trim padding bases off the allele
            while( iii < refBases.length && refBases[iii] == altBases[iii] ) { iii++; }
        }
        final ArrayList<Allele> mergedAlleles = new ArrayList<Allele>();
        mergedAlleles.add( Allele.create( ArrayUtils.subarray(refBases, iii, refBases.length), true ) );
        mergedAlleles.add( Allele.create( ArrayUtils.subarray(altBases, iii, altBases.length), false ) );
        return new VariantContextBuilder("merged", thisVC.getChr(), thisVC.getStart() + iii, nextVC.getEnd(), mergedAlleles).make();
    }

    @Requires({"x11 >= 0", "x12 >= 0", "x21 >= 0", "x22 >= 0"})
    @Ensures({"result >= 0.0", "result <= 1.0"})
    protected static double calculateR2LD( final int x11, final int x12, final int x21, final int x22 ) {
        final int total = x11 + x12 + x21 + x22;
        final double pa1b1 = ((double) x11) / ((double) total);
        final double pa1b2 = ((double) x12) / ((double) total);
        final double pa2b1 = ((double) x21) / ((double) total);
        final double pa1 = pa1b1 + pa1b2;
        final double pb1 = pa1b1 + pa2b1;
        return ((pa1b1 - pa1*pb1) * (pa1b1 - pa1*pb1)) / ( pa1 * (1.0 - pa1) * pb1 * (1.0 - pb1) );
    }

    @Requires({"haplotypes.size() >= eventsAtThisLoc.size() + 1"})
    @Ensures({"result.size() == eventsAtThisLoc.size() + 1"})
    protected static ArrayList<ArrayList<Haplotype>> createAlleleMapper( final int loc, final ArrayList<VariantContext> eventsAtThisLoc, final ArrayList<Haplotype> haplotypes ) {
        final ArrayList<ArrayList<Haplotype>> alleleMapper = new ArrayList<ArrayList<Haplotype>>();
        final ArrayList<Haplotype> refList = new ArrayList<Haplotype>();
        for( final Haplotype h : haplotypes ) {
            if( h.getEventMap().get(loc) == null ) { // no event at this location so this is a reference-supporting haplotype
                refList.add(h);
            } else {
                boolean foundInEventList = false;
                for( final VariantContext vcAtThisLoc : eventsAtThisLoc ) {
                    if( h.getEventMap().get(loc).hasSameAllelesAs(vcAtThisLoc) ) {
                        foundInEventList = true;
                    }
                }
                if( !foundInEventList ) { // event at this location isn't one of the genotype-able options (during GGA) so this is a reference-supporting haplotype
                    refList.add(h);
                }
            }
        }
        alleleMapper.add(refList);
        for( final VariantContext vcAtThisLoc : eventsAtThisLoc ) {
            final ArrayList<Haplotype> list = new ArrayList<Haplotype>();
            for( final Haplotype h : haplotypes ) {
                if( h.getEventMap().get(loc) != null && h.getEventMap().get(loc).hasSameAllelesAs(vcAtThisLoc) ) {
                    list.add(h);
                }
            }
            alleleMapper.add(list);
        }
        return alleleMapper;
    }

    @Ensures({"result.size() == haplotypeAllelesForSample.size()"})
    protected static List<Allele> findEventAllelesInSample( final List<Allele> eventAlleles, final List<Allele> haplotypeAlleles, final List<Allele> haplotypeAllelesForSample, final ArrayList<ArrayList<Haplotype>> alleleMapper, final ArrayList<Haplotype> haplotypes ) {
        if( haplotypeAllelesForSample.contains(Allele.NO_CALL) ) { return noCall; }
        final ArrayList<Allele> eventAllelesForSample = new ArrayList<Allele>();
        for( final Allele a : haplotypeAllelesForSample ) {
            final Haplotype haplotype = haplotypes.get(haplotypeAlleles.indexOf(a));
            for( int iii = 0; iii < alleleMapper.size(); iii++ ) {
                final ArrayList<Haplotype> mappedHaplotypes = alleleMapper.get(iii);
                if( mappedHaplotypes.contains(haplotype) ) {
                    eventAllelesForSample.add(eventAlleles.get(iii));
                    break;
                }
            }
        }
        return eventAllelesForSample;
    }

    protected static boolean containsVCWithMatchingAlleles( final List<VariantContext> list, final VariantContext vcToTest ) {
        for( final VariantContext vc : list ) {
            if( vc.hasSameAllelesAs(vcToTest) ) {
                return true;
            }
        }
        return false;
    }

    protected static HashMap<Integer,VariantContext> generateVCsFromAlignment( final int alignmentStartHapwrtRef, final Cigar cigar, final byte[] ref, final byte[] alignment, final GenomeLoc refLoc, final String sourceNameToAdd, final int MNP_LOOK_AHEAD ) {
        return generateVCsFromAlignment(null, alignmentStartHapwrtRef, cigar, ref, alignment, refLoc, sourceNameToAdd, MNP_LOOK_AHEAD); // BUGBUG: needed for compatibility with HaplotypeResolver code
    }

    protected static HashMap<Integer,VariantContext> generateVCsFromAlignment( final Haplotype haplotype, final int alignmentStartHapwrtRef, final Cigar cigar, final byte[] ref, final byte[] alignment, final GenomeLoc refLoc, final String sourceNameToAdd, final int MNP_LOOK_AHEAD ) {
        final HashMap<Integer,VariantContext> vcs = new HashMap<Integer,VariantContext>();

        int refPos = alignmentStartHapwrtRef;
        if( refPos < 0 ) { return null; } // Protection against SW failures
        int alignmentPos = 0;

        for( final CigarElement ce : cigar.getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                    final byte[] insertionBases = Arrays.copyOfRange( alignment, alignmentPos, alignmentPos + elementLength );
                    boolean allN = true;
                    for( final byte b : insertionBases ) {
                        if( b != (byte) 'N' ) {
                            allN = false;
                            break;
                        }
                    }
                    if( !allN ) {
                        final ArrayList<Allele> insertionAlleles = new ArrayList<Allele>();
                        final int insertionStart = refLoc.getStart() + refPos - 1;
                        if( haplotype != null && (haplotype.leftBreakPoint + alignmentStartHapwrtRef + refLoc.getStart() - 1 == insertionStart + elementLength + 1 || haplotype.rightBreakPoint + alignmentStartHapwrtRef + refLoc.getStart() - 1 == insertionStart + elementLength + 1) ) {
                            insertionAlleles.add( Allele.create(ref[refPos-1], true) );
                            insertionAlleles.add( SYMBOLIC_UNASSEMBLED_EVENT_ALLELE );
                            vcs.put(insertionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), insertionStart, insertionStart, insertionAlleles).make());
                        } else {
                            insertionAlleles.add( Allele.create(Allele.NULL_ALLELE_STRING, true) );
                            insertionAlleles.add( Allele.create(insertionBases, false) );
                            vcs.put(insertionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), insertionStart, insertionStart, insertionAlleles).referenceBaseForIndel(ref[refPos-1]).make());
                        }

                    }
                    alignmentPos += elementLength;
                    break;
                case S:
                    alignmentPos += elementLength;
                    break;
                case D:
                    final byte[] deletionBases = Arrays.copyOfRange( ref, refPos, refPos + elementLength );
                    final ArrayList<Allele> deletionAlleles = new ArrayList<Allele>();
                    final int deletionStart = refLoc.getStart() + refPos - 1;
                    // BUGBUG: how often does this symbolic deletion allele case happen?
                    //if( haplotype != null && ( (haplotype.leftBreakPoint + alignmentStartHapwrtRef + refLoc.getStart() + elementLength - 1 >= deletionStart && haplotype.leftBreakPoint + alignmentStartHapwrtRef + refLoc.getStart() + elementLength - 1 < deletionStart + elementLength)
                    //        || (haplotype.rightBreakPoint + alignmentStartHapwrtRef + refLoc.getStart() + elementLength - 1 >= deletionStart && haplotype.rightBreakPoint + alignmentStartHapwrtRef + refLoc.getStart() + elementLength - 1 < deletionStart + elementLength) ) ) {
                    //    deletionAlleles.add( Allele.create(ref[refPos-1], true) );
                    //    deletionAlleles.add( SYMBOLIC_UNASSEMBLED_EVENT_ALLELE );
                    //    vcs.put(deletionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), deletionStart, deletionStart, deletionAlleles).make());
                    //} else {
                        deletionAlleles.add( Allele.create(deletionBases, true) );
                        deletionAlleles.add( Allele.create(Allele.NULL_ALLELE_STRING, false) );
                        vcs.put(deletionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), deletionStart, deletionStart + elementLength, deletionAlleles).referenceBaseForIndel(ref[refPos-1]).make());
                    //}
                    refPos += elementLength;
                    break;
                case M:
                    int numSinceMismatch = -1;
                    int stopOfMismatch = -1;
                    int startOfMismatch = -1;
                    int refPosStartOfMismatch = -1;
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        if( ref[refPos] != alignment[alignmentPos] && alignment[alignmentPos] != ((byte) 'N') ) {
                            // SNP or start of possible MNP
                            if( stopOfMismatch == -1 ) {
                                startOfMismatch = alignmentPos;
                                stopOfMismatch = alignmentPos;
                                numSinceMismatch = 0;
                                refPosStartOfMismatch = refPos;
                            } else {
                                stopOfMismatch = alignmentPos;
                            }
                        }
                        if( stopOfMismatch != -1) {
                            numSinceMismatch++;
                        }
                        if( numSinceMismatch > MNP_LOOK_AHEAD || (iii == elementLength - 1 && stopOfMismatch != -1) ) {
                            final byte[] refBases = Arrays.copyOfRange( ref, refPosStartOfMismatch, refPosStartOfMismatch + (stopOfMismatch - startOfMismatch) + 1 );
                            final byte[] mismatchBases = Arrays.copyOfRange( alignment, startOfMismatch, stopOfMismatch + 1 );
                            final ArrayList<Allele> snpAlleles = new ArrayList<Allele>();
                            snpAlleles.add( Allele.create( refBases, true ) );
                            snpAlleles.add( Allele.create( mismatchBases, false ) );
                            final int snpStart = refLoc.getStart() + refPosStartOfMismatch;
                            vcs.put(snpStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), snpStart, snpStart + (stopOfMismatch - startOfMismatch), snpAlleles).make());
                            numSinceMismatch = -1;
                            stopOfMismatch = -1;
                            startOfMismatch = -1;
                            refPosStartOfMismatch = -1;
                        }
                        refPos++;
                        alignmentPos++;
                    }
                    break;
                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator created during SW alignment: " + ce.getOperator() );
            }
        }
        return vcs;
    }
}