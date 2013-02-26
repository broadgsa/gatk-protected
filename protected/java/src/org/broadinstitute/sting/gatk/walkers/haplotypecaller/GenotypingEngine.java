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

package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.variant.variantcontext.*;

import java.io.PrintStream;
import java.util.*;

public class GenotypingEngine {

    private final boolean DEBUG;
    private final boolean USE_FILTERED_READ_MAP_FOR_ANNOTATIONS;
    private final static List<Allele> noCall = new ArrayList<Allele>(); // used to noCall all genotypes until the exact model is applied
    private final static Allele SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele.create("<UNASSEMBLED_EVENT>", false);
    private final VariantAnnotatorEngine annotationEngine;

    public GenotypingEngine( final boolean DEBUG, final VariantAnnotatorEngine annotationEngine, final boolean USE_FILTERED_READ_MAP_FOR_ANNOTATIONS ) {
        this.DEBUG = DEBUG;
        this.annotationEngine = annotationEngine;
        this.USE_FILTERED_READ_MAP_FOR_ANNOTATIONS = USE_FILTERED_READ_MAP_FOR_ANNOTATIONS;
        noCall.add(Allele.NO_CALL);
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * @param UG_engine                              UG Engine with basic input parameters
     * @param haplotypes                             Haplotypes to assign likelihoods to
     * @param samples                                Samples to genotype
     * @param haplotypeReadMap                       Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList
     * @param ref                                    Reference bytes at active region
     * @param refLoc                                 Corresponding active region genome location
     * @param activeRegionWindow                     Active window
     * @param genomeLocParser                        GenomeLocParser
     * @param activeAllelesToGenotype                Alleles to genotype
     * @return                                       List of VC's with genotyped events
     */
    @Requires({"refLoc.containsP(activeRegionWindow)", "haplotypes.size() > 0"})
    @Ensures("result != null")
    // TODO - can this be refactored? this is hard to follow!
    public List<VariantContext> assignGenotypeLikelihoods( final UnifiedGenotyperEngine UG_engine,
                                                           final List<Haplotype> haplotypes,
                                                           final List<String> samples,
                                                           final Map<String, PerReadAlleleLikelihoodMap> haplotypeReadMap,
                                                           final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList,
                                                           final byte[] ref,
                                                           final GenomeLoc refLoc,
                                                           final GenomeLoc activeRegionWindow,
                                                           final GenomeLocParser genomeLocParser,
                                                           final List<VariantContext> activeAllelesToGenotype ) {
        // sanity check input arguments
        if (UG_engine == null)
            throw new IllegalArgumentException("UG_Engine input can't be null, got "+UG_engine);
        if (haplotypes == null || haplotypes.isEmpty())
            throw new IllegalArgumentException("haplotypes input should be non-empty and non-null, got "+haplotypes);
        if (samples == null || samples.isEmpty())
            throw new IllegalArgumentException("samples input must be non-empty and non-null, got "+samples);
        if (haplotypeReadMap == null || haplotypeReadMap.isEmpty())
            throw new IllegalArgumentException("haplotypeReadMap input should be non-empty and non-null, got "+haplotypeReadMap);
        if (ref == null || ref.length == 0 )
            throw new IllegalArgumentException("ref bytes input should be non-empty and non-null, got "+ref);
        if (refLoc == null || refLoc.getStop()-refLoc.getStart()+1 != ref.length)
            throw new IllegalArgumentException(" refLoc must be non-null and length must match ref bytes, got "+refLoc);
        if (activeRegionWindow == null )
            throw new IllegalArgumentException("activeRegionWindow must be non-null, got "+activeRegionWindow);
        if (activeAllelesToGenotype == null )
            throw new IllegalArgumentException("activeAllelesToGenotype must be non-null, got "+activeAllelesToGenotype);
        if (genomeLocParser == null )
            throw new IllegalArgumentException("genomeLocParser must be non-null, got "+genomeLocParser);

        final List<VariantContext> returnCalls = new ArrayList<VariantContext>();
        final boolean in_GGA_mode = !activeAllelesToGenotype.isEmpty();

        // Using the cigar from each called haplotype figure out what events need to be written out in a VCF file
        final TreeSet<Integer> startPosKeySet = new TreeSet<Integer>();
        int count = 0;
        if( DEBUG ) { System.out.println("=== Best Haplotypes ==="); }
        for( final Haplotype h : haplotypes ) {
            // Walk along the alignment and turn any difference from the reference into an event
            h.setEventMap( generateVCsFromAlignment( h, h.getAlignmentStartHapwrtRef(), h.getCigar(), ref, h.getBases(), refLoc, "HC" + count++ ) );
            if( !in_GGA_mode ) { startPosKeySet.addAll(h.getEventMap().keySet()); }
            if( DEBUG ) {
                System.out.println( h.toString() );
                System.out.println( "> Cigar = " + h.getCigar() );
                System.out.println( ">> Events = " + h.getEventMap());
            }
        }

        cleanUpSymbolicUnassembledEvents( haplotypes );
        if( !in_GGA_mode && samples.size() >= 10 ) { // if not in GGA mode and have at least 10 samples try to create MNP and complex events by looking at LD structure
            mergeConsecutiveEventsBasedOnLD( haplotypes, samples, haplotypeReadMap, startPosKeySet, ref, refLoc );
            cleanUpSymbolicUnassembledEvents( haplotypes ); // the newly created merged events could be overlapping the unassembled events
        }
        if( in_GGA_mode ) {
            for( final VariantContext compVC : activeAllelesToGenotype ) {
                startPosKeySet.add( compVC.getStart() );
            }
        }

        // Walk along each position in the key set and create each event to be outputted
        for( final int loc : startPosKeySet ) {
            if( loc >= activeRegionWindow.getStart() && loc <= activeRegionWindow.getStop() ) { // genotyping an event inside this active region
                final List<VariantContext> eventsAtThisLoc = new ArrayList<VariantContext>(); // the overlapping events to merge into a common reference view
                final List<String> priorityList = new ArrayList<String>(); // used to merge overlapping events into common reference view

                if( !in_GGA_mode ) {
                    for( final Haplotype h : haplotypes ) {
                        final Map<Integer,VariantContext> eventMap = h.getEventMap();
                        final VariantContext vc = eventMap.get(loc);
                        if( vc != null && !containsVCWithMatchingAlleles(eventsAtThisLoc, vc) ) {
                            eventsAtThisLoc.add(vc);
                            priorityList.add(vc.getSource());
                        }
                    }
                } else { // we are in GGA mode!
                    int compCount = 0;
                    for( final VariantContext compVC : activeAllelesToGenotype ) {
                        if( compVC.getStart() == loc ) {
                            int alleleCount = 0;
                            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                                List<Allele> alleleSet = new ArrayList<Allele>(2);
                                alleleSet.add(compVC.getReference());
                                alleleSet.add(compAltAllele);
                                final String vcSourceName = "Comp" + compCount + "Allele" + alleleCount;
                                // check if this event is already in the list of events due to a repeat in the input alleles track
                                final VariantContext candidateEventToAdd = new VariantContextBuilder(compVC).alleles(alleleSet).source(vcSourceName).make();
                                boolean alreadyExists = false;
                                for( final VariantContext eventToTest : eventsAtThisLoc ) {
                                    if( eventToTest.hasSameAllelesAs(candidateEventToAdd) ) {
                                        alreadyExists = true;
                                    }
                                }
                                if( !alreadyExists ) {
                                    priorityList.add(vcSourceName);
                                    eventsAtThisLoc.add(candidateEventToAdd);
                                }
                                alleleCount++;
                            }
                        }
                        compCount++;
                    }
                }

                if( eventsAtThisLoc.isEmpty() ) { continue; }

                // Create the event mapping object which maps the original haplotype events to the events present at just this locus
                final Map<Event, List<Haplotype>> eventMapper = createEventMapper(loc, eventsAtThisLoc, haplotypes);

                // Sanity check the priority list for mistakes
                validatePriorityList( priorityList, eventsAtThisLoc );

                // Merge the event to find a common reference representation
                final VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(eventsAtThisLoc, priorityList, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);
                if( mergedVC == null ) { continue; }

                if( eventsAtThisLoc.size() != mergedVC.getAlternateAlleles().size() ) {
                    throw new ReviewedStingException("Record size mismatch! Something went wrong in the merging of alleles.");
                }
                final Map<VariantContext, Allele> mergeMap = new LinkedHashMap<VariantContext, Allele>();
                mergeMap.put(null, mergedVC.getReference()); // the reference event (null) --> the reference allele
                for(int iii = 0; iii < mergedVC.getAlternateAlleles().size(); iii++) {
                    mergeMap.put(eventsAtThisLoc.get(iii), mergedVC.getAlternateAllele(iii)); // BUGBUG: This is assuming that the order of alleles is the same as the priority list given to simpleMerge function
                }

                final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(mergeMap, eventMapper);

                if( DEBUG ) {
                    System.out.println("Genotyping event at " + loc + " with alleles = " + mergedVC.getAlleles());
                    //System.out.println("Event/haplotype allele mapping = " + alleleMapper);
                }

                final Map<String, PerReadAlleleLikelihoodMap> alleleReadMap = convertHaplotypeReadMapToAlleleReadMap( haplotypeReadMap, alleleMapper, UG_engine.getUAC().CONTAMINATION_FRACTION, UG_engine.getUAC().contaminationLog );

                final GenotypesContext genotypes = calculateGLsForThisEvent( samples, alleleReadMap, mergedVC );
                final VariantContext call = UG_engine.calculateGenotypes(new VariantContextBuilder(mergedVC).genotypes(genotypes).make(), UG_engine.getUAC().GLmodel);
                if( call != null ) {
                    final Map<String, PerReadAlleleLikelihoodMap> alleleReadMap_annotations = ( USE_FILTERED_READ_MAP_FOR_ANNOTATIONS ? alleleReadMap :
                            convertHaplotypeReadMapToAlleleReadMap( haplotypeReadMap, alleleMapper, 0.0, UG_engine.getUAC().contaminationLog ) );
                    final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap = filterToOnlyOverlappingReads( genomeLocParser, alleleReadMap_annotations, perSampleFilteredReadList, call );
                    VariantContext annotatedCall = annotationEngine.annotateContext(stratifiedReadMap, call);

                    if( annotatedCall.getAlleles().size() != mergedVC.getAlleles().size() ) { // some alleles were removed so reverseTrimming might be necessary!
                        annotatedCall = GATKVariantContextUtils.reverseTrimAlleles(annotatedCall);
                    }

                    returnCalls.add( annotatedCall );
                }
            }
        }
        return returnCalls;
    }

    /**
     * For a particular event described in inputVC, form PL vector for each sample by looking into allele read map and filling likelihood matrix for each allele
     * @param samples                List of samples to genotype
     * @param alleleReadMap          Allele map describing mapping from reads to alleles and corresponding likelihoods
     * @param mergedVC               Input VC with event to genotype
     * @return                       GenotypesContext object wrapping genotype objects with PLs
     */
    @Requires({"samples != null","alleleReadMap!= null", "mergedVC != null"})
    @Ensures("result != null")
    private GenotypesContext calculateGLsForThisEvent( final List<String> samples, final Map<String, PerReadAlleleLikelihoodMap> alleleReadMap, final VariantContext mergedVC ) {
        final GenotypesContext genotypes = GenotypesContext.create(samples.size());
        // Grab the genotype likelihoods from the appropriate places in the haplotype likelihood matrix -- calculation performed independently per sample
        for( final String sample : samples ) {
            final int numHaplotypes = mergedVC.getAlleles().size();
            final double[] genotypeLikelihoods = new double[numHaplotypes * (numHaplotypes+1) / 2];
            final double[][] haplotypeLikelihoodMatrix = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(sample, alleleReadMap, mergedVC.getAlleles());
            int glIndex = 0;
            for( int iii = 0; iii < numHaplotypes; iii++ ) {
                for( int jjj = 0; jjj <= iii; jjj++ ) {
                    genotypeLikelihoods[glIndex++] = haplotypeLikelihoodMatrix[iii][jjj]; // for example: AA,AB,BB,AC,BC,CC
                }
            }
            genotypes.add( new GenotypeBuilder(sample).alleles(noCall).PL(genotypeLikelihoods).make() );
        }
        return genotypes;
    }

    private void validatePriorityList( final List<String> priorityList, final List<VariantContext> eventsAtThisLoc ) {
        for( final VariantContext vc : eventsAtThisLoc ) {
            if( !priorityList.contains(vc.getSource()) ) {
                throw new ReviewedStingException("Event found on haplotype that wasn't added to priority list. Something went wrong in the merging of alleles.");
            }
        }
        for( final String name : priorityList ) {
            boolean found = false;
            for( final VariantContext vc : eventsAtThisLoc ) {
                if(vc.getSource().equals(name)) { found = true; break; }
            }
            if( !found ) {
                throw new ReviewedStingException("Event added to priority list but wasn't found on any haplotype. Something went wrong in the merging of alleles.");
            }
        }
    }

    private static Map<String, PerReadAlleleLikelihoodMap> filterToOnlyOverlappingReads( final GenomeLocParser parser,
                                                                                         final Map<String, PerReadAlleleLikelihoodMap> perSampleReadMap,
                                                                                         final Map<String, List<GATKSAMRecord>> perSampleFilteredReadList,
                                                                                         final VariantContext call ) {

        final Map<String, PerReadAlleleLikelihoodMap> returnMap = new LinkedHashMap<String, PerReadAlleleLikelihoodMap>();
        final GenomeLoc callLoc = parser.createGenomeLoc(call);
        for( final Map.Entry<String, PerReadAlleleLikelihoodMap> sample : perSampleReadMap.entrySet() ) {
            final PerReadAlleleLikelihoodMap likelihoodMap = new PerReadAlleleLikelihoodMap();

            for( final Map.Entry<GATKSAMRecord,Map<Allele,Double>> mapEntry : sample.getValue().getLikelihoodReadMap().entrySet() ) {
                // only count the read if it overlaps the event, otherwise it is not added to the output read list at all
                if( callLoc.overlapsP(parser.createGenomeLoc(mapEntry.getKey())) ) { // BUGBUG: This uses alignment start and stop, NOT soft start and soft end...
                    for( final Map.Entry<Allele,Double> alleleDoubleEntry : mapEntry.getValue().entrySet() ) {
                        likelihoodMap.add(mapEntry.getKey(), alleleDoubleEntry.getKey(), alleleDoubleEntry.getValue());
                    }
                }
            }

            // add all filtered reads to the NO_CALL list because they weren't given any likelihoods
            for( final GATKSAMRecord read : perSampleFilteredReadList.get(sample.getKey()) ) {
                // only count the read if it overlaps the event, otherwise it is not added to the output read list at all
                if( callLoc.overlapsP(parser.createGenomeLoc(read)) ) {
                    for( final Allele allele : call.getAlleles() ) {
                        likelihoodMap.add(read, allele, 0.0);
                    }
                }
            }

            returnMap.put(sample.getKey(), likelihoodMap);
        }
        return returnMap;
    }

    /**
     * Removes symbolic events from list of haplotypes
     * @param haplotypes       Input/output list of haplotypes, before/after removal
     */
    // TODO - split into input haplotypes and output haplotypes as not to share I/O arguments
    @Requires("haplotypes != null")
    protected static void cleanUpSymbolicUnassembledEvents( final List<Haplotype> haplotypes ) {
        final List<Haplotype> haplotypesToRemove = new ArrayList<Haplotype>();
        for( final Haplotype h : haplotypes ) {
            for( final VariantContext vc : h.getEventMap().values() ) {
                if( vc.isSymbolic() ) {
                    for( final Haplotype h2 : haplotypes ) {
                        for( final VariantContext vc2 : h2.getEventMap().values() ) {
                            if( vc.getStart() == vc2.getStart() && (vc2.isIndel() || vc2.isMNP()) ) { // unfortunately symbolic alleles can't currently be combined with non-point events
                                haplotypesToRemove.add(h);
                                break;
                            }
                        }
                    }
                }
            }
        }
        haplotypes.removeAll(haplotypesToRemove);
    }

    // BUGBUG: ugh, too complicated
    protected Map<String, PerReadAlleleLikelihoodMap> convertHaplotypeReadMapToAlleleReadMap( final Map<String, PerReadAlleleLikelihoodMap> haplotypeReadMap,
                                                                                              final Map<Allele, List<Haplotype>> alleleMapper,
                                                                                              final double downsamplingFraction,
                                                                                              final PrintStream downsamplingLog ) {

        final Map<String, PerReadAlleleLikelihoodMap> alleleReadMap = new LinkedHashMap<String, PerReadAlleleLikelihoodMap>();
        for( final Map.Entry<String, PerReadAlleleLikelihoodMap> haplotypeReadMapEntry : haplotypeReadMap.entrySet() ) { // for each sample
            final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
            for( final Map.Entry<Allele, List<Haplotype>> alleleMapperEntry : alleleMapper.entrySet() ) { // for each output allele
                final List<Haplotype> mappedHaplotypes = alleleMapperEntry.getValue();
                for( final Map.Entry<GATKSAMRecord, Map<Allele,Double>> readEntry : haplotypeReadMapEntry.getValue().getLikelihoodReadMap().entrySet() ) { // for each read
                    double maxLikelihood = Double.NEGATIVE_INFINITY;
                    for( final Map.Entry<Allele,Double> alleleDoubleEntry : readEntry.getValue().entrySet() ) { // for each input allele
                        if( mappedHaplotypes.contains( new Haplotype(alleleDoubleEntry.getKey())) ) { // exact match of haplotype base string
                            maxLikelihood = Math.max( maxLikelihood, alleleDoubleEntry.getValue() );
                        }
                    }
                    perReadAlleleLikelihoodMap.add(readEntry.getKey(), alleleMapperEntry.getKey(), maxLikelihood);
                }
            }
            perReadAlleleLikelihoodMap.performPerAlleleDownsampling(downsamplingFraction, downsamplingLog); // perform contamination downsampling
            alleleReadMap.put(haplotypeReadMapEntry.getKey(), perReadAlleleLikelihoodMap);
        }

        return alleleReadMap;
    }

    /**
     * TODO - comment me, clean me, refactor me!
     * @param haplotypes
     * @param samples
     * @param haplotypeReadMap
     * @param startPosKeySet
     * @param ref
     * @param refLoc
     */
    protected void mergeConsecutiveEventsBasedOnLD( final List<Haplotype> haplotypes,
                                                    final List<String> samples,
                                                    final Map<String, PerReadAlleleLikelihoodMap> haplotypeReadMap,
                                                    final TreeSet<Integer> startPosKeySet,
                                                    final byte[] ref,
                                                    final GenomeLoc refLoc ) {

        final int MAX_SIZE_TO_COMBINE = 15;
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
                    double x11 = Double.NEGATIVE_INFINITY;
                    double x12 = Double.NEGATIVE_INFINITY;
                    double x21 = Double.NEGATIVE_INFINITY;
                    double x22 = Double.NEGATIVE_INFINITY;

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
                        for( final String sample : samples ) {
                            final double haplotypeLikelihood = LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods( Collections.singleton(sample), haplotypeReadMap, Collections.singletonList(Allele.create(h, true)) )[0][0];
                            if( thisHapVC == null ) {
                                if( nextHapVC == null ) { x11 = MathUtils.approximateLog10SumLog10(x11, haplotypeLikelihood); }
                                else { x12 = MathUtils.approximateLog10SumLog10(x12, haplotypeLikelihood); }
                            } else {
                                if( nextHapVC == null ) { x21 = MathUtils.approximateLog10SumLog10(x21, haplotypeLikelihood); }
                                else { x22 = MathUtils.approximateLog10SumLog10(x22, haplotypeLikelihood); }
                            }
                        }
                    }
                    if( thisVC == null || nextVC == null ) {
                        continue;
                    }
                    if( isBiallelic ) {
                        final double R2 = calculateR2LD( Math.pow(10.0, x11), Math.pow(10.0, x12), Math.pow(10.0, x21), Math.pow(10.0, x22) );
                        if( DEBUG ) {
                            System.out.println("Found consecutive biallelic events with R^2 = " + String.format("%.4f", R2));
                            System.out.println("-- " + thisVC);
                            System.out.println("-- " + nextVC);
                        }
                        if( R2 > MERGE_EVENTS_R2_THRESHOLD ) {

                            final VariantContext mergedVC = createMergedVariantContext(thisVC, nextVC, ref, refLoc);

                            // remove the old event from the eventMap on every haplotype and the start pos key set, replace with merged event
                            for( final Haplotype h : haplotypes ) {
                                final Map<Integer, VariantContext> eventMap = h.getEventMap();
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
                                final Map<Integer, VariantContext> eventMap = h.getEventMap();
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
        byte[] refBases = new byte[]{};
        byte[] altBases = new byte[]{};
        refBases = ArrayUtils.addAll(refBases, thisVC.getReference().getBases());
        altBases = ArrayUtils.addAll(altBases, thisVC.getAlternateAllele(0).getBases());
        int locus;
        for( locus = thisStart + refBases.length; locus < nextStart; locus++ ) {
            final byte refByte = ref[locus - refLoc.getStart()];
            refBases = ArrayUtils.add(refBases, refByte);
            altBases = ArrayUtils.add(altBases, refByte);
        }
        refBases = ArrayUtils.addAll(refBases, ArrayUtils.subarray(nextVC.getReference().getBases(), locus > nextStart ? 1 : 0, nextVC.getReference().getBases().length)); // special case of deletion including the padding base of consecutive indel
        altBases = ArrayUtils.addAll(altBases, nextVC.getAlternateAllele(0).getBases());

        int iii = 0;
        if( refBases.length == altBases.length ) { // insertion + deletion of same length creates an MNP --> trim common prefix bases off the beginning of the allele
            while( iii < refBases.length && refBases[iii] == altBases[iii] ) { iii++; }
        }
        final List<Allele> mergedAlleles = new ArrayList<Allele>();
        mergedAlleles.add( Allele.create( ArrayUtils.subarray(refBases, iii, refBases.length), true ) );
        mergedAlleles.add( Allele.create( ArrayUtils.subarray(altBases, iii, altBases.length), false ) );
        return new VariantContextBuilder("merged", thisVC.getChr(), thisVC.getStart() + iii, nextVC.getEnd(), mergedAlleles).make();
    }

    protected static double calculateR2LD( final double x11, final double x12, final double x21, final double x22 ) {
        final double total = x11 + x12 + x21 + x22;
        final double pa1b1 = x11 / total;
        final double pa1b2 = x12 / total;
        final double pa2b1 = x21 / total;
        final double pa1 = pa1b1 + pa1b2;
        final double pb1 = pa1b1 + pa2b1;
        return ((pa1b1 - pa1*pb1) * (pa1b1 - pa1*pb1)) / ( pa1 * (1.0 - pa1) * pb1 * (1.0 - pb1) );
    }

    protected static Map<Allele, List<Haplotype>> createAlleleMapper( final Map<VariantContext, Allele> mergeMap, final Map<Event, List<Haplotype>> eventMap ) {
        final Map<Allele, List<Haplotype>> alleleMapper = new LinkedHashMap<Allele, List<Haplotype>>();
        for( final Map.Entry<VariantContext, Allele> entry : mergeMap.entrySet() ) {
            alleleMapper.put(entry.getValue(), eventMap.get(new Event(entry.getKey())));
        }
        return alleleMapper;
    }

    @Requires({"haplotypes.size() >= eventsAtThisLoc.size() + 1"})
    @Ensures({"result.size() == eventsAtThisLoc.size() + 1"})
    protected static Map<Event, List<Haplotype>> createEventMapper( final int loc, final List<VariantContext> eventsAtThisLoc, final List<Haplotype> haplotypes ) {

        final Map<Event, List<Haplotype>> eventMapper = new LinkedHashMap<Event, List<Haplotype>>(eventsAtThisLoc.size()+1);
        VariantContext refVC = eventsAtThisLoc.get(0); // the genome loc is the only safe thing to pull out of this VC because ref/alt pairs might change reference basis
        eventMapper.put(new Event(null), new ArrayList<Haplotype>());
        for( final VariantContext vc : eventsAtThisLoc ) {
            eventMapper.put(new Event(vc), new ArrayList<Haplotype>());
        }

        final List<Haplotype> undeterminedHaplotypes = new ArrayList<Haplotype>(haplotypes.size());
        for( final Haplotype h : haplotypes ) {
            if( h.isArtificialHaplotype() && loc == h.getArtificialAllelePosition() ) {
                final List<Allele> alleles = new ArrayList<Allele>(2);
                alleles.add(h.getArtificialRefAllele());
                alleles.add(h.getArtificialAltAllele());
                final Event artificialVC = new Event( (new VariantContextBuilder()).source("artificialHaplotype")
                                                                                            .alleles(alleles)
                                                                                            .loc(refVC.getChr(), refVC.getStart(), refVC.getStart() + h.getArtificialRefAllele().length() - 1).make() );
                if( eventMapper.containsKey(artificialVC) ) {
                    eventMapper.get(artificialVC).add(h);
                }
            } else if( h.getEventMap().get(loc) == null ) { // no event at this location so let's investigate later
                undeterminedHaplotypes.add(h);
            } else {
                boolean haplotypeIsDetermined = false;
                for( final VariantContext vcAtThisLoc : eventsAtThisLoc ) {
                    if( h.getEventMap().get(loc).hasSameAllelesAs(vcAtThisLoc) ) {
                        eventMapper.get(new Event(vcAtThisLoc)).add(h);
                        haplotypeIsDetermined = true;
                        break;
                    }
                }

                if( !haplotypeIsDetermined )
                    undeterminedHaplotypes.add(h);
            }
        }

        for( final Haplotype h : undeterminedHaplotypes ) {
            Event matchingEvent = new Event(null);
            for( final Map.Entry<Event, List<Haplotype>> eventToTest : eventMapper.entrySet() ) {
                // don't test against the reference allele
                if( eventToTest.getKey().equals(new Event(null)) )
                    continue;

                final Haplotype artificialHaplotype = eventToTest.getValue().get(0);
                if( isSubSetOf(artificialHaplotype.getEventMap(), h.getEventMap(), true) ) {
                    matchingEvent = eventToTest.getKey();
                    break;
                }
            }

            eventMapper.get(matchingEvent).add(h);
        }

        return eventMapper;
    }

    protected static boolean isSubSetOf(final Map<Integer, VariantContext> subset, final Map<Integer, VariantContext> superset, final boolean resolveSupersetToSubset) {

        for ( final Map.Entry<Integer, VariantContext> fromSubset : subset.entrySet() ) {
            final VariantContext fromSuperset = superset.get(fromSubset.getKey());
            if ( fromSuperset == null )
                return false;

            List<Allele> supersetAlleles = fromSuperset.getAlternateAlleles();
            if ( resolveSupersetToSubset )
                supersetAlleles = resolveAlternateAlleles(fromSubset.getValue().getReference(), fromSuperset.getReference(), supersetAlleles);

            if ( !supersetAlleles.contains(fromSubset.getValue().getAlternateAllele(0)) )
                return false;
        }

        return true;
    }

    private static List<Allele> resolveAlternateAlleles(final Allele targetReference, final Allele actualReference, final List<Allele> currentAlleles) {
        if ( targetReference.length() <= actualReference.length() )
            return currentAlleles;

        final List<Allele> newAlleles = new ArrayList<Allele>(currentAlleles.size());
        final byte[] extraBases = Arrays.copyOfRange(targetReference.getBases(), actualReference.length(), targetReference.length());
        for ( final Allele a : currentAlleles ) {
            newAlleles.add(Allele.extend(a, extraBases));
        }
        return newAlleles;
    }

    @Ensures({"result.size() == haplotypeAllelesForSample.size()"})
    protected static List<Allele> findEventAllelesInSample( final List<Allele> eventAlleles, final List<Allele> haplotypeAlleles, final List<Allele> haplotypeAllelesForSample, final List<List<Haplotype>> alleleMapper, final List<Haplotype> haplotypes ) {
        if( haplotypeAllelesForSample.contains(Allele.NO_CALL) ) { return noCall; }
        final List<Allele> eventAllelesForSample = new ArrayList<Allele>();
        for( final Allele a : haplotypeAllelesForSample ) {
            final Haplotype haplotype = haplotypes.get(haplotypeAlleles.indexOf(a));
            for( int iii = 0; iii < alleleMapper.size(); iii++ ) {
                final List<Haplotype> mappedHaplotypes = alleleMapper.get(iii);
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

    protected static Map<Integer,VariantContext> generateVCsFromAlignment( final Haplotype haplotype, final int alignmentStartHapwrtRef, final Cigar cigar, final byte[] ref, final byte[] alignment, final GenomeLoc refLoc, final String sourceNameToAdd ) {
        final Map<Integer,VariantContext> vcs = new LinkedHashMap<Integer,VariantContext>();

        int refPos = alignmentStartHapwrtRef;
        if( refPos < 0 ) { return null; } // Protection against SW failures
        int alignmentPos = 0;

        for( int cigarIndex = 0; cigarIndex < cigar.numCigarElements(); cigarIndex++ ) {
            final CigarElement ce = cigar.getCigarElement(cigarIndex);
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    final List<Allele> insertionAlleles = new ArrayList<Allele>();
                    final int insertionStart = refLoc.getStart() + refPos - 1;
                    final byte refByte = ref[refPos-1];
                    if( BaseUtils.isRegularBase(refByte) ) {
                        insertionAlleles.add( Allele.create(refByte, true) );
                    }
                    if( cigarIndex == 0 || cigarIndex == cigar.getCigarElements().size() - 1 ) { // if the insertion isn't completely resolved in the haplotype then make it a symbolic allele
                        insertionAlleles.add( SYMBOLIC_UNASSEMBLED_EVENT_ALLELE );
                    } else {
                        byte[] insertionBases = new byte[]{};
                        insertionBases = ArrayUtils.add(insertionBases, ref[refPos-1]); // add the padding base
                        insertionBases = ArrayUtils.addAll(insertionBases, Arrays.copyOfRange( alignment, alignmentPos, alignmentPos + elementLength ));
                        if( BaseUtils.isAllRegularBases(insertionBases) ) {
                            insertionAlleles.add( Allele.create(insertionBases, false) );
                        }
                    }
                    if( insertionAlleles.size() == 2 ) { // found a proper ref and alt allele
                        vcs.put(insertionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), insertionStart, insertionStart, insertionAlleles).make());
                    }
                    alignmentPos += elementLength;
                    break;
                }
                case S:
                {
                    alignmentPos += elementLength;
                    break;
                }
                case D:
                {
                    final byte[] deletionBases = Arrays.copyOfRange( ref, refPos - 1, refPos + elementLength );  // add padding base
                    final List<Allele> deletionAlleles = new ArrayList<Allele>();
                    final int deletionStart = refLoc.getStart() + refPos - 1;
                    final byte refByte = ref[refPos-1];
                    if( BaseUtils.isRegularBase(refByte) && BaseUtils.isAllRegularBases(deletionBases) ) {
                        deletionAlleles.add( Allele.create(deletionBases, true) );
                        deletionAlleles.add( Allele.create(refByte, false) );
                        vcs.put(deletionStart, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), deletionStart, deletionStart + elementLength, deletionAlleles).make());
                    }
                    refPos += elementLength;
                    break;
                }
                case M:
                case EQ:
                case X:
                {
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        final byte refByte = ref[refPos];
                        final byte altByte = alignment[alignmentPos];
                        if( refByte != altByte ) { // SNP!
                            if( BaseUtils.isRegularBase(refByte) && BaseUtils.isRegularBase(altByte) ) {
                                final List<Allele> snpAlleles = new ArrayList<Allele>();
                                snpAlleles.add( Allele.create( refByte, true ) );
                                snpAlleles.add( Allele.create( altByte, false ) );
                                vcs.put(refLoc.getStart() + refPos, new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), refLoc.getStart() + refPos, refLoc.getStart() + refPos, snpAlleles).make());
                            }
                        }
                        refPos++;
                        alignmentPos++;
                    }
                    break;
                }
                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator created during SW alignment: " + ce.getOperator() );
            }
        }
        return vcs;
    }

    private static class Event {
        public VariantContext vc;

        public Event( final VariantContext vc ) {
            this.vc = vc;
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof Event && ((((Event) obj).vc == null && vc == null) || (((Event) obj).vc != null && vc != null && ((Event) obj).vc.hasSameAllelesAs(vc))) ;
        }

        @Override
        public int hashCode() {
            return (vc == null ? -1 : vc.getAlleles().hashCode());
        }
    }
}