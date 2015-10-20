package org.broadinstitute.gatk.tools.walkers.phasing;

import org.broadinstitute.gatk.engine.refdata.*;
import org.broadinstitute.gatk.engine.contexts.*;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.samples.Gender;

import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.commandline.Tags;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.SampleUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;


import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFConstants;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;



import java.util.HashMap;
import java.util.ArrayList;
import java.lang.String;
import java.util.Collection;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.Map;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.lang.System;
import java.math.BigDecimal;
import java.math.RoundingMode;

import java.io.PrintStream;



/**
 * Computes the most likely genotype configuration and phases trios and parent/child pairs
 *
 * <p>
 * PhaseByTransmission is a GATK tool that 1) computes the most likely genotype configuration and phases trios and parent/child pairs given their genotype likelihoods and a mutation prior and 2) phases
 * all sites were parent/child transmission can be inferred unambiguously. It reports the genotype configuration (and hence phasing) probability. Additionally, if ReadBackPhasing was used prior to PhaseByTransmission, on the same vcf file (i.e.: HP tags are present in the genotype data), parental (haplotype) assignment of detected de novo mutations can be computed
 * Ambiguous sites are:
 * <ul>
 *     <li>Sites where all individuals are heterozygous</li>
 *     <li>Sites where there is a Mendelian violation</li>
 * </ul>
 *
 * <h2>Input</h2>
 * <p>
 * <ul>
 *     <li>A VCF variant set containing trio(s) and/or parent/child pair(s).</li>
 *     <li>A PED pedigree file containing the description of the individuals relationships and sex information if chromosome X is being called.</li>
 * </ul>
 * </p>
 *
 * <h2>Options</h2>
 * <p>
 *     <ul>
 *         <li>MendelianViolationsFile: An optional argument for reporting. If a file is specified, all sites that remain in mendelian violation after being assigned the most likely genotype
 *         configuration will be reported there. Information reported: chromosome, position, filter, allele count in VCF, family, transmission probability,
 *         and each individual genotype, depth, allelic depth and likelihoods.</li>
 *         <li>DeNovoPrior: Mutation prio; default is 1e-8</li>
 *     </ul>
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An VCF with genotypes recalibrated as most likely under the familial constraint and phased by descent where non ambiguous..
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T PhaseByTransmission \
 *   -V input.vcf \
 *   -ped input.ped \
 *   -o output.vcf
 * </pre>
 *
 */
//@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class} )
public class PhaseByTransmission extends RodWalker<HashMap<Byte,Double>, HashMap<Byte,Double>> 
{

	@ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(shortName = "mvf",required = false,fullName = "MendelianViolationsFile", doc="File to output the mendelian violation details.")
    private PrintStream mvFile = null;

    @Argument(shortName = "prior",required = false,fullName = "DeNovoPrior", doc="Prior for de novo mutations. Default: 1e-8")
    private double deNovoPrior=1e-8;
   
    //@Argument(shortName = "Xchr",required = false,fullName = "X-chromozome", doc="indicate whether the input variants are from the X chromozome")
    //private boolean Xchr = false;

    @Argument(shortName = "useAF",required = false,fullName = "UseAlleleFrequency", doc="Use of an AF-based prior on the founders. Accepted values: NO | AFDOS | GT | EXTERNAL. AFDOS - allele frequency dosage - computes the AF prior from the allele dosages of all the founder samples in the input; GT computes the AF prior from the best-guess genotypes of all founders found in input; default: NO - no allele frequency prior used")
    private String useAF = "NO";

    @Argument(shortName = "afFile", required = false, fullName = "alleleFrequencyFile", doc = "A vcf file containing the minor allele frequency in its INFO field, for each site")
    private RodBinding<VariantContext> afRod = null;// = new RodBinding<VariantContext>(VariantContext.class, "afRod", "afFile", "vcf", new Tags());
    
    @Argument(shortName = "afTag",required = false,fullName = "alleleFrequencyTag", doc="name of the allele frequency tag to be searched for in the externally supplied allele vcf file")
    private String af_tag = VCFConstants.ALLELE_FREQUENCY_KEY;
    
    @Argument(shortName = "af_cap",required = false,fullName = "alleleFrequencyCap", doc="Prevents the allele frequency prior from taking extreme values when the minor allele frequency is low and few samples are supplied; Very high AF prior can introduce false positive DNM calls. By default, the minor allele frequency is thus capped at 10^-3 i.e.: heterozigosity level. Any value can be speciffied in the (0,1) interval; if other values are speciffied the AF cap option is disabled")
    private double af_cap = 0.0001;
    
    @Argument (shortName = "assignPO", required = false, fullName = "assignParentalOrigin", doc="attempt to assign parental origin to detected mutations; PBT cannot phase de novos or HET/HET/HET genotype combinations so this argument is  only aplicable if the supplied input VCF has been run through ReadBackedPhasing, which fills in this gap; We can therefore combine trio-aware phasing with NGS data-driven phasing to identify the parental haplotype on which the de novo mutation occured")
    private int assignPO = 0;
    
    @Argument (shortName = "POET", required = false, fullName = "ParentalOriginEvidenceThreshold", doc = "This parameter sets the desired strictness when assigning a muation to one of the parental haplotypes; it represents the proportion of the total adjacent phasing sites that is to indicate for one of the parental haplotypes; default value is currently 1, i.e.: all adjacent phasing sites must be consistent with one of the parental haplotypes")
    private double POEvidenceThreshold = 1.0;
    
    //@Argument(shortName = "contrast",required = false,fullName = "UseContrastiveModel", doc="When set, EM is used to compute P(AF=0) and P(AF=X) and these 2 hypotheses are contrasted.")
    //private boolean contrast = false;

    @Hidden
    @Argument(shortName = "updatePLs",required = false,fullName = "UpdatePLs", doc="Updates the PLs to take the trio information (and AC information if useAF option is used).")
    private boolean updatePLs = false;

    @Hidden
    @Argument(shortName = "gt",required = false,fullName = "GenotypeAssignment", doc="Selects how the genotypes are assigned/updated. TRIOCONF => Most likely trio configuration; INDIV => Most likely individual genotype from trio-adjusted PLs. Can only be used in combination with updatePLs")
    private String gtAssignment = "TRIOCONF";
    
    @Hidden
    @Argument(shortName = "multipleOffspring",required = false,fullName = "multipleOffspringAsTrios", doc=" if a family has multiple offspring, build an independent trio for each one? - in case a family has multiple offspring, we may evaluate the joint probability of parents' genotypes and each individual offspring, as opposed to the joint probability of all family genotypes. This argument implements the former option. The PBT model is consistent and correct in evaluating the probability of a DNM in each parents-offspring trios, but, as a result of the same parents' genotypes being evaluated under different trio combinations (for the different offspring), the parents' genotypes in the output vcf will be the most likely genotypes from the last trio configuration evaluated (i.e.: not the most likely genotypes given all the available information - all offspring). In these unlikely cases, the DNMs are correctly represented in the mvFile, but may not be in the outputted vcf")
    private boolean multiple_offspring = false;

    //@Argument(shortName = "EMtest",required = false,fullName = "EMtest", doc="TRUE => AF=0 AND AF=x cases are considered in case of an MV; FALSE => Only AF=0 is considered.")
    //private boolean emTest = false;

    @Argument(shortName = "fatherAlleleFirst",required = false,fullName = "FatherAlleleFirst", doc="Ouputs the father allele as the first allele in phased child genotype. i.e. father|mother rather than mother|father.")
    private boolean fatherAlleleFirst=false;
    
    @Argument(shortName="hapIntervals", required=false, fullName="haploidMaleIntervals", doc="Haploid regions of the genome in males. Male diploid genotypes and likelihoods in these regions will be converted to haploid values. Overrides the site_male_ploidy argument")
	private String x_file = null; // locations to be considerred haploid for males
    
    @Argument(shortName="fh", required=false, fullName="fake_hets", doc="output dummy, very big PL values for het genotypes as well -- relevant only for X-like inheritance pattern")
	private int fake_hets = 0;
    
    @Hidden
    @Argument(shortName="allX", required=false, fullName="allInputIsX", doc="consider all input positions to be haploid for males? (i.e.: X-chr)")
	private boolean allX = false; // locations to be considerred
    
    @Hidden
    @Argument(shortName="noTP", required=false, fullName="noTransmissionPosterior", doc="do not output the transmission posterior for a the genotyped trio configuration (i.e.: if DNM identification is not the main purpose)")
	private boolean noTP = false; // locations to be considerred
    
    @Hidden
    @Argument(shortName = "ufad",required = false,fullName = "useFatherAsNeeded", doc="does not use the father in the likelihood computation if not needed; i.e.: for calling X chromozome in families with male children; useful for dnms")
    private boolean ufad=true;
    
    @Hidden
    @Argument (shortName = "op", required = false, fullName="onlyPhase", doc = "only phase found trios; i.e.: do not evaluate genotype calls")
    private boolean op = false;

	@Hidden
	@Argument (shortName = "pep", required = false, fullName="phaseEverythingPossible", doc = "attempt to trio-phase sites that are normally ignored by PBT such as multi-allelic sites")
	private boolean phaseEverything = false;

    @Output
    protected VariantContextWriter vcfWriter = null;

	private int site_male_ploidy = 2;
	private int no_af_found = 0;
	private int no_af_dnm_call = 0;
	
	//-----Parameter constants - TODO: turn them to enums
    private final String NO_AF = "NO";
    private final String AFDOS =  "AFDOS";
    private final String GT = "GT";
    private final String EXTERNAL = "EXTERNAL";

    private final String TRIOCONF = "TRIOCONF";
    private final String INDIV = "INDIV";
	
	private final String TRANSMISSION_PROBABILITY_TAG_NAME = "TP";
	private final String RBP_TAG = "HP";
	private final String RBP_CANCELLED = "PI";
	private final String SOURCE_NAME = "PhaseByTransmission";
	private final int RBP_CANCELLED_VAL = 1;
	//-----Parameter constants
	
	public final short NO_TRANSMISSION_PROB = -1;
	
	private final Allele dummy_ref  = Allele.create("A", true);
	private final Allele dummy_alt  = Allele.create("C", false);
	//private final String dummy_name = new String("DUMMY_SAMPLE");
	
	//Cache some log values
    private final double log10_1_3 = Math.log10(1.0/3.0);
    private final double log10_1_2 = Math.log10(1.0/2.0);
    private final double log10_2 = Math.log10(2.0);
	
	
	private ArrayList<Sample> trios = new ArrayList<Sample>();
    private HashSet<String> founderIds = new HashSet<String>();
    
    private ArrayList<GenomeLoc> haploid_intervals = null;
    private ArrayList<String> haploidIndiv = null;
    
    // hash adjacent variants as we traverse the input vcf in order to assign parental origin of DNMs, where possible
    private ArrayList<VariantContext> variants_q = null;
    
    
    private EnumMap<GenotypeType,Double> neutralAF = new EnumMap<GenotypeType, Double>(GenotypeType.class);
    
    private EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap< GenotypeType, ArrayList<TrioPhase> >>> transmissionMatrix;
	
    private HashMap<String,TrioRbpContext> rbp_queue = new HashMap<String,TrioRbpContext>();
    
    //Random number generator
    private Random rand = new Random(1);
    
    
    //Metrics counters hash keys
    private final Byte NUM_TRIO_GENOTYPES_CALLED = 0;
    private final Byte NUM_TRIO_GENOTYPES_NOCALL = 1;
    private final Byte NUM_TRIO_GENOTYPES_PHASED = 2;
    private final Byte NUM_TRIO_HET_HET_HET = 3;
    private final Byte NUM_TRIO_VIOLATIONS = 4;
    private final Byte NUM_TRIO_DOUBLE_VIOLATIONS = 10;
    private final Byte NUM_PAIR_GENOTYPES_CALLED = 5;
    private final Byte NUM_PAIR_GENOTYPES_NOCALL = 6;
    private final Byte NUM_PAIR_GENOTYPES_PHASED = 7;
    private final Byte NUM_PAIR_HET_HET = 8;
    private final Byte NUM_PAIR_VIOLATIONS = 9;
    private final Byte NUM_GENOTYPES_MODIFIED = 11;
    private final Byte NUM_MUTATIONS = 12;
    private final Byte NUM_INCOMPLETE_TRIO_PHASED = 13;
    private final Byte NUM_INCOMPLETE_TRIO_VIOLATIONS = 14;
    private final Byte NUM_TRIO_SINGLES_PHASED = 15;
    private final Byte NUM_PAIR_SINGLES_PHASED = 16;
	
	
	private enum FamilyMember
	{
		MOTHER,FATHER,CHILD;
	}
	
	/**
	 * Minimalist datastructure to hold the cached, generic information about genotypes from trio-combinations
	 */
	private class DummyGT
	{
		protected ArrayList<Allele> alleles;
		protected boolean phased = false;
		
		public DummyGT(ArrayList<Allele> alleles, boolean phased)
		{
			this.alleles = alleles;
			this.phased = phased;
		}
		
	}
	
	/**
	 * Minimalist data structure to keep track of meta-data needed to assign parental origin to DNMs (requires RBP be run before PBT)
	 */
	private class TrioRbpContext
	{
		Sample child = null;
		ArrayList<Set<String>> filters = new ArrayList<Set<String>>();
		ArrayList<String> ADs = new ArrayList<String>();
		ArrayList<GenomeLoc> gtLoci = new ArrayList<GenomeLoc> ();
		ArrayList<TrioConfiguration> trioGTs = new ArrayList<TrioConfiguration>();
		// for each locus: 0=no missing parent; 1=father genotype is missing; 2=mother genotype is missing
		ArrayList<Integer> missingParents = new ArrayList<Integer>();
		//ArrayList<Genotype> childGTs = new ArrayList<Genotype>();
		long phaseRefPos = -1;
		boolean hasDnms = false;
	}
	
	/**
	 *	class to build, store and access cached information about a trio-genotype configuration. It is ploidy aware. currently only ploidy 1 or 2; X-like ploidy
	 */
	private class TrioPhase
	{

		
		private Gender childSex; // only relevant for haploid sites
		private int mvCount = 0;
		private double transmissionProb = 0; // LOG-space
		
		
		private EnumMap<FamilyMember,DummyGT> trioPhasedGenotypes;
		
		// method to populate current configuration's genotypes; perform any possible phasing given one single genotype
		private void addSinglePhasedGenotype(FamilyMember member, GenotypeType g_type, Gender sex)
		{
			if (site_male_ploidy == 1 && sex == Gender.MALE)
				trioPhasedGenotypes.put( member, 
						new DummyGT(getAlleles(g_type,sex), (g_type == GenotypeType.HOM_REF || g_type == GenotypeType.HOM_VAR)) );
			else
				trioPhasedGenotypes.put( member, 
						new DummyGT(getAlleles(g_type,sex), (g_type == GenotypeType.HOM_REF || g_type == GenotypeType.HOM_VAR || g_type == GenotypeType.HET)) );
		}
		
		// method to populate current configuration's genotypes; perform any possible phasing given the genotype of the child and that of one parent
		private void addPairPhasedGenotypes(FamilyMember parent, GenotypeType parent_g_type, GenotypeType child_g_type)
		{
			// HET(parent) HET(child) cannot be phased
			if (parent_g_type == GenotypeType.HET && child_g_type == GenotypeType.HET)
				{
					trioPhasedGenotypes.put(parent, 			new DummyGT(getAlleles(parent_g_type, (parent == FamilyMember.FATHER)? Gender.MALE:Gender.FEMALE), false) );
					trioPhasedGenotypes.put(FamilyMember.CHILD, new DummyGT(getAlleles(child_g_type,childSex), false) );
					return;
				}
			
			ArrayList<Allele> parentAlleles = getAlleles(parent_g_type, (parent == FamilyMember.FATHER)? Gender.MALE:Gender.FEMALE);
			ArrayList<Allele> childAlleles  = getAlleles(child_g_type, childSex);
			
			// check returned alleles to see if we have what to phase
			if (childAlleles == null)
			{
				if (parentAlleles == null)
					return;
				else
					addSinglePhasedGenotype(parent,parent_g_type, (parent == FamilyMember.FATHER)? Gender.MALE:Gender.FEMALE);
				return;
			}
			if (parentAlleles == null)
			{
				addSinglePhasedGenotype(FamilyMember.CHILD, child_g_type, childSex);
				return;
			}
			
			//phase
			int found_transmitted_allele = childAlleles.indexOf(parentAlleles.get(0));
			if (found_transmitted_allele > -1)
			{
				childAlleles.add(0,childAlleles.remove(found_transmitted_allele));
				trioPhasedGenotypes.put(parent, 			new DummyGT(parentAlleles,true));
				trioPhasedGenotypes.put(FamilyMember.CHILD, new DummyGT(childAlleles,true));
			}
			else
				if (parentAlleles.size() > 1)
					if ( (found_transmitted_allele = childAlleles.indexOf(parentAlleles.get(1))) > -1 )
					{
						parentAlleles.add(0, parentAlleles.remove(1));
						childAlleles.add(0, childAlleles.remove(found_transmitted_allele));
						trioPhasedGenotypes.put(parent, 			new DummyGT(parentAlleles,true));
						trioPhasedGenotypes.put(FamilyMember.CHILD, new DummyGT(childAlleles,true));
					}
					else // Mendelian Violation (diploid parent) => no phasing possible in pair
					{
						trioPhasedGenotypes.put(parent, 			new DummyGT(parentAlleles, false));
						trioPhasedGenotypes.put(FamilyMember.CHILD, new DummyGT(childAlleles, false));
					}
				else // Mendelian Violation (haploid parent) => no phasing possible in pair
				{
					trioPhasedGenotypes.put(parent, 			new DummyGT(parentAlleles, false));
					trioPhasedGenotypes.put(FamilyMember.CHILD, new DummyGT(childAlleles, false));
				}	
		}
		
		// puts transmitted allele first in parents; perform any phasing possible given all the genotypes of a trio-family
		private void addFamilyPhasedGenotypes(GenotypeType mother_g_type, GenotypeType father_g_type, GenotypeType child_g_type)
		{
			ArrayList<Allele> fatherAlleles = getAlleles(father_g_type, Gender.MALE);
			ArrayList<Allele> motherAlleles = getAlleles(mother_g_type, Gender.FEMALE);
			ArrayList<Allele> childAlleles  = getAlleles(child_g_type,  childSex);
			
			
			if (childAlleles == null || fatherAlleles == null || motherAlleles == null) 
				return;
			
			if (mother_g_type == GenotypeType.HET && father_g_type == GenotypeType.HET && child_g_type == GenotypeType.HET)
			{
				trioPhasedGenotypes.put(FamilyMember.MOTHER, new DummyGT(motherAlleles, false));
				trioPhasedGenotypes.put(FamilyMember.FATHER, new DummyGT(fatherAlleles, false));
				trioPhasedGenotypes.put(FamilyMember.CHILD,  new DummyGT(childAlleles, 	false));
				return;
			}
			
			int child_first_allele_idx = -1;
			int child_second_allele_idx = -1;
			
			// as of currently, the ONLY reason why the child would have only one allele at this point is because we are in a male-haploid region with X-like inheritance pattern
			if (childAlleles.size() == 1)
			{
				child_second_allele_idx = motherAlleles.indexOf(childAlleles.get(0));
				if (child_second_allele_idx > -1)
				{
					motherAlleles.add(0,motherAlleles.remove(child_second_allele_idx));
					trioPhasedGenotypes.put(FamilyMember.MOTHER, new DummyGT(motherAlleles, true));
					trioPhasedGenotypes.put(FamilyMember.FATHER, new DummyGT(fatherAlleles, true));
					trioPhasedGenotypes.put(FamilyMember.CHILD,  new DummyGT(childAlleles, 	true));
				}
				else // Mendelian Violation in X-like inheritance pattern => no phasing
				{
					trioPhasedGenotypes.put(FamilyMember.MOTHER, new DummyGT(motherAlleles, false));
					trioPhasedGenotypes.put(FamilyMember.FATHER, new DummyGT(fatherAlleles, false));
					trioPhasedGenotypes.put(FamilyMember.CHILD,  new DummyGT(childAlleles, 	false));
				}
			}
			else // child is diploid
			{
				
				child_first_allele_idx  = fatherAlleles.indexOf(childAlleles.get(0));
				child_second_allele_idx = motherAlleles.indexOf(childAlleles.get(1));
				// child's alleles can be either first from father - second from mother, or the other way around
				if (child_first_allele_idx > -1 && child_second_allele_idx > -1)
				{
					if (fatherAlleleFirst == false)
					{
						childAlleles.add(0, childAlleles.remove(1));
					}
					fatherAlleles.add(0,fatherAlleles.remove(child_first_allele_idx));
					motherAlleles.add(0,motherAlleles.remove(child_second_allele_idx));
					trioPhasedGenotypes.put(FamilyMember.MOTHER, new DummyGT(motherAlleles, true));
					trioPhasedGenotypes.put(FamilyMember.FATHER, new DummyGT(fatherAlleles, true));
					trioPhasedGenotypes.put(FamilyMember.CHILD,  new DummyGT(childAlleles, 	true));
				}
				else // check the other way around
				{
					child_first_allele_idx  = fatherAlleles.indexOf(childAlleles.get(1));
					child_second_allele_idx = motherAlleles.indexOf(childAlleles.get(0));
					// child's returned alleles can be either first from father - second from mother, or the other way around (i.e.: consistent with Mendelian transmission)
					if (child_first_allele_idx > -1 && child_second_allele_idx > -1)
					{
						if (fatherAlleleFirst == true)
						{
							childAlleles.add(0, childAlleles.remove(1));
						}
						fatherAlleles.add(0, fatherAlleles.remove(child_first_allele_idx));
						motherAlleles.add(0, motherAlleles.remove(child_second_allele_idx));
						trioPhasedGenotypes.put(FamilyMember.MOTHER, new DummyGT(motherAlleles, true));
						trioPhasedGenotypes.put(FamilyMember.FATHER, new DummyGT(fatherAlleles, true));
						trioPhasedGenotypes.put(FamilyMember.CHILD,  new DummyGT(childAlleles, 	true));
					}
					else // Mendelian Violation in diploid child => no phasing
					{
						trioPhasedGenotypes.put(FamilyMember.MOTHER, new DummyGT(motherAlleles, false));
						trioPhasedGenotypes.put(FamilyMember.FATHER, new DummyGT(fatherAlleles, false));
						trioPhasedGenotypes.put(FamilyMember.CHILD,  new DummyGT(childAlleles, 	false));
					}
				}
			}
		
		}
		
		// compute trio configuration prior for Mendelian consistent trio-genotypes; only for bi-allelic sites
		private double getNonDNMCombinationPrior(GenotypeType mother, GenotypeType father, GenotypeType child)
        {
            //If this prior is not used, then just return 1
            //if(nomtp) return 1.0;          
            if (site_male_ploidy == 2)
            {
            	if(mother == GenotypeType.HET)
            	{
            		if(father == GenotypeType.HET && !(child == GenotypeType.HET))
            		{
            			return 0.25;
            		}
            		return 0.5;
            	}
            	else 
            		if(father == GenotypeType.HET)
            		{
            			return 0.5;
            		}
            }
            else
            {
            	if (childSex == Gender.MALE)
            		if (mother == GenotypeType.HET)
            			return 0.5;
            		else
            			return 1;
            	else
            		if (childSex == Gender.FEMALE)
            			if (mother == GenotypeType.HET)
            				return 0.75;
            			else
            				return 1;
            }
            
            return 1.0;
        }
		
		/**
		 * Constructor for the class; builds and stores trio-combination information and phasing
		 * @param mother the symbolic type of the mother's genotype
		 * @param father the symbolic type of the father's genotype
		 * @param child the symbolic type of the child's genotype
		 * @param deNovoPrior mutation prior
		 * @param childSex the sex of the offspring; relevant for haploid regions
		 */
		public TrioPhase(GenotypeType mother, GenotypeType father, GenotypeType child, double deNovoPrior, Gender childSex)
		{
			this.childSex = childSex;
			this.mvCount = getConfigurationMVCount(mother,father,child,childSex);
			trioPhasedGenotypes = new EnumMap<FamilyMember,DummyGT>(FamilyMember.class);
			// set transmission prior -- LOG-space
			switch (mvCount)
			{
			case 0: 
				transmissionProb = Math.log10(getNonDNMCombinationPrior(mother,father,child));
				break;
			case 1: 
				transmissionProb = Math.log10(deNovoPrior);
				break;
			case 2:
				if (childSex == Gender.MALE && site_male_ploidy == 1)
					throw new UserException.BadArgumentValue("mvCount", "---> GTs probably malformed, found 2 DNMs in male offspring on the X chromosome (at some locus)");
				transmissionProb = 2 * Math.log10(deNovoPrior);
				break;
			default:
				throw new UserException.BadArgumentValue("mvCount", "---> GTs probably malformed, found an inconsistent number of mendelian violations");
			
			}
			
			// produce the possible phasing for this combination
			if (!isPhasable(child, childSex))
			{
				addSinglePhasedGenotype(FamilyMember.CHILD, child, childSex);
				addSinglePhasedGenotype(FamilyMember.MOTHER, mother, Gender.FEMALE);
				addSinglePhasedGenotype(FamilyMember.FATHER, father, Gender.MALE);
				
			}
			else
			{
				if (!isPhasable(mother, Gender.FEMALE))
				{
					addSinglePhasedGenotype(FamilyMember.MOTHER, mother, Gender.FEMALE);
					if (!isPhasable(father, Gender.MALE))
					{
						addSinglePhasedGenotype(FamilyMember.CHILD, child, childSex);
						addSinglePhasedGenotype(FamilyMember.FATHER, father, Gender.MALE);
					}
					else
						addPairPhasedGenotypes(FamilyMember.FATHER, father, child);
				}
				else
				{
					if (!isPhasable(father, Gender.MALE))
					{
						addSinglePhasedGenotype(FamilyMember.FATHER, father, Gender.MALE);
						addPairPhasedGenotypes(FamilyMember.MOTHER, mother, child);
					}
					else
					{
						//if (mother == GenotypeType.HET && father == GenotypeType.HET && child == GenotypeType.HET)
						//{
						//	addSinglePhasedGenotype(FamilyMember.CHILD, child, childSex);
						//	addSinglePhasedGenotype(FamilyMember.MOTHER, mother, Gender.FEMALE);
						//	addSinglePhasedGenotype(FamilyMember.FATHER, father, Gender.MALE);
						//}
						//else
							addFamilyPhasedGenotypes(mother,father,child);
					}
				}
			}
			
		
		}
		
		public DummyGT getDummyGT(FamilyMember dude){ return trioPhasedGenotypes.get(dude); }		
		public int getMVCount() { return mvCount; }
		public double getTransmissionProb() { return transmissionProb; }
		
	}
	
	//Create the transmission matrices for ploidy 1 AND 2
    private void buildMatrices(double deNovoPrior)
    {
        transmissionMatrix = new EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap< GenotypeType, ArrayList<TrioPhase> >>>(GenotypeType.class);
        
        	int _site_male_ploidy = site_male_ploidy;
        	for(GenotypeType mother : GenotypeType.values())
        	{
           		transmissionMatrix.put(mother,new EnumMap<GenotypeType,EnumMap<GenotypeType,ArrayList<TrioPhase>>>(GenotypeType.class));
           		for(GenotypeType father : GenotypeType.values())
           		{
               		transmissionMatrix.get(mother).put(father,new EnumMap<GenotypeType,ArrayList<TrioPhase>>(GenotypeType.class));
               		for(GenotypeType child : GenotypeType.values())
               		{
                   		ArrayList<TrioPhase> sexDiff = new ArrayList<TrioPhase>(3);
                   		site_male_ploidy = 1;
                   		sexDiff.add(new TrioPhase(mother,father,child,deNovoPrior,Gender.FEMALE));
                   		sexDiff.add(new TrioPhase(mother,father,child,deNovoPrior,Gender.MALE));
                   		
                   		site_male_ploidy = 2;
                   		sexDiff.add(new TrioPhase(mother,father,child,deNovoPrior,Gender.FEMALE));
               			transmissionMatrix.get(mother).get(father).put(child,sexDiff);
               		}
           		}
        	}
        	
        	site_male_ploidy = _site_male_ploidy;
    }

    // parse input families from ped file
    private void setTrios()
    {

        Map<String,Set<Sample>> families = this.getSampleDB().getFamilies();
        Set<Sample> family;
        ArrayList<Sample> parents;
        
        for(String familyID : families.keySet())
        {
            family = families.get(familyID);
            
            if(family.size() <= 2 || (family.size()>3 && !multiple_offspring) )
            {
                logger.info(String.format("Caution: Family %s has %d members; At the moment Phase By Transmission only supports trios and parent/child pairs. Family skipped.",familyID,family.size()));
                continue;
            }
            else
            {
                for(Sample familyMember : family)
                {
                    parents = familyMember.getParents();
                    if(parents.size()>0)
                    {
                   
                        if(family.containsAll(parents))
                        {
                            this.trios.add(familyMember);
                            this.founderIds.add(familyMember.getMaternalID());
                            this.founderIds.add(familyMember.getPaternalID());
                            this.rbp_queue.put(familyMember.getID(), null);
                        }
                        else
                        {
                            if (family.contains(familyMember.getFather()))
                            {
                            	this.trios.add(familyMember);
                                this.founderIds.add(familyMember.getPaternalID());
                            }
                            else 
                            	if (family.contains(familyMember.getMother()))
                                {
                                	this.trios.add(familyMember);
                                    this.founderIds.add(familyMember.getMaternalID());
                                } 
                            	else 
                            		logger.info(String.format("Caution: Family %s skipped as it is not a trio nor a parent/child pair; At the moment Phase By Transmission only supports trios and parent/child pairs. Family skipped.",familyID));
                          
                        }
                    }
                }
            }
        }
    }

	/**
     * Parse the familial relationship specification, build the transmission matrices and initialize VCF writer
     */
    public void initialize() 
	{	
    	if (af_cap <= 0.0 || af_cap >= 1)
    		af_cap = -1;
    	
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());
        if (afRod != null)
        	rodNames.add(afRod.getName());
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        //Get the trios from the families passed as ped
        setTrios();
        if(trios.size() < 1)
            throw new UserException.BadInput("No PED file passed or no trios found in PED file. Aborted.");

        //check for supplied X-like inheritance intervals and select samples that may require genotype adjustments
        if ( x_file != null && !x_file.isEmpty() )
        	haploid_intervals = HaploidWriter.readLocations(x_file, this.getMasterSequenceDictionary());
        if (haploid_intervals != null || allX)
        	haploidIndiv = HaploidWriter.subsetMaleSamples(this.getSampleDB().getSamples());
        // in case ALL input variants are to be considerred from an X-like inheritance pattern
        if (haploid_intervals == null && allX)
        	site_male_ploidy = 1;

        //Check arguments passed
        if( !(useAF.equals(AFDOS) || useAF.equals(NO_AF) || useAF.equals(GT) || useAF.equals(EXTERNAL)) )
            throw new UserException.BadArgumentValue("useAlleleFrequency",String.format("Value '%s' is not an acceptable value for argument useAlleleFrequency. Possible values are: NO | AFDOS | GT | EXTERNAL",useAF));

        if( !(gtAssignment.equals(TRIOCONF) || gtAssignment.equals(INDIV)) )
            throw new UserException.BadArgumentValue("GenotypeAssignment",String.format("Value '%s' is not an acceptable value for argument GenotypeAssignment. Possible values are: TRIOCONF | INDIV",gtAssignment));
        else if(gtAssignment.equals(INDIV) && !updatePLs)
            throw new UserException.BadArgumentValue("GenotypeAssignment","Value 'INDIV' is only possible in conjunction with the 'updatePLs' option");

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(GATKVCFUtils.getHeaderFields(this.getToolkit()));
        if (!noTP)
        	headerLines.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred score of the posterior probability of the genotype configuration and phased."));
        headerLines.add(new VCFFormatHeaderLine(RBP_CANCELLED, 1, VCFHeaderLineType.Integer, "Tag signaling whether existing RBP information (in HP tags) has been made inconsistent by applying PBT"));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));

        buildMatrices(deNovoPrior);

        //Initialize the neutral AF matrix
        for(GenotypeType genotype : GenotypeType.values())
        {
            neutralAF.put(genotype,0.0);
        }
        
        // header of mvFile
        if(mvFile != null)
		{
            mvFile.println("#CHROM\tPOS\tFILTER\tAC\tTP\t#MV\tFAMILY\tDNM_ORIGIN\tPAT_HAPS\tMAT_HAPS\tMOTHER_ID\tMOTHER_GT\tMOTHER_AL\tMOTHER_DP\tMOTHER_AD\tMOTHER_PL\tFATHER_ID\tFATHER_GT\tFATHER_AL\tFATHER_DP\tFATHER_AD\tFATHER_PL\tCHILD_ID\tCHILD_GT\tCHILD_AL\tCHILD_DP\tCHILD_AD\tCHILD_PL");
            //mvFile.println("#CHROM\tPOS\tFILTER\tAC\tTP\t#MV\tFAMILY\tMOTHER_GT\tMOTHER_AL\tMOTHER_DP\tMOTHER_RAD\tMOTHER_AAD\tMOTHER_HRPL\tMOTHER_HETPL\tMOTHER_HAPL\tFATHER_GT\tFATHER_AL\tFATHER_DP\tFATHER_RAD\tFATHER_AAD\tFATHER_HRPL\tFATHER_HETPL\tFATHER_HAPL\tCHILD_GT\tCHILD_AL\tCHILD_DP\tCHILD_RAD\tCHILD_AAD\tCHILD_HRPL\tCHILD_HETPL\tCHILD_HAPL");
            //mvFile.println("#CHROM\tPOS\tFILTER\tAC\tFAMILY\tTP\tMOTHER_GT\tMOTHER_DP\tMOTHER_RAD\tMOTHER_AAD\tMOTHER_HRPL\tMOTHER_HETPL\tMOTHER_HAPL\tFATHER_GT\tFATHER_DP\tFATHER_RAD\tFATHER_AAD\tFATHER_HRPL\tFATHER_HETPL\tFATHER_HAPL\tCHILD_GT\tCHILD_DP\tCHILD_RAD\tCHILD_AAD\tCHILD_HRPL\tCHILD_HETPL\tCHILD_HAPL");
        }
    }

	/**
     * Phases the genotypes of the given trio. If one of the parents is null, it is considered a parent/child pair.
     * @param ref: Reference allele
     * @param alt: Alternative allele
     * @param mother: Mother's genotype
     * @param father: Father's genotype
     * @param child: Child's genotype
     * @return
    */
    private TrioConfiguration phaseTrioGenotypes(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child, EnumMap<GenotypeType,Double> noDnmGtPriors, Gender childSex) 
    {    	

        //If the child wasn't called, just return the phased genotypes if homozygous, without using the familial information
        if(child==null || !child.isCalled() )
		{
            return new UnrelatedTrioConfiguration(ref,alt,mother,father,child,childSex);
        }
        if (op)
        {
        	TrioConfiguration op_conf = new TrioConfiguration(mother.getType(), father.getType(), child.getType(), (int)0, (short)0.5, childSex);
        	op_conf.populateGenotypes(ref,alt,mother,father,child);
        	return op_conf;
        }
        //Check whether it is  a pair or trio
        FamilyConfigurations configurations;
        //Create the configurations container for either a parent/child pair or a trio
        if(mother == null || !mother.isCalled())
        {
            //If both parents are not called, phase the child if homozygous without considering the familial information
            if(father == null || !father.isCalled())
            {
                return new UnrelatedTrioConfiguration(ref,alt,mother,father,child,childSex);
            }
            configurations = new PairConfigurations(mother,FamilyMember.MOTHER, childSex, rand);
        }
        else 
        	if(father == null || !father.isCalled())
        	{
        		configurations = new PairConfigurations(father,FamilyMember.FATHER, childSex, rand);
        	}
        	else
        	{
        		configurations = new TrioConfigurations(rand,childSex);
        	}

        //Genotype likelihoods as EnumMap. In case of a missing parent, equal likelihoods are assigned to each of its genotypes.
        EnumMap<GenotypeType,Double> motherLikelihoods = getLikelihoodsAsMapSafeNull(mother,Gender.FEMALE);
        EnumMap<GenotypeType,Double> fatherLikelihoods = getLikelihoodsAsMapSafeNull(father,Gender.MALE);
        EnumMap<GenotypeType,Double> childLikelihoods = getLikelihoodsAsMapSafeNull(child,childSex);

        //Compute the likelihood for each trio configurations and store them
        TrioConfiguration bestConfiguration;
        int mvCount;
        
        double configurationLikelihood;
        TrioPhase currentConfiguration;

        for(Map.Entry<GenotypeType,Double> childGenotype : childLikelihoods.entrySet())
        {
            for(Map.Entry<GenotypeType,Double> motherGenotype : motherLikelihoods.entrySet())
            {
                for(Map.Entry<GenotypeType,Double> fatherGenotype : fatherLikelihoods.entrySet())
                {                    	
                	if (site_male_ploidy == 1) 
                    	if (childSex == Gender.MALE)
                    		currentConfiguration = transmissionMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey()).get(1);
                    	else 
                    		currentConfiguration = transmissionMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey()).get(0);
                    else
                    	currentConfiguration = transmissionMatrix.get(motherGenotype.getKey()).get(fatherGenotype.getKey()).get(childGenotype.getKey()).get(2);
                    
                    mvCount = currentConfiguration.getMVCount();
                    //System.out.printf("site_male_ploidy:%d\tcurr_conf_father_gt:%s\tfather GT: %s\n", site_male_ploidy, currentConfiguration.getDummyGT(FamilyMember.FATHER).alleles.size() == 1 ? currentConfiguration.getDummyGT(FamilyMember.FATHER).alleles.get(0).getBaseString() : currentConfiguration.getDummyGT(FamilyMember.FATHER).alleles.get(0).getBaseString() + currentConfiguration.getDummyGT(FamilyMember.FATHER).alleles.get(1).getBaseString(), father.getAlleles().size() == 1 ? father.getAlleles().get(0).getBaseString() : father.getAlleles().get(0).getBaseString() + father.getAlleles().get(1).getBaseString() );
                    
                    if (ufad && child.getPloidy() == 1 && site_male_ploidy == 1)
                    	configurationLikelihood =  currentConfiguration.getTransmissionProb() + noDnmGtPriors.get(motherGenotype.getKey()) + motherGenotype.getValue() + childGenotype.getValue();
                    else
                    	configurationLikelihood =  currentConfiguration.getTransmissionProb() + noDnmGtPriors.get(motherGenotype.getKey()) + motherGenotype.getValue() + father.getPloidy() * 0.5 * noDnmGtPriors.get(fatherGenotype.getKey()) + fatherGenotype.getValue() + childGenotype.getValue();
                    
                    configurations.addConfiguration(motherGenotype.getKey(), fatherGenotype.getKey(), childGenotype.getKey(), mvCount, configurationLikelihood);
                }
            }
        }

        if(updatePLs)
        {
            GenotypeLikelihoods childPLs = configurations.getChildPLs();
            GenotypeLikelihoods motherPLs = configurations.getMotherPLs();
            GenotypeLikelihoods fatherPLs = configurations.getFatherPLs();

            if(gtAssignment.equals(TRIOCONF))
            {
                bestConfiguration = configurations.getMostLikelyConfiguration(father);
            }
            else
            {
                bestConfiguration = configurations.getMostLikelyConfigurationFromPLs(motherPLs, fatherPLs, childPLs);
            }
            bestConfiguration.populateGenotypes(ref,alt,mother,father,child,motherPLs,fatherPLs,childPLs);
        }
        else
        {
            bestConfiguration = configurations.getMostLikelyConfiguration(father);
            bestConfiguration.populateGenotypes(ref,alt,mother,father,child);
        }

        return bestConfiguration;
    }
     
    /**
     * 
     * Stores information regarding some trio-combination, after evaluation; has multiple data access patterns
     *
     */
    private class TrioConfiguration 
    {
        private GenotypeType phasedMotherGT;
        private GenotypeType phasedFatherGT;
        private GenotypeType phasedChildGT;
        private short posterior;
        private int mvCount;
        private Genotype phasedMother;
        private Genotype phasedFather;
        private Genotype phasedChild;
        //private BigDecimal mvBurden;
        private Gender childSex;

        public TrioConfiguration(GenotypeType phasedMotherGT, GenotypeType phasedFatherGT, GenotypeType phasedChildGT, int numMVs, short posterior, Gender childSex)
        {
            this.phasedMotherGT = phasedMotherGT;
            this.phasedFatherGT = phasedFatherGT;
            this.phasedChildGT = phasedChildGT;
            this.posterior = posterior;
            this.mvCount = numMVs;
            //this.mvBurden = mvBurden;
            this.childSex = childSex;
        }
        
        /**
         * build the actual genotypes of a trio-combination by aggregating the cached phased information, original Genotype and pbt results
         * @param ref reference allele at locus
         * @param alt alternative allele at locus
         * @param motherGenotype initial mother's GT
         * @param fatherGenotype initial father's GT
         * @param childGenotype initial child's GT
         */
        public void populateGenotypes(Allele ref, Allele alt, Genotype motherGenotype, Genotype fatherGenotype, Genotype childGenotype)
        {
        	TrioPhase phasedTrioGenotypes;
        	//System.out.printf("site_male_ploidy: %d\n", site_male_ploidy);
            if (site_male_ploidy == 1)
            	if (childSex == Gender.MALE)
            		phasedTrioGenotypes = transmissionMatrix.get(phasedMotherGT).get(phasedFatherGT).get(phasedChildGT).get(1);
            	else
            		phasedTrioGenotypes = transmissionMatrix.get(phasedMotherGT).get(phasedFatherGT).get(phasedChildGT).get(0);
            else
            	phasedTrioGenotypes = transmissionMatrix.get(phasedMotherGT).get(phasedFatherGT).get(phasedChildGT).get(2);
            
            //System.out.printf("##### --> mother hashed GT: %s %s\t father hashed GT: %s %s\t child hashed GT: %s %s\n",phasedTrioGenotypes.getDummyGT(FamilyMember.MOTHER).alleles.get(0).getBaseString(), phasedTrioGenotypes.getDummyGT(FamilyMember.MOTHER).alleles.size() > 1 ? phasedTrioGenotypes.getDummyGT(FamilyMember.MOTHER).alleles.get(1).getBaseString() : "", phasedTrioGenotypes.getDummyGT(FamilyMember.FATHER).alleles.get(0).getBaseString(), phasedTrioGenotypes.getDummyGT(FamilyMember.FATHER).alleles.size() > 1 ? phasedTrioGenotypes.getDummyGT(FamilyMember.FATHER).alleles.get(1).getBaseString() : "", phasedTrioGenotypes.getDummyGT(FamilyMember.CHILD).alleles.get(0).getBaseString(), phasedTrioGenotypes.getDummyGT(FamilyMember.CHILD).alleles.size() > 1 ? phasedTrioGenotypes.getDummyGT(FamilyMember.CHILD).alleles.get(1).getBaseString() : "" );
            //System.out.printf("\t\tfather(%s) GT before call: %s%s\n", fatherGenotype.getSampleName(), fatherGenotype.getAlleles().get(0).getBaseString(), fatherGenotype.getAlleles().size() > 1 ? fatherGenotype.getAlleles().get(1).getBaseString(): "");
            //System.out.printf("\t\tmother(%s) GT before call: %s%s\n", motherGenotype.getSampleName(), motherGenotype.getAlleles().get(0).getBaseString(), motherGenotype.getAlleles().size() > 1 ? motherGenotype.getAlleles().get(1).getBaseString(): "");
            phasedMother = getPhasedGenotype(ref, alt,motherGenotype, phasedTrioGenotypes.getDummyGT(FamilyMember.MOTHER), posterior, null, Gender.FEMALE);
            phasedFather = getPhasedGenotype(ref,alt,fatherGenotype, phasedTrioGenotypes.getDummyGT(FamilyMember.FATHER), posterior, null, Gender.MALE);
            phasedChild  = getPhasedGenotype(ref,alt,childGenotype, phasedTrioGenotypes.getDummyGT(FamilyMember.CHILD), posterior, null, childSex);
        }
        
        /**
         * build the actual genotypes of a trio-combination by aggregating the cached phased information, original Genotype and pbt results; update individual genotype PLs as well
         * @param ref reference allele at locus
         * @param alt alternative allele at locus
         * @param motherGenotype initial mother's GT
         * @param fatherGenotype initial father's GT
         * @param childGenotype initial child's GT
         */
        public void populateGenotypes(Allele ref, Allele alt, Genotype motherGenotype, Genotype fatherGenotype, Genotype childGenotype, GenotypeLikelihoods motherPLs, GenotypeLikelihoods fatherPLs, GenotypeLikelihoods childPLs)
        {
        	TrioPhase phasedTrioGenotypes;
        	if (site_male_ploidy == 1)
        		if (childSex == Gender.MALE)
            		phasedTrioGenotypes = transmissionMatrix.get(phasedMotherGT).get(phasedFatherGT).get(phasedChildGT).get(1);
             	else
            		phasedTrioGenotypes = transmissionMatrix.get(phasedMotherGT).get(phasedFatherGT).get(phasedChildGT).get(0);
            else
            	phasedTrioGenotypes = transmissionMatrix.get(phasedMotherGT).get(phasedFatherGT).get(phasedChildGT).get(2);
            
        	phasedMother = getPhasedGenotype(ref, alt,motherGenotype, phasedTrioGenotypes.getDummyGT(FamilyMember.MOTHER), posterior, motherPLs.getAsPLs(), Gender.FEMALE);
            phasedFather = getPhasedGenotype(ref,alt,fatherGenotype, phasedTrioGenotypes.getDummyGT(FamilyMember.FATHER), posterior, fatherPLs.getAsPLs(), Gender.MALE);
            phasedChild  = getPhasedGenotype(ref,alt,childGenotype, phasedTrioGenotypes.getDummyGT(FamilyMember.CHILD), posterior, childPLs.getAsPLs(), childSex);
        }
        
        // build an actual genotype by aggregating the cached phased information, original Genotype and pbt results; update the PLs field value as necessary
		private Genotype getPhasedGenotype(Allele ref, Allele alt, Genotype genotype, DummyGT phasedGenotype, int transmissionPost, int []pls, Gender sex)
		{
			if (genotype == null || !isPhasable(genotype.getType(), sex))
				return genotype;
			
			// System.out.printf("**** --> hashed GT: %s %s\t GT: %s %s\t %s\n",phasedGenotype.alleles.get(0).getBaseString(), phasedGenotype.alleles.size() > 1 ? phasedGenotype.alleles.get(1).getBaseString() : "", genotype.getAlleles().get(0).getBaseString(), genotype.getAlleles().size() > 1 ? genotype.getAlleles().get(1).getBaseString() : "", genotype.getSampleName());
			GenotypeBuilder gb = new GenotypeBuilder();
			ArrayList<Allele> alleles = null;
			alleles = new ArrayList<Allele>(2);
			if (phasedGenotype != null)
				{
					//alleles = new ArrayList<Allele>(2);
					if (phasedGenotype.alleles.get(0).isReference())
						alleles.add(ref);
					else
						alleles.add(alt);
					if (phasedGenotype.alleles.size() == 2 && genotype.getPloidy() == 2)
						if (phasedGenotype.alleles.get(1).isReference())
							alleles.add(ref);
						else
							alleles.add(alt);
					else
						if ( phasedGenotype.alleles.size() != genotype.getPloidy() )
						//throw error; print the conflicting allele (sets)
							throw new UserException(String.format("Error: ploidy missmatch between sample genotype(:%s) and corresponding phased hash entry(:%s) in sample %s\t%s\n", (genotype.getAlleles().size() > 1)? genotype.getAlleles().get(0).getBaseString() + genotype.getAlleles().get(1).getBaseString()  : genotype.getAlleles().get(0).getBaseString(), (phasedGenotype.alleles.size() == 2)?phasedGenotype.alleles.get(0).getBaseString() + phasedGenotype.alleles.get(1).getBaseString() : phasedGenotype.alleles.get(0).getBaseString(), genotype.getSampleName(), sex.name()));
				}
			Map<String,Object> new_attributes = new HashMap<String, Object>();
			Map<String,Object> initial_attributes = genotype.getExtendedAttributes();
			
			// ----------- Compatibility with existing RBP information
			if ( initial_attributes.containsKey(RBP_TAG) )
				{
					String rb_phasing = initial_attributes.remove(RBP_TAG).toString();
					if (rb_phasing.compareTo(".") != 0)
					{
						ArrayList<Allele> initial_als = new ArrayList<>(2);
						initial_als.addAll(genotype.getAlleles());
						// if the alleles of the GT are changed by PBT, the existing RBP haplotype info  becomes unreliable ( changing a genotype HET-to-HOM at this position, influences at least the quality of any downstream HET genotype that is phased w.r.t. this genotype or the same genotype as this genotype)
						if ( ( genotype.getPloidy() == 1 && initial_als.get(0) != alleles.get(0) ) || 
								( genotype.getPloidy() == 2 && !( initial_als.get(0) == alleles.get(0) && initial_als.get(1) == alleles.get(1) || 
																  initial_als.get(1) == alleles.get(0) && initial_als.get(0) == alleles.get(1) ) ) )
						{
							// leave the RBP info as it was, just another tag indicating that the GT was changed
							// TODO:check if this a VALID thing to do and also if there are better options 
							new_attributes.put(RBP_TAG, rb_phasing);
							new_attributes.put(RBP_CANCELLED, RBP_CANCELLED_VAL); 
						}
						else
						{
							if ( genotype.getPloidy() == 2 && initial_als.get(1) == alleles.get(0) && initial_als.get(0) == alleles.get(1)) // alleles got inverted by pbt phasing => invert HP tags
								{
									String tokens[] = rb_phasing.split(",");
									if (tokens.length == 2)
									{
										rb_phasing = tokens[1] + "," +  tokens[0];
										new_attributes.put(RBP_TAG, rb_phasing);
									}
									else
										throw new UserException(String.format("Error: phasing tag of diploid genotype malformed: %s\n", rb_phasing));
								}
							else
								new_attributes.put(RBP_TAG, rb_phasing);
						}
					}
					else // RBP phasing is "."....put it right back in
						new_attributes.put(RBP_TAG, rb_phasing);
						
				}
			// ----------- Compatibility with existing RBP information
			
			new_attributes.putAll(initial_attributes);
			if (!noTP)
				new_attributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, transmissionPost);
			
			gb.alleles(alleles);
			gb.phased(phasedGenotype.phased);
			gb.attributes(new_attributes);
			gb.DP(genotype.getDP());
			gb.AD(genotype.getAD());
			gb.name(genotype.getSampleName());
			
			int GQ;
			if (updatePLs && pls != null)
			{
				if (pls.length == 3)
					GQ = Math.max(Math.max(Math.min(pls[0], pls[1]), Math.min(pls[0],pls[2])), Math.min(pls[1],pls[2]));
				else GQ = Math.max(pls[0],pls[1]);
				gb.PL(pls);
				gb.GQ(GQ);
					
			}
			else
			{
				gb.GQ(genotype.getGQ());
				gb.PL(genotype.getPL());
			}
			
			return gb.make();
		}
		        
		public Allele getMutantAllele()
		{
			if ( getConfigurationMVCount(phasedMotherGT, phasedFatherGT, phasedChildGT, childSex) != 1 )
				return null;
			Allele _ret = null;
			List<Allele> child_als = phasedChild.getAlleles();
			List<Allele> father_als = phasedFather.getAlleles();
			List<Allele> mother_als = phasedMother.getAlleles();
			if ( child_als.size() == 1 ) // child haploid => male child, X-chromosome => can only inherit from the mother
			{
				if ( !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1))  )
					_ret = child_als.get(0);
			}
			else // child diploid
				if ( father_als.size() == 1 ) // father is haploid ( => X-chromosome )
					if ( !child_als.get(0).basesMatch(father_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
						_ret = child_als.get(0);
					else
						if ( !child_als.get(1).basesMatch(father_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
							_ret = child_als.get(1);
						else // next cases are for dnms that in practice are usually discarded (non HOM_REF parents) by post-calling filtering (i.e.: outside PBT); but still... 
							if ( child_als.get(0).basesMatch(father_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
								_ret = child_als.get(1);
							else
							{
								if ( child_als.get(1).basesMatch(father_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
									_ret = child_als.get(0);
							}
				else // father is diploid
					if ( !child_als.get(0).basesMatch(father_als.get(0)) && !child_als.get(0).basesMatch(father_als.get(1)) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
						_ret = child_als.get(0);
					else
						if ( !child_als.get(1).basesMatch(father_als.get(0)) && !child_als.get(1).basesMatch(father_als.get(1)) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
							_ret = child_als.get(1);
						else // next cases are for dnms that in practice are usually discarded (non HOM_REF parents) by post-calling filtering (i.e.: outside PBT); but still... 
							if ( ( child_als.get(0).basesMatch(father_als.get(0)) || child_als.get(0).basesMatch(father_als.get(1)) ) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
								_ret = child_als.get(1);
							else
								if ( ( child_als.get(1).basesMatch(father_als.get(0)) || child_als.get(1).basesMatch(father_als.get(1)) ) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
									_ret = child_als.get(0);
			
			return _ret;					
		}
		
		public int getMutantAlleleIdx()
		{
			if ( getConfigurationMVCount(phasedMotherGT, phasedFatherGT, phasedChildGT, childSex) != 1 )
				return -1;
			int _ret = -1;
			List<Allele> child_als = phasedChild.getAlleles();
			List<Allele> father_als = phasedFather.getAlleles();
			List<Allele> mother_als = phasedMother.getAlleles();
			if ( child_als.size() == 1 ) // child haploid => male child, X-chromosome => can only inherit from the mother
			{
				if ( !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1))  )
					_ret = 0;
			}
			else // child diploid
				if ( father_als.size() == 1 ) // father is haploid ( => X-chromosome )
					if ( !child_als.get(0).basesMatch(father_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
						_ret = 0;
					else
						if ( !child_als.get(1).basesMatch(father_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
							_ret = 1;
						else // next cases are for dnms that in practice are usually discarded (non HOM_REF parents) by post-calling filtering (i.e.: outside PBT); but still... 
							if ( child_als.get(0).basesMatch(father_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
								_ret = 1;
							else
							{
								if ( child_als.get(1).basesMatch(father_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
									_ret = 0;
							}
				else // father is diploid
					if ( !child_als.get(0).basesMatch(father_als.get(0)) && !child_als.get(0).basesMatch(father_als.get(1)) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
						_ret = 0;
					else
						if ( !child_als.get(1).basesMatch(father_als.get(0)) && !child_als.get(1).basesMatch(father_als.get(1)) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
							_ret = 1;
						else // next cases are for dnms that in practice are usually discarded (non HOM_REF parents) by post-calling filtering (i.e.: outside PBT); but still... 
							if ( ( child_als.get(0).basesMatch(father_als.get(0)) || child_als.get(0).basesMatch(father_als.get(1)) ) && !child_als.get(1).basesMatch(mother_als.get(0)) && !child_als.get(1).basesMatch(mother_als.get(1)) )
								_ret = 1;
							else
								if ( ( child_als.get(1).basesMatch(father_als.get(0)) || child_als.get(1).basesMatch(father_als.get(1)) ) && !child_als.get(0).basesMatch(mother_als.get(0)) && !child_als.get(0).basesMatch(mother_als.get(1)) )
									_ret = 0;
			
			return _ret;					
		}
		
		// just interface/wrapper --- remember to delete
		public void buildFamilyMemberGT(Allele ref, Allele alt, Genotype genotype, DummyGT phasedGT, FamilyMember individual, int transmissionPost, int []pls)
		{
			switch (individual)
			{
			case MOTHER:
				phasedMother = getPhasedGenotype(ref,alt,genotype,phasedGT, transmissionPost, pls, Gender.FEMALE);
			case FATHER:
				phasedFather = getPhasedGenotype(ref,alt,genotype, phasedGT, transmissionPost, pls, Gender.MALE);
			case CHILD:
				phasedChild  = getPhasedGenotype(ref,alt,genotype, phasedGT, transmissionPost, pls, childSex);
			
			default: 
				return;	
			}
		}
		
		// just interface/wrapper --- remember to delete
		public void buildFamilyMemberGT(Allele ref, Allele alt, Genotype genotype, DummyGT phasedGT, FamilyMember individual, int transmissionPost)
		{
			switch (individual)
			{
			case MOTHER:
				phasedMother = getPhasedGenotype(ref,alt,genotype, phasedGT, transmissionPost, null, Gender.FEMALE);
			case FATHER:
				phasedFather = getPhasedGenotype(ref,alt,genotype, phasedGT, transmissionPost, null, Gender.MALE);
			case CHILD:
				phasedChild  = getPhasedGenotype(ref,alt,genotype, phasedGT, transmissionPost, null, childSex);

			default: 
				return;	
			}
		}

        public Genotype getPhasedMother() {return phasedMother;}
        public Genotype getPhasedFather() {return phasedFather;}
        public Genotype getPhasedChild() {return phasedChild;}
        public short getPosterior() {return posterior;}
        public int getMVCount() {return mvCount;}
        //public BigDecimal getMVBurden() {return mvBurden;}
        public boolean isMV() {return mvCount>0;}
        public Gender getChildGender() {return childSex;}
        
    }

    
    private class UnrelatedTrioConfiguration extends TrioConfiguration 
	{

        public UnrelatedTrioConfiguration(Allele ref, Allele alt, Genotype mother, Genotype father, Genotype child,Gender childSex){
            super((mother != null)?mother.getType():GenotypeType.UNAVAILABLE,(father != null)?father.getType():GenotypeType.UNAVAILABLE,(child != null)?child.getType():GenotypeType.UNAVAILABLE,0,NO_TRANSMISSION_PROB,childSex);
            populateGenotypes(ref,alt,mother,father,child);
        }

    }
    
    /**
     * 
     * Base class for all the different types of family configurations that can be processed by PBT currently. Used in storing and managing the trio-combinations when evaluating the most likely one
     *
     */
    private abstract class FamilyConfigurations 
    {
        
        protected EnumMap<GenotypeType,EnumMap<GenotypeType,EnumMap<GenotypeType, Pair<Integer,Double>>>> configurations;
        protected BigDecimal norm = new BigDecimal(0.0);
        //protected BigDecimal totalMV = new BigDecimal(0.0);
        protected Random rand;
        protected Gender childSex;

        //Store only the informative genotypes as not to loop over unused genotypes        	
        protected GenotypeType[] informativeGenotypes = {GenotypeType.HOM_REF,GenotypeType.HET,GenotypeType.HOM_VAR};
        protected GenotypeType[] informativeGenotypesHap = {GenotypeType.HOM_REF,GenotypeType.HOM_VAR};

        public FamilyConfigurations(Random rand,Gender childSex)
        {

            //Create an empty configuration matrix
            configurations = new EnumMap<GenotypeType, EnumMap<GenotypeType, EnumMap<GenotypeType, Pair<Integer,Double>>>>(GenotypeType.class);
            for(GenotypeType firstParent : GenotypeType.values()) 
            {
                configurations.put(firstParent, new EnumMap<GenotypeType, EnumMap<GenotypeType, Pair<Integer, Double>>>(GenotypeType.class));
                for(GenotypeType secondParent : GenotypeType.values()) 
                {
                    configurations.get(firstParent).put(secondParent, new EnumMap<GenotypeType, Pair<Integer,Double>>(GenotypeType.class));
                }
            }

            this.rand = rand;
            this.childSex = childSex;

        }
        
        public abstract TrioConfiguration getMostLikelyConfiguration( Genotype father );

        public abstract TrioConfiguration getConfiguration(GenotypeType mother, GenotypeType father, GenotypeType child);

        public abstract TrioConfiguration getMostLikelyConfigurationFromPLs(GenotypeLikelihoods motherPLs, GenotypeLikelihoods fatherPLs, GenotypeLikelihoods childPLs);
        
        public void addConfiguration(GenotypeType mother, GenotypeType father, GenotypeType child, int numMVs, double configurationLikelihood)
        {
            configurations.get(mother).get(father).put(child,new Pair<Integer, Double>(numMVs,configurationLikelihood));
            
            BigDecimal linearConfigurationLikelihood =new BigDecimal(Math.pow(10.0,configurationLikelihood));
            norm = norm.add(linearConfigurationLikelihood);
        }
        
        public GenotypeLikelihoods getMotherPLs()
        {
            return getIndividualPLs(FamilyMember.MOTHER);
        }

        public GenotypeLikelihoods getFatherPLs()
        {
            return getIndividualPLs(FamilyMember.FATHER);
        }

        public GenotypeLikelihoods getChildPLs()
        {
            return getIndividualPLs(FamilyMember.CHILD);
        }

        //Returns the combination posterior probability from the linear scale likelihood
        protected short getPosteriorPhred(BigDecimal likelihood, BigDecimal norm)
        {
            likelihood =  likelihood.divide(norm,Math.max(norm.scale(),likelihood.scale()),RoundingMode.CEILING);

            //Division is using the biggest scale; this is a tricky issue as the norm and combination prob can be very similar
            likelihood =  new BigDecimal(1.0).subtract(likelihood);
            short posterior = (short)(-10*Math.log10(likelihood.doubleValue()));

            //Take care of overflow
            if(posterior<0)
                return Short.MAX_VALUE;
            else
                return posterior;
        }

        //Returns the combination posterior probability from the log scale likelihood
        protected short getPosteriorPhred(double logLikelihood, BigDecimal norm)
        {
            return getPosteriorPhred(new BigDecimal(Math.pow(10.0,logLikelihood)),norm);
        }

        /**
         * marginalizes from the joint PLs distribution in a trio for each/any family individual; ploidy aware 
         * @param individual
         * @return GenotypeLikelihoods object containing the appropriate, recomputed PLs
         * 
         */
        
        protected GenotypeLikelihoods getIndividualPLs(FamilyMember individual)
        {
            ArrayList<Double> configurationsLikelihoods = new ArrayList<Double>(9);
            ArrayList<Double> rawLikelihoods = new ArrayList<Double>(3);
            
            // ---> !!!: We assume that the trio-combination storing structure's nesting ("configurations") is mother-father-child (which is currently the case)
            switch(individual)
            {
            case MOTHER:
            	for(GenotypeType individualGT : informativeGenotypes)
            	{
            		//Sum up the configuration likelihoods for a given genotype
            		for(GenotypeType fatherGT : (site_male_ploidy == 1)? informativeGenotypesHap : informativeGenotypes)
            			for(GenotypeType childGT : (site_male_ploidy == 1 && childSex == Gender.MALE)? informativeGenotypesHap : informativeGenotypes)
            				if (configurations.containsKey(individualGT))
            					if (configurations.get(individualGT).containsKey(fatherGT))
            						if (configurations.get(individualGT).get(fatherGT).containsKey(childGT))
            							configurationsLikelihoods.add(configurations.get(individualGT).get(fatherGT).get(childGT).getSecond());
            		
            		if (configurationsLikelihoods.size() != 0)
                    {
            			rawLikelihoods.add(MathUtils.log10sumLog10(getPrimitiveArray(configurationsLikelihoods)));
                    }
            		
            		configurationsLikelihoods.clear();
            	}
            	break;
            case FATHER:
            	for(GenotypeType individualGT : (site_male_ploidy == 1)? informativeGenotypesHap : informativeGenotypes)
            	{
            		//Sum up the configuration likelihoods for a given genotype
            		for(GenotypeType motherGT : informativeGenotypes)
            			for(GenotypeType childGT : (site_male_ploidy == 1 && childSex == Gender.MALE)? informativeGenotypesHap : informativeGenotypes)
            				if (configurations.containsKey(motherGT))
            					if (configurations.get(motherGT).containsKey(individualGT))
            						if (configurations.get(motherGT).get(individualGT).containsKey(childGT))
            							configurationsLikelihoods.add(configurations.get(motherGT).get(individualGT).get(childGT).getSecond());
            		
            		if (configurationsLikelihoods.size() != 0)
                    {
            			rawLikelihoods.add(MathUtils.log10sumLog10(getPrimitiveArray(configurationsLikelihoods)));
                    }
            		
            		configurationsLikelihoods.clear();
            	}
            	break;
            case CHILD:
            	for(GenotypeType individualGT : (site_male_ploidy == 1 && childSex == Gender.MALE)? informativeGenotypesHap : informativeGenotypes)
            	{
            		//Sum up the configuration likelihoods for a given genotype
            		for(GenotypeType motherGT : informativeGenotypes)
            			for(GenotypeType fatherGT : (site_male_ploidy == 1)? informativeGenotypesHap : informativeGenotypes)
            				if (configurations.containsKey(motherGT))
            					if (configurations.get(motherGT).containsKey(fatherGT))
            						if (configurations.get(motherGT).get(fatherGT).containsKey(individualGT))
            							configurationsLikelihoods.add(configurations.get(motherGT).get(fatherGT).get(individualGT).getSecond());
            		
            		if (configurationsLikelihoods.size() != 0)
                    {
            			rawLikelihoods.add(MathUtils.log10sumLog10(getPrimitiveArray(configurationsLikelihoods)));
                    }
            		
            		configurationsLikelihoods.clear();
            	}
            	break;
            }
            
            Double []dummy1 = new Double [0];
            Double []dummy;
            dummy = rawLikelihoods.toArray(dummy1);
            double []rawLikelihoods2 =new double[dummy.length];
            for (int i = 0; i < rawLikelihoods2.length; i++)
            	rawLikelihoods2[i] = dummy[i].doubleValue();
            return GenotypeLikelihoods.fromLog10Likelihoods(MathUtils.normalizeFromLog10(rawLikelihoods2,true,true));
        }
        
    }
    
	private class TrioConfigurations extends FamilyConfigurations 
	{

        public TrioConfigurations(Random rand, Gender childSex)
        {
            super(rand,childSex);
        }

        public TrioConfiguration getConfiguration(GenotypeType firstParent, GenotypeType secondParent, GenotypeType child)
        {
            return new TrioConfiguration(firstParent,secondParent,child, configurations.get(firstParent).get(secondParent).get(child).getFirst(), 
									getPosteriorPhred(configurations.get(firstParent).get(secondParent).get(child).getSecond(), norm), childSex);
        }

        /**
         * @param father: the most likely individual genotype for the father; needed for the -ufad option for X-like inheritance pattern only
         * @return a TrioConfiguration storing the most likely trio-genotype
         */
        public TrioConfiguration getMostLikelyConfiguration(Genotype father)
        {

            //Keeps all most likely configurations
            double maxLikelihood = Double.NEGATIVE_INFINITY;
            ArrayList<GenotypeType> bestMother = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestFather = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestChild = new ArrayList<GenotypeType>(2);

            //Find the most likely configuration(s)
            for(GenotypeType motherGT : informativeGenotypes)
			{
                for(GenotypeType fatherGT : (site_male_ploidy == 1)? informativeGenotypesHap : informativeGenotypes)
				{
                    for(GenotypeType childGT : (site_male_ploidy == 1 && childSex == Gender.MALE)? informativeGenotypesHap : informativeGenotypes)
					{
                        if ( configurations.get(motherGT).get(fatherGT).get(childGT) != null )
                        {
                        	if(configurations.get(motherGT).get(fatherGT).get(childGT).getSecond() > maxLikelihood)
                        	{
                        		maxLikelihood = configurations.get(motherGT).get(fatherGT).get(childGT).getSecond();
                        		bestMother.clear();
                        		bestFather.clear();
                        		bestChild.clear();
                        	}
                        	if(configurations.get(motherGT).get(fatherGT).get(childGT).getSecond() == maxLikelihood)
                        	{
                        		bestMother.add(motherGT);
                        		bestFather.add(fatherGT);
                        		bestChild.add(childGT);
                        	}
                        }
                    }
                }
            }
            
            int configuration_index;
            if (ufad && childSex == Gender.MALE && site_male_ploidy == 1)
            	configuration_index = bestFather.indexOf(father.getType());
            else
            //Report the most likely configuration
            //Pick one at random if multiple ones are reported
            	configuration_index = bestMother.size()>1 ? rand.nextInt(bestMother.size()-1) : 0;
            return new TrioConfiguration(bestMother.get(configuration_index), bestFather.get(configuration_index), bestChild.get(configuration_index), 
									configurations.get(bestMother.get(configuration_index)).get(bestFather.get(configuration_index)).get(bestChild.get(configuration_index)).getFirst(), getPosteriorPhred(maxLikelihood,norm), childSex);

        }

        /**
         * Find the most likely genotype configuration given the recalibrated trio-aware individual PLs. The -ufad option for X-like inheritance pattern is NOT enforced in this case --- unsure whether this can ever differ from the best TRIO-GENOTYPE assignment
         * @param motherPLs: recalibrated PLs of the mother
         * @param fatherPLs: recalibrated PLs of the father
         * @param childPLs: recalibrated PLs of the father
         * @return a TrioConfiguration storing the most likely trio-genotype
         */
        public TrioConfiguration getMostLikelyConfigurationFromPLs(GenotypeLikelihoods motherPLs, GenotypeLikelihoods fatherPLs, GenotypeLikelihoods childPLs)
		{

            // Get the most likely genotype(s) individually
            ArrayList<GenotypeType> bestMother = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestFather = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestChild = new ArrayList<GenotypeType>(2);
            //TODO: Verify somehow (or better just test) if .getAsMap() method of GenotypeLikelihoods works correctly for ploidy 1
            for(Map.Entry<GenotypeType,Double> firstParent : getLikelihoodsAsMapSafeNull(motherPLs).entrySet())
			{
                if(firstParent.getValue()==0)
                    bestMother.add(firstParent.getKey());
            }

            for(Map.Entry<GenotypeType,Double> child : getLikelihoodsAsMapSafeNull(childPLs).entrySet())
			{
                if(child.getValue()==0)
                    bestChild.add(child.getKey());
            }
            //If a parent is missing, just assign any genotype as it will be ignored.
            for(Map.Entry<GenotypeType,Double> secondParent : getLikelihoodsAsMapSafeNull(fatherPLs).entrySet())
			{
                if(secondParent.getValue()==0)
                    bestFather.add(secondParent.getKey());
            }

            //Get the corresponding most likely configuration(s)
            ArrayList<TrioConfiguration> bestConfigurations = new ArrayList<TrioConfiguration>(2);
            TrioConfiguration current;
            short maxPosterior = 0;
            for(GenotypeType firstParent : bestMother)
			{
                for(GenotypeType secondParent : bestFather)
				{
                    for(GenotypeType child : bestChild)
					{
                        current = getConfiguration(firstParent, secondParent, child);
                        if(current.getPosterior() > maxPosterior)
						{
                            bestConfigurations.clear();
                            maxPosterior = current.getPosterior();
                        }
                        if(current.getPosterior()== maxPosterior)
                            bestConfigurations.add(current);
                    }
                }
            }
            int configuration_index = bestConfigurations.size()>1 ? rand.nextInt(bestConfigurations.size()-1) : 0;
            return bestConfigurations.get(configuration_index);
        }


    }
	
	//
	private class PairConfigurations extends FamilyConfigurations 
	{

        private Genotype missingParentGT;
        private FamilyMember missingParent;

        public PairConfigurations(Genotype missingParentGT, FamilyMember missingParent, Gender childSex, Random rand)
        {
            super(rand,childSex);
            this.missingParent = missingParent;
            this.missingParentGT = missingParentGT;
        }

        public void addConfiguration(GenotypeType mother, GenotypeType father, GenotypeType child, int numMVs, double configurationLikelihood)
		{
        	configurations.get(mother).get(father).put(child, new Pair<Integer,Double>(numMVs, configurationLikelihood));
        	
            norm = norm.add(new BigDecimal(Math.pow(10.0, configurationLikelihood)));
        }

        public TrioConfiguration getConfiguration(GenotypeType mother, GenotypeType father, GenotypeType child)
		{
            //Get the configuration likelihood based the only 2 individuals available
            int MVcount = 0;
            GenotypeType availableParent = (missingParent == FamilyMember.MOTHER)? father : mother;
            BigDecimal configurationLikelihood = new BigDecimal(0.0);
            //if one parent is missing, average him out
            for(GenotypeType missingParentGenotype : (missingParent == FamilyMember.FATHER && site_male_ploidy == 1) ? informativeGenotypesHap : informativeGenotypes)
			{
                        
            	configurationLikelihood = 
            			missingParent == FamilyMember.MOTHER ?
            			configurationLikelihood.add(new BigDecimal(Math.pow(10.0, configurations.get(missingParentGenotype).get(availableParent).get(child).getSecond()))) : 
            			configurationLikelihood.add(new BigDecimal(Math.pow(10.0, configurations.get(availableParent).get(missingParentGenotype).get(child).getSecond())));
            	
            	/*
            	configurationLikelihood = configurationLikelihood.add(new BigDecimal(Math.pow(10.0, configurations
                        																					.get(missingParent == FamilyMember.MOTHER ? missingParentGenotype : availableParent)
                        																					.get(missingParent == FamilyMember.FATHER ? missingParentGenotype : availableParent)
                        																					.get(child).getSecond())));
                */
                MVcount +=  configurations.get(availableParent).get(missingParentGenotype).get(child).getFirst();
            }

            return new TrioConfiguration(mother, father, child, MVcount/3, getPosteriorPhred(configurationLikelihood,norm),childSex);

        }

        public GenotypeLikelihoods getMotherPLs()
        {
            //Handle missing case
            if(missingParent == FamilyMember.MOTHER)
            	if (missingParentGT == null)
            		return null;
            	else
            		return missingParentGT.getLikelihoods();
            return getIndividualPLs(FamilyMember.MOTHER);
        }

        public GenotypeLikelihoods getFatherPLs()
        {
            //Handle missing case
            if(missingParent == FamilyMember.FATHER)
            	if (missingParentGT != null)
            		return missingParentGT.getLikelihoods();
            	else
            		return null;
            //Otherwise, return the "mother" genotype as the available genotype
            //was always placed in the first map. Unnatural-looking but OK. --- made it look all natural again!
            return getIndividualPLs(FamilyMember.FATHER);
        }

        public TrioConfiguration getMostLikelyConfiguration(Genotype father)
		{

            //Keeps all most likely configurations
            BigDecimal configurationLikelihood;
            BigDecimal maxLikelihood = new BigDecimal(0.0);
            ArrayList<GenotypeType> bestAvailableParent = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestChild = new ArrayList<GenotypeType>(2);
            ArrayList<Integer> bestMVs = new ArrayList<Integer>(2);
            int MVcount;

            //Find the most likely configuration(s)
            for(GenotypeType availableParent : (site_male_ploidy == 1 && missingParent == FamilyMember.MOTHER) ? informativeGenotypesHap : informativeGenotypes)
			{
                for(GenotypeType child : (site_male_ploidy == 1 && childSex == Gender.MALE) ? informativeGenotypesHap : informativeGenotypes)
				{
                    configurationLikelihood = new BigDecimal(0.0);
                    MVcount = 0;
                    for(GenotypeType absentParent : (site_male_ploidy == 1 && missingParent == FamilyMember.FATHER) ? informativeGenotypesHap : informativeGenotypes)
					{
                    	if ( (missingParent == FamilyMember.MOTHER ?
                    			configurations.get(absentParent).get(availableParent).get(child) : 
                    			configurations.get(availableParent).get(absentParent).get(child)) != null )
                    	{
                    		configurationLikelihood = configurationLikelihood.add(new BigDecimal(Math.pow(10.0, 
                    				missingParent == FamilyMember.MOTHER ?
                                		configurations.get(absentParent).get(availableParent).get(child).getSecond() : 
                                		configurations.get(availableParent).get(absentParent).get(child).getSecond() )));
                    		
                    		MVcount +=  missingParent == FamilyMember.MOTHER ?
                            		configurations.get(absentParent).get(availableParent).get(child).getFirst() : 
                            		configurations.get(availableParent).get(absentParent).get(child).getFirst();
                    	}
                    }

                    if ( configurationLikelihood.compareTo(maxLikelihood) > 0 )
					{
                        maxLikelihood = configurationLikelihood;
                        bestAvailableParent.clear();
                        bestChild.clear();
                    }
                    if ( configurationLikelihood.compareTo(maxLikelihood) == 0 )
					{
                        bestAvailableParent.add(availableParent);
                        bestChild.add(child);
                        bestMVs.add( (site_male_ploidy == 1 && missingParent == FamilyMember.FATHER) ? MVcount/2 : MVcount/3);
                    }
                }
            }
            
            int configuration_index;
          
       
            //Report the most likely configuration
            //Pick one at random if multiple ones are reported
            configuration_index = (bestAvailableParent.size() > 1) ? rand.nextInt(bestAvailableParent.size()-1) : 0;

            //Return the genotypes in familial order
            if(missingParent == FamilyMember.MOTHER)
            {
                if (missingParentGT != null)
                	return new TrioConfiguration(missingParentGT.getType(), bestAvailableParent.get(configuration_index), bestChild.get(configuration_index), bestMVs.get(configuration_index), getPosteriorPhred(maxLikelihood,norm),childSex);
            	return new TrioConfiguration(GenotypeType.UNAVAILABLE, bestAvailableParent.get(configuration_index), bestChild.get(configuration_index), bestMVs.get(configuration_index), getPosteriorPhred(maxLikelihood,norm),childSex);
            }
            if (missingParentGT != null)
            	return new TrioConfiguration(bestAvailableParent.get(configuration_index), missingParentGT.getType(), bestChild.get(configuration_index), bestMVs.get(configuration_index), getPosteriorPhred(maxLikelihood,norm),childSex);
            return new TrioConfiguration(bestAvailableParent.get(configuration_index), GenotypeType.UNAVAILABLE, bestChild.get(configuration_index), bestMVs.get(configuration_index), getPosteriorPhred(maxLikelihood,norm),childSex);
        }

        //Find the most likely genotype configuration given the trio-aware PLs
        public TrioConfiguration getMostLikelyConfigurationFromPLs(GenotypeLikelihoods motherPLs, GenotypeLikelihoods fatherPLs, GenotypeLikelihoods childPLs)
		{

            //Get the most likely genotype(s) individually
            ArrayList<GenotypeType> bestMother = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestFather = new ArrayList<GenotypeType>(2);
            ArrayList<GenotypeType> bestChild = new ArrayList<GenotypeType>(2);

            if(missingParent == FamilyMember.MOTHER)
			{
                bestMother.add(GenotypeType.NO_CALL);
                for(Map.Entry<GenotypeType,Double> fatherGQ : fatherPLs.getAsMap(false).entrySet())
				{
                    if(fatherGQ.getValue()==0)
                        bestFather.add(fatherGQ.getKey());
                }
            } 
			else
			{
                bestFather.add(GenotypeType.NO_CALL);
                for(Map.Entry<GenotypeType,Double> motherGQ : motherPLs.getAsMap(false).entrySet())
				{
                    if(motherGQ.getValue()==0)
                        bestMother.add(motherGQ.getKey());
                }
            }

            for(Map.Entry<GenotypeType,Double> childGQ : childPLs.getAsMap(false).entrySet())
			{
                if(childGQ.getValue()==0)
                    bestChild.add(childGQ.getKey());
            }

            //Get the corresponding most likely configuration(s)
            ArrayList<TrioConfiguration> bestConfigurations = new ArrayList<TrioConfiguration>(2);
            TrioConfiguration current;
            short maxPosterior = 0;
            for(GenotypeType motherGT : bestMother)
			{
                for(GenotypeType fatherGT : bestFather)
				{
                    for(GenotypeType childGT : bestChild)
					{
                        current = getConfiguration(motherGT, fatherGT, childGT);
                        if(current.getPosterior() > maxPosterior)
						{
                            bestConfigurations.clear();
                            maxPosterior = current.getPosterior();
                        }
                        if(current.getPosterior() == maxPosterior)
                            bestConfigurations.add(current);
                    }
                }
            }

            int configuration_index = (bestConfigurations.size() > 1)? rand.nextInt(bestConfigurations.size()-1) : 0;
            return bestConfigurations.get(configuration_index);

        }

    }
	
	
	
	/**
     * Process all individuals in the vcf; where trios or pairs can be formed, as speciffied by the ped file, evaluate the joint probability of these individuals' genotypes and assign new genotypes accordingly. Also phase where possible 
     *
     * @param tracker  the reference meta-data tracker
     * @param ref      the reference context
     * @param context  the alignment context
     * @return information for output statistics
     */
    @Override
    public HashMap<Byte,Double> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
    {	
        //Check that there is a bi-allelic VC
        if (tracker == null)
            return parseGenotypesAtLocus(null, null, null);

        VariantContext vc_raw = tracker.getFirstValue(variantCollection.variants, context.getLocation());
        
        if(vc_raw == null) // || !vc_raw.isBiallelic())
            return parseGenotypesAtLocus(null, null, null);
        
        //if (vc_raw.getGenotype("gonl-124c").hasAnyAttribute(ReadBackedPhasing.HP_KEY))
        //{
        //		System.out.printf("TEST: RBP tag for child: %s\n", vc_raw.getGenotype("gonl-124c").getAnyAttribute(ReadBackedPhasing.HP_KEY).toString());
        //}
        //else
        //{
        //	System.out.println("TEST: no RBP tag");
        //}
        
        //determine site ploidy and transform GTs to haploid if necessary
        VariantContext vc;
        if ( haploid_intervals != null )
        {
        	boolean haploid_site = false;
            for (GenomeLoc i: haploid_intervals)
            	if (i.containsP(ref.getLocus()))
            		{ haploid_site = true; break;}
            
            if ( haploid_site )
            {
            	site_male_ploidy = 1;
            	vc = HaploidWriter.processVC(vc_raw, haploidIndiv, haploid_intervals, 0, -1000, false); // transform initially diploid GTs to haploid
            }
            else 
            {
            	site_male_ploidy = 2;
            	vc = new VariantContextBuilder(vc_raw).make(); // just keep the original variant context if not
            	//vc.getAlternateAllele(0).getBaseString();
            }
        }
        else // if no specific haploid intervals were supplied, check whether all of the input is to be treated as haploid
        	if (allX) 
        	{
        		site_male_ploidy = 1;
            	vc = HaploidWriter.processVC(vc_raw, haploidIndiv, haploid_intervals, 0, -1000, false); // transform initially diploid GTs to haploid
        	}
        	else
        	{
        		site_male_ploidy = 2;
            	vc = new VariantContextBuilder(vc_raw).make(); // just keep the original variant context if not
        	}
        
        //Prepare the AF priors;
        VariantContext vcAF = null;
        if (afRod != null)
        	vcAF = tracker.getFirstValue(afRod,context.getLocation());
        EnumMap<GenotypeType,Double> noDnmGtPriors = buildAFPrior(vc,vcAF);
        
        
        if ( !vc_raw.isBiallelic() && phaseEverything )
		{
			phaseNonBiallelicSites(vc);
			return parseGenotypesAtLocus(null, null, null);
		}

        return parseGenotypesAtLocus(vc, noDnmGtPriors, ref.getLocus());
        
    }
    
	
	/**
     * Initializes the reporting counters.
     *
     * @return All counters initialized to 0
     */
    @Override
    public HashMap<Byte,Double> reduceInit() {
        HashMap<Byte,Double> metricsCounters = new HashMap<Byte, Double>(10);
        metricsCounters.put(NUM_TRIO_GENOTYPES_CALLED,0.0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_NOCALL,0.0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_PHASED,0.0);
        metricsCounters.put(NUM_TRIO_HET_HET_HET,0.0);
        metricsCounters.put(NUM_TRIO_VIOLATIONS,0.0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_CALLED,0.0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_NOCALL,0.0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_PHASED,0.0);
        metricsCounters.put(NUM_PAIR_HET_HET,0.0);
        metricsCounters.put(NUM_PAIR_VIOLATIONS,0.0);
        metricsCounters.put(NUM_TRIO_DOUBLE_VIOLATIONS,0.0);
        metricsCounters.put(NUM_GENOTYPES_MODIFIED,0.0);
        metricsCounters.put(NUM_MUTATIONS,0.0);
        metricsCounters.put(NUM_INCOMPLETE_TRIO_PHASED,0.0);
        metricsCounters.put(NUM_INCOMPLETE_TRIO_VIOLATIONS,0.0);
        metricsCounters.put(NUM_TRIO_SINGLES_PHASED,0.0);
        metricsCounters.put(NUM_PAIR_SINGLES_PHASED,0.0);

        return metricsCounters;
    }
	
	/**
     * Adds the value of the site phased to the reporting counters.
     *
     * @param value Site values
     * @param sum   accumulator for the reporting counters
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public HashMap<Byte,Double> reduce(HashMap<Byte,Double> value, HashMap<Byte,Double> sum) {
        sum.put(NUM_TRIO_GENOTYPES_CALLED,value.get(NUM_TRIO_GENOTYPES_CALLED)+sum.get(NUM_TRIO_GENOTYPES_CALLED));
        sum.put(NUM_TRIO_GENOTYPES_NOCALL,value.get(NUM_TRIO_GENOTYPES_NOCALL)+sum.get(NUM_TRIO_GENOTYPES_NOCALL));
        sum.put(NUM_TRIO_GENOTYPES_PHASED,value.get(NUM_TRIO_GENOTYPES_PHASED)+sum.get(NUM_TRIO_GENOTYPES_PHASED));
        sum.put(NUM_TRIO_HET_HET_HET,value.get(NUM_TRIO_HET_HET_HET)+sum.get(NUM_TRIO_HET_HET_HET));
        sum.put(NUM_TRIO_VIOLATIONS,value.get(NUM_TRIO_VIOLATIONS)+sum.get(NUM_TRIO_VIOLATIONS));
        sum.put(NUM_PAIR_GENOTYPES_CALLED,value.get(NUM_PAIR_GENOTYPES_CALLED)+sum.get(NUM_PAIR_GENOTYPES_CALLED));
        sum.put(NUM_PAIR_GENOTYPES_NOCALL,value.get(NUM_PAIR_GENOTYPES_NOCALL)+sum.get(NUM_PAIR_GENOTYPES_NOCALL));
        sum.put(NUM_PAIR_GENOTYPES_PHASED,value.get(NUM_PAIR_GENOTYPES_PHASED)+sum.get(NUM_PAIR_GENOTYPES_PHASED));
        sum.put(NUM_PAIR_HET_HET,value.get(NUM_PAIR_HET_HET)+sum.get(NUM_PAIR_HET_HET));
        sum.put(NUM_PAIR_VIOLATIONS,value.get(NUM_PAIR_VIOLATIONS)+sum.get(NUM_PAIR_VIOLATIONS));
        sum.put(NUM_TRIO_DOUBLE_VIOLATIONS,value.get(NUM_TRIO_DOUBLE_VIOLATIONS)+sum.get(NUM_TRIO_DOUBLE_VIOLATIONS));
        sum.put(NUM_GENOTYPES_MODIFIED,value.get(NUM_GENOTYPES_MODIFIED)+sum.get(NUM_GENOTYPES_MODIFIED));
        sum.put(NUM_MUTATIONS,value.get(NUM_MUTATIONS)+sum.get(NUM_MUTATIONS));
        sum.put(NUM_INCOMPLETE_TRIO_PHASED,value.get(NUM_INCOMPLETE_TRIO_PHASED)+sum.get(NUM_INCOMPLETE_TRIO_PHASED));
        sum.put(NUM_INCOMPLETE_TRIO_VIOLATIONS,value.get(NUM_INCOMPLETE_TRIO_VIOLATIONS)+sum.get(NUM_INCOMPLETE_TRIO_VIOLATIONS));
        sum.put(NUM_TRIO_SINGLES_PHASED,value.get(NUM_TRIO_SINGLES_PHASED)+sum.get(NUM_TRIO_SINGLES_PHASED));
        sum.put(NUM_PAIR_SINGLES_PHASED,value.get(NUM_PAIR_SINGLES_PHASED)+sum.get(NUM_PAIR_SINGLES_PHASED));
        return sum;
    }
	
	/**
     * Reports statistics on the phasing by transmission process.
     * @param result Accumulator with all counters.
     */
    @Override
    public void onTraversalDone(HashMap<Byte,Double> result) 
	{
        //Only report trios if any are found.
        if(result.get(NUM_TRIO_GENOTYPES_CALLED) > 0 || result.get(NUM_TRIO_GENOTYPES_NOCALL) > 0)
		{
            logger.info("Number of complete trio-genotypes: " + result.get(NUM_TRIO_GENOTYPES_CALLED).intValue());
            logger.info("Number of trio-genotypes containing no call(s): " + result.get(NUM_TRIO_GENOTYPES_NOCALL).intValue());
            logger.info("Number of genotypes phased in complete trios: " + result.get(NUM_TRIO_GENOTYPES_PHASED).intValue());
            logger.info("Number of resulting Het/Het/Het trios: " + result.get(NUM_TRIO_HET_HET_HET).intValue());
            logger.info("Number of remaining single mendelian violations in complete trios: " + result.get(NUM_TRIO_VIOLATIONS).intValue());
            logger.info("Number of remaining double mendelian violations in complete trios: " + result.get(NUM_TRIO_DOUBLE_VIOLATIONS).intValue());
            logger.info("Number of incomplete (2 individuals called) genotypes phased in trios: " + result.get(NUM_INCOMPLETE_TRIO_PHASED).intValue());
            logger.info("Number of single individual homozygous genotypes phased in trios: " + result.get(NUM_TRIO_SINGLES_PHASED).intValue());
            logger.info("Number of remaining mendelian violations in incomplete (2 individuals called) trios: " + result.get(NUM_TRIO_DOUBLE_VIOLATIONS).intValue());
            logger.info("Estimated number of mutations in complete trios: " + result.get(NUM_MUTATIONS));
        }
        else
		{
            logger.info("No trios with genotypes on given sites found.");
        }
        //Report pairs if any found
        if(result.get(NUM_PAIR_GENOTYPES_CALLED) > 0 || result.get(NUM_PAIR_GENOTYPES_NOCALL) > 0)
		{
            logger.info("Number of complete pair-genotypes: " + result.get(NUM_PAIR_GENOTYPES_CALLED).intValue());
            logger.info("Number of pair-genotypes containing no call(s): " + result.get(NUM_PAIR_GENOTYPES_NOCALL).intValue());
            logger.info("Number of genotypes phased in pairs: " + result.get(NUM_PAIR_GENOTYPES_PHASED).intValue());
            logger.info("Number of resulting Het/Het pairs: " + result.get(NUM_PAIR_HET_HET).intValue());
            logger.info("Number of remaining mendelian violations in pairs: " + result.get(NUM_PAIR_VIOLATIONS).intValue());
            logger.info("Number of single individual homozygous genotypes phased in pairs: " + result.get(NUM_PAIR_SINGLES_PHASED).intValue());
        }
        else
		{
            logger.info("No parent/child pairs with genotypes on given sites found.");
        }
        logger.info("Number of genotypes updated: " + result.get(NUM_GENOTYPES_MODIFIED).intValue());
        logger.info("Number of sites where no AF tag was found in the supplied VCF: " + no_af_found);
        logger.info("Number of DNM-called sites where no AF tag was found in the supplied VCF: " + no_af_dnm_call);
        
        
		
    }
	
	/**
	 * Returns an array of symbolic alleles for internal, temporary, use
	 * @param g_type: the sample genotype to build alleles for
	 * @param sex: the sex of the sample having g_type; used to determine ploidy
	 * @return an ArrayList with the inferred alleles
    */
	protected ArrayList<Allele> getAlleles(GenotypeType g_type, Gender sex)
	{
		ArrayList<Allele> ret = new ArrayList<Allele>(2);
		ret.clear();
		switch (g_type)
		{
		case HOM_REF:
			if (site_male_ploidy == 1 && sex == Gender.MALE )
				ret.add(dummy_ref);
			else
			{
				ret.add(dummy_ref);
				ret.add(dummy_ref);
			}
			break;
		case HOM_VAR:
			if (site_male_ploidy == 1 && sex == Gender.MALE )
				ret.add(dummy_alt);
			else
			{
				ret.add(dummy_alt);
				ret.add(dummy_alt);
			}
			break;
		case HET:
			if (site_male_ploidy == 1 && sex == Gender.MALE )
				return null;
			else
			{
				ret.add(dummy_ref);
				ret.add(dummy_alt);
			}
			break;
		default:
			return null;
		}
		
		//return null if no alleles could be added (i.e.: GT = UNAVAILABLE, NO_CALL, MIXED - we cannot use these GTs further anyway as they have no PLs )
		return ret.size() > 0? ret:null;
	}
	
	// determine based on genotype value and ploidy whether a GT can be phased
	protected boolean isPhasable(GenotypeType g_type, Gender sex)
	{
		return ( (site_male_ploidy == 1 && sex == Gender.MALE && (g_type == GenotypeType.HOM_REF || g_type == GenotypeType.HOM_VAR)) ||
				 ((site_male_ploidy == 2 || sex == Gender.FEMALE) && (g_type == GenotypeType.HOM_REF || g_type == GenotypeType.HOM_VAR || g_type == GenotypeType.HET)) );
	}
    
	//Given 3 symbolic GT types compute the number of mendelian violations present
	private int getConfigurationMVCount(GenotypeType mother, GenotypeType father, GenotypeType child, Gender sex)
    {

        int mv_count = 0;
        //get alleles for comparison
        ArrayList<Allele> als_father = getAlleles(father, Gender.MALE);
        ArrayList<Allele> als_mother = getAlleles(mother, Gender.FEMALE);
        ArrayList<Allele> als_child  = getAlleles(child, sex);
        
        // if child has a non-usable genotype (i.e.: UNAVAILABLE,MIXED,NO_CALL) => can't count
        if (als_child == null)
        	return 0;
        // if both parents have non-usable genotypes => can't count
        if (als_father == null && als_mother == null)
        	return 0;
        
        // count MVs
        if (site_male_ploidy == 2 || sex == Gender.FEMALE)
        	// check if both child alleles can be found, in either parent
        	mv_count = Math.min((als_father == null? 0:als_father.indexOf(als_child.get(0)) >= 0? 0:1) + (als_mother == null? 0:als_mother.indexOf(als_child.get(1)) >= 0? 0:1),
        			(als_father == null? 0:als_father.indexOf(als_child.get(1)) >= 0? 0:1) + (als_mother == null? 0:als_mother.indexOf(als_child.get(0)) >= 0? 0:1));
        	
        else
        	// male X-inheritance: check if child allele can be found in mother
        	mv_count = als_mother == null? 0:als_mother.indexOf(als_child.get(0)) >= 0? 0:1;
        
        return mv_count;

    }
    
	// retrieve the PL values as a map(i.e.: w.r.t GT type); wrapper function to treat missing GTs and ploidy 1
	private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(Genotype genotype, Gender sex)
    {
        if(genotype == null || !genotype.isCalled())
        {
        	if (site_male_ploidy == 1 && sex == Gender.MALE)
        	{
        		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        		likelihoods.put(GenotypeType.HOM_REF,log10_1_2);
        		likelihoods.put(GenotypeType.HOM_VAR,log10_1_2);
        		return likelihoods;
        		
        	}
        	else
        	{
        		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        		likelihoods.put(GenotypeType.HOM_REF,log10_1_3);
        		likelihoods.put(GenotypeType.HET,log10_1_3);
        		likelihoods.put(GenotypeType.HOM_VAR,log10_1_3);
        		return likelihoods;
        	}
        }
		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
		
		// check for ploidy 1 manually -- code for Genotype class not visible => not sure the getAsMap works well for ploidy 1 
		double []pls = genotype.getLikelihoods().getAsVector();
		if (pls.length == 2)
		{
			likelihoods.put(GenotypeType.HOM_REF, pls[0]);
			likelihoods.put(GenotypeType.HOM_VAR, pls[1]);
			return likelihoods;
		}
		else
			return genotype.getLikelihoods().getAsMap(false);
    }
	
	
	// retrieve the PL values as a map(i.e.: w.r.t GT type); wrapper function to treat missing GTs and ploidy 1
	private EnumMap<GenotypeType,Double> getLikelihoodsAsMapSafeNull(GenotypeLikelihoods gl)
    {
        if(gl == null || gl.getAsPLs().length != 2 || gl.getAsPLs().length != 3 )
        {
        	int ploidy = gl.getAsPLs().length == 2? 1:2; 
        	if (ploidy == 1)
        	{
        		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        		likelihoods.put(GenotypeType.HOM_REF,log10_1_2);
        		likelihoods.put(GenotypeType.HOM_VAR,log10_1_2);
        		return likelihoods;
        		
        	}
        	else
        	{
        		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        		likelihoods.put(GenotypeType.HOM_REF,log10_1_3);
        		likelihoods.put(GenotypeType.HET,log10_1_3);
        		likelihoods.put(GenotypeType.HOM_VAR,log10_1_3);
        		return likelihoods;
        	}
        }
		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
		
		// check for ploidy 1 manually -- code for Genotype class not visible => not sure the getAsMap works well for ploidy 1 
		double []pls = gl.getAsVector();
		if (pls.length == 2)
		{
			likelihoods.put(GenotypeType.HOM_REF, pls[0]);
			likelihoods.put(GenotypeType.HOM_VAR, pls[1]);
			return likelihoods;
		}
		else
			return gl.getAsMap(false);
    }
	
	// retrieve the PL values as a map(i.e.: w.r.t GT type) in linear space; wrapper function to treat missing GTs and ploidy 1
	private EnumMap<GenotypeType,Double> getLinearLikelihoodsAsMapSafeNull(Genotype genotype, Gender sex)
    {
        if(genotype == null || !genotype.isCalled())
        {
        	if (site_male_ploidy == 1 && sex == Gender.MALE)
        	{
        		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        		likelihoods.put(GenotypeType.HOM_REF,0.5);
        		likelihoods.put(GenotypeType.HOM_VAR,0.5);
        		return likelihoods;
        		
        	}
        	else
        	{
        		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        		likelihoods.put(GenotypeType.HOM_REF,1.0/3.0);
        		likelihoods.put(GenotypeType.HET,1.0/3.0);
        		likelihoods.put(GenotypeType.HOM_VAR,1.0/3.0);
        		return likelihoods;
        	}
        }
		EnumMap<GenotypeType,Double> likelihoods = new EnumMap<GenotypeType, Double>(GenotypeType.class);
		double []pls = genotype.getLikelihoods().getAsVector();
		if (pls.length == 2)
		{
			double []pls_normalized = MathUtils.normalizeFromLog10(pls);
			likelihoods.put(GenotypeType.HOM_REF, pls_normalized[0]);
			likelihoods.put(GenotypeType.HOM_VAR, pls_normalized[1]);
			return likelihoods;
		}
		else
			return genotype.getLikelihoods().getAsMap(true);
    }

	// looks ugly(i.e.: 200 lines of ONE if statement); does the job; mostly for internal use for now; not sure how relevant it would be for a general context
	private void phaseNonBiallelicSites (VariantContext vc)
	{
		if (vc == null)
			return;

		VariantContextBuilder builder = new VariantContextBuilder(vc);

		GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());
		for (Sample sample : trios)
		{
			Genotype mother = vc.getGenotype(sample.getMaternalID());
			Genotype father = vc.getGenotype(sample.getPaternalID());
			Genotype child = vc.getGenotype(sample.getID());

			if (mother == null || father == null || child == null)
				continue;


			List<Allele> siteAls = vc.getAlleles();
			List<Allele> motherAls = mother.getAlleles();
			List<Allele> fatherAls = father.getAlleles();
			List<Allele> childAls  = child.getAlleles();
			if (child.sameGenotype(mother, true) && child.sameGenotype(father, true) && child.isHet()) // the only unphasable combination
				continue;

			Genotype new_kid, new_mother, new_father;
			if ( childAls.get(0).equals(motherAls.get(0)) )
				if ( childAls.get(1).equals(fatherAls.get(0)) )
				{
					if (fatherAlleleFirst)
					{
						new_kid = flipAndPhaseGenotype(child,true,true);
						new_mother = flipAndPhaseGenotype(mother,false,true);
						new_father = flipAndPhaseGenotype(father,false,true);
					}
					else
					{
						new_kid = flipAndPhaseGenotype(child,false,true);
						new_mother = flipAndPhaseGenotype(mother,false,true);
						new_father = flipAndPhaseGenotype(father,false,true);
					}
				}
				else
				{
					if (childAls.get(1).equals(fatherAls.get(1)))
					{
						if (fatherAlleleFirst)
						{
							new_kid = flipAndPhaseGenotype(child,true,true);
							new_mother = flipAndPhaseGenotype(mother,false,true);
							new_father = flipAndPhaseGenotype(father,true,true);
						}
						else
						{
							new_kid = flipAndPhaseGenotype(child,false,true);
							new_mother = flipAndPhaseGenotype(mother,false,true);
							new_father = flipAndPhaseGenotype(father,true,true);
						}
					}
					else
						continue; // MV or otherwise inconsistent combination of GTs
				}
			else
				if ( childAls.get(1).equals(motherAls.get(0)) )
					if ( childAls.get(0).equals(fatherAls.get(0)) )
					{
						if (fatherAlleleFirst)
						{
							new_kid = flipAndPhaseGenotype(child,false,true);
							new_mother = flipAndPhaseGenotype(mother,false,true);
							new_father = flipAndPhaseGenotype(father,false,true);
						}
						else
						{
							new_kid = flipAndPhaseGenotype(child,true,true);
							new_mother = flipAndPhaseGenotype(mother,false,true);
							new_father = flipAndPhaseGenotype(father,false,true);
						}
					}
					else
					{
						if (childAls.get(0).equals(fatherAls.get(1)))
						{
							if (fatherAlleleFirst)
							{
								new_kid = flipAndPhaseGenotype(child,false,true);
								new_mother = flipAndPhaseGenotype(mother,false,true);
								new_father = flipAndPhaseGenotype(father,true,true);
							}
							else
							{
								new_kid = flipAndPhaseGenotype(child,true,true);
								new_mother = flipAndPhaseGenotype(mother,false,true);
								new_father = flipAndPhaseGenotype(father,true,true);
							}
						}
						else
							continue; // MV or otherwise inconsistent combination of GTs
					}
				else
					if ( childAls.get(0).equals(motherAls.get(1)) )
						if ( childAls.get(1).equals(fatherAls.get(0)) )
						{
							if (fatherAlleleFirst)
							{
								new_kid = flipAndPhaseGenotype(child,true,true);
								new_mother = flipAndPhaseGenotype(mother,true,true);
								new_father = flipAndPhaseGenotype(father,false,true);
							}
							else
							{
								new_kid = flipAndPhaseGenotype(child,false,true);
								new_mother = flipAndPhaseGenotype(mother,true,true);
								new_father = flipAndPhaseGenotype(father,false,true);
							}
						}
						else
						{
							if (childAls.get(1).equals(fatherAls.get(1)))
							{
								if (fatherAlleleFirst)
								{
									new_kid = flipAndPhaseGenotype(child,true,true);
									new_mother = flipAndPhaseGenotype(mother,true,true);
									new_father = flipAndPhaseGenotype(father,true,true);
								}
								else
								{
									new_kid = flipAndPhaseGenotype(child,false,true);
									new_mother = flipAndPhaseGenotype(mother,true,true);
									new_father = flipAndPhaseGenotype(father,true,true);
								}
							}
							else
								continue; // MV or otherwise inconsistent combination of GTs
						}
					else
						if ( childAls.get(1).equals(motherAls.get(1)) )
							if ( childAls.get(0).equals(fatherAls.get(0)) )
							{
								if (fatherAlleleFirst)
								{
									new_kid = flipAndPhaseGenotype(child,false,true);
									new_mother = flipAndPhaseGenotype(mother,true,true);
									new_father = flipAndPhaseGenotype(father,false,true);
								}
								else
								{
									new_kid = flipAndPhaseGenotype(child,true,true);
									new_mother = flipAndPhaseGenotype(mother,true,true);
									new_father = flipAndPhaseGenotype(father,false,true);
								}
							}
							else
							{
								if (childAls.get(0).equals(fatherAls.get(1)))
								{
									if (fatherAlleleFirst)
									{
										new_kid = flipAndPhaseGenotype(child,false,true);
										new_mother = flipAndPhaseGenotype(mother,true,true);
										new_father = flipAndPhaseGenotype(father,true,true);
									}
									else
									{
										new_kid = flipAndPhaseGenotype(child,true,true);
										new_mother = flipAndPhaseGenotype(mother,true,true);
										new_father = flipAndPhaseGenotype(father,true,true);
									}
								}
								else
									continue; // MV or otherwise inconsistent combination of GTs
							}
						else
							continue; // MV or otherwise inconsistent combination of GTs
			genotypesContext.replace(new_mother);
			genotypesContext.replace(new_father);
			genotypesContext.replace(new_kid);
		}
		builder.genotypes(genotypesContext);
		vcfWriter.add(builder.make());

	}

	private Genotype flipAndPhaseGenotype(Genotype gt, boolean flip, boolean phase)
	{
		GenotypeBuilder gb = new GenotypeBuilder();

		List<Allele> old_alleles = gt.getAlleles();
		List<Allele> new_alleles = new ArrayList<Allele>(2);
		if ( flip )
		{
			new_alleles.add(old_alleles.get(1));
			new_alleles.add(old_alleles.get(0));
		}
		else
		{
			new_alleles.add(old_alleles.get(0));
			new_alleles.add(old_alleles.get(1));
		}
		gb.alleles(new_alleles);
		if ( phase )
			gb.phased(true);
		else
			gb.phased(false);

		gb.AD(gt.getAD());
		gb.DP(gt.getDP());
		gb.GQ(gt.getGQ());
		gb.PL(gt.getPL());
		gb.name(gt.getSampleName());

		HashMap<String,Object> new_attributes = new HashMap<>();
		new_attributes.putAll(gt.getExtendedAttributes());
		if (!noTP)
			new_attributes.put(TRANSMISSION_PROBABILITY_TAG_NAME, new String("."));
		gb.attributes(new_attributes);

		return gb.make();


	}

	// update output statistics
	private void updatePairMetricsCounters(TrioConfiguration trioGenotypes, HashMap<Byte,Double> counters, EnumMap<GenotypeType,Double> noDnmGtPriors)
	{
        Genotype parent;
        if (trioGenotypes.isMV())
        	if (useAF == EXTERNAL && noDnmGtPriors.equals(neutralAF))
     		   no_af_dnm_call += 1;
        if(trioGenotypes.getPhasedMother() == null)
            parent = trioGenotypes.getPhasedFather();
        else
            parent = trioGenotypes.getPhasedMother();

        //Increment metrics counters
        if(parent.isCalled() && trioGenotypes.getPhasedChild().isCalled()){
           counters.put(NUM_PAIR_GENOTYPES_CALLED,counters.get(NUM_PAIR_GENOTYPES_CALLED)+1);
            //Note that both parent and child should be either phased or unphased
            if(parent.isPhased())
               counters.put(NUM_PAIR_GENOTYPES_PHASED,counters.get(NUM_PAIR_GENOTYPES_PHASED)+1);
           if(trioGenotypes.getPhasedChild().isPhased())
               counters.put(NUM_PAIR_GENOTYPES_PHASED,counters.get(NUM_PAIR_GENOTYPES_PHASED)+1);

            counters.put(NUM_PAIR_VIOLATIONS,counters.get(NUM_PAIR_VIOLATIONS)+trioGenotypes.getMVCount());
           if(parent.isHet() && trioGenotypes.getPhasedChild().isHet())
               counters.put(NUM_PAIR_HET_HET,counters.get(NUM_PAIR_HET_HET)+1);
        }else{
            counters.put(NUM_PAIR_GENOTYPES_NOCALL,counters.get(NUM_PAIR_GENOTYPES_NOCALL)+1);
            if(parent.isPhased())
                counters.put(NUM_PAIR_SINGLES_PHASED,counters.get(NUM_PAIR_SINGLES_PHASED)+1);
            if(trioGenotypes.getPhasedChild().isPhased())
                counters.put(NUM_PAIR_SINGLES_PHASED,counters.get(NUM_PAIR_SINGLES_PHASED)+1);
        }

    }
    
	//update output statistics
    private void updateTrioMetricsCounters(TrioConfiguration trioGenotypes, HashMap<Byte,Double> counters, EnumMap<GenotypeType,Double> noDnmGtPriors)
    {

        //Increment metrics counters
        if(trioGenotypes.getPhasedMother().isCalled() && trioGenotypes.getPhasedFather().isCalled() && trioGenotypes.getPhasedChild().isCalled()){
           counters.put(NUM_TRIO_GENOTYPES_CALLED,counters.get(NUM_TRIO_GENOTYPES_CALLED)+1);
           counters.put(NUM_MUTATIONS,counters.get(NUM_MUTATIONS)+trioGenotypes.getMVCount());
            //Note that normally they should either be all phased or not
           if(trioGenotypes.getPhasedMother().isPhased())
               counters.put(NUM_TRIO_GENOTYPES_PHASED,counters.get(NUM_TRIO_GENOTYPES_PHASED)+1);
           if(trioGenotypes.getPhasedFather().isPhased())
               counters.put(NUM_TRIO_GENOTYPES_PHASED,counters.get(NUM_TRIO_GENOTYPES_PHASED)+1);
           if(trioGenotypes.getPhasedChild().isPhased())
               counters.put(NUM_TRIO_GENOTYPES_PHASED,counters.get(NUM_TRIO_GENOTYPES_PHASED)+1);

           if(trioGenotypes.isMV()){
        	   if (useAF == EXTERNAL && noDnmGtPriors.equals(neutralAF))
        		   no_af_dnm_call += 1;
               if(trioGenotypes.getMVCount() > 1)
                    counters.put(NUM_TRIO_DOUBLE_VIOLATIONS,counters.get(NUM_TRIO_DOUBLE_VIOLATIONS)+1);
               else
                   counters.put(NUM_TRIO_VIOLATIONS,counters.get(NUM_TRIO_VIOLATIONS)+1);
           }
           else if(trioGenotypes.getPhasedMother().isHet() && trioGenotypes.getPhasedFather().isHet() && trioGenotypes.getPhasedChild().isHet())
               counters.put(NUM_TRIO_HET_HET_HET,counters.get(NUM_TRIO_HET_HET_HET)+1);

        }
        else{
            counters.put(NUM_TRIO_GENOTYPES_NOCALL,counters.get(NUM_TRIO_GENOTYPES_NOCALL)+1);
            if(trioGenotypes.getPhasedChild().isCalled() && (trioGenotypes.getPhasedFather().isCalled() || trioGenotypes.getPhasedMother().isCalled())){
                if(trioGenotypes.getPhasedChild().isPhased())
                    counters.put(NUM_INCOMPLETE_TRIO_PHASED,counters.get(NUM_INCOMPLETE_TRIO_PHASED)+1);
                if(trioGenotypes.getPhasedFather().isPhased())
                    counters.put(NUM_INCOMPLETE_TRIO_PHASED,counters.get(NUM_INCOMPLETE_TRIO_PHASED)+1);
                if(trioGenotypes.getPhasedMother().isPhased())
                    counters.put(NUM_INCOMPLETE_TRIO_PHASED,counters.get(NUM_INCOMPLETE_TRIO_PHASED)+1);
                if(trioGenotypes.isMV())
                    counters.put(NUM_INCOMPLETE_TRIO_VIOLATIONS,counters.get(NUM_INCOMPLETE_TRIO_VIOLATIONS)+1);
            }
            else if(trioGenotypes.getPhasedMother().isPhased()){
                    counters.put(NUM_TRIO_SINGLES_PHASED,counters.get(NUM_TRIO_SINGLES_PHASED)+1);
            }
            else if(trioGenotypes.getPhasedFather().isPhased()){
                    counters.put(NUM_TRIO_SINGLES_PHASED,counters.get(NUM_TRIO_SINGLES_PHASED)+1);
            }
            else if(trioGenotypes.getPhasedChild().isPhased()){
                    counters.put(NUM_TRIO_SINGLES_PHASED,counters.get(NUM_TRIO_SINGLES_PHASED)+1);
            }

        }
    }
	
	/**
	 * Effectively runs PBT on each trio at a locus
	 * @param vc: VariantContext with all the samples to be considered
	 * @param noDnmGtPriors: AF prior computed from the supplied population
	 * @return the output statistics computed for this locus
	 */
    private HashMap<Byte,Double> parseGenotypesAtLocus(VariantContext vc, EnumMap<GenotypeType,Double> noDnmGtPriors, GenomeLoc locus)
	{
		HashMap<Byte,Double> metricsCounters = new HashMap<Byte, Double>(10);
        metricsCounters.put(NUM_TRIO_GENOTYPES_CALLED,0.0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_NOCALL,0.0);
        metricsCounters.put(NUM_TRIO_GENOTYPES_PHASED,0.0);
        metricsCounters.put(NUM_TRIO_HET_HET_HET,0.0);
        metricsCounters.put(NUM_TRIO_VIOLATIONS,0.0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_CALLED,0.0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_NOCALL,0.0);
        metricsCounters.put(NUM_PAIR_GENOTYPES_PHASED,0.0);
        metricsCounters.put(NUM_PAIR_HET_HET,0.0);
        metricsCounters.put(NUM_PAIR_VIOLATIONS,0.0);
        metricsCounters.put(NUM_TRIO_DOUBLE_VIOLATIONS,0.0);
        metricsCounters.put(NUM_GENOTYPES_MODIFIED,0.0);
        metricsCounters.put(NUM_MUTATIONS,0.0);
        metricsCounters.put(NUM_INCOMPLETE_TRIO_PHASED,0.0);
        metricsCounters.put(NUM_INCOMPLETE_TRIO_VIOLATIONS,0.0);
        metricsCounters.put(NUM_TRIO_SINGLES_PHASED,0.0);
        metricsCounters.put(NUM_PAIR_SINGLES_PHASED,0.0);

        if (vc == null)
        	return metricsCounters;
				
		String mvfLine = "";
		VariantContextBuilder builder = new VariantContextBuilder(vc);

        GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());
		for (Sample sample : trios) 
		{
			mvfLine="";
            Genotype mother = vc.getGenotype(sample.getMaternalID());
            Genotype father = vc.getGenotype(sample.getPaternalID());
            Genotype child = vc.getGenotype(sample.getID());

            //Treat only trios and parent/child pairs
            if(mother == null && father == null || child == null)
                continue;
            
            TrioConfiguration trioGenotypes;
            trioGenotypes = phaseTrioGenotypes(vc.getReference(), vc.getAltAlleleWithHighestAlleleCount(), mother, father, child, noDnmGtPriors, sample.getGender());
            
            
            Genotype phasedMother = trioGenotypes.getPhasedMother();
            Genotype phasedFather = trioGenotypes.getPhasedFather();
            Genotype phasedChild = trioGenotypes.getPhasedChild();
            
            
            /*
            if (child.hasAnyAttribute(ReadBackedPhasing.HP_KEY))
            {
            		System.out.printf("TEST: RBP tag for child: %s\n", child.getAnyAttribute(ReadBackedPhasing.HP_KEY).toString());
            }
            else
            {
            	System.out.println("TEST: no RBP tag");
            }
            */
            if ( assignPO == 0 )
            {
            	String parentalAssignment = "unknown";
            	if ( mother != null && father != null )
            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY), makeADString(phasedMother), phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
            	else if ( father == null && mother != null )
            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedMother),phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
            	else // ( mother == null && father != null )
            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t.\t.\t.\t.\t.\t\t.%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedFather.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,phasedChild.getGenotypeString(),sample.getID(),phasedChild.getType(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
            	
            	//Report violation if set so
            	if( mother != null && father != null && mother.isCalled() && father.isCalled() && trioGenotypes.isMV() && mvFile != null )
                    mvFile.println(mvfLine);
            }
            else
            {
	            if ( vc.getAttributes().containsKey("PhasingInconsistent") )
	            {
	            	String parentalAssignment = "PhasingUnreliable";
	            	if ( mother != null && father != null )
	            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY), makeADString(phasedMother), phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
	            	else if ( father == null && mother != null )
	            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedMother),phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
	            	else // ( mother == null && father != null )
	            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t.\t.\t.\t.\t.\t\t.%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedFather.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,phasedChild.getGenotypeString(),sample.getID(),phasedChild.getType(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
	            	
	            	//Report violation if set so
	            	if( mother != null && father != null && mother.isCalled() && father.isCalled() && trioGenotypes.isMV() && mvFile != null )
	                    mvFile.println(mvfLine);
	            }
	            else
	            {
		            int POaction = handlePOInformation(trioGenotypes, sample, locus, mother, father, vc.getFilters(), vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY).toString());
		            if ( POaction == -1 )
		            {
		            	String parentalAssignment = "unknown";
		            	if ( mother != null && father != null )
		            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY), makeADString(phasedMother), phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
		            	else if ( father == null && mother != null )
		            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedMother),phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
		            	else // ( mother == null && father != null )
		            		mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t-1\t-1\t.\t.\t.\t.\t.\t\t.%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedFather.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),parentalAssignment,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,phasedChild.getGenotypeString(),sample.getID(),phasedChild.getType(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
		            	
		            	//Report violation if set so
		            	if( mother != null && father != null && mother.isCalled() && father.isCalled() && trioGenotypes.isMV() && mvFile != null )
		                    mvFile.println(mvfLine);
		                
		            }
		            else if ( POaction == 1 )
		            {
		            	mvfLine = buildMvfLinesPo(sample, vc); // may actually be more lines
		            	if ( mvFile != null && mvfLine != null )
		                    mvFile.println(mvfLine);
		            }
	            }
            }
            //if( mother != null && father != null && mother.isCalled() && father.isCalled() && trioGenotypes.isMV() && mvFile != null )
            //    mvFile.println(mvfLine);
            
            //Fill the genotype map with the new genotypes and increment metrics counters
            //genotypesContext.replace(phasedChild);
            //System.out.printf("mother DP:%s\tmother AD:%s, %d,%d\n", mother.getAnyAttribute("DP"), mother.getAttributeAsString("AD","mata"), mother.getAD()[0],mother.getAD()[1]);
            //System.out.printf("mother DP:%s\tmother AD:%s, %d,%d\n\n", phasedMother.getAnyAttribute("DP"), phasedMother.getAttributeAsString("AD", "mata"), phasedMother.getAD()[0], phasedMother.getAD()[1]);
            if(mother != null)
            {
                genotypesContext.replace(phasedMother);
                if(father != null)
                {
                	//genotypesContext.replace(phasedFather);
                    updateTrioMetricsCounters(trioGenotypes,metricsCounters,noDnmGtPriors);
                    //int []adm = phasedMother.getAD();
                    //if (adm.length == 2)
                    //	mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t% s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY),adm[0], adm[1],phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),phasedFather.getAnyAttribute("AD"),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),phasedChild.getAnyAttribute("AD"),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
                    //else
                    	
                    //mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY), makeADString(phasedMother), phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
                    
                    if( !(phasedMother.getType()==mother.getType() && phasedFather.getType()==father.getType() && phasedChild.getType()==child.getType()) )
                        metricsCounters.put(NUM_GENOTYPES_MODIFIED,metricsCounters.get(NUM_GENOTYPES_MODIFIED)+1);
                    if (fake_hets != 0 && phasedFather.getPloidy() == 1)
                    	genotypesContext.replace(haploidToFakeDiploid(phasedFather, -800));
                    else
                    	genotypesContext.replace(phasedFather);
                }
                else
                {
                    if( !vc.isFiltered() )
                        updatePairMetricsCounters(trioGenotypes,metricsCounters,noDnmGtPriors);
                    if( !(phasedMother.getType()==mother.getType() && phasedChild.getType()==child.getType()) )
                        metricsCounters.put(NUM_GENOTYPES_MODIFIED,metricsCounters.get(NUM_GENOTYPES_MODIFIED)+1);
                    //mvfLine = String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedMother),phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
                }
            }
            else
			{
                //genotypesContext.replace(phasedFather);
                 if( !vc.isFiltered() )
                    updatePairMetricsCounters(trioGenotypes,metricsCounters,noDnmGtPriors);
                if(!(phasedFather.getType()==father.getType() && phasedChild.getType()==child.getType()))
                    metricsCounters.put(NUM_GENOTYPES_MODIFIED,metricsCounters.get(NUM_GENOTYPES_MODIFIED)+1);
                //mvfLine =   String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t\t.%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",vc.getChr(),vc.getStart(),vc.getFilters(),vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY),phasedFather.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),trioGenotypes.mvCount,sample.getFamilyID(),sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,phasedChild.getGenotypeString(),sample.getID(),phasedChild.getType(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
                if (fake_hets != 0 && phasedFather.getPloidy() == 1)
                	genotypesContext.replace(haploidToFakeDiploid(phasedFather, -800));
                else
                	genotypesContext.replace(phasedFather);
            }
            if (fake_hets != 0 && phasedChild.getPloidy() == 1)
            	genotypesContext.replace(haploidToFakeDiploid(phasedChild, -800));
            else
            	genotypesContext.replace(phasedChild);
            

            //Report violation if set so
            //if(mother != null && father != null && mother.isCalled() && father.isCalled() && trioGenotypes.isMV() && mvFile != null)
            //    mvFile.println(mvfLine);
            }
		/*
		if (fake_hets != 0)
		{
			genotypesContext = haploidToFakeDiploid(genotypesContext, -800);
		}
		*/
		builder.genotypes(genotypesContext);
        vcfWriter.add(builder.make());

        return metricsCounters;
	}

    
    private GenotypesContext haploidToFakeDiploid(GenotypesContext _genotypesContext, int fake_hets_val)
    {
    	GenotypesContext ret_gts = GenotypesContext.copy(_genotypesContext);
    	
    	for (Sample sample : trios)
    	{
    		//Genotype mother = vc.getGenotype(sample.getMaternalID());
            Genotype father =  _genotypesContext.get(sample.getPaternalID());
            Genotype child 	=  _genotypesContext.get(sample.getID());
            
            
            if (father.getPloidy() == 1)
            	ret_gts.replace(haploidToFakeDiploid(father, fake_hets_val));
            if (child.getPloidy() == 1)
            	ret_gts.replace(haploidToFakeDiploid(child, fake_hets_val));	
    	}
    	return ret_gts;
    		
    }
    /**
     * Transforms a haploid called GT to its diploid version by putting dummy values for all het-informative fields:
     * 	- transforms the alleles and PLs of the supplied GT (i.e.: does not change GQ value)
     * 	- currently designed only for PBT use (i.e.: behaviour hard-coded for bi-allelic GTs)
     * @param gt: haploid GT to modify
     * @return diploid GT with an "irrelevantly" large PL value for the het value
     */
    private Genotype haploidToFakeDiploid(Genotype gt, int fake_hets_val)
    {
    	//if (!(gt != null && gt.getPloidy() == 1) ) // to use the short-circuiting (I think..)
    	if ( gt == null || gt.getPloidy() != 1 )  
    		return gt;
    	
    	GenotypeBuilder JB = new GenotypeBuilder(gt);
    	
    	// double the allele found
    	List<Allele> _als = gt.getAlleles();
    	_als.add(_als.get(0));
    	JB.alleles(_als);
    	
    	// add the fake het PL
    	int []pls = gt.getPL();
    	int []_pls = new int[3];
    	if (pls == null)
    		JB.noPL();
    	else
    		if (pls.length == 2)
    		{
    			_pls[0] = pls[0]; 
    			_pls[1] = -10 * fake_hets_val; 
    			_pls[2] = pls[1];
    			JB.PL(_pls);
    		}
    		else
    			if (pls.length == 0)
    				JB.noPL();
    			else
    				throw new UserException(String.format("Error: expected haploid bi-allelic genotype, PL field found to have: %d values\n", pls.length ));
    				
    				
    	return JB.make();    	
    }
	
    private String buildMvfLinesPo(Sample _sample, VariantContext vc)
    {
    	//System.out.printf( "DBG: analyzing RBP Haplotype %s of length %d\n", rbp_queue.get(_sample.getID()).trioGTs.get(0) != null ? rbp_queue.get(_sample.getID()).trioGTs.get(0).getPhasedChild().getAnyAttribute(RBP_TAG).toString() : "null", rbp_queue.get(_sample.getID()).trioGTs.size());
    	
    	TrioRbpContext oneHapPlus = rbp_queue.get(_sample.getID());
    	//TrioRbpContext oneHap = 
    	//if ( !oneHapPlus.hasDnms )
    	//	{ System.out.println("buildMvfLines; retrun 1"); return null;}
    	if ( oneHapPlus.trioGTs.size() == 0 )
    	//{ System.out.println("buildMvfLines; retrun 2"); return null;}
    		return null;
    	ArrayList<Set<String>> _filters = new ArrayList<Set<String>>();
		ArrayList<String> _ADs = new ArrayList<String>();
		ArrayList<GenomeLoc> _gtLoci = new ArrayList<GenomeLoc> ();
		ArrayList<TrioConfiguration> _trioGTs = new ArrayList<TrioConfiguration>();
		// for each locus: 0=no missing parent; 1=father genotype is missing; 2=mother genotype is missing
		ArrayList<Integer> _missingParents = new ArrayList<Integer>();
		long _phaseRefPos = Integer.parseInt(oneHapPlus.trioGTs.get(0).getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[0].split("-")[0]);
		boolean _hasDnms = oneHapPlus.hasDnms;
    	
    	for (;;)
    	{
    		if ( oneHapPlus.trioGTs.size() == 0 )
    			break;
    		TrioConfiguration itTG = oneHapPlus.trioGTs.remove(0);
    		//very non-safe line and subsequent operations but everything was checked before
    		long itPhasingRef = Integer.parseInt(itTG.getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[0].split("-")[0]);
    		if ( itPhasingRef != _phaseRefPos )
    		{
    			oneHapPlus.trioGTs.add(0,itTG);
    			oneHapPlus.phaseRefPos = itPhasingRef;
    			if (itTG.isMV())
    				oneHapPlus.hasDnms = true;
    			else
    				oneHapPlus.hasDnms = false;
    			break;
    		}
    		_trioGTs.add(itTG);
    		_missingParents.add(oneHapPlus.missingParents.remove(0));
    		_gtLoci.add(oneHapPlus.gtLoci.remove(0));
    		_ADs.add(oneHapPlus.ADs.remove(0));
    		_filters.add(oneHapPlus.filters.remove(0));
    	}
    	rbp_queue.put(_sample.getID(), oneHapPlus);
    	if (_trioGTs.size() == 0)
    	//{ System.out.println("buildMvfLines; retrun 3"); return null;}
    		return null;
    	
    	// return value
    	String _ret = new String("");
    	
    	ArrayList<Integer> dnmsIdx = new ArrayList<Integer>();
    	for (int i = 0; i < _trioGTs.size(); i++)
    		if ( _trioGTs.get(i).isMV() )
    			dnmsIdx.add(new Integer(i));
    	if ( dnmsIdx.size() == 0 )
    	//{ System.out.println("buildMvfLines; retrun 4"); return null;}
    		return null;
    		
    	int []maternal_evidences = new int[dnmsIdx.size()];
    	int []paternal_evidences = new int[dnmsIdx.size()];
    	
    	int father_haplotype = fatherAlleleFirst == true ? 0 : 1;
    	boolean phaseInconsistent = false;
    	int _first_non_dnm_idx = -1;
    	for (int i = 0; i < _trioGTs.size(); i++)
    	{
    		if ( !_trioGTs.get(i).isMV() )
    		{
    			_first_non_dnm_idx = i;
    			break;
    		}
    	}
    	for (int i = 0; i < _trioGTs.size(); i++)
    	{
    		if ( _first_non_dnm_idx == -1 ) // => only dnms in this phased hap => cannot assign PO to any of them
    			break;
    		TrioConfiguration itTG = _trioGTs.get(i);
    		if ( itTG.isMV() )
    			continue;
    		if ( _trioGTs.get(_first_non_dnm_idx).getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[0].split("-")[1].compareTo(itTG.getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[0].split("-")[1]) != 0 )
    		{
    			//if ( dnmsIdx.size() != 0 && _gtLoci.get(dnmsIdx.get(0)).getStart() == 79279844)
    			//{
    			//	System.out.printf("DBG investigation: the \"ref\" RBP info %s\t%s\n\t..and the evaluation RBP info: %s\t%s\n",_trioGTs.get(_first_non_dnm_idx).getPhasedChild().getAnyAttribute(RBP_TAG).toString(),_trioGTs.get(_first_non_dnm_idx).getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[0].split("-")[1],itTG.getPhasedChild().getAnyAttribute(RBP_TAG).toString(),itTG.getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[0].split("-")[1] );
    			//}
    			phaseInconsistent = true;
    			//break;
    		}
    			
    		//Genotype _childEvidence = itTG.getPhasedChild();
    		int _childEvidencePaternal = Integer.parseInt(itTG.getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[father_haplotype].split("-")[1]);
    		for ( int j = 0; j < dnmsIdx.size(); j++)
    		{
    			int mutant_allele = _trioGTs.get(dnmsIdx.get(j).intValue()).getMutantAlleleIdx();
    			if (mutant_allele < 0)
    				continue;
    			int mutant_hap = Integer.parseInt(_trioGTs.get(dnmsIdx.get(j)).getPhasedChild().getAnyAttribute(RBP_TAG).toString().split(",")[mutant_allele].split("-")[1]);
    			if ( mutant_hap == _childEvidencePaternal )
    				paternal_evidences[j] += 1;
    			else
    				maternal_evidences[j] += 1;
    			
    		}
    	}
    	
    	for (int i = 0; i < dnmsIdx.size(); i++)
    	{
    		Genotype mother = vc.getGenotype(_sample.getMaternalID());
            Genotype father = vc.getGenotype(_sample.getPaternalID());
            Genotype phasedMother = _trioGTs.get(dnmsIdx.get(i)).getPhasedMother();
            Genotype phasedFather = _trioGTs.get(dnmsIdx.get(i)).getPhasedFather();
            Genotype phasedChild  = _trioGTs.get(dnmsIdx.get(i)).getPhasedChild();
            
    		double paternal_evidence_ratio = (paternal_evidences[i] + maternal_evidences[i]) != 0 ? paternal_evidences[i]/(paternal_evidences[i] + maternal_evidences[i]) : 0;
    		String parentalAssignment = phaseInconsistent ? "phaseInconsistent" : (maternal_evidences[i] == 0 && paternal_evidences[i] == 0) ? "unknown" : paternal_evidence_ratio >= POEvidenceThreshold ? "paternal" : 1.0 - (paternal_evidence_ratio) >= POEvidenceThreshold ? "maternal" : "undecided"; 
    		
        	if ( mother != null && father != null )
        		_ret += String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",_gtLoci.get(dnmsIdx.get(i)).getContig(),_gtLoci.get(dnmsIdx.get(i)).getStart(),_filters.get(dnmsIdx.get(i)),_ADs.get(dnmsIdx.get(i)),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),_trioGTs.get(dnmsIdx.get(i)).mvCount,_sample.getFamilyID(),parentalAssignment,paternal_evidences[i],maternal_evidences[i],_sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY), makeADString(phasedMother), phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,_sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,_sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
        	else if ( father == null && mother != null )
        		_ret += String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",_gtLoci.get(dnmsIdx.get(i)).getContig(),_gtLoci.get(dnmsIdx.get(i)).getStart(),_filters.get(dnmsIdx.get(i)),_ADs.get(dnmsIdx.get(i)),phasedMother.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),_trioGTs.get(dnmsIdx.get(i)).mvCount,_sample.getFamilyID(),parentalAssignment,paternal_evidences[i],maternal_evidences[i],_sample.getMaternalID(),phasedMother.getType(),phasedMother.getGenotypeString(),phasedMother.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedMother),phasedMother.getLikelihoods() != null? phasedMother.getLikelihoods().toString():null,_sample.getID(),phasedChild.getType(),phasedChild.getGenotypeString(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
        	else // ( mother == null && father != null )
        		_ret += String.format("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t.\t.\t.\t.\t.\t\t.%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",_gtLoci.get(dnmsIdx.get(i)).getContig(),_gtLoci.get(dnmsIdx.get(i)).getStart(),_filters.get(dnmsIdx.get(i)),_ADs.get(dnmsIdx.get(i)),phasedFather.getAnyAttribute(TRANSMISSION_PROBABILITY_TAG_NAME),_trioGTs.get(dnmsIdx.get(i)).mvCount,_sample.getFamilyID(),parentalAssignment,paternal_evidences[i],maternal_evidences[i],_sample.getPaternalID(),phasedFather.getType(),phasedFather.getGenotypeString(),phasedFather.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedFather),phasedFather.getLikelihoods() != null? phasedFather.getLikelihoods().toString():null,phasedChild.getGenotypeString(),_sample.getID(),phasedChild.getType(),phasedChild.getAnyAttribute(VCFConstants.DEPTH_KEY),makeADString(phasedChild),phasedChild.getLikelihoods() != null? phasedChild.getLikelihoods().toString():null);
        	
    	}
    	return _ret.substring(0, _ret.length() - 1);
    }
    
    private int handlePOInformation(TrioConfiguration _trioGenotypes, Sample child, GenomeLoc _locus, Genotype _mother, Genotype _father, Set<String> _filters, String _ad)
    {
    	Genotype phasedMother = _trioGenotypes.getPhasedMother();
    	Genotype phasedFather = _trioGenotypes.getPhasedFather();
    	Genotype phasedChild  = _trioGenotypes.getPhasedChild();
    	
    	if ( !phasedChild.isHet() || ( !phasedChild.isPhased() && !_trioGenotypes.isMV() ) )
    		return -1;
    	TrioRbpContext currTrio = rbp_queue.get(child.getID());
    	//System.out.printf("DBG: just RBP tag: %s\n",phasedChild.getAnyAttribute(RBP_TAG) != null ? phasedChild.getAnyAttribute(RBP_TAG).toString() : "null");
    	// start of a new RBP-haplotype in a child
    	if (currTrio == null)
    	{
    		
    		if ( phasedChild.hasAnyAttribute(RBP_CANCELLED) )
    			return -1;
    		if ( !phasedChild.hasAnyAttribute(RBP_TAG) )
    			//if ( !phasedChild.isPhased() )
    				return -1;
    		if ( phasedChild.getAnyAttribute(RBP_TAG).toString().compareTo(".") == 0 )
    			return -1;
    		currTrio = new TrioRbpContext();
    		currTrio.child = child;
    		currTrio.trioGTs.add(_trioGenotypes);
    		currTrio.missingParents.add( new Integer(_mother==null ? 2 : _father==null ? 1 : 0) );
    		currTrio.filters.add(_filters);
    		currTrio.ADs.add(_ad);
    		//currTrio.childGTs.add(phasedChild);
    		currTrio.gtLoci.add(_locus);
    		
    		String _rbp_val = phasedChild.getAnyAttribute(RBP_TAG).toString();
        	String []_rbp_toks = _rbp_val.split(",");
        	if ( _rbp_toks.length != 2 )
        		throw new UserException(String.format("Error: RBP phasing tag of diploid genotype malformed at position: %s - %d, for child: %s --> %s", _locus.getContig(), _locus.getStart(), phasedChild.getSampleName(), _rbp_val) );
        	
        	//if ( Integer.parseInt(_rbp_toks[0].split("-")[0]) != _locus.getStart() )
        	//	throw new UserException(String.format("Error: Start of a new RBP haplotype appears to be phased to a DIFFERENT position than the first variant in the haplotype; phasing reference:%d\t first position in hap:%s",Integer.parseInt(_rbp_toks[0].split("-")[0]), _locus.getStart()) );
        	
        	currTrio.phaseRefPos = Integer.parseInt(_rbp_toks[0].split("-")[0]); //_locus.getStart();
    		if(_trioGenotypes.isMV())
    			currTrio.hasDnms = true;
    		
    		
    		rbp_queue.put(child.getID(), currTrio);
    		return 0;
    	}
    	// handle a het when an RBP-haplotype is already in progress 
    	if ( phasedChild.hasAnyAttribute(RBP_CANCELLED) )
			return -1;
    	if ( !phasedChild.hasAnyAttribute(RBP_TAG) )
    		return -1;
    	if ( phasedChild.getAnyAttribute(RBP_TAG).toString().compareTo(".") == 0 )
			return -1;
    	
    	String _rbp_val = phasedChild.getAnyAttribute(RBP_TAG).toString();
    	String []_rbp_toks = _rbp_val.split(",");
    	if ( _rbp_toks.length != 2 )
    		throw new UserException(String.format("Error: RBP phasing tag of diploid genotype malformed at position:%d, for child:%s", _locus.getStart(), phasedChild.getSampleName()) );
    	
    	long _refLocusForCurrKid = Integer.parseInt(_rbp_toks[0].split("-")[0]);
    	
    	// if the current genotye is the start of a new haplotype ...
    	if ( _refLocusForCurrKid != currTrio.phaseRefPos )
    	{
    		//checkAndWriteMutationsChild(child);
    		TrioRbpContext newTrioHap = rbp_queue.get(child.getID());
    		newTrioHap.child = child;
    		
    		newTrioHap.trioGTs.add(_trioGenotypes); //newTrioHap.childGTs.add(phasedChild);
    		newTrioHap.missingParents.add( new Integer(_mother==null ? 2 : _father==null ? 1 : 0) );
    		newTrioHap.filters.add(_filters);
    		newTrioHap.ADs.add(_ad);
    		newTrioHap.gtLoci.add(_locus);
    		
    		// NOT a valid check anymore since the first position in a haplotype can be a het-het-het combination and we cannot use them in the parental assignment
    		//if ( _refLocusForCurrKid != _locus.getStart() )
        	//	throw new UserException(String.format("Error: Start of a new RBP haplotype appears to be phased to a DIFFERENT position than the first variant in the haplotype; phasing reference:%d\t first position in hap:%s", _refLocusForCurrKid, _locus.getStart()) );
        	//newTrioHap.phaseRefPos = _refLocusForCurrKid;
        	//if ( _trioGenotypes.isMV() )
        	//	newTrioHap.hasDnms = true;
    		
        	rbp_queue.put(child.getID(), newTrioHap);
    		return 1;
    	}
    	// ... otherwise the current genotype continues the haplotype currently being processed ...
    	currTrio.trioGTs.add( _trioGenotypes );
    	currTrio.missingParents.add( new Integer(_mother==null ? 2 : _father==null ? 1 : 0) );
    	currTrio.filters.add(_filters);
		currTrio.ADs.add(_ad);
		//currTrio.gtLoci.add(_locus);
    	//currTrio.childGTs.add(phasedChild);
    	currTrio.gtLoci.add(_locus);
    	if ( !currTrio.hasDnms && _trioGenotypes.isMV() )
    		currTrio.hasDnms = true;
    	
    	rbp_queue.put(child.getID(), currTrio);
    	return 0;
    }
    
    
    private EnumMap<GenotypeType,Double> buildAFPrior(VariantContext vc, VariantContext vcAF)
    {
    	EnumMap<GenotypeType,Double> noDnmGtPriors = new EnumMap<GenotypeType, Double> (neutralAF);
    	
    	
        if(!useAF.equals(NO_AF))
        {
            double refAF = -1.0;
            double altAF = -1.0;
            double calledChromCount = (double)vc.getCalledChrCount(founderIds);
            
            
            if (useAF.equals(EXTERNAL))
            {
            	if (vcAF == null)
            	{
            		
            		no_af_found += 1;
            		if (af_cap > 0.0 && af_cap < 1.0)
                    {
            			// no point printing this, as it seems to happen quite often
            			//System.out.println("AF rod didn't bind - setting Af to the allele frequency cap");
                    	
            			refAF = 1 - af_cap;
                    	altAF = af_cap;
                    	
                    }
            		else
            		{
            			//System.out.println("AF rod didn't bind - af cap disabled - using neutral AF");
            			return noDnmGtPriors;
            		}
            	}
            	else            	
            		altAF = vcAF.getAttributeAsDouble(af_tag, -1.0);
            	//TODO: make this default nicer
            	if ( altAF == -1.0 )
            		//throw new UserException.BadArgumentValue("useAlleleFrequency",String.format("No allele frequency tag was found in the supplied ROD file "));
            		altAF = 0.000000000000001; // 1e-15
            	else
            	{
            		refAF = 1 - altAF;
            		//System.out.printf("---> AF[tag] at some site is: %f", altAF);
            	}
            	
            	if (af_cap > 0.0 && af_cap < 1.0)
                {
                	if (refAF < af_cap && refAF >= 0)
                		refAF = af_cap;
                	if (refAF > 1 - af_cap)
                		refAF = 1 - af_cap;
                
                	if (altAF < af_cap && altAF >= 0)
                		altAF = af_cap;
                	if (altAF > 1 - af_cap)
                		altAF = 1 - af_cap;
                }
            }
            
            //Compute the AF from the allelic dosage if set to
            if(useAF.equals(AFDOS))
            {
            	
                //Get the AF from the founders
                refAF = getCalledAlleleDosage(vc.getReference(),founderIds, vc)/calledChromCount;
                
                altAF = getCalledAlleleDosage(vc.getAltAlleleWithHighestAlleleCount(),founderIds, vc)/calledChromCount; 
                if (af_cap > 0.0 && af_cap < 1.0)
                {
                	if (refAF < af_cap && refAF >= 0)
                		refAF = af_cap;
                	if (refAF > 1 - af_cap)
                		refAF = 1 - af_cap;
                
                	if (altAF < af_cap && altAF >= 0)
                		altAF = af_cap;
                	if (altAF > 1 - af_cap)
                		altAF = 1 - af_cap;
                }
                
            }
            
            
            if (useAF.equals(GT))
            {
            	
            	int refAC = vc.getCalledChrCount(vc.getReference(),founderIds);
            	int altAC = vc.getCalledChrCount(vc.getAltAlleleWithHighestAlleleCount(),founderIds);
            	//refAF = -1.0;
            	//altAF = -1.0;
            	if (vc.getCalledChrCount(founderIds) > 0)
            	{
            		refAF = refAC / vc.getCalledChrCount(founderIds);
            		altAF = altAC / vc.getCalledChrCount(founderIds);
            	}
            	
            	if (af_cap > 0.0 && af_cap < 1.0)
            	{
            		if (refAF < af_cap && refAF >= 0)
            			refAF = af_cap;
            		if (refAF > 1 - af_cap)
            			refAF = 1 - af_cap;
            		if (altAF < af_cap && refAF >= 0)
                		altAF = af_cap;
                	if (altAF > 1 - af_cap)
                		altAF = 1 - af_cap;
            	}            	
            }
            
            //TODO: Find a better way to get around the numerical stability
            //For now, if one of the frequencies is 0, just disable the frequency; if af_cap is used (and it should) this will not be a problem
            if(refAF > 0.0 && altAF > 0.0)
            {
                double logRefAF = Math.log10(refAF);
                double logAltAF = Math.log10(altAF);

                noDnmGtPriors.put(GenotypeType.HOM_REF, 2 * logRefAF);
                noDnmGtPriors.put(GenotypeType.HET, log10_2 + logAltAF + logRefAF);
                noDnmGtPriors.put(GenotypeType.HOM_VAR, 2 * logAltAF);
            }
        }

        return noDnmGtPriors;
    }
    
    
    /**
     * Returns the total dosage of the allele a in the genotypes
     * NOTE: !!! Only works for bi-allelic sites - ploidy 1 or 2 !!!
     * @param a - the allele
     * @param sampleIds - a set of samples to take into account. If empty, all samples are included
     * @return - The total allele dosage. -1 if site was multi-allelic or if no sample had genotype called with ploidy == 2
     */
    public double getCalledAlleleDosage(Allele a, Set<String> sampleIds, VariantContext vc)
    {

    	
        if(vc.getAlleles().size() > 2)
            return -1;

        double n = 0;
        double ad;
        boolean genotypeOKFound = false;

        GenotypesContext genotypes = (sampleIds != null && !sampleIds.isEmpty() ) ? vc.getGenotypes( sampleIds ) : vc.getGenotypes();

        for ( final Genotype g : genotypes ) 
        {
            ad = getAlleleDosage(a,g);
            if(ad > -1)
            {
                n += ad;
                genotypeOKFound = true;
            }
        }

        return genotypeOKFound ? n : -1;
    }

    /**
     * Returns the allele dosage as:
     * 2*P(Hom)+P(het) if PLs are available; or 2*P(HOM) for haploid GTs,
     * or a discrete value if not (2 for HomRef, 1 for Het, 0 for HomVar)
     * Only works on sites containing the Ref allele.
     * Returns -1 in case of missing alleles
     * @param a - Allele to return the dosage of
     * @return the allelic dosage for the reference allele
     */
    public double getAlleleDosage(Allele a, Genotype g)
    {

        double refDosage = getRefAlleleDosage(g);
        if(!a.isReference() && refDosage != -1)
        {
                return 2-refDosage;
        }
        return refDosage;
    }
    
    
    /**
     * Returns the reference allele dosage as:
     * 2*P(HomRef)+P(het) if PLs are available, or 2*P(HOM) for haploid GTs,
     * or a discrete value if not (2 for HomRef, 1 for Het, 0 for HomVar)
     * Only works on sites containing the Ref allele.
     * Returns -1 in case of non diploid or missing alleles
     * @return the allelic dosage for the reference allele
     */
    public double getRefAlleleDosage(Genotype g)
    {

         if ( !g.isCalled() )
            return -1.0;

        if(g.hasLikelihoods())
        {
            //EnumMap<GenotypeType,Double> gtProbs = getLinearLikelihoodsAsMapSafeNull(g);
            if (g.getPloidy() == 2)
            {
            	EnumMap<GenotypeType,Double> gtProbs = getLinearLikelihoodsAsMapSafeNull(g,Gender.FEMALE);
            	return 2*gtProbs.get(GenotypeType.HOM_REF)+gtProbs.get(GenotypeType.HET);
            }
            if (g.getPloidy() == 1)
            {
            	EnumMap<GenotypeType,Double> gtProbs = getLinearLikelihoodsAsMapSafeNull(g,Gender.MALE);
            	return gtProbs.get(GenotypeType.HOM_REF);
            }
            return -1.0;
        }
        else
        {
            switch( g.getType() )
            {
                case HOM_REF:
                    if (g.getPloidy() == 1)
                    	return 1.0;
                    return 2.0;
                case HET:
                	if (g.getPloidy() == 2)
                		return 1.0;
                	return -1;
                case HOM_VAR:
                    return 0.0;
                default:
                    return -1.0;
            }
        }
    }
    
    // because java is lame
    public double[] getPrimitiveArray(ArrayList<Double> in)
    {
    	if (in == null) return null;
    	double temp[] = new double[in.size()];
    	//if (in.size() != 0)
        //{
        	
        	int ii = 0, i;        	
        	for (i = 0; i < in.size(); i++)// (Double[]) in.toArray() )
        	{
        		temp[ii] = in.get(i).doubleValue();
        		ii++;
        	}
        //}
    	return temp;
    }
    
    // make a string out of the AD field of a genotype; for outputting it to the mvf file
    // we only care about biallelic 
    private String makeADString(Genotype _gt)
    {
    	if (_gt == null) return "";
    	int []ad = _gt.getAD();
    	if (ad == null) return "";
    	if (ad.length == 2)
    		return "" + ad[0] + "," + ad[1];
    	else
    		return "-1";
    }
    
    private boolean isGenotypeRBPed(Genotype _gt, int consistency)
    {
    	if ( !_gt.hasExtendedAttribute(RBP_TAG) )
    		return false;
    	
    	String rbp = String.valueOf(_gt.getExtendedAttribute(RBP_TAG));
    	if ( rbp.compareTo(".") == 0 || rbp.compareTo("null") == 0 )
    		return false;
    	if ( rbp.split(",").length != 2 )
    		return false;
    	
    	if ( consistency == 1 && _gt.hasExtendedAttribute(RBP_CANCELLED) && String.valueOf(_gt.getExtendedAttribute(RBP_CANCELLED)).compareTo(String.valueOf(RBP_CANCELLED_VAL)) == 0 )
    		return false;
    	
    	return true;
    		
    }
}











