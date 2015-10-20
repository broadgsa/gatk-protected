package org.broadinstitute.gatk.tools.walkers.phasing;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.samples.Gender;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
//import org.broadinstitute.sting.utils.variantcontext.Genotype;


import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import htsjdk.samtools.SAMSequenceDictionary;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.math.BigDecimal;

import java.lang.Void;
import static java.lang.Math.pow;
import static java.lang.Math.log10;

/**
 * Traverses a VCF file and ensures that positions within specified intervals are called as haploid (for all male samples). It corrects diploid initial genotypes based on the information stored in the PL field. 
 * I supports multi-allelic sites. It ignores genotypes for which the PL field is not appropriately built, missing genotypes, already haploid genotypes. In case a call for the most likely haploid genotype cannot be made from the diploid
 * genotype PL values, a haploid NO_CALL genotype is outputted
 *
 * <p>
 * Current corrections are:
 * <ul>
 *     <li>genotype alleles</li>
 *     <li>genotype quality</li>
 *     <li>genotype phred scaled likelihoods</li>
 * </ul>
 * </p>
 * <h2>Input</h2>
 * <p>
 * <ul>
 *     <li>A VCF containing variants in haploid regions.</li>
 *     <li>A PED pedigree file containing the sex of the samples .</li>
 *     <li>An interval list file specifying the intervals that are to be treated as haploid</li>
 * </ul>
 * </p>
 *
 * <h2>Options</h2>
 * <p>
 *     <ul>
 *         <li>fake_hets: correct the genotype call and PL values for haploid samples, but leave the output in diploid representation; i.e.: HOM_REF = AA instead of A as alleles and heterozygous genotypes are assigned a prohibitively large PL value (=8000)</li>
 *     </ul>
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A VCF with where all genotypes for male samples, that fall within specified regions, are now haploid
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T HaploidWriter \
 *   -V input.vcf \
 *   -ped input.ped \
 *   -locs intervals.list \
 *   -o output.vcf
 * </pre>
 *
 */
public class HaploidWriter extends RodWalker<Void,Void>
{

	@ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();
	@Argument(shortName="locs", required=true, fullName="locations", doc="genomic regions that are to be transformed to haploid calls")
	private String filename; // locations to be considerred
	@Argument(shortName="fh", required=false, fullName="fake_hets", doc="output dummy, very big PL values for het genotypes as well")
	private int fake_hets = 0; 
	@Argument(shortName="oa", required=false, fullName="only_alleles", doc="whether we wish to modify only the allele content of the haploid genotypes; i.e: for vcfs holding phasing information")
	private boolean only_alleles = false; 
	
	@Output
    protected VariantContextWriter vcfWriter = null;
	
    private final String SOURCE_NAME = "HaploidWriter";
	
	private ArrayList<String> haploidIndiv = new ArrayList<String>();
	private ArrayList<GenomeLoc> intervals = new ArrayList<GenomeLoc>();
    //boolean flag =true;
    // constant to assign heterozygous PL values of to-be-transformed genotypes; only if fake_hets is set.
    private final int fake_hets_val = -800;
    // counter for number of genotypes that could not be transformed due to bad diploid PLs; i.e.: 2 PL vales = 0
	private int bad_likes = 0;
	// counter for number of genotypes that could not be transformed due to bad diploid alleles; i.e.: HET; only if only_alleles is set
	private int bad_genotypes = 0;
	
		
    public void initialize() 
    {
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantCollection.variants.getName());
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        
        //Get the trios from the families passed as ped
        haploidIndiv = subsetMaleSamples(this.getSampleDB().getSamples());
        intervals = readLocations(filename, this.getMasterSequenceDictionary());
        //if(trios.size()<1)
        //    throw new UserException.BadInput("No PED file passed or no trios found in PED file. Aborted.");

        //Set<Sample> families = this.getSampleDB().getSamples();
        
        //for (Sample s: families)
        //{
        //	if (s.getGender() == Gender.MALE && !haploidIndiv.contains(s.getID()))
        //		haploidIndiv.add(s.getID());
        		
        //}
       
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.addAll(GATKVCFUtils.getHeaderFields(this.getToolkit()));
        //headerLines.add(new VCFFormatHeaderLine(TRANSMISSION_PROBABILITY_TAG_NAME, 1, VCFHeaderLineType.Integer, "Phred score of the posterior probability of the genotype configuration and phased."));
        headerLines.add(new VCFHeaderLine("source", SOURCE_NAME));
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
        

    }
    
    public static ArrayList<GenomeLoc> readLocations(String filename, SAMSequenceDictionary dict)
    {
    	ArrayList<GenomeLoc> intervals = new ArrayList<GenomeLoc>();
    	
    	try
    	{
    		File f = new File(filename);
			Scanner scan = new Scanner(f);
			while(scan.hasNextLine())
			{
				String nextLoc = scan.nextLine();
				String []words1 = nextLoc.split(":");
				//if only contig is specified transform it all
				if ( words1.length == 1 )
				{
					GenomeLocParser builder = new GenomeLocParser(dict);
					if ( !builder.contigIsInDictionary(words1[0]) )
					{
						System.out.println("Invalid contig region");
						continue;
					}
					else
					{
						//int indx = builder.getContigIndex(words1[0]);
						GenomeLoc newRegion = builder.createOverEntireContig(words1[0]);
						intervals.add(newRegion);
					}
					continue;
				}
				if (words1.length != 2)
				{
					System.out.println("Location badly formatted");
					continue;
				}
				
				
				String contig = words1[0];
				String []words2 = words1[1].split("-");
				if ( words2.length == 1 )
				{
					GenomeLocParser builder = new GenomeLocParser(dict);
					if ( !builder.contigIsInDictionary(contig) )
					{
						System.out.println("Invalid contig region");
						continue;
					}
					else
					{
						int start = Integer.parseInt(words2[0]);
						int stop = builder.getContigInfo(words1[0]).getSequenceLength();
						int indx = builder.getContigIndex(contig);
						GenomeLoc newRegion = builder.createGenomeLoc(words1[0],indx,start,stop);
						intervals.add(newRegion);
					}
										
					continue;
				}
				
				if (words2.length != 2)
				{
					System.out.println("Location badly formatted");
					continue;
				}
				
				int start = Integer.parseInt(words2[0]);
				int stop = Integer.parseInt(words2[1]);
				
				if ( start > stop )
				{
					System.out.println("Invalid loacation");
					continue;
				}
				
				GenomeLocParser builder = new GenomeLocParser(dict);
				if ( !builder.contigIsInDictionary(contig) )
				{
					System.out.println("Invalid contig region");
					continue;
				}
				else
				{
					int indx = builder.getContigIndex(contig);
					GenomeLoc newRegion = builder.createGenomeLoc(contig,indx,start,stop);
					intervals.add(newRegion);
				}
			}
    	}
    	catch(Exception e)
    	{
    		System.out.println("Error: could not open locations file");
    		e.printStackTrace();
    	}
    	
    	return intervals;    	
    }

    public static ArrayList<String> subsetMaleSamples(Set<Sample> families)
    {

        //Set<Sample> families = this.getSampleDB().getSamples();
        ArrayList<String> haploidIndiv = new ArrayList<String>();
    	
        for (Sample s: families)
        {
        	if (s.getGender() == Gender.MALE && !haploidIndiv.contains(s.getID()))
        		haploidIndiv.add(s.getID());
        	//Sample pops = s.getFather();
        	if ( !haploidIndiv.contains(s.getPaternalID()) )
        		haploidIndiv.add(s.getPaternalID());
        		
        }
        
        return haploidIndiv;
    }
	
    private static Genotype onlyConsiderAlleles(Genotype temp)
    {
    	GenotypeBuilder JB = new GenotypeBuilder(temp);
		List<Allele> oldAls = temp.getAlleles();
		ArrayList<Allele> newAlleles = new ArrayList<Allele>(1);
		
		if ( oldAls.size() == 1 )
			newAlleles.add(oldAls.get(0));
		else
			if ( oldAls.size() == 2 )
			{
				if ( oldAls.get(0).equals(oldAls.get(1),true))
					newAlleles.add(oldAls.get(0));
				// if the 2 diploid alleles are not identical no further discrimination can be made
				else
				{
					//bad_genotypes++;
					newAlleles.add(oldAls.get(0));
				}
			}
			// more than 2 alleles = bad
			else
			{
				//bad_genotypes++;
				newAlleles.add(oldAls.get(0));
			}
		JB.alleles(newAlleles);
		Genotype replacement = JB.make();
		
		return replacement;
    }
    
    public static Genotype transformMultiAllelicGT(Genotype temp, int numAls, Allele ref, List<Allele> alts, int fake_hets, int fake_hets_val)
    {
    	int []oldLikes = temp.getPL();
    	
    	GenotypeBuilder JB = new GenotypeBuilder(temp);
    	
    	ArrayList<Allele> newAlleles = new ArrayList<Allele>(2);
    	
		int len = (((numAls-1) * numAls) / 2) + numAls;
		double []probs = new double[len];
		// return to probability space
		for (int j = 0; j < len; j++)
		{
			probs[j] = pow(10, (double)oldLikes[j]/(-10));
			if (probs[j] == 0)
				probs[j] = Double.MIN_VALUE;
		}
		// compute the probabiltiy of the 0-forced likelihood
		for (int j = 0; j < len; j++)
		{
			if ( probs[j] == 1)
			{
				//double sum = 0.0;
				for (int j1 = 0; j1 < len; j1++)
					if (j1 != j)
						probs[j] -= probs[j1];
				break;
			}
		}
		
		//find the new, haploid, genotype from the likelihood values
		double max = 0.0;
		int maxInd = 0;
		int alIndx = 0;
		int gqIndx = 0;
		// only compare likelihoods corresponding to homozygous initial genotypes
		for (int j = 0; j < numAls; j++)
		{
			if (probs[(j*(j-1))/2 + j] > max)
			{
				max = probs[(j*(j-1))/2 + j];
				maxInd = (j*(j-1))/2 + j;
				gqIndx = alIndx;
				alIndx = j;
			}
		}
		
		// renormalize genotype posteriors, considering only hamozygous GTs
		double newNorm = 0.0;
		//boolean stroke = false;
		for (int j = 0; j < numAls; j++ )
		{
			
			newNorm += probs[(j*(j+1))/2 + j];
			
			// if the most likely homozygous genotype is not unique, i.e.: there are 2 PLs having the same most likely value, a decision cannot be made and haploid genotype is set to NO_CALL
			for (int j1 = 0; j1 < j; j1++)
			{
				
				if (oldLikes[(j*(j+1))/2 + j] == oldLikes[(j1*(j1+1))/2 + j1] && oldLikes[(j*(j+1))/2 + j] == (int)(log10(max) * (-10)) )
				{
					newAlleles.add(Allele.NO_CALL);
        			
        			JB.alleles(newAlleles);
                	JB.GQ(-1);
                	JB.noPL();
                	Genotype replacement = JB.make();
                	return replacement;
				}
			}
		}
		
		
		for ( int j = 0; j < numAls; j++ )
		{
			
			probs[(j*(j+1))/2 + j] /= newNorm;
		}
		
		// force called GT posterior to 1, i.e.: PL = 0
		probs[maxInd] = 1.0;
		
		if (maxInd == 0)
			newAlleles.add(ref);
		else
			newAlleles.add(alts.get(alIndx-1));
		
		// if heterozygous PLs are to be kept, they are assigned irrelevantly large values 
		if (fake_hets == 1)
		{
			double []Likes = new double[oldLikes.length];
			for(int j = 0; j < numAls; j++)
				for (int k = 0; k < numAls; k++)
				{
					if (k == j)
						Likes[k*(k+1)/2 + j] = log10(probs[k*(k+1)/2 + j]);
					else
						Likes[k*(k+1)/2 + j] = fake_hets_val;//log10(Double.MIN_VALUE);
				}
			
			JB.GQ( (int) (Likes[gqIndx] * (-10)) );
			JB.alleles(newAlleles);
        	JB.PL(Likes);
        	Genotype replacement = JB.make();
        	
        	return replacement;
		}
		// new PLs are computed and assigned
		else
		{
			double []Likes = new double[numAls];
			for(int j = 0; j < numAls; j++)
			{
				Likes[j] = log10(probs[(j*(j-1))/2 + j]);
			}
			
			JB.GQ( (int) (Likes[gqIndx] * (-10)) );
			JB.alleles(newAlleles);
        	JB.PL(Likes);
        	Genotype replacement = JB.make();
        	
        	return replacement;
		}
		        			
	}
      
    public static Genotype transformBiAllelicGT(Genotype temp, Allele ref, Allele alt, int fake_hets, int fake_hets_val)
    {
    	int []oldLikes = temp.getPL();
    	if (oldLikes.length != 3)
    		return temp;
    	
    	GenotypeBuilder JB = new GenotypeBuilder(temp);
    	
    	ArrayList<Allele> newAlleles = new ArrayList<Allele>(2);
    	
    	double[] probs = new double[3];
    	for (int j=0; j<3;j++)
    	{
    		probs[j] = pow(10, (double)oldLikes[j]/-10.0);
    		// precision limitations when transforming from PL-log-space to probability-space
    		if (probs[j] == 0)
    			probs[j] = Double.MIN_VALUE;
    	}
    	// retrieve the forced-to-0 PL
    	if (probs[0] == 1.0)
    		probs[0] = 1 - probs[1] - probs[2];
    	else
    		if (probs[1] == 1.0)
    			probs[1] = 1 - probs[0] - probs[2];
    		else
    			if (probs[2] == 1.0)
    				probs[2] = 1 - probs[1] - probs[0];
    	
    	// renormalise using only homozygous genotypes
    	double newNorm = probs[0] + probs[2];
    	probs[0] = probs[0]/newNorm;
    	probs[2] = probs[2]/newNorm;
    	if (probs[0] == 0)
			probs[0] = Double.MIN_VALUE;
    	if (probs[2] == 0)
			probs[2] = Double.MIN_VALUE;
    	
    	// if both homozygous genotypes are equally likely from PL values => no decision can be made => NO_CALL haploid genotype
    	if (oldLikes[0] == oldLikes[2] && temp.isHet())
    	{
    		//probs[0] = probs[2];
    		newAlleles.add(Allele.NO_CALL);
    		
    		JB.alleles(newAlleles);
        	JB.GQ(-1);
        	JB.noPL();
        	
        	            	
        	Genotype replacement = JB.make();
        	return replacement;
    		
    	}
    	else
    	{
    		if (probs[0] > probs[2])
    		{
    			probs[0] = 1.0;
    			newAlleles.add(ref);
    			//try{newAlleles.add(temp.getAllele(0));}
    			//catch (Exception e){
    			//	int a =1;}
    		}
    		else
    		{
    			probs[2] = 1.0;
    			newAlleles.add(alt);
    		}
    		// if heterozygous PLs are to be kept, assign them the predefined value
    		if (fake_hets == 1)
    		{
    			double[] Likes = new double[3];
    			Likes[0] = log10(probs[0]);// * (-10.0);
    			Likes[1] = fake_hets_val;//log10(Double.MIN_VALUE);
    			Likes[2] = log10(probs[2]);// * (-10.0);
    			
    			
    			// Populate the modiffied Genotype Object
    			JB.alleles(newAlleles);
    			JB.PL(Likes);
    			if (Math.abs(Likes[0]) == 0)
    				JB.GQ((int) (Likes[2] * (-10)));
    			else
    				JB.GQ((int) (Likes[0] * (-10)));
    			
    			Genotype replacement = JB.make();
            	return replacement;
    		}
    		else
    		{
    			double[] Likes = new double[2];

    			// the PL() method of GenotypeBuilder normalizes the input again before it assigns it so don't normalize
    			Likes[0] = log10(probs[0]);// * (-10.0);
    			Likes[1] = log10(probs[2]);// * (-10.0);
    			// Populate the modiffied Genotype Object
    			JB.alleles(newAlleles);
    			JB.PL(Likes);
    			if (Math.abs(Likes[0]) == 0)
    				JB.GQ((int) (Likes[1] * (-10)));
    			else
    				JB.GQ((int) (Likes[0] * (-10)));
        	           	           	
    			Genotype replacement = JB.make();
            	//EnumMap<GenotypeType,Double> fatherLikelihoods2 = replacement.getLikelihoods().getAsMap(false);

    			return replacement;
    		}
    	}
    }
    
    public static VariantContext processVC(VariantContext vc, ArrayList<String> haploidIndiv, ArrayList<GenomeLoc> intervals, int fake_hets, int fake_hets_val, boolean only_alleles)
    {
    	VariantContextBuilder builder = new VariantContextBuilder(vc);
        // get all genotypes at this position
        GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());
        for (int i = 0; i < genotypesContext.size(); i++)
        {
        	
        	Genotype temp = genotypesContext.get(i);
        	
        	// skip genotypes that are not to be made haploid as specified by the ped file
        	if (!haploidIndiv.contains(temp.getSampleName()))
        	{
        		genotypesContext.replace(temp);
        		continue;
        	}
        	
        	// if only alleles are available/to be used for the genotypes to be transformed, take care of it fast; i.e.: phasing data
        	if (only_alleles)
        	{
        		Genotype replacement = onlyConsiderAlleles(temp);
        		genotypesContext.replace(replacement);
        		continue;
        	}
        	       	
        	// this makes likelihoods unusable in determining the most likely haploid GT
        	//if (temp.getGQ() == 0)
        	//	bad_likes++;
        	
        	//newAlleles.add(temp.getAlleles().get(0));        	        	
        	int[] oldLikes = null;
        	oldLikes = temp.getPL();
        	
        	if (oldLikes == null)
        	{
        		// Initialize the modified Genotype by the original Genotype then compute differences
            	GenotypeBuilder JB = new GenotypeBuilder(temp);
            	ArrayList<Allele> newAlleles = new ArrayList<Allele>(2);
        		List<Allele> nonStandardAls = temp.getAlleles();
        		if (nonStandardAls.get(0).equals(nonStandardAls.get(1),true))
        			newAlleles.add(nonStandardAls.get(0));
        		else
        		{
        			newAlleles.add(nonStandardAls.get(0));
        			newAlleles.add(nonStandardAls.get(1));
        		}
            	JB.alleles(newAlleles);
            	JB.GQ(-1);
            	JB.noPL();
            	Genotype replacement = JB.make();
            	genotypesContext.replace(replacement);
        		
            	
        		//System.out.println("No likes accompanying the Genotype\n");
        		continue;
        	}
        	
        	// multi-allelic case
        	if (oldLikes.length > 3)
        	{
        		int numAls = vc.getNAlleles();
        		
        		// check if the number of produced likelihoods is consistent with the number of alleles found
        		if ( oldLikes.length == (((numAls-1) * numAls) / 2) + numAls )
        		{        			
        			Genotype replacement = transformMultiAllelicGT(temp,numAls, vc.getReference(), vc.getAlternateAlleles(), fake_hets, fake_hets_val);
        			genotypesContext.replace(replacement);
        		}
        		// if PLs are not well specified leave genotype unchanged
        		else
        		{
        			genotypesContext.replace(temp);
        		}
        		//System.out.println("Either not diploid or multi-allelic site\n");
        		
        	}
        	else
        	{
        		// Most common, bi-allelic case
        		Genotype replacement = transformBiAllelicGT(temp, vc.getReference(), vc.getAlternateAlleles().get(0), fake_hets, fake_hets_val);
        		genotypesContext.replace(replacement);
        	}
        }
        
        builder.genotypes(genotypesContext);
        return builder.make();
    }
    
    
    @Override
	public Void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) 
	{
		if (tracker == null)
			return null;
		VariantContext vc = tracker.getFirstValue(variantCollection.variants, context.getLocation());
        
        if (vc == null)
        	return null;
        
        // check if current variant position should he changed, write it like it is if not
        boolean transform = false;
        for (GenomeLoc i: intervals)
        	if (i.containsP(ref.getLocus()))
        		transform = true;
        if (!transform)
        {
        	VariantContextBuilder builder = new VariantContextBuilder(vc);
        	vcfWriter.add(builder.make());
        	return null;
        }
        
        VariantContextBuilder builder = new VariantContextBuilder(vc);
        // get all genotypes at this position
        GenotypesContext genotypesContext = GenotypesContext.copy(vc.getGenotypes());
        for (int i = 0; i < genotypesContext.size(); i++)
        {
        	
        	Genotype temp = genotypesContext.get(i);
        	
        	// skip genotypes that are not to be made haploid as specified by the ped file
        	if (!haploidIndiv.contains(temp.getSampleName()))
        	{
        		genotypesContext.replace(temp);
        		continue;
        	}
        	
        	// if only alleles are available/to be used for the genotypes to be transformed, take care of it fast; i.e.: phasing data
        	if (only_alleles)
        	{
        		Genotype replacement = onlyConsiderAlleles(temp);
        		genotypesContext.replace(replacement);
        		continue;
        	}
        	       	
        	// this makes likelihoods unusable in determining the most likely haploid GT
        	if (temp.getGQ() == 0)
        		bad_likes++;
        	
        	//newAlleles.add(temp.getAlleles().get(0));        	        	
        	int[] oldLikes = null;
        	oldLikes = temp.getPL();
        	
        	if (oldLikes == null)
        	{
        		// Initialize the modified Genotype by the original Genotype then compute differences
            	GenotypeBuilder JB = new GenotypeBuilder(temp);
            	ArrayList<Allele> newAlleles = new ArrayList<Allele>(2);
        		List<Allele> nonStandardAls = temp.getAlleles();
        		if (nonStandardAls.get(0).equals(nonStandardAls.get(1),true))
        			newAlleles.add(nonStandardAls.get(0));
        		else
        		{
        			newAlleles.add(nonStandardAls.get(0));
        			newAlleles.add(nonStandardAls.get(1));
        		}
            	JB.alleles(newAlleles);
            	JB.GQ(-1);
            	JB.noPL();
            	Genotype replacement = JB.make();
            	genotypesContext.replace(replacement);
        		
            	
        		//System.out.println("No likes accompanying the Genotype\n");
        		continue;
        	}
        	
        	// multi-allelic case
        	if (oldLikes.length != 3)
        	{
        		int numAls = vc.getNAlleles();
        		
        		// check if the number of produced likelihoods is consistent with the number of alleles found
        		if ( oldLikes.length == (((numAls-1) * numAls) / 2) + numAls )
        		{        			
        			Genotype replacement = transformMultiAllelicGT(temp,numAls, vc.getReference(), vc.getAlternateAlleles(), fake_hets, fake_hets_val);
        			genotypesContext.replace(replacement);
        		}
        		// if PLs are not well specified leave genotype unchanged
        		else
        		{
        			genotypesContext.replace(temp);
        		}
        		//System.out.println("Either not diploid or multi-allelic site\n");
        		
        	}
        	else
        	{
        		// Most common, bi-allelic case
        		Genotype replacement = transformBiAllelicGT(temp, vc.getReference(), vc.getAlternateAlleles().get(0), fake_hets, fake_hets_val);
        		genotypesContext.replace(replacement);
        	}
        }
        
        builder.genotypes(genotypesContext);
        vcfWriter.add(builder.make());       
        
        
        return null;
	}
	
	@Override
	public Void reduceInit()
	{
		return null;
	}
	
	@Override
	public Void reduce(Void ala, Void bala)
	{
		return null;
	}
	
	@Override
	public void onTraversalDone(Void main)
	{
		System.out.printf("\t bad likes: %d\n",bad_likes);
		System.out.printf("\t haploid genotypes imputed as het: %d\n",bad_genotypes);
		
	}
	



}
