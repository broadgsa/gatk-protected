/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Allele balance across all samples
 *
 * <p>The allele balance is the fraction of ref bases over ref + alt bases.</p>
 *
 * <h3>Caveats</h3>
 * <p>Note that this annotation will only work properly for biallelic samples that are called as heterozygous.</p>
 */
public class AlleleBalance extends InfoFieldAnnotation {


    char[] BASES = {'A','C','G','T'};
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        //if ( stratifiedContexts.size() == 0 )
        //    return null;

        if ( !vc.isBiallelic() )
            return null;
        final GenotypesContext genotypes = vc.getGenotypes();
        if ( !vc.hasGenotypes() )
            return null;

        double ratioHom = 0.0;
        double ratioHet = 0.0;
        double weightHom = 0.0;
        double weightHet = 0.0;
        double overallNonDiploid = 0.0;
        for ( Genotype genotype : genotypes ) {
            
            if ( vc.isSNP() ) {
            
            	final int[] counts = genotype.getAD();
            	// If AD was not calculated, we can't continue
            	if(counts == null)
            		continue;
            	
            	final int n_allele = counts.length;
            	int count_sum = 0;
            	for(int i=0; i<n_allele; i++){
            		count_sum += counts[i];
            	}
                double pTrue = 1.0 - Math.pow(10.0,-genotype.getGQ() / (double) 10 );
                if ( genotype.isHet() ) {
                	
                	final int otherCount = count_sum - (counts[0] + counts[1]);
                    // sanity check
                    if ( counts[0] + counts[1] == 0 )
                        continue;

                    // weight the allele balance by genotype quality so that e.g. mis-called homs don't affect the ratio too much
                    ratioHet += pTrue * ((double)counts[0] / (double)(counts[0] + counts[1]));
                    weightHet += pTrue;
                    overallNonDiploid += ( (double) otherCount )/((double) count_sum*genotypes.size());
                } else if ( genotype.isHom() ) {
                	final int alleleIdx = genotype.isHomRef() ?  0 : 1 ;
                    final int alleleCount = counts[alleleIdx];
                    int bestOtherCount = 0;
                    for(int i=0; i<n_allele; i++){
                    	if( i == alleleIdx )
                    		continue;
                    	if( counts[i] > bestOtherCount )
                    		bestOtherCount = counts[i];
                    }
                    
                    final int otherCount = count_sum - alleleCount;
                    ratioHom += pTrue*( (double) alleleCount)/((double) (alleleCount+bestOtherCount));
                    weightHom += pTrue;
                    overallNonDiploid += ((double ) otherCount)/((double) count_sum*genotypes.size());
                }
                // Allele Balance for indels was not being computed correctly (since there was no allele matching).  Instead of
                // prolonging the life of imperfect code, I've decided to delete it.  If someone else wants to try again from
                // scratch, be my guest - but make sure it's done correctly!  [EB]
            }
        }

        // make sure we had a het genotype

        Map<String, Object> map = new HashMap<String, Object>();
        if ( weightHet > 0.0 ) {
            map.put("ABHet",ratioHet/weightHet);
        }

        if ( weightHom > 0.0 ) {
            map.put("ABHom",ratioHom/weightHom);
        }

        if ( overallNonDiploid > 0.0 ) {
            map.put("OND",overallNonDiploid);
        }
        return map;
    }


    public List<String> getKeyNames() { return Arrays.asList("ABHet","ABHom","OND"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ABHet", 1, VCFHeaderLineType.Float, "Allele Balance for hets (ref/(ref+alt))"),
            new VCFInfoHeaderLine("ABHom", 1, VCFHeaderLineType.Float, "Allele Balance for homs (A/(A+O))"),
            new VCFInfoHeaderLine("OND", 1, VCFHeaderLineType.Float, "Overall non-diploid ratio (alleles/(alleles+non-alleles))")); }
}
