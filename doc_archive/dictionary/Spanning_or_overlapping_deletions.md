## Spanning or overlapping deletions

http://gatkforums.broadinstitute.org/gatk/discussion/6926/spanning-or-overlapping-deletions

<p>We use the term <strong>spanning deletion</strong> or <strong>overlapping deletion</strong> to refer to a deletion that spans a position of interest.  </p>
<p>The presence of a spanning deletion affects how we can represent genotypes at any site(s) that it spans for those samples that carry the deletion, whether in heterozygous or homozygous variant form. Page 8, item 5 of the <a href="https://samtools.github.io/hts-specs/VCFv4.3.pdf">VCF v4.3 specification</a> reserves the <code>*</code> allele to reference overlapping deletions. This is not to be confused with the bracketed asterisk <code>&lt;*&gt;</code> used to denote symbolic alternate alleles.</p>
<hr />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/3e/6389220c22db2a69857811121c3dd1.png" width="400" align="right" border="10"/> 
<p>Here we illustrate with four human samples. Bob and Lian each have a heterozygous <code>A</code> to <code>T</code> single polymorphism at position 20, our position of interest. Kyra has a 9 bp deletion from position 15 to 23 on both homologous chromosomes that extends across position 20. Lian and Omar each are heterozygous for the same 9 bp deletion. Omar and Bob's other allele is the reference <code>A</code>.</p>
<p><strong>What are the genotypes for each individual at position 20?</strong> For Bob, the reference A and variant T alleles are clearly present for a genotype of <code>A/T</code>.</p>
<p>What about Lian? Lian has a variant T allele plus a 9 bp deletion overlapping position 20. To notate the deletion as we do single nucleotide deletions is technically inaccurate. We need a placeholder notation to signify absent sequence that extends beyond the position of interest and that is listed for an earlier position, in our case position 14. The solution is to use a star or asterisk <code>*</code> at position 20 to refer to the spanning deletion. Using this convention, Lian's genotype is <code>T/*</code>.</p>
<p>At the sample-level, Kyra and Omar would not have records for position 20. However, we are comparing multiple samples and so we indicate the spanning deletion at position 20 with <code>*</code>. Omar's genotype is <code>A/*</code> and Kyra's is <code>*/*</code>. </p>
<hr />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/f6/72037063701d84343ea6469fe64d2e.png" height="180" align="right" border="10"/>
<p><strong>In the VCF</strong>, depending on the format used by tools, positions equivalent to our example position 20 may or may not be listed. If listed, such as in the first example VCF shown, the spanning deletion is noted with the asterisk <code>*</code> under the <code>ALT</code> column. The spanning deletion is then referred to in the genotype <code>GT</code> for Kyra, Lian and Omar. Alternatively, a VCF may altogether avoid referencing the spanning deletion by listing the variant with the spanning deletion together with the deletion. This is shown in the second example VCF at position 14. </p>