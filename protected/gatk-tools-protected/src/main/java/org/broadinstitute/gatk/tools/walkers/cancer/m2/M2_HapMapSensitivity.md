# CRSP HapMap Sensitivity Evaluation

###Current M2 Performance
(gsa-unstable 9/1/15, commit:a08903d)

| Mixture | type | sensitvity |
|------|----------------------|
| 5-plex |SNP|0.9691274|
| 5-plex |INDEL|0.87466127|
| 10-plex |SNP|0.97179496|
| 10-plex |INDEL|0.8888889|
| 20-plex |SNP|0.9537307|
| 20-plex |INDEL|0.83281654|


###Run Procedure
Run the script separately for each HapMap mixture bam:

inputDir=/dsde/working/mutect/laura/hapmapSensitivity/inputs/
Queue_Jar=<Queue jar of interest>

```
java -jar $Queue_Jar -S Qscript_HapMapPlex.scala \
-intervals $inputDir/agilent_5plex_intervalFiles.list \
-tumors $inputDir/agilent_5plex_bams.list \
-truthVCF $inputDir/agilent_5plex_truth_intervals.vcf \
-snpCounts $inputDir/agilent_5plex_truth_intervals.snpCounts.list \
-indelCounts $inputDir/agilent_5plex_truth_intervals.indelCounts.list \
-o <output.5plex.sensitivity.report> \
-qsub -jobQueue gsa -jobResReq virtual_free=5G -sc 50
```

The HapMap bams get run as tumors without normals because we're not interested in specificity here, so we don't need the normals to filter out noise

###Inputs
Bam lists:
5- and 10-plex have 3 replicates, 20-plex has 9

Interval files:
If we're only interested in sensitivity, then we only need to run the caller around known true positive sites, which we take from the truth VCFs
This workaround repeats the truth filename for the number of bams -- in theory each could have a separate truth VCF, but they are the same titration mixture

SNP/INDEL counts:
This is the number of events in the truth VCFs so we can find the sensitivity across all samples
TODO: this could be generalized

###Outputs
Each run outputs its own SNP and INDEL sensitivity combined across all samples:
```
Sensitivity across all samples:
SNPs: 0.95156
INDELs: 0.7328859
```

Note that these are not filtered for depth as described in the CRSP documentation

###Resources
Truth file preparation for 5-plex:
Start with /cga/tcga-gsc/benchmark/data/crsp-truth/1kg_5plex_wgs_hc_calls.codingIndelSnp.db135.recode.vcf
Select out allele fraction greater than 20% using "vc.isBiallelic() ? AF >= 0.2 : vc.hasGenotypes() && vc.getCalledChrCount(vc.getAltAlleleWithHighestAlleleCount())/(1.0*vc.getCalledChrCount()) >= 0.2"

Similar for 10-plex source:
/cga/tcga-gsc/benchmark/data/crsp-truth/1kg_10plex_wgs_hc_calls.codingIndelSnp.db135.recode.vcf

And 20-plex source:
/cga/tcga-gsc/benchmark/data/crsp-truth/1kg_20plex_wgs_hc_calls.codingIndelSnp.db135.recode.vcf

both also using AF filter of 0.2