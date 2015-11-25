# CRSP ICE NA12878 Specificity Evaluation

In order to evaluate the specificity of M2, we sequenced replicates of NA12878 using ICE (Illumina Content Exomes) all call all pairwise combinations as tumor-normals.  By definition, everything called is a false positive.

The target territory is ```/dsde/working/mutect/crsp_nn/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.no_empty.interval_list```

All scripts referenced here are relative to the current working directory of ```
/dsde/working/mutect/crsp_nn```

### Current M2 Performance

(gsa-unstable 7/13/15, commit:9e93a70)

| type | # of false positives |
|------|----------------------|
|SNP|99|
|INDEL|15|


TODO: write a simple tool to do this more easily

To calculate per pair-counts, run:
```
# for SNPs
for vcf in *.bam.vcf 
do     
    cat $vcf | grep PASS | awk '{ if ( length($4) + length($5) == 2) print $0 }' | wc -l 
done

# for INDELs
for vcf in *.bam.vcf 
do     
    cat $vcf | grep PASS | awk '{ if ( length($4) + length($5) != 2) print $0 }' | wc -l 
done
```

### Current M1 and Indelocator Performance
For comparison, the M1 & Indelocator calls have been made on this same data set in the Firehose workspace ```CRSP_ICE_NA1878_Production_Analysis``` in the pair set ```NA12878_Replicate_Pairs``` which contains 4 samples and 12 pairwise combinations.

| type | # of false positives |
|------|----------------------|
|SNP|181|
|INDEL|106|

These results can be obtained (from a LSF / CGA node running the FuSE daemon)

```
SNP:
cat /local/cga-fh/cga/CRSP_ICE_NA1878_Production_Analysis/Pair_Set/NA12878_Replicate_Pairs/Pair/*/jobs/capture/mut/calls/latest/*.call_stats.txt | grep KEEP | wc -l

INDEL (need to restrict to target territory):
reuse BEDTools
cat /dsde/working/mutect/crsp_nn/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.no_empty.interval_list | grep -v "@" | awk '{ print $1 "\t" $2-1 "\t" $3 }' > ice.bed
cat /local/cga-fh/cga/CRSP_ICE_NA1878_Production_Analysis/Pair_Set/NA12878_Replicate_Pairs/Pair/*/jobs/capture/indel/maflite/latest/*.full.maf | grep KEEP | cut -f2-4 | awk '{ print $1 "\t" $2-1 "\t" $3 }' > /tmp/indels.bed
bedtools intersect -wa -a /tmp/ice.bed -b /tmp/indels.bed | wc -l
```


### How To Run
The SCALA script for running M2 can be found in the gsa-unstable repository under ```private/gatk-tools-private/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2```

First, chose the appropriate settings (runnable as environment variables here)
```
QUEUE_JAR=<your-queue-jar>
OUT_VCF=<your-output-vcf>
GSA_UNSTABLE_HOME=<path-to-your-gsa-unstable-checkout>
TEMPDIR=/broad/hptmp/$USER
```

and then run the following Queue command
```
java \
  -Djava.io.tmpdir=$TEMPDIR \
  -jar $QUEUE_JAR \
  -S $GSA_UNSTABLE_HOME/private/gatk-tools-private/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2/run_M2_ICE_NN.scala \
  -sc 50 \
  --job_queue gsa -qsub -jobResReq virtual_free=5G -startFromScratch  \
  --allbams /humgen/gsa-hpprojects/NA12878Collection/bams/crsp_ice_validation//NA12878.intra.flowcell.replicate.bam_list \
  -o <your-output-prefix>
  -run
```
