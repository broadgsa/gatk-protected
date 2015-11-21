# Dream Challenge Evaluation

In order to evaluate the performance of M2, we use two sets of data from the SMC DREAM Challenge.  Specifically challenges #3 and #4.

All scripts referenced here are relative to the current working directory of ```
/dsde/working/mutect/dream_smc```

### Current Performance (Unmasked)
From the output of the evaluation method 

(gsa-unstable 7/13/15, commit:9e93a70)

|set | subset | type | sensitivity | specificity | accuracy |
|----|--------|------|-------------|-------------|----------|
|SMC 3|chr21|SNP|0.935897435897|0.935897435897|0.935897435897|
|SMC 3|chr21|INDEL|0.904255319149|0.977011494253|0.940633406701|
|SMC 3|wgs|SNP|0.930532709098|0.955188985583|0.94286084734|
|SMC 3|wgs|INDEL|0.902139907396|0.970516962843|0.93632843512|
|SMC 4|chr21|SNP|0.769607843137|0.969135802469|0.869371822803|
|SMC 4|chr21|INDEL|0.771241830065|0.991596638655|0.88141923436|
|SMC 4|wgs|SNP|0.764507007622|0.975374480433|0.869940744028|
|SMC 4|wgs|INDEL|0.768634634353|0.989389679877|0.879012157115|


 
### How To Run
The SCALA script for running M2 can be found in the gsa-unstable repository under ```private/gatk-tools-private/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2```

First, chose the appropriate settings (runnable as environment variables here)
```
QUEUE_JAR=<your-queue-jar>
OUT_VCF=<your-output-vcf>
GSA_UNSTABLE_HOME=<path-to-your-gsa-unstable-checkout>

# for Dream 3
NORMAL_BAM=/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set3.normal.bam
TUMOR_BAM=/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set3.tumor.bam

# for Dream 4
NORMAL_BAM=/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set4.normal.bam
TUMOR_BAM=/dsde/working/mutect/dream_smc/bams/synthetic.challenge.set4.tumor.bam

# for WGS
INTERVALS=/dsde/working/mutect/dream_smc/bams/wgs_calling_regions.v1.interval_list

# for chromosome 21 only
INTERVALS=/dsde/working/mutect/ts/c21_wgs_calling_regions.v1.interval_list

TEMPDIR=/broad/hptmp/kcibul/mutect
```

and then run the following Queue command
```
java \
  -Djava.io.tmpdir=$TEMPDIR \
  -jar $QUEUE_JAR \
  -S $GSA_UNSTABLE_HOME/private/gatk-tools-private/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2/run_M2_dream.scala \
  --job_queue gsa -qsub -jobResReq virtual_free=5G -startFromScratch \
  -sc 200 \
  -normal $NORMAL_BAM \
  -tumor $TUMOR_BAM \
  -L $INTERVALS \
  -o $OUT_VCF \
  -run
```

### How To Evaluate

Run the following
```
/dsde/working/mutect/dream_smc/dream_eval.pl [3|4] [wgs|21] [SNV|INDEL] input.vcf
```
where 
    - [3|4] the dream challenge round
    - [wgs|21] evaluate the whole genome, or just a subset (chromosome 21)
    - [SNV|INDEL] evaulate SNV (SNPs) or INDELS

