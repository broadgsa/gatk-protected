# M2 Panel of Normals

In order to reduce false positives, we use a panel of "normal" (ie non-cancer) samples to filter out both germline events as well as systematic noise.  The form of the panel is a VCF file, which is produced via a Queue scripts

You must supply:
- the reference (defaults to hg19)
- the intervals to evaluate
- the list of BAM files 

### How To Run
The Queue script for producing the PON can be found in the gsa-unstable repository under ```private/gatk-tools-private/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2```

First, chose the appropriate settings (runnable as environment variables here)
```
QUEUE_JAR=<your-queue-jar>
GSA_UNSTABLE_HOME=<path-to-your-gsa-unstable-checkout>
TEMPDIR=/broad/hptmp/$USER
```

and then run the following Queue command
```
java \
  -Djava.io.tmpdir=$TEMPDIR \
  -jar $QUEUE_JAR \
  -S $GSA_UNSTABLE_HOME/private/gatk-tools-private/src/main/java/org/broadinstitute/gatk/tools/walkers/cancer/m2/run_M2_ICE_NN.scala \
  --job_queue gsa -qsub -jobResReq virtual_free=5G -startFromScratch  \
  -sc 50 \
  --allbams <your-list-of-bams> \
  --intervals <your-intervals> \
  --outputprefix <your-output-prefix> \
  --start_from_scratch --keep_intermediate_outputs \
  -run
```

This will produce many VCFs (1 per sample), plus \<your-output-prefix\>.genotypes.vcf and \<your-output-prefix\>.vcf which are the panel of normals VCF both with and without sample-genotype information.  Typically the latter is the one used as input to M2, although either will work.