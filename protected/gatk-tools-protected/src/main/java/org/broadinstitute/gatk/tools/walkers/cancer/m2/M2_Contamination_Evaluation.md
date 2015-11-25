# M2 Contamination Correction Evaluation

In order to evaluate the efficacy of the contamination correction in M2 (adapted from HaplotypeCaller), we created synthetic data consisting of the four intra-run CRSP NA12878 replicates, each contaminated with 1-5% of the HCC-1143 normal. 

### Creating Artificial Contamination Data
``` QUEUE_JAR=<your-queue-jar>
    GSA_UNSTABLE_HOME=<path-to-your-gsa-unstable-checkout>
    BASELINE_BAM=/crsp/picard_aggregation/000007820918/SM-612V3/current/SM-612V3.bam
    CONTAMINANT_BAM=/seq/tier3b/picard_aggregation/C970/HCC1143_BL/v1/HCC1143_BL.bam
    ```

```
java -jar $QUEUE_JAR \
-S $GSA_UNSTABLE_HOME/private/gatk-queue-extensions-internal/src/main/qscripts/org/broadinstitute/gatk/queue/qscripts/dev/CreateSyntheticContaminationScript.scala \
-b2 CONTAMINANT_BAM -b1 BASELINE_BAM \
-f 0.01 -f 0.02 -f 0.03 -f 0.04 -f 0.05
```
     
Repeat for the other three NA12878 replicates.

    BASELINE_BAM=/crsp/picard_aggregation/000007820818/SM-612V4/current/SM-612V4.bam
    
    BASELINE_BAM=/crsp/picard_aggregation/000007820718/SM-612V5/current/SM-612V5.bam
    
    BASELINE_BAM=/crsp/picard_aggregation/000007820618/SM-612V6/current/SM-612V6.bam

Use ContEst to get the contamination estimate in the data to be passed into M2.  (Note that for these data, the ContEst estimate is on the order of 1% higher than the value used to generate the contaminated data.)

    TEMPDIR=/broad/hptmp/$USER
    BAM1=HCC1143_BL.small.0.04.contaminated.with.SM-612V3.small.0.96.bam
    BAM2=/crsp/picard_aggregation/000007820818/SM-612V4/current/SM-612V4.bam
    OUT_TXT=ContEst_0.04HCC1143inNA12878.txt
    
    java -Djava.io.tmpdir=$TEMPDIR \
    -Xmx512m -jar /xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/ContaminationAnalysis/broadinstitute.org/cancer.genome.analysis/00262/107//Queue-1.4-437-g6b8a9e1-svn-35362.jar \
    -S /xchip/tcga/gdac_prod/applications/process_mgmt/firehose_task_registry/cga/ContaminationAnalysis/broadinstitute.org/cancer.genome.analysis/00262/107//ContaminationPipeline.scala \
    -reference /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta -interval /xchip/cga/reference/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list \
    -out $OUT_TXT \
    -bam $BAM1 -nbam $BAM2 \
    -array /xchip/gcc_data/results2/production/TEENY/BySample/hg19.vcf.txt/TEENY_p_TCGA_302_303_N_GenomeWideSNP_6_F04_1344608.hg19.vcf.txt.store/verstore.00000.TEENY_p_TCGA_302_303_N_GenomeWideSNP_6_F04_1344608.hg19.vcf.txt \
    -pop /xchip/cga/reference/hg19/hg19_population_stratified_af_hapmap_3.3.fixed.vcf -faf true -run -array_interval /xchip/cga/reference/hg19/SNP6.hg19.interval_list
    
    BAM1=HCC1143_BL.small.0.03.contaminated.with.SM-612V3.small.0.97.bam
    BAM2=/crsp/picard_aggregation/000007820818/SM-612V4/current/SM-612V4.bam
    OUT_TXT=ContEst_0.03HCC1143inNA12878.txt
    
    BAM1=HCC1143_BL.small.0.02.contaminated.with.SM-612V3.small.0.98.bam
    BAM2=/crsp/picard_aggregation/000007820818/SM-612V4/current/SM-612V4.bam
    OUT_TXT=ContEst_0.02HCC1143inNA12878.txt
    
    BAM1=HCC1143_BL.small.0.01.contaminated.with.SM-612V3.small.0.99.bam
    BAM2=/crsp/picard_aggregation/000007820818/SM-612V4/current/SM-612V4.bam
    OUT_TXT=ContEst_0.01HCC1143inNA12878.txt
    
   And so on for the other replicates.
   
   ContEst estimates for the four replicates at the five contamination levels are as follows:
   |Sample|Input Contamination Level|ContEst Estimate|
   |------|-------------------------|----------------|
   |SM-612V3|0.01|0.016|
   |SM-612V3|0.02|0.030|
   |SM-612V3|0.03|0.042|
   |SM-612V3|0.04|0.055|
   |SM-612V3|0.05|0.067|
   |SM-612V4|0.01|0.015|
   |SM-612V4|0.02|0.028|
   |SM-612V4|0.03|0.039|
   |SM-612V4|0.04|0.051|
   |SM-612V4|0.05|0.063|
   |SM-612V5|0.01|0.016|
   |SM-612V5|0.02|0.030|
   |SM-612V5|0.03|0.042|
   |SM-612V5|0.04|0.054|
   |SM-612V5|0.05|0.066|
   |SM-612V6|0.01|0.015|
   |SM-612V6|0.02|0.027|
   |SM-612V6|0.03|0.040|
   |SM-612V6|0.04|0.051|
   |SM-612V6|0.05|0.062|
   
###Prepare the inputs for the normal-normal calling script
Create a list of all contamination levels for each replicate

    ls -1 HCC*contam*V3*.bam > HCC1143withNA12878_3.bams.list
    ls -1 HCC*contam*V4*.bam > HCC1143withNA12878_4.bams.list
    ls -1 HCC*contam*V5*.bam > HCC1143withNA12878_5.bams.list
    ls -1 HCC*contam*V6*.bam > HCC1143withNA12878_6.bams.list
    
Create a list of the other, uncontaminated normals to call against

    ls -1 /humgen/gsa-hpprojects/NA12878Collection/bams/crsp_ice_validation/SM-612V[^37D].bam  > /dsde/working/mutect/laura/contamination/NA12878_not3.list
    ls -1 /humgen/gsa-hpprojects/NA12878Collection/bams/crsp_ice_validation/SM-612V[^47D].bam  > /dsde/working/mutect/laura/contamination/NA12878_not4.list
    ls -1 /humgen/gsa-hpprojects/NA12878Collection/bams/crsp_ice_validation/SM-612V[^57D].bam  > /dsde/working/mutect/laura/contamination/NA12878_not5.list
    ls -1 /humgen/gsa-hpprojects/NA12878Collection/bams/crsp_ice_validation/SM-612V[^67D].bam  > /dsde/working/mutect/laura/contamination/NA12878_not6.list
     
###Run the Caller
Run M2 on contaminated bams versus other all other replicates. Run one loop for each contaminated replicate, passing in contamination estimates as given above
   
    java -jar $QUEUE_JAR -S /dsde/working/mutect/laura/contamination/Qscript_M2_normalNormalLoop.scala -normal /dsde/working/mutect/laura/contamination/NA12878_not3.list -tumor /dsde/working/mutect/laura/contamination/HCC1143withNA12878_3.bams.list -o M2_NA12878run3_ -f 0.016 -f 0.03 -f 0.042 -f 0.055 -f 0.067
    
    java -jar $QUEUE_JAR -S /dsde/working/mutect/laura/contamination/Qscript_M2_normalNormalLoop.scala -normal /dsde/working/mutect/laura/contamination/NA12878_not4.list -tumor /dsde/working/mutect/laura/contamination/HCC1143withNA12878_4.bams.list -o M2_NA12878run4_ -f 0.015 -f 0.028 -f 0.039 -f 0.051 -f 0.063
    
    java -jar $QUEUE_JAR -S /dsde/working/mutect/laura/contamination/Qscript_M2_normalNormalLoop.scala -normal /dsde/working/mutect/laura/contamination/NA12878_not5.list -tumor /dsde/working/mutect/laura/contamination/HCC1143withNA12878_5.bams.list -o M2_NA12878run5_ -f 0.016 -f 0.030 -f 0.042 -f 0.054 -f 0.066
    
    java -jar $QUEUE_JAR -S /dsde/working/mutect/laura/contamination/Qscript_M2_normalNormalLoop.scala -normal /dsde/working/mutect/laura/contamination/NA12878_not6.list -tumor /dsde/working/mutect/laura/contamination/HCC1143withNA12878_6.bams.list -o M2_NA12878run6_ -f 0.015 -f 0.027 -f 0.040 -f 0.051 -f 0.062
        
###Count the False Positives        
Pull out passing SNPs not in PON for each contamination level:

    for vcf in M2_NA12878run[0-9]_HCC1143_BL.small.0.01.contaminated.with.SM-612V*.bam.vcf 
    do     
        bedtools intersect -a $vcf -b ICE.corrected.bed | grep PASS | awk '{ if ( length($4) + length($5) == 2) print $0 }' | wc -l
    done
    
    for vcf in M2_NA12878run[0-9]_HCC1143_BL.small.0.02.contaminated.with.SM-612V*.bam.vcf 
    do     
        bedtools intersect -a $vcf -b ICE.corrected.bed | grep PASS | awk '{ if ( length($4) + length($5) == 2) print $0 }' | wc -l
    done
    
    for vcf in M2_NA12878run[0-9]_HCC1143_BL.small.0.03.contaminated.with.SM-612V*.bam.vcf 
    do     
        bedtools intersect -a $vcf -b ICE.corrected.bed | grep PASS | awk '{ if ( length($4) + length($5) == 2) print $0 }' | wc -l
    done
    
    for vcf in M2_NA12878run[0-9]_HCC1143_BL.small.0.04.contaminated.with.SM-612V*.bam.vcf 
    do     
        bedtools intersect -a $vcf -b ICE.corrected.bed | grep PASS | awk '{ if ( length($4) + length($5) == 2) print $0 }' | wc -l
    done
    
    for vcf in M2_NA12878run[0-9]_HCC1143_BL.small.0.05.contaminated.with.SM-612V*.bam.vcf 
    do     
        bedtools intersect -a $vcf -b ICE.corrected.bed | grep PASS | awk '{ if ( length($4) + length($5) == 2) print $0 }' | wc -l
    done
    
(I pasted the results from the terminal into Excel because it's just so easy.) 
   
###Comparison Without Downsampling
To run normal-normal contaminated calling without downsampling, the above /dsde/working/mutect/laura/contamination/Qscript_M2_normalNormalLoop.scala commands can be used, passing in -f 0 for each contamination level instead, e.g.:

    java -jar $QUEUE_JAR -S /dsde/working/mutect/laura/contamination/Qscript_M2_normalNormalLoop.scala -normal /dsde/working/mutect/laura/contamination/NA12878_not3.list -tumor /dsde/working/mutect/laura/contamination/HCC1143withNA12878_3.bams.list -o M2_NA12878run3_noContam_ -f 0.0 -f 0.0 -f 0.0 -f 0.0 -f 0.0 

###Comparison to M1
To run normal-normal contaminated calling using M1, run the above Queue commands using a Queue jar containing MuTect and passing in /dsde/working/mutect/laura/contamination/Qscript_M1_normalNormalLoop.scala instead of Qscript_M2_normalNormalLoop.scala, e.g.:

    java -jar $QUEUE_JAR_WITH_M1 -S /dsde/working/mutect/laura/contamination/Qscript_M1_normalNormalLoop.scala -normal /dsde/working/mutect/laura/contamination/NA12878_not3.list -tumor /dsde/working/mutect/laura/contamination/HCC1143withNA12878_3.bams.list -o M1_NA12878run3_ -f 0.016 -f 0.03 -f 0.042 -f 0.055 -f 0.067

(The MuTect-containing Queue jar can be built from the gsa-unstable branch ldg_MuTect1.)

###Latest Results
|M2 SNPs no correction|M2 SNPs with correction|M1 SNPs no correction|M1 SNPs with correction|M2 INDELs no correction|M2 INDELs with correction|
|---------------------|-----------------------|---------------------|-----------------------|-----------------------|-------------------------|
|0%|93|93|181|181|25|25|
|1%|938|258|854|317|68|30|
|2%|2550|464|1941|385|92|21|
|3%|4171|596|3061|515|134|18|
|4%|5513|707|4002|589|162|21|
|5%|6475|794|4854|624|188|29|