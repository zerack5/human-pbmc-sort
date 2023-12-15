#PBS -N testcount06_1
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1200:00:00
#PBS -o 1_8806testresult.logo
#PBS -e 1_8806testerror.loge
#PBS -V
#PBS -S /bin/bash

cd /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/CleanData

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_002 1_002_clean_R1.fq.gz 1_002_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_003 1_003_clean_R1.fq.gz 1_003_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_005 1_005_clean_R1.fq.gz 1_005_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_006 1_006_clean_R1.fq.gz 1_006_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_008 1_008_clean_R1.fq.gz 1_008_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_009 1_009_clean_R1.fq.gz 1_009_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_010 1_010_clean_R1.fq.gz 1_010_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_011 1_011_clean_R1.fq.gz 1_011_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_012 1_012_clean_R1.fq.gz 1_012_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_013 1_013_clean_R1.fq.gz 1_013_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_014 1_014_clean_R1.fq.gz 1_014_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_015 1_015_clean_R1.fq.gz 1_015_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_1_016 1_016_clean_R1.fq.gz 1_016_clean_R2.fq.gz
