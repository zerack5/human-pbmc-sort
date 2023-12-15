#PBS -N testcount06_2
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1200:00:00
#PBS -o 2_8806testresult.logo
#PBS -e 2_8806testerror.loge
#PBS -V
#PBS -S /bin/bash

cd /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/CleanData

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_002 2_002_clean_R1.fq.gz 2_002_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_003 2_003_clean_R1.fq.gz 2_003_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_004 2_004_clean_R1.fq.gz 2_004_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_005 2_005_clean_R1.fq.gz 2_005_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_006 2_006_clean_R1.fq.gz 2_006_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_008 2_008_clean_R1.fq.gz 2_008_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_009 2_009_clean_R1.fq.gz 2_009_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_010 2_010_clean_R1.fq.gz 2_010_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_011 2_011_clean_R1.fq.gz 2_011_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_012 2_012_clean_R1.fq.gz 2_012_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_013 2_013_clean_R1.fq.gz 2_013_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_014 2_014_clean_R1.fq.gz 2_014_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_015 2_015_clean_R1.fq.gz 2_015_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_2_016 2_016_clean_R1.fq.gz 2_016_clean_R2.fq.gz
