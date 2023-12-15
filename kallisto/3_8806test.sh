#PBS -N testcount06_3
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1200:00:00
#PBS -o 3_8806testresult.logo
#PBS -e 3_8806testerror.loge
#PBS -V
#PBS -S /bin/bash

cd /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/CleanData

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_002 3_002_clean_R1.fq.gz 3_002_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_003 3_003_clean_R1.fq.gz 3_003_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_004 3_004_clean_R1.fq.gz 3_004_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_006 3_006_clean_R1.fq.gz 3_006_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_008 3_008_clean_R1.fq.gz 3_008_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_010 3_010_clean_R1.fq.gz 3_010_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_011 3_011_clean_R1.fq.gz 3_011_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_012 3_012_clean_R1.fq.gz 3_012_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_013 3_013_clean_R1.fq.gz 3_013_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_015 3_015_clean_R1.fq.gz 3_015_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL2022121946/06_3_016 3_016_clean_R1.fq.gz 3_016_clean_R2.fq.gz
