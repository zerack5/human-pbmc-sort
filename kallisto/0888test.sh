#PBS -N testcount08
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1200:00:00
#PBS -o 8808testresult.logo
#PBS -e 8808testerror.loge
#PBS -V
#PBS -S /bin/bash

cd /lustre/wangjiucunlab/zhangqing/il17a/TPL202301508/CleanData

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o  /lustre/wangjiucunlab/zhangqing/il17a/TPL202301508/08_1_004 1_004_clean_R1.fq.gz 1_004_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o  /lustre/wangjiucunlab/zhangqing/il17a/TPL202301508/08_3_009 3_009_clean_R1.fq.gz 3_009_clean_R2.fq.gz
