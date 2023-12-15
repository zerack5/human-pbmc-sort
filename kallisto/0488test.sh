#PBS -N tquantrun04
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1200:00:00
#PBS -o 8804testresult.logo
#PBS -e 8804testerror.loge
#PBS -V
#PBS -S /bin/bash

cd /lustre/wangjiucunlab/zhangqing/il17a/TPL202301504/Clean_cutData

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL202301504/04_4_014 4_014_clean_R1.fq.gz 4_014_clean_R2.fq.gz

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto quant -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx -o /lustre/wangjiucunlab/zhangqing/il17a/TPL202301504/04_3_005 3_005_clean_R1.fq.gz 3_005_clean_R2.fq.gz
