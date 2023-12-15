#PBS -N test88indexrun1
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1200:00:00
#PBS -o 88testresult.logo
#PBS -e 88testerror.loge
#PBS -V
#PBS -S /bin/bash

cd /lustre/wangjiucunlab/zhangqing/il17a

/lustre/wangjiucunlab/zhangqing/il17a/kallisto/kallisto index -i /lustre/wangjiucunlab/zhangqing/il17a/88hg38allcdna.idx /lustre/wangjiucunlab/zhangqing/il17a/Homo_sapiens.88.GRCh38.cdna.all.fa.gz

