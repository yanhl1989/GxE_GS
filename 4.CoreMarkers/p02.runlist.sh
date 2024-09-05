#!/bin/bash
#SBATCH -N 6
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=10
#SBATCH -p normal
#SBATCH -t 60-5
#SBATCH --mem=100G

runlist=$1
module load apps/singularity/3.6.1
#module load compiler/devtoolset/7.3.1
#module load compiler/rocm/2.9
#module load mpi/hpcx/2.4.1/gcc-7.3.1

# trap "exec 1000>&-;exec 1000<&-;exit 0" 2 #定义Ctrl+C是取消管道符1000

# 创建管道文件，绑定管道号1000到文件，删除文件（防止写入读出只能同步进行）
mkfifo testfifo
exec 1000<>testfifo
rm -fr testfifo

# 通过向管道符1000写入空行来设定后台最大进程数
for ((n=1;n<=$SLURM_JOB_NUM_NODES;n++))
do
    echo >&1000
done

# 记录开始时间
start=`date "+%s"`

if [ -e p02.log ] ; then
    rm p02.log
fi

# 总任务数100，每个任务先从读取管道1000中读取1个空行，读取后管道中就减少一个，管道中没有空行时任务暂停，直到有空行可读，任务完成后向管道写入空行
cat $runlist | while read line
do
    read -u1000
    {
        #echo `date` "Now runing with $line">>p02.log
	srun -N1 -n1 -c10  -o ${line}.out -e ${line}.error sh $line
	#echo `date` "Now Job $line is done">>Job.log
        echo >&1000
    }&
done
#sleep 120m
#sh p03.combin.sh
# 等待所有任务完成后向下进行
wait

sleep 2d 

# 记录结束时间
end=`date "+%s"`

echo "TIME: `expr $end - $start `"

# 删除管道符1000
exec 1000>&-
exec 1000<&-

