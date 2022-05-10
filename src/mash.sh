ref=$1
query=$2
out=$3



for i in `ls $ref/*msh`
do
mash dist $i $query/*msh -s 1000 >>$out
echo $i $query
date
done

