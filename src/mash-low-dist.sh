ref=$1
query=$2
out=$3



for i in `ls $ref/*msh`
do
mash dist -d 0.1 $i -s 10000 $query/*msh >>$out
echo $i $query
date
done

