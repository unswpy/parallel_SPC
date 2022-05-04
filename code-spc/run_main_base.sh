datasets="pe fl tw"
for i in $datasets
do
    sh run_base_script.sh $i.txt 1 > result_${i}_base_3y.txt
done

