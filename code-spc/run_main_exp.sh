datasets="fb gw wi go db be yt pe fl tw"
for i in $datasets
do
    sh run_base_script.sh $i.txt 1 > result_${i}_base_3y.txt
done
for i in $datasets
do
   sh run_script.sh $i.txt 32 > result_${i}_3y_32.txt 
done

