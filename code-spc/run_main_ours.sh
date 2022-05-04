datasets="fb gw wi go db be yt pe fl tw"
for i in $datasets
do
   sh run_script.sh $i.txt 32 > result_${i}_3y_32.txt 
done

