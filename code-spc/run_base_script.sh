#! /bin/bash -e
#make
#make u_index
echo "compile completed, u_index and u_query are generated.\n"

echo "data 1 SPC index constructing baseline"
#./u_index -g graph/1.txt -l label/1.txt -b n -o degree -t $1

./u_index_base -g graph/$1 -l label/$1 -o sigpath -s y -e y -i y -t $2
#./u_index_base -g graph/1.txt -l label/1.txt -o degree -s y -e y -i y -t $1
echo $1
echo "data 1 SPC querying baseline"
./u_query_base -l label/$1 -q query/$1 -a answer/$1 -t o
echo "-----------------------------------------------------"

#echo "data 1 SCCnt index constructing"
#./u_index -g graph/1_bi.txt -l label/1_bi.txt -b y -o bidegree
#echo "data 1 SCCnt querying"
#./u_query -l label/1_bi.txt -q query/ -a answer/ -t b
#echo "-----------------------------------------------------"

#echo "data 1 SCCnt (removed edges) index constructing"
#./u_index -g graph/1_bi_lack.txt -l label/1_bi_lack.txt -b y -o bidegree
#echo "data 1 SCCnt edge adding and querying without minimality"
#./u_query -l label/1_bi_lack.txt -i graph/1_bi_inc.txt -q query/ -a answer/ -t i
#echo "-----------------------------------------------------"

#echo "data 1 SCCnt (removed edges) index constructing"
#./u_index -g graph/1_bi_lack.txt -l label/1_bi_lack.txt -b y -o bidegree
#echo "data 1 SCCnt edge adding and querying with minimality"
#./u_query -l label/1_bi_lack.txt -i graph/1_bi_inc.txt -m m -q query/ -a answer/ -t i
#echo "-----------------------------------------------------"

echo "All done."

