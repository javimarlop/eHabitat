cat hriRes_* |grep -v p|grep 0|grep  -v - > hriRes_tmp.csv
cat hriRes2_* |grep -v p|grep 0|grep  -v - > hriRes2_tmp.csv
cat hriInRes_* |grep -v p|grep 0|grep  -v - > hriInRes_tmp.csv
cat hriInRes2_* |grep -v p|grep 0|grep  -v - > hriInRes2_tmp.csv

cat header.txt hriRes_tmp.csv > hriRes.csv
cat header.txt hriRes2_tmp.csv > hriRes2.csv
cat header.txt hriInRes_tmp.csv > hriInRes.csv
cat header.txt hriInRes2_tmp.csv > hriInRes2.csv

rm *tmp*

