KMAX=17
KCOUNT=2
while [ $KCOUNT -lt $KMAX ]
do
for i in $(ls -rL --sort=size *.el)
do
graph="${i//@}"
./inc_topk_o3_random $graph $KCOUNT 0 20000
done
KCOUNT=$[2*$KCOUNT]
done
