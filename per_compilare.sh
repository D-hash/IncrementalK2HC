g++ -O3 -fprefetch-loop-arrays main_random.cpp incremental_topk.h incremental_topk.cpp progressBar.h -o inc_topk_o3_random -lnetworkit
g++ -O3 -fprefetch-loop-arrays main.cpp incremental_topk.h incremental_topk.cpp progressBar.h -o inc_topk_o3_rem_add -lnetworkit
g++ -O3 -fprefetch-loop-arrays main_temporal.cpp incremental_topk.h incremental_topk.cpp progressBar.h -o inc_o3_temporal -lnetworkit
