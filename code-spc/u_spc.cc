#include "u_spc.h"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>
//#include <omp.h>
#include "macros.h"
#include <string.h>
#include <fstream>
#include <vector>

namespace spc {
// construct an index for the graph "const_graph"
void USPCQuery::ValidWithBfsQueries(int nq, std::vector<std::vector<int>>& adj, bool newCount)//nq is the number of queries
{
    std::cout << "enter the valid bfs queries function" << std::endl;
    for(int i = 0; i < nq; i++)
    {
       if(i % 1000 == 0) 
       std::cout << i << std::endl;
       uint32_t s = rand() % n_;
       uint32_t t = rand() % n_;
       std::pair<int, uint32_t> r1 = SpcBfs(adj, s, t);
       //std::cout << "end of spc bfs function" << std::endl;
       //uint32_t r2 = Count(s, t);
       uint32_t r2 = 0;
       if(newCount)
       {
            r2 = Count_New(s,t);//Count(s,t);//Count_New(s, t);
       }
       else
       {
           r2 = Count(s, t);
       }
       if(r1.second != r2)
       {
           std::cout << "error count for " << s << " : " << t << " r1 bfs : " << r1.second << " index count is " << r2 << std::endl;
       }
    }
}

std::pair<int, uint32_t> USPCQuery::SpcBfs(std::vector<std::vector<int>>& adj , uint32_t u, uint32_t v)
{    //first element of result is the distance, and the second is the shortest path counting
    //std::cout << u << " : " << v << " : " << n_ << std::endl;
    std::pair<int, uint32_t> result;
    std::queue<uint32_t> Q;
    int dis[n_+1];
    int count[n_+1];
    memset(dis, -1, sizeof dis);
    memset(count, 0, sizeof count);
    Q.push(u);
    count[u] = 1;
    dis[u] = 0;
    int d = -1;
    bool earlyStop = false;
    while(!Q.empty())
    {
        if(earlyStop) break;
        uint32_t vtx = Q.front(); Q.pop();
        //std::cout << "start of " << vtx << std::endl;
        //std::cout << adj[vtx].size() << "is the size of " << vtx << std::endl;
        for(const uint32_t j: adj[vtx])
        {
            //std::cout << j << " : ";//std::endl;
            if(dis[j] == -1)
            {
                dis[j] = dis[vtx] + 1, Q.push(j);
                if(d > 0 && dis[j] > d){
                    earlyStop = true;
                    break;
                }
                count[j] += count[vtx];
            }
            else if(dis[j] == dis[vtx] + 1)
            {
                count[j] += count[vtx];
            }
            if(j ==  v){
                d = dis[v];
               // earlyStop = true;
               // break;
            }
        }
    }
    result.first = dis[v];
    result.second = count[v];
    return result;
}

void USPCIndex::BuildIndexParallel(const Graph& const_graph, int num_threads)
{
  std::cout << "test" << std::endl;
  const auto init_start = std::chrono::steady_clock::now();
  int queryCount = 0;
  bool debug = false;
  bool hybrid = false;
  bool print_index = true;
  ASSERT(dL_.empty() && cL_.empty() && G_.empty());
  G_ = const_graph;
  omp_set_num_threads(num_threads);
  printf("%d threads in build function \n", num_threads);
  n_ = G_.size();
  reduced_.resize(n_, false);
  local_.resize(n_, false);
  eqm_.resize(n_, 1);
  // reduction
  if (opt_shell_) ShellReduction();
  if (opt_equiv_) EquivReduction();
  // use the significant-path-based ordering
  if (OrderScheme::kSigPath == os_) {
    SigPathIndex();
    return;
  }
  // initialization
  std::cout << "end of init 1" << std::endl;
  dL_.resize(n_);
  cL_.resize(n_);
  std::cout << "end of init 2" << std::endl;
  order_.resize(n_);
  rank_.resize(n_);
  std::cout << "end of init 3" << std::endl;
  // ordering
  (this->*of_[os_])(G_);
  OrderRank();
  if (opt_local_) {
    for (uint32_t u = 0; u < n_; ++u) {
      if (G_[u].size() == 0) continue;
      uint32_t cnt = 0;
      for (const uint32_t v : G_[u]) {
        if (rank_[v] < rank_[u]) ++cnt;
      }
      if (G_[u].size() == cnt) local_[u] = true;
    }
  }
  // some auxiliary structures
  std::vector<uint32_t> dLu;//(n_, UINT32_MAX);
  std::vector<std::vector<uint32_t>> dI;//(n_, UINT32_MAX);
  int LMnum = 0;
  dI.resize(LMnum);
  if(debug) std::cout << "before init DI" << std::endl;
  for(int i = 0; i < LMnum; i++)
  {
      dI[i].resize(n_);
      for(int j = 0; j < n_; j++)
      {
          dI[i][j] = UINT32_MAX;
      }
  }
  if(debug) std::cout << "end init DI" << std::endl;
  std::vector<uint32_t> D(n_, UINT32_MAX);
  std::vector<uint32_t> C(n_, 0);
  // hub pushing
  int np = 0;
  int maxD = 20;//the diameter of the input graph.
  
  const auto start = std::chrono::steady_clock::now();
  std::vector<std::unordered_map<uint32_t, uint32_t>> cands;//<vertex, count> pair, the distance is the d
  std::vector<std::unordered_map<uint32_t, uint32_t>> tempCands;
  const auto end = std::chrono::steady_clock::now();
  auto dif = end - start;
  auto dif_ui = end - start;
  auto dif_erase = end - start;
  auto dif_IC = end - start;
  auto queryTime = end - start;
  auto oneTime = end - start;
  cands.resize(n_);
  tempCands.resize(n_);
  if(debug)
  {
    for(int i = 0; i < n_; i++)
    {
      std::cout << i << " rank is " << order_[i] << std::endl;
    }
  }
  const auto init_end = std::chrono::steady_clock::now();
  auto initTime = init_end - init_start;

  std::cout << "2time for init ms" << std::chrono::duration<double, std::milli>(initTime).count() << std::endl;
  int cur_i = 0;
  int inc = 100*num_threads;
  std::cout << "n is " << n_ << std::endl;
  bool whileEnd= false;
  while(cur_i < n_ &&(!whileEnd))
  {
    if(!hybrid) whileEnd = true;
    if(hybrid)
    {
    for(int c = 0; c < n_;c++)
    {
        cands[c].clear();
        tempCands[c].clear();
    }
    }
    for(int d = 1; d <= maxD ; d++)
   {
        //std::cout << d << " is the current d" << std::endl;
        bool earlyStop = false;
        if(d != 1)
        {
            earlyStop = true;
            for(size_t i = 0; i < n_; ++i)
            {
                if(!cands[i].empty())
                {
                    earlyStop = false;
                    break;
                }
            }
        }
    if(earlyStop)
    {
        std::cout << "d is " << d << " at the early stop" << std::endl;
        break;
    }

  if(debug)
  {
    for(int i = 0; i < n_; i++)
    {
        for(auto iter = dL_[i].begin(); iter != dL_[i].end(); iter++)
        {
            std::cout << i << " v d c" << LEExtractV(*iter) << " : " <<LEExtractD(*iter) << " : " << LEExtractC(*iter) << std::endl;
        }
    }
    if(debug)std::cout << " cL_" << std::endl;
    for(int i = 0; i < n_; i++)
    {
        for(auto iter = cL_[i].begin(); iter != cL_[i].end(); iter++)
        {
            std::cout << i << " v d c" << LEExtractV(*iter) << " : " <<LEExtractD(*iter) << " : " << LEExtractC(*iter) << std::endl;
        }
    }
  }
  
  #pragma omp parallel for private(dLu)
  for(size_t i = 0 ; i < n_ ;i++)
  {
  //for (size_t j = 0; j < inc; ++j)
    //size_t i = j + cur_i;
    //std::cout << "id : " << i << std::endl;
    if(i >= n_ && hybrid) continue;
    //i = j;
    dLu.resize(n_);
    //dI.resize(n_);
    //memset(dLu, UINT32_MAX, sizeof dLu);
    //#pragma omp parallel for
    for(int ei = 0; ei < n_; ei++)
    {
        dLu[ei] = UINT32_MAX;
        //dI[ei] = UINT32_MAX;
    }
    np = omp_get_num_threads();
   // printf("%d")
    const uint32_t u = order_[i];
    //std::cout << u << std::endl;
    //tempCands[u].clear();
    
    if(d == 1)
    {//special case for d==1
        //if(i < cur_i || i >= cur_i + inc)continue;
        if(debug) std::cout << "d == 1 i : u : " << i << " \t " << u << std::endl;
        const auto one_time_beg = std::chrono::steady_clock::now();
        //#pragma omp parallel for
        for(int wi = 0; wi < G_[u].size();wi++)
        //for(const uint32_t w : G_[u])
        {
            const uint32_t w = G_[u][wi];
            if(debug)std::cout << "w: "<< w << "\t" << rank_[w] << std::endl;
            //std::cout << dL_[w].size() << " : " << dL_[u].size() << " : " << w << " : " << u <<std::endl;
            if(rank_[w] >= rank_[u]) continue;
 
            if(rank_[w] < cur_i && hybrid)continue;
            if(rank_[w] >= cur_i+inc && hybrid) continue;
            const auto begQT = std::chrono::steady_clock::now();
            //if(Distance(dL_[w], dL_[u], w, u) == 1) continue;
            int r = rank_[w];
            if(r < LMnum)
            {
                dI[r][u] = 1;// 1 means it could domiate all the later
            }

            queryCount++;
            const auto endQT = std::chrono::steady_clock::now();
            queryTime += endQT - begQT;

            if(debug) std::cout << u << " : " << w << std::endl;
            {
                cands[u].insert(std::make_pair(w,eqm_[w]));
                //if(dL_[u].size() == 0) dL_[u].push_back(LEMerge(w,1,1));
               // else
                if(Distance(dL_[w], dL_[u], w, u) == 1)
                {
                    int low = 0;
                    int len = cL_[u].size();
                    int high = len - 1;
                    while(low <= high)
                    {
                        uint32_t mid = (low + high) / 2;
                        if(rank_[LEExtractV(dL_[u][mid])] < rank_[w])
                        {
                         low = mid + 1;
                        }
                        else
                        {
                            high = mid - 1;
                        }
                    
                }
                cL_[u].insert(cL_[u].begin()+low, LEMerge(w,1,eqm_[w]));// assume there is at most one edge between two nodes and the weight is always one
                }
                else
                {
                int low = 0;
                int len = dL_[u].size();
                int high = len - 1;
                while(low <= high)
                {
                    uint32_t mid = (low + high) / 2;
                    if(rank_[LEExtractV(dL_[u][mid])] < rank_[w])
                    {
                        low = mid + 1;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                    
                }
                //if(high < 0) high = 0;
                if(debug) std::cout << u << " : " << " insert " << low << " : "  << w << std::endl;
                dL_[u].insert(dL_[u].begin() + low, LEMerge(w,1,eqm_[w]));// assume there is at most one edge between two nodes and the weight is always one
                //tempCands[i].insert(std::make_pair(w,1));
                }
            }
        }
        const auto one_time_end = std::chrono::steady_clock::now();
        oneTime += one_time_end - one_time_beg;
    }
    else
    {
        
        //for (const auto e : dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);
        //#pragma omp parallel for
        // for(const uint32_t w : G_[u])
        const auto begIC = std::chrono::steady_clock::now();
        for(int ci = 0; ci < G_[u].size(); ci ++)
        {
            const uint32_t w = G_[u][ci];
            if(rank_[w] < cur_i && hybrid)continue;
            //if(rank_[w] <= rank_[u]) continue;
            //else
            {
               for(auto iter = cands[w].begin(); iter != cands[w].end(); iter++)
               {
                    if(debug)std::cout << u << " : " << iter->first << std::endl;
                    if(rank_[iter->first] >= rank_[u]) continue;
                    if(rank_[iter->first] < cur_i && hybrid)continue;
                   /* const auto begQT = std::chrono::steady_clock::now();
                    uint32_t tempD = Distance(dL_[u], dL_[iter->first], u, iter->first);
                    queryCount++;
                    const auto endQT = std::chrono::steady_clock::now();
                    queryTime += endQT - begQT;
                    if(debug)std::cout << u << " : " << iter->first << " dis is " << tempD << std::endl;
                    if(tempD < d)
                    {
                        continue;
                    }
                    else
                    {*/
                        if(debug)
                        std::cout << "insert the temp cands" << std::endl;
                       //dL_[w].push_back(LEMerge(iter->first,d,iter->second));
                       auto temp_iter = tempCands[u].find(iter->first);
                       if(temp_iter == tempCands[u].end())
                       {
                            tempCands[u].insert(std::make_pair(iter->first, iter->second * eqm_[w]));
                       }
                       else
                       {
                            temp_iter->second += iter->second* eqm_[w];
                       }
                    //}
               }
            }
        }
        const auto endIC = std::chrono::steady_clock::now();
        dif_IC += endIC - begIC;
       // for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;
       
        for (const auto e : dL_[u]) {
            dLu[LEExtractV(e)] = LEExtractD(e);
            //dI[LEExtractV(e)] = LEExtractD(e);
        }

        for(auto iter = tempCands[u].begin(); iter != tempCands[u].end(); )
        {
            const auto begQT = std::chrono::steady_clock::now();
            uint32_t tempD = UINT32_MAX;
            int r = rank_[iter->first];
            //if(r >= rank_[u]) continue;
            if(r < cur_i && hybrid) continue;
            if(r < LMnum)
            {
            if(dI[r][u] == UINT32_MAX)
            {
                tempD = PrunedDistance(dLu, dL_[iter->first],d, u, iter->first);//Distance(dL_[u], dL_[iter->first], u, iter->first);
                //tempD = Distance(dL_[u], dL_[iter->first], u, iter->first);
                dI[r][u] = (tempD == 0) ? 1 : tempD;
            }
            else{
                queryCount --;
                tempD = dI[r][u];
            }
            }
            else
            {

                tempD = PrunedDistance(dLu, dL_[iter->first],d, u, iter->first);//Distance(dL_[u], dL_[iter->first], u, iter->first);
            }
            queryCount++;
            const auto endQT = std::chrono::steady_clock::now();
            queryTime += endQT - begQT;

            if( tempD ==0)//== d)
            {
                int low = 0;
                int len = cL_[u].size();
                int high = len - 1;
                while(low <= high)
                {
                    uint32_t mid = (low + high) / 2;
                    if(rank_[LEExtractV(cL_[u][mid])] < rank_[iter->first])
                    {
                        low = mid + 1;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                    
                }
                //if(high < 0) high = 0;
                cL_[u].insert(cL_[u].begin()+low, LEMerge(iter->first,d,iter->second));
                int r = rank_[iter->first];
                if(r < LMnum)
                {
                    dI[r][u] = 1;// 1 means it could domiate all the later
                }
            }
            else if(tempD == -1 ) //> d)
            {
                int low = 0;
                int len = dL_[u].size();
                int high = len - 1;
                if(debug) std::cout << "insert dL_ low, high len" << low << " : " << high << " : " << len << " "<< std::endl;
                while(low <= high)
                {
                    uint32_t mid = (low + high) / 2;
                    if(rank_[LEExtractV(dL_[u][mid])] < rank_[iter->first])
                    {
                        low = mid + 1;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                    
                }
                //if(high < 0) high = 0;
                if(debug) std::cout << "insert at " << high << std::endl;
                dL_[u].insert(dL_[u].begin() + low,LEMerge(iter->first,d,iter->second));
                 int r = rank_[iter->first];
                if(r < LMnum)
                {
                    dI[r][u] = 1;// 1 means it could domiate all the later
                }
            }
            else
            {
            //queryTime += endQT - begQT;
                const auto beg2 = std::chrono::steady_clock::now();
                iter = tempCands[u].erase(iter);
                const auto end2 = std::chrono::steady_clock::now();
                dif_erase += end2 - beg2;
                continue;
            }
            iter++;
        }
        //for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;
       }
  }
  


    const auto start1 = std::chrono::steady_clock::now();
    #pragma omp barrier
    if(d != 1)
    {
        #pragma omp parallel for
        for(size_t i = 0; i < n_; i++)
        {
            cands[i].clear();
        }
        #pragma omp barrier
        std::vector<bool> used(n_, false);
        #pragma omp parallel for
        for(size_t i = 0; i < n_; i++)
        {
            if(used[order_[i]])
            {
                if(debug)std::cout << i  << " : " << order_[i] << "is used " << std::endl;
            }

            used[order_[i]] = true;
            const uint32_t u = order_[i];
            for(auto iter = tempCands[u].begin(); iter != tempCands[u].end(); iter++)
            {
                if(debug)
                {
                    std::cout << u << ":" << iter->first << " : " << iter->second << std::endl;
                }
                cands[u].insert(std::make_pair(iter->first, iter->second));
            }
        }
        #pragma omp barrier
        #pragma omp parallel for
        for(size_t i = 0; i < n_; i++)
        {
            tempCands[i].clear();
        }
    }
    else
    {
        if(debug)
        {
            for(int i = 0; i < n_; i++)
            {
                for(auto iter = cands[i].begin(); iter != cands[i].end(); iter++)
                {
                    std::cout << i  << " : "<< iter->first << " : " << iter->second << std::endl;
                }
            }
        }
    }

    const auto end1 = std::chrono::steady_clock::now();
    dif += end1 - start1;
    }
    if(hybrid)cur_i += inc;
  }
  std::cout << "time for update the in and out vectors ms" << std::chrono::duration<double, std::milli>(dif).count() << std::endl;
  std::cout << "time for d=1 ms" << std::chrono::duration<double, std::milli>(oneTime).count() << std::endl;
  std::cout << "time for queryTime ms " << std::chrono::duration<double, std::milli>(queryTime).count() << std::endl;
  std::cout << "time for erase  ms " << std::chrono::duration<double, std::milli>(dif_erase).count() << std::endl;
  std::cout << "time for insert candidate  ms " << std::chrono::duration<double, std::milli>(dif_IC).count() << std::endl;
  std::cout << np << " threads in total" << std::endl;
  std::cout << queryCount << " is the number of Distance calls" << std::endl;
}

void USPCIndex::BuildIndexParallel_vector(const Graph& const_graph, int num_threads)
{
    std::cout << "test" << std::endl;
    const auto init_start = std::chrono::steady_clock::now();
    int queryCount = 0;
    bool debug = false;
    bool hybrid = false;
    bool print_index = true;
    ASSERT(dL_.empty() && cL_.empty() && G_.empty());
    G_ = const_graph;
    omp_set_num_threads(num_threads);
    printf("%d threads in build function \n", num_threads);
    n_ = G_.size();
    reduced_.resize(n_, false);
    local_.resize(n_, false);
    eqm_.resize(n_, 1);
    // reduction
    if (opt_shell_) ShellReduction();
    if (opt_equiv_) EquivReduction();
    // use the significant-path-based ordering
    if (OrderScheme::kSigPath == os_) {
        SigPathIndex();
        return;
    }
    // initialization
    std::cout << "end of init 1" << std::endl;
    dL_.resize(n_);
    cL_.resize(n_);
    std::cout << "end of init 2" << std::endl;
    order_.resize(n_);
    rank_.resize(n_);
    std::cout << "end of init 3" << std::endl;
    // ordering
    (this->*of_[os_])(G_);
    OrderRank();
    if (opt_local_) {
        for (uint32_t u = 0; u < n_; ++u) {
            if (G_[u].size() == 0) continue;
            uint32_t cnt = 0;
            for (const uint32_t v : G_[u]) {
                if (rank_[v] < rank_[u]) ++cnt;
            }
            if (G_[u].size() == cnt) local_[u] = true;
        }
    }
    // some auxiliary structures
    std::vector<uint32_t> dLu(n_, UINT32_MAX);
    //std::unordered_map<uint32_t, uint32_t> dLu;
    //std::map<uint32_t, uint32_t> dLu;
    //std::map<uint32_t, uint32_t> dLu;

    std::vector<std::vector<uint32_t>> dI;//(n_, UINT32_MAX);
    int LMnum = 0;
    dI.resize(LMnum);
    if (debug) std::cout << "before init DI" << std::endl;
    for (int i = 0; i < LMnum; i++)
    {
        dI[i].resize(n_);
        for (int j = 0; j < n_; j++)
        {
            dI[i][j] = UINT32_MAX;
        }
    }
    if (debug) std::cout << "end init DI" << std::endl;
    std::vector<uint32_t> D(n_, UINT32_MAX);
    std::vector<uint32_t> C(n_, 0);
    // hub pushing
    int np = 0;
    int maxD = 100;//the diameter of the input graph.

    const auto start = std::chrono::steady_clock::now();
    std::vector<std::unordered_map<uint32_t, uint32_t>> cands;//<vertex, count> pair, the distance is the d
    std::vector<std::unordered_map<uint32_t, uint32_t>> tempCands;
    const auto end = std::chrono::steady_clock::now();
    auto dif = end - start;
    auto dif_ui = end - start;
    auto dif_erase = end - start;
    auto dif_IC = end - start;
    auto queryTime = end - start;
    auto oneTime = end - start;
    cands.resize(n_);
    tempCands.resize(n_);
    if (debug)
    {
        for (int i = 0; i < n_; i++)
        {
            std::cout << i << " rank is " << order_[i] << std::endl;
        }
    }
    const auto init_end = std::chrono::steady_clock::now();
    auto initTime = init_end - init_start;

    std::cout << "2time for init ms" << std::chrono::duration<double, std::milli>(initTime).count() << std::endl;
    int cur_i = 0;
    int inc = 100 * num_threads;
    std::cout << "n is " << n_ << std::endl;
    bool whileEnd = false;
    std::vector<std::vector<LabelEntry>> dL_temp, cL_temp;
    //std::vector<LabelV> dL_temp, cL_temp;
    std::vector<std::vector<std::vector<LabelEntry>>> dL_pre, cL_pre;
    //dL_temp.resize(n_);
    //cL_temp.resize(n_);
    dL_temp.resize(n_);
    cL_temp.resize(n_);
    dL_pre.resize(maxD);
    cL_pre.resize(maxD);
    //#pragma omp declare reduction(MyMerge: std::vector<LabelV>: omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    for (int i = 0; i < maxD; i++)
    {
        dL_pre[i].resize(n_);
        cL_pre[i].resize(n_);
    }

    while (cur_i < n_ && (!whileEnd))
    {
        if (!hybrid) whileEnd = true;
        if (hybrid)
        {
            for (int c = 0; c < n_; c++)
            {
                cands[c].clear();
                tempCands[c].clear();
            }
        }
        for (int d = 1; d <= maxD; d++)
        {
       //     std::cout << d << " is the current d" << std::endl;
            bool earlyStop = false;
            if (d != 1)
            {
                earlyStop = true;
                for (size_t i = 0; i < n_; ++i)
                {
                    if ((!dL_pre[d - 1][i].empty()) || (!cL_pre[d - 1][i].empty()))
                    {
                        earlyStop = false;
                        break;
                    }
                }
            }
            if (debug)
            {
                for (int i = 0; i < n_; i++)
                {
                    for (auto iter = dL_[i].begin(); iter != dL_[i].end(); iter++)
                    {
                        std::cout << i << " v d c" << LEExtractV(*iter) << " : " << LEExtractD(*iter) << " : " << LEExtractC(*iter) << std::endl;
                    }
                }
                if (debug)std::cout << " cL_" << std::endl;
                for (int i = 0; i < n_; i++)
                {
                    for (auto iter = cL_[i].begin(); iter != cL_[i].end(); iter++)
                    {
                        std::cout << i << " v d c" << LEExtractV(*iter) << " : " << LEExtractD(*iter) << " : " << LEExtractC(*iter) << std::endl;
                    }
                }
            }
            if (earlyStop)
            {
                std::cout << "d is " << d << " at the early stop" << std::endl;
                break;
            }


#pragma omp parallel for firstprivate(dLu) schedule(dynamic, 1)
            //#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < n_; i++)
            {
                const uint32_t u = order_[i];
                const uint32_t o_mul = eqm_[u];
                eqm_[u] = 1;
                dL_temp[i].clear();
                cL_temp[i].clear();
                //for (size_t j = 0; j < inc; ++j) 
                  //size_t i = j + cur_i;
                  //std::cout << "id : " << i << std::endl;
                if (i >= n_ && hybrid) continue;
                //i = j;
                //std::vector<uint32_t> dLuT(n_, UINT32_MAX);
                //swap(dLu, dLuT);

                //dLu.clear();


                np = omp_get_num_threads();
                // printf("%d")
                //std::cout << u << std::endl;
                //tempCands[u].clear();
                if (!local_[u])
                {
                    for (const auto e : dL_[u])
                    {
                        dLu[LEExtractV(e)] = LEExtractD(e);
                        //dI[LEExtractV(e)] = LEExtractD(e);
                    }
                }
                if (d == 1)
                {//special case for d==1
                    //if(i < cur_i || i >= cur_i + inc)continue;
                    if (debug) std::cout << "d == 1 i : u : " << i << " \t " << u << std::endl;
                    const auto one_time_beg = std::chrono::steady_clock::now();
                    //#pragma omp parallel for
                    for (int wi = 0; wi < G_[u].size(); wi++)
                        //for(const uint32_t w : G_[u])
                    {
                        const uint32_t w = G_[u][wi];
                        if (debug)std::cout << "w: " << w << "\t" << rank_[w] << std::endl;
                        //std::cout << dL_[w].size() << " : " << dL_[u].size() << " : " << w << " : " << u <<std::endl; 
                        if (rank_[w] >= rank_[u]) continue;

                        if (rank_[w] < cur_i && hybrid)continue;
                        if (rank_[w] >= cur_i + inc && hybrid) continue;
                        //const auto begQT = std::chrono::steady_clock::now();
                        //if(Distance(dL_[w], dL_[u], w, u) == 1) continue;
                        int r = rank_[w];
                        if (r < LMnum)
                        {
                            dI[r][u] = 1;// 1 means it could domiate all the later
                        }

                        queryCount++;
                        //const auto endQT = std::chrono::steady_clock::now();
                       // queryTime += endQT - begQT;

                        if (debug) std::cout << u << " : " << w << std::endl;
                        {
                            //cands[u].insert(std::make_pair(w,eqm_[w]));
                            if (Distance(dL_[w], dL_[u], w, u) == 1)
                            {
                                /*    int low = 0;
                                    int len = cL_[u].size();
                                    int high = len - 1;
                                    while(low <= high)
                                    {
                                        uint32_t mid = (low + high) / 2;
                                        if(rank_[LEExtractV(dL_[u][mid])] < rank_[w])
                                        {
                                         low = mid + 1;
                                        }
                                        else
                                        {
                                            high = mid - 1;
                                        }

                                    }
                                cL_[u].insert(cL_[u].begin()+low, LEMerge(w,1,eqm_[w]));*/
                                cL_temp[i].push_back(LEMerge(w, 1, 1));

                                // assume there is at most one edge between two nodes and the weight is always one
                            }
                            else
                            {
                                /*int low = 0;
                                int len = dL_[u].size();
                                int high = len - 1;
                                while(low <= high)
                                {
                                    uint32_t mid = (low + high) / 2;
                                    if(rank_[LEExtractV(dL_[u][mid])] < rank_[w])
                                    {
                                        low = mid + 1;
                                    }
                                    else
                                    {
                                        high = mid - 1;
                                    }

                                }
                                //if(high < 0) high = 0;
                                if(debug) std::cout << u << " : " << " insert " << low << " : "  << w << std::endl;
                                dL_[u].insert(dL_[u].begin() + low, LEMerge(w,1,eqm_[w]));*/
                                dL_temp[i].push_back(LEMerge(w, 1, 1));
                                // assume there is at most one edge between two nodes and the weight is always one
                            }
                        }
                    }
                    if (debug) std::cout << "end of d == 1" << " u == " << u << std::endl;
                    const auto one_time_end = std::chrono::steady_clock::now();
                    oneTime += one_time_end - one_time_beg;
                }
                else
                {
                    if (debug) std::cout << "test d == " << d << std::endl;
                    //for (const auto e : dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);
                    //#pragma omp parallel for
                    // for(const uint32_t w : G_[u])
                    //const auto begIC = std::chrono::steady_clock::now();
                    std::vector<LabelEntry> inCands;
                    //#pragma omp parallel for
                    for (int ci = 0; ci < G_[u].size(); ci++)
                    {
                        const uint32_t w = G_[u][ci];
                        if (rank_[w] < cur_i && hybrid)continue;
                        //if(rank_[w] <= rank_[u]) continue;
                        //else
                        {
                            for (auto iter : dL_pre[d - 1][w])
                            {
                                inCands.push_back(LEMerge(LEExtractV(iter), LEExtractD(iter), LEExtractC(iter) * eqm_[w]));
                            }
                            for (auto iter : cL_pre[d - 1][w])
                            {
                                inCands.push_back(LEMerge(LEExtractV(iter), LEExtractD(iter), LEExtractC(iter) * eqm_[w]));
                            }
                        }
                    }
                    //qsort(inCands.begin(), inCands.size(),);
                    sort(inCands.begin(), inCands.end());
                    std::vector<LabelEntry> ansd;
                    int ansdI = 0;


                    //if (debug) std::cout << u << " is the current vertex " << std::endl;
                    //#pragma omp parallel for
                    for (int vi = 0; vi < inCands.size(); vi++)
                    {
                        uint32_t v1 = LEExtractV(inCands[vi]);
                        if (debug) std::cout << "v ==" << v1 << std::endl;
                        if (vi == 0)
                        {
                            ansd.push_back(inCands[vi]);
                        }
                        else {
                            uint32_t v2 = LEExtractV(ansd[ansdI]);
                            if (v1 != v2)
                            {
                                ansd.push_back(inCands[vi]);
                                ansdI++;
                            }
                            else
                            {
                                uint32_t d1 = LEExtractD(ansd[ansdI]);
                                uint32_t d2 = LEExtractD(inCands[vi]);
                                uint32_t c1 = LEExtractC(ansd[ansdI]);
                                uint32_t c2 = LEExtractC(inCands[vi]);
                                if (d2 > d1) continue;
                                else if (d2 < d1)
                                {
                                    ansd[ansdI] = LEMerge(v1, d2, c2);
                                }
                                else {
                                    c1 += LEExtractC(inCands[vi]);
                                    ansd[ansdI] = LEMerge(v1, d1, c1);
                                }
                            }

                        }

                    }
                    swap(inCands, ansd);
                    ansd.clear();

                    /*std::vector<std::vector<LabelEntry>> dL_t, cL_t;

                    dL_t.resize(num_threads+1);
                    cL_t.resize(num_threads+1);*/
                    //#pragma omp parallel for //reduction(+= : dL_temp[i]) reduction(+=: cL_temp[i])
                    for (int vi = 0; vi < inCands.size(); vi++)//auto iter : inCands)
                    {
                        //int pid = omp_get_thread_num();
                        if (!local_[u])
                        {
                            auto iter = inCands[vi];
                            uint32_t vt = LEExtractV(iter);
                            uint32_t ct = LEExtractC(iter);
                            if (rank_[vt] >= rank_[u]) continue;
                            if (rank_[vt] < cur_i && hybrid)continue;
                            // const auto begQT = std::chrono::steady_clock::now();
                            uint32_t tempD = UINT32_MAX;
                            int r = rank_[vt];
                            //if(r >= rank_[u]) continue;
                            if (r < cur_i && hybrid) continue;
                            if (r < LMnum)
                            {
                                if (dI[r][u] == UINT32_MAX)
                                {
                                    tempD = PrunedDistance(dLu, dL_[vt], d, u, vt);//Distance(dL_[u], dL_[iter->first], u, iter->first);
                                    //tempD = Distance(dL_[u], dL_[iter->first], u, iter->first);
                                    dI[r][u] = (tempD == 0) ? 1 : tempD;
                                }
                                else {
                                    queryCount--;
                                    tempD = dI[r][u];
                                }
                            }
                            else
                            {

                                tempD = PrunedDistance(dLu, dL_[vt], d, u, vt);//Distance(dL_[u], dL_[iter->first], u, iter->first);
                            }
                            queryCount++;
                            // const auto endQT = std::chrono::steady_clock::now();
                            // queryTime += endQT - begQT;

                            if (tempD == 0)//== d)
                            {
                                cL_temp[i].push_back(LEMerge(vt, d, ct));
                                int r = rank_[vt];
                                if (r < LMnum)
                                {
                                    dI[r][u] = 1;// 1 means it could domiate all the later
                                }
                            }
                            else if (tempD == -1) //> d)
                            {
                                dL_temp[i].push_back(LEMerge(vt, d, ct));
                                int r = rank_[vt];
                                if (r < LMnum)
                                {
                                    dI[r][u] = 1;// 1 means it could domiate all the later
                                }
                            }
                        }
                        else {
                            auto iter = inCands[vi];
                            uint32_t vt = LEExtractV(iter);
                            uint32_t ct = LEExtractC(iter);
                            dL_temp[i].push_back(LEMerge(vt, d, ct));
                        }
                    }
                    /*for (int vi = 0; vi < num_threads; vi++)//auto iter : inCands)
                    {
                        cL_temp[i].insert(cL_temp[i].end(),cL_t[vi].begin(), cL_t[vi].end());
                        dL_temp[i].insert(dL_temp[i].end(),dL_t[vi].begin(), dL_t[vi].end());
                    }*/
                }




                //const auto endIC = std::chrono::steady_clock::now();
                //dif_IC += endIC - begIC;
                // for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;



                 /*for(auto iter = tempCands[u].begin(); iter != tempCands[u].end(); )
                 {
                     const auto begQT = std::chrono::steady_clock::now();
                     uint32_t tempD = UINT32_MAX;
                     int r = rank_[iter->first];
                     //if(r >= rank_[u]) continue;
                     if(r < cur_i && hybrid) continue;
                     if(r < LMnum)
                     {
                     if(dI[r][u] == UINT32_MAX)
                     {
                         tempD = PrunedDistance(dLu, dL_[iter->first],d, u, iter->first);//Distance(dL_[u], dL_[iter->first], u, iter->first);
                         //tempD = Distance(dL_[u], dL_[iter->first], u, iter->first);
                         dI[r][u] = (tempD == 0) ? 1 : tempD;
                     }
                     else{
                         queryCount --;
                         tempD = dI[r][u];
                     }
                     }
                     else
                     {
                         tempD = PrunedDistance(dLu, dL_[iter->first],d, u, iter->first);//Distance(dL_[u], dL_[iter->first], u, iter->first);
                     }
                     queryCount++;
                     const auto endQT = std::chrono::steady_clock::now();
                     queryTime += endQT - begQT;
                     if( tempD ==0)//== d)
                     {
                        cL_temp.push_back(LEMerge(iter->first,d,iter->second));
                         int r = rank_[iter->first];
                         if(r < LMnum)
                         {
                             dI[r][u] = 1;// 1 means it could domiate all the later
                         }
                     }
                     else if(tempD == -1 ) //> d)
                     {
                        dL_temp.push_back(LEMerge(iter->first,d,iter->second));
                          int r = rank_[iter->first];
                         if(r < LMnum)
                         {
                             dI[r][u] = 1;// 1 means it could domiate all the later
                         }
                     }
                     else
                     {
                     //queryTime += endQT - begQT;
                         const auto beg2 = std::chrono::steady_clock::now();
                         iter = tempCands[u].erase(iter);
                         const auto end2 = std::chrono::steady_clock::now();
                         dif_erase += end2 - beg2;
                         continue;
                     }
                     iter++;
                 }*/
                 //for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;
                if (debug) std::cout << "enter the sort" << std::endl;
                sort(dL_temp[i].begin(), dL_temp[i].end());
                std::vector<LabelEntry> ansd;
                int ansdI = 0;
                if (debug) std::cout << u << " is the current vertex " << std::endl;
                //#pragma omp parallel for
                for (int vi = 0; vi < dL_temp[i].size(); vi++)
                {
                    uint32_t v1 = LEExtractV(dL_temp[i][vi]);
                    if (debug) std::cout << "v ==" << v1 << std::endl;
                    if (vi == 0)
                    {
                        ansd.push_back(dL_temp[i][vi]);
                    }
                    else {
                        uint32_t v2 = LEExtractV(ansd[ansdI]);
                        if (v1 != v2)
                        {
                            ansd.push_back(dL_temp[i][vi]);
                            ansdI++;
                        }
                        else
                        {
                            uint32_t d1 = LEExtractD(ansd[ansdI]);
                            uint32_t d2 = LEExtractD(dL_temp[i][vi]);
                            uint32_t c1 = LEExtractC(ansd[ansdI]);
                            uint32_t c2 = LEExtractC(dL_temp[i][vi]);
                            if (d2 > d1) continue;
                            else if (d2 < d1)
                            {
                                ansd[ansdI] = LEMerge(v1, d2, c2);
                            }
                            else {
                                c1 += LEExtractC(dL_temp[i][vi]);
                                ansd[ansdI] = LEMerge(v1, d1, c1);
                            }
                        }

                    }

                }
                swap(dL_temp[i], ansd);
                ansd.clear();

                //remove duplicate elements
                if (debug) std::cout << "enter the sort cL" << std::endl;
                sort(cL_temp[i].begin(), cL_temp[i].end());
                std::vector<LabelEntry> ans;
                int ansI = 0;
                //#pragma omp parallel for
                for (int vi = 0; vi < cL_temp[i].size(); vi++)
                {
                    if (vi == 0)
                    {
                        ans.push_back(cL_temp[i][vi]);
                    }
                    else {
                        uint32_t v1 = LEExtractV(cL_temp[i][vi]);
                        uint32_t v2 = LEExtractV(ans[ansI]);
                        if (v1 != v2)
                        {
                            ans.push_back(cL_temp[i][vi]);
                            ansI++;
                        }
                        else
                        {
                            uint32_t d1 = LEExtractD(ans[ansI]);
                            uint32_t c1 = LEExtractC(ans[ansI]);
                            c1 += LEExtractC(cL_temp[i][vi]);
                            ans[ansI] = LEMerge(v1, d1, c1);
                        }

                    }

                }
                swap(cL_temp[i], ans);
                if (!local_[u])
                {
                    for (const auto e : dL_[u])
                    {
                        dLu[LEExtractV(e)] = UINT32_MAX;
                        //dI[LEExtractV(e)] = LEExtractD(e);
                    }
                }
                eqm_[u] = o_mul;
            }



            const auto start1 = std::chrono::steady_clock::now();
            //#pragma omp barrier
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < n_; i++)
            {
                const uint32_t u = order_[i];

                if (debug) std::cout << "enter the insert dL" << std::endl;
                if (!local_[u])
                {
                    for (auto e : dL_temp[i])
                    {
                        if (debug) std::cout << "dL_temp v d c " << LEExtractV(e) << " \t " << LEExtractD(e) << "\t" << LEExtractC(e) << std::endl;
                        int tv = LEExtractV(e);
                        int low = 0;
                        int len = dL_[u].size();
                        int high = len - 1;
                        if (debug) std::cout << "insert dL_ low, high len" << low << " : " << high << " : " << len << " " << std::endl;
                        while (low <= high)
                        {
                            uint32_t mid = (low + high) / 2;
                            if (rank_[LEExtractV(dL_[u][mid])] < rank_[tv])
                            {
                                low = mid + 1;
                            }
                            else
                            {
                                high = mid - 1;
                            }

                        }
                        //if(high < 0) high = 0;
                        if (debug) std::cout << "insert at " << high << std::endl;
                        dL_[u].insert(dL_[u].begin() + low, e);

                    }
                    for (auto e : cL_temp[i])
                    {
                        int tv = LEExtractV(e);
                        int low = 0;
                        int len = cL_[u].size();
                        int high = len - 1;
                        if (debug) std::cout << "insert dL_ low, high len" << low << " : " << high << " : " << len << " " << std::endl;
                        while (low <= high)
                        {
                            uint32_t mid = (low + high) / 2;
                            if (rank_[LEExtractV(cL_[u][mid])] < rank_[tv])
                            {
                                low = mid + 1;
                            }
                            else
                            {
                                high = mid - 1;
                            }

                        }
                        //if(high < 0) high = 0;
                        if (debug) std::cout << "insert at " << high << std::endl;
                        cL_[u].insert(cL_[u].begin() + low, e);

                    }
                }
                swap(dL_temp[i], dL_pre[d][u]);
                swap(cL_temp[i], cL_pre[d][u]);
                if (debug) std::cout << "end the sort cL" << std::endl;
            }
            const auto end1 = std::chrono::steady_clock::now();
            dif += end1 - start1;
        }
        if (hybrid)cur_i += inc;
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n_; i++)
    {
        if (local_[i]) continue;
        dL_[i].push_back(LEMerge(i, 0, 1));
        //cL_[i].push_back(LEMerge(i, 0, 1));
    }
    std::cout << "time for update the in and out vectors ms" << std::chrono::duration<double, std::milli>(dif).count() << std::endl;
    std::cout << "time for d=1 ms" << std::chrono::duration<double, std::milli>(oneTime).count() << std::endl;
    std::cout << "time for queryTime ms " << std::chrono::duration<double, std::milli>(queryTime).count() << std::endl;
    std::cout << "time for erase  ms " << std::chrono::duration<double, std::milli>(dif_erase).count() << std::endl;
    std::cout << "time for insert candidate  ms " << std::chrono::duration<double, std::milli>(dif_IC).count() << std::endl;
    std::cout << np << " threads in total" << std::endl;
    std::cout << queryCount << " is the number of Distance calls" << std::endl;
}

void USPCIndex::BuildIndex(const Graph& const_graph) {
  ASSERT(dL_.empty() && cL_.empty() && G_.empty());
  G_ = const_graph;
  int queryCount = 0;
  //int temp = 2;
  //omp_set_num_threads(2);
  n_ = G_.size();
  reduced_.resize(n_, false);
  local_.resize(n_, false);
  eqm_.resize(n_, 1);
  // reduction
  if (opt_shell_) ShellReduction();
  if (opt_equiv_) EquivReduction();
  // use the significant-path-based ordering
  if (OrderScheme::kSigPath == os_) {
    SigPathIndex();
    return;
  }
  // initialization
  dL_.resize(n_);
  cL_.resize(n_);
  order_.resize(n_); // order[i]: i-th rank -> node
  rank_.resize(n_); // rank[i]: i-th node -> rank
  // ordering
  (this->*of_[os_])(G_);
  OrderRank();
  if (opt_local_) {
    for (uint32_t u = 0; u < n_; ++u) {
      if (G_[u].size() == 0) continue;
      uint32_t cnt = 0;
      for (const uint32_t v : G_[u]) {
        if (rank_[v] < rank_[u]) ++cnt;
      }
      if (G_[u].size() == cnt) local_[u] = true;
    }
  }
  // some auxiliary structures
  std::vector<uint32_t> dLu(n_, UINT32_MAX);
  std::vector<uint32_t> D(n_, UINT32_MAX);
  std::vector<uint32_t> C(n_, 0);
  // hub pushing
  int np = 0;
  #pragma omp parallel for num_threads(1)
  for (size_t i = 0; i < n_; ++i) {
    np = omp_get_num_threads();
    const uint32_t u = order_[i];
    const uint32_t o_mul = eqm_[u];
    eqm_[u] = 1;
    // for fast distance computation
    for (const auto e : dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);
    // bfs
    std::vector<uint32_t> reset({u});
    std::queue<uint32_t> Q({u});
    D[u] = 0; C[u] = 1;
    while (!Q.empty()) {
      const uint32_t v = Q.front(); Q.pop();
      const uint32_t dSoFar = Distance(dLu, dL_[v]);
      queryCount ++;
      if (D[v] > dSoFar) continue;
      // add a corresponding entry
      NormalD(D[v]); NormalC(C[v]);
      if (!local_[v]) {
        (D[v] < dSoFar? dL_[v] : cL_[v]).push_back(LEMerge(u, D[v], C[v]));
      }
      // correct C[v]
      if (likely(kUBC / eqm_[v] > C[v])) C[v] *= eqm_[v];
      else C[v] = kUBC;
      for (const uint32_t w : G_[v]) {
        if (rank_[w] <= rank_[u]) continue;
        if (UINT32_MAX == D[w]) {
          D[w] = D[v] + 1;
          C[w] = C[v];
          Q.push(w);
          reset.push_back(w);
        } else if (D[w] == D[v] + 1) {
          if (likely(kUBC - C[v] >= C[w])) C[w] += C[v];
          else C[w] = kUBC;
        }
      }
    }
    // clear
    eqm_[u] = o_mul;
    //#pragma omp parallel for
    for (const uint32_t v : reset) {
      D[v] = UINT32_MAX; C[v] = 0;
    }
    //#pragma omp parallel for
    for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;
  }
  printf("%d threads in total in the build function", np);
  std::cout << queryCount << " is the number of Distance calls" << std::endl;
}

void USPCIndex::SigPathIndex() {
  std::cout << "significant-path-based ordering" << std::endl;
  // initialization
  dL_.resize(n_);
  cL_.resize(n_);
  // some auxiliary structures for BFS
  std::vector<uint32_t> dLu(n_, UINT32_MAX);
  std::vector<uint32_t> D(n_, UINT32_MAX);
  std::vector<uint32_t> C(n_, 0);
  // breadth-first tree
  std::vector<uint32_t> parent(n_, UINT32_MAX);
  std::vector<std::vector<uint32_t>> children(n_);
  std::vector<uint32_t> des(n_, 1);
  // sort the vertices by degree
  std::vector<uint32_t> deg_order(n_);
  std::iota(deg_order.begin(), deg_order.end(), 0);
  std::sort(deg_order.begin(), deg_order.end(),
            [this](const uint32_t v1, const uint32_t v2) {
              return G_[v1].size() > G_[v2].size();
            });
  // pruned bfs's
  uint32_t u = UINT32_MAX;
  uint32_t imp = 0;
  std::vector<bool> chosen(n_, false);
  for (size_t i = 0; i < n_; ++i) {
    if (UINT32_MAX == u) {
      while (chosen[deg_order[imp]]) ++imp;
      u = deg_order[imp];
    }
    ASSERT(!chosen[u]);
    chosen[u] = true;
    order_.push_back(u);
    // temporily set the multiplicity to 1
    const uint32_t o_mul = eqm_[u];
    eqm_[u] = 1;
    // for fast distance computation
    for (const auto e : dL_[u]) dLu[LEExtractV(e)] = LEExtractD(e);
    // bfs
    std::vector<uint32_t> reset({u});
    std::queue<uint32_t> Q({u});
    std::vector<uint32_t> R;
    D[u] = 0; C[u] = 1;
    while (!Q.empty()) {
      const uint32_t v = Q.front(); Q.pop();
      R.push_back(v);
      const uint32_t dSoFar = Distance(dLu, dL_[v]);
      if (D[v] > dSoFar) continue;
      // add a corresponding entry
      NormalD(D[v]); NormalC(C[v]);
      (D[v] < dSoFar? dL_[v] : cL_[v]).push_back(LEMerge(u, D[v], C[v]));
      // correct C[v]
      if (likely(kUBC / eqm_[v] > C[v])) C[v] *= eqm_[v];
      else C[v] = kUBC;
      for (const uint32_t w : G_[v]) {
        if (chosen[w]) continue;
        if (UINT32_MAX == D[w]) {
          // breadth-first tree update
          parent[w] = v;
          children[v].push_back(w);
          // update other information
          D[w] = D[v] + 1;
          C[w] = C[v];
          Q.push(w);
          reset.push_back(w);
        } else if (D[w] == D[v] + 1) {
          if (likely(kUBC - C[v] >= C[w])) C[w] += C[v];
          else C[w] = kUBC;
        }
      }
    }

    // identify the significant path
    while (!R.empty()) {
      const uint32_t v = R.back(); R.pop_back();
      if (UINT32_MAX == parent[v]) continue;
      des[parent[v]] += des[v];
    }
    std::vector<uint32_t> sig({u});
    uint32_t cur = u;
    while (!children[cur].empty()) {
      uint32_t largest_des = 0;
      uint32_t next = 0;
      for (const uint32_t cand : children[cur]) {
        if (des[cand] > largest_des) {
          next = cand;
          largest_des = des[cand];
        }
      }
      cur = next;
      sig.push_back(cur);
    }
    // clean first
    eqm_[u] = o_mul;
    for (const auto e : dL_[u]) dLu[LEExtractV(e)] = UINT32_MAX;

    if (opt_local_ && G_[u].size() > 0) {
      uint32_t cnt = 0;
      for (const uint32_t v : G_[u]) {
        if (chosen[v]) ++cnt;
      }
      if (G_[u].size() == cnt) {
        local_[u] = true;
        std::vector<LabelEntry>().swap(dL_[u]);
        std::vector<LabelEntry>().swap(cL_[u]);
      }
    }
    // identify the significant vertex
    u = UINT32_MAX;
    uint64_t largest_score = 0;
    for (size_t i = 1; i < sig.size(); ++i) {
      const uint32_t v = sig[i];
      if (chosen[v]) continue;
      const uint64_t score = static_cast<uint64_t>(G_[v].size()) *
                             (des[parent[v]] - des[v]);
      if (score > largest_score) {
        largest_score = score;
        u = v;
      }
    }
    // clear
    for (const uint32_t v : reset) {
      D[v] = UINT32_MAX; C[v] = 0;
      parent[v] = UINT32_MAX;
      children[v].clear();
      des[v] = 1;
    }
  }
  rank_ = std::vector<uint32_t>(n_);
  OrderRank();
}

uint64_t USPCQuery::Count(uint32_t v1, uint32_t v2) const {
  //ASSERT(v1 != v2);
  if(v1 == v2) return 1;// 1 for the self loops
  if (opt_shell_) {
    v1 = shr_[v1];
    v2 = shr_[v2];
    if (v1 == v2) return 1;
  }
  if (opt_equiv_) {
    v1 = eqr_[v1];
    v2 = eqr_[v2];
    if (v1 == v2) return eqc_[v1];
  }
  // count the # of shortest paths
  uint32_t sp_d = UINT32_MAX;
  uint64_t sp_c = 0;
  if (!opt_local_) {
    size_t p1 = 0, p2 = 0;
    while (p1 < cL_[v1].size() && p2 < cL_[v2].size()) {
      const uint32_t w1 = LEExtractV(cL_[v1][p1]);
      const uint32_t w2 = LEExtractV(cL_[v2][p2]);
      if (rank_[w1] < rank_[w2]) ++p1;
      else if (rank_[w1] > rank_[w2]) ++p2;
      else {
        const uint32_t d = LEExtractD(cL_[v1][p1]) +
            LEExtractD(cL_[v2][p2]);
        if (d < sp_d) {
          sp_d = d;
          sp_c = static_cast<uint64_t>(LEExtractC(cL_[v1][p1])) *
              LEExtractC(cL_[v2][p2]);
          if (w1 != v1 && w1 != v2) sp_c *= eqm_[w1];
        } else if (d == sp_d) {
          uint64_t c = static_cast<uint64_t>(LEExtractC(cL_[v1][p1])) *
              LEExtractC(cL_[v2][p2]);
          if (w1 != v1 && w1 != v2) c *= eqm_[w1];
          sp_c += c;
        }
        ++p1; ++p2;
      }
    }
  } else {
    if (local_[v1]) std::swap(v1, v2);
    std::vector<uint32_t> dirty;
    std::vector<uint32_t> R1 = {v1}, R2 = {v2};
    if (local_[v1]) R1 = G_[v1];
    if (local_[v2]) R2 = G_[v2];

    if (local_[v2]) {
      for (const uint32_t r1 : R1) {
        for (const auto e : dL_[r1]) {
          const uint32_t w = LEExtractV(e);
          const uint32_t d = LEExtractD(e);
          if (UINT32_MAX == dH_[w]) dirty.push_back(w);
          if (d < dH_[w]) dH_[w] = d;
        }
      }
      uint32_t s2 = 0;
      for (const uint32_t r2 : R2) {
        uint32_t min_d = UINT32_MAX;
        for (const auto e : dL_[r2]) {
          const uint32_t w = LEExtractV(e);
          if (UINT32_MAX == dH_[w]) continue;
          const uint32_t d = LEExtractD(e);
          if (dH_[w] + d < min_d) {
            min_d = dH_[w] + d;
          }
        }
        if (min_d < sp_d) {
          sp_d = min_d;
          R2[0] = r2;
          s2 = 1;
        } else if (min_d == sp_d) {
          R2[s2++] = r2;
        }
      }
      R2.resize(s2);
      for (const uint32_t w : dirty) dH_[w] = UINT32_MAX;
      dirty.clear();
    }

    if (local_[v1]) {
      for (const uint32_t r2 : R2) {
        for (const auto e : dL_[r2]) {
          const uint32_t w = LEExtractV(e);
          const uint32_t d = LEExtractD(e);
          if (UINT32_MAX == dH_[w]) dirty.push_back(w);
          if (d < dH_[w]) dH_[w] = d;
        }
      }
      uint32_t s1 = 0;
      for (const uint32_t r1 : R1) {
        uint32_t min_d = UINT32_MAX;
        for (const auto e : dL_[r1]) {
          const uint32_t w = LEExtractV(e);
          if (UINT32_MAX == dH_[w]) continue;
          const uint32_t d = LEExtractD(e);
          if (dH_[w] + d < min_d) min_d = dH_[w] + d;
        }
        if (min_d == sp_d) R1[s1++] = r1;
      }
      R1.resize(s1);
      for (const uint32_t w : dirty) dH_[w] = UINT32_MAX;
      dirty.clear();
    }

    for (const uint32_t r1 : R1) {
      for (const auto e : dL_[r1]) {
        const uint32_t w = LEExtractV(e);
        const uint32_t d = LEExtractD(e);
        uint64_t c = LEExtractC(e);
        if (UINT32_MAX == dH_[w]) dirty.push_back(w);
        if (d < dH_[w]) {
          dH_[w] = d;
          if (w != r1 && r1 != v1) c *= eqm_[r1];
          cH_[w] = c;
        } else if (d == dH_[w]) {
          if (w != r1 && r1 != v1) c *= eqm_[r1];
          cH_[w] += c;
        }
      }
      for (const auto e : cL_[r1]) {
        const uint32_t w = LEExtractV(e);
        const uint32_t d = LEExtractD(e);
        uint64_t c = LEExtractC(e);
        if (UINT32_MAX == dH_[w]) dirty.push_back(w);
        if (d < dH_[w]) {
          dH_[w] = d;
          if (w != r1 && r1 != v1) c *= eqm_[r1];
          cH_[w] = c;
        } else if (d == dH_[w]) {
          if (w != r1 && r1 != v1) c *= eqm_[r1];
          cH_[w] += c;
        }
      }
    }

    for (const uint32_t r2 : R2) {
      for (const auto e : dL_[r2]) {
        const uint32_t w = LEExtractV(e);
        const uint32_t d = LEExtractD(e);
        if (UINT32_MAX == dH_[w]) continue;
        if (dH_[w] + d < sp_d) {
          sp_d = dH_[w] + d;
          uint64_t c = LEExtractC(e);
          if (w != r2 && r2 != v2) c *= eqm_[r2];
          c *= cH_[w];
          if (w != v1 && w != v2) c *= eqm_[w];
          sp_c = c;
        } else if (dH_[w] + d == sp_d) {
          uint64_t c = LEExtractC(e);
          if (w != r2 && r2 != v2) c *= eqm_[r2];
          c *= cH_[w];
          if (w != v1 && w != v2) c *= eqm_[w];
          sp_c += c;
        }
      }
      for (const auto e : cL_[r2]) {
        const uint32_t w = LEExtractV(e);
        const uint32_t d = LEExtractD(e);
        if (UINT32_MAX == dH_[w]) continue;
        if (dH_[w] + d < sp_d) {
          sp_d = dH_[w] + d;
          uint64_t c = LEExtractC(e);
          if (w != r2 && r2 != v2) c *= eqm_[r2];
          c *= cH_[w];
          if (w != v1 && w != v2) c *= eqm_[w];
          sp_c = c;
        } else if (dH_[w] + d == sp_d) {
          uint64_t c = LEExtractC(e);
          if (w != r2 && r2 != v2) c *= eqm_[r2];
          c *= cH_[w];
          if (w != v1 && w != v2) c *= eqm_[w];
          sp_c += c;
        }
      }
    }
    for (const uint32_t v : dirty) dH_[v] = UINT32_MAX, cH_[v] = 0;
  }
  return sp_c;
}
uint32_t USPCQuery::Count_New(uint32_t u, uint32_t v) const 
{
  uint32_t v1 = u;
  uint32_t v2 = v;
  bool debug = false;
  uint32_t spd = UINT32_MAX;
  uint32_t count = 0;
  if(u == v) return 1;
  if (opt_shell_) {
    u = shr_[u];
    v = shr_[v];
    if (u == v) return 1;
  }
  if (opt_equiv_) {
    u = eqr_[u];
    v = eqr_[v];
   if(u == v) return eqc_[u]; 
  }
  
  size_t p1 = 0, p2 = 0;
  //std::cout << "before enter count new" << std::endl;
  //std::cout << "enter count new" << std::endl;
  auto dLu = cL_[u];
  auto dLv = cL_[v];
  if( (dLu.size() == 0 &&  dLv.size() > 0) ||  (dLu.size() > 0 &&  dLv.size() == 0) )
  {
  if(debug)std::cout << "enter spd query" << std::endl;
      p1 = 0;
      p2 = 0;
      while(p1 < dLu.size())
      {
        const uint32_t w1 = LEExtractV(dLu[p1]);
        const uint32_t d1 = LEExtractD(dLu[p1]);
        if( w1 == v || w1 == u )
        {
            if(d1 < spd)
            {
                spd = d1;
                count = LEExtractC(dLu[p1]);
            }
            else if(d1 == spd)
            {
                count += LEExtractC(dLu[p1]);
            }
        }
        ++p1;
      }
      while(p2 < dLv.size())
      {
        const uint32_t w2 = LEExtractV(dLv[p2]);
        const uint32_t d2 = LEExtractD(dLv[p2]);
        if(w2 == u || w2 == v)
        {
            if(d2 < spd)
            {
                spd = d2;
                count = LEExtractC(dLv[p2]);
            }
            else if(d2 == spd)
            {
                count += LEExtractC(dLv[p2]);
            }
        }
        ++p2;
      }
  }
  else
  {
  if(debug)std::cout << "enter spd query double " << std::endl;
  while(p1 < dLu.size() && p2 < dLv.size())
  {
    const uint32_t w1 = LEExtractV(dLu[p1]);
    const uint32_t w2 = LEExtractV(dLv[p2]);
    const uint32_t d1 = LEExtractD(dLu[p1]);
    const uint32_t d2 = LEExtractD(dLv[p2]);
    const uint32_t c1 = LEExtractC(dLu[p1]);
    const uint32_t c2 = LEExtractC(dLv[p2]);
    //std::cout << "v d c :" << w1 << " : " << d1 << " : " << c1 << " the second " << w2 << " : " << d2 << " : " << c2 << std::endl;
    if(w1 == v || w1 == u)
    {
        if(d1 < spd)
        {
            spd = d1;
            count = c1;
        }
        else if(spd == d1)
        {
            count += c1;
        }
        p1++;
        continue;
    }
    if( w2 == u|| w2 == v)
    {
        if(d2 < spd)
        {
            spd = d2;
            count = c2;
        }
        else if(d2 == spd)
        {
            count += c2;
        }
        p2++;
        continue;
    }
    if(rank_[w1] < rank_[w2]) ++p1;
    else if(rank_[w1] > rank_[w2]) ++p2;
    else
    {
        const uint32_t d = d1 + d2;//LEExtractD(dLu[p1]) + LEExtractD(dLv[p2]);
        if(d < spd)
        {
            spd = d;
            count = c1 * c2;
        }
        else if(d == spd)
        {
            count += c1 * c2;
        }

        ++p1; ++p2;
    }
  }

  if(debug)std::cout << "before enter the single entry query, the distance for" << u << " : " << v << " is " << spd << std::endl;
  while(p1 < dLu.size())
  {
    if(debug)std::cout << "enter dLu" << u << " : " << dLu.size() << " : " << p1 << std::endl;
    const uint32_t w1 = LEExtractV(dLu[p1]);
    const uint32_t d1 = LEExtractD(dLu[p1]);
    const uint32_t c1 = LEExtractC(dLu[p1]);
    if( w1 == v || w1 == u)
    {
        if(d1 < spd)
        {
            spd = d1;
            count = c1;
        }
        else if(d1 == spd)
        {
            count += c1;
        }
    }
    ++p1;
  }
  while(p2 < dLv.size())
  {
    if(debug)std::cout << "enter dLv" << v << " : " << dLv.size() << " : " << p2 << std::endl;
    const uint32_t w2 = LEExtractV(dLv[p2]);
    const uint32_t d2 = LEExtractD(dLv[p2]);
    const uint32_t c2 = LEExtractC(dLv[p2]);
    if(w2 == u || w2 == v)
    {
        if(d2 < spd)
        {
            spd = d2;
            count = c2;
        }
        else if(d2 == spd)
        {
            count += c2;
        }
    }
    ++p2;
  }
  }
  //add the shortest path count with the cL_ label entries
  /*if(spd = UINT32_MAX) return 0;// there is not shortest path from u to v
  else
  {
    p1 = 0, p2 = 0;
    while(p1 < cL_[u].size() && p2 < cL_[v].size())
    {
        std::cout << "enter the count query with cL_, p1 and p2 together" << std::endl;
        const uint32_t w1 = LEExtractV(cL_[u][p1]);
        const uint32_t w2 = LEExtractV(cL_[v][p2]);
        if(rank_[w1] < rank_[w2]) ++p1;
        else if (rank_[w1] > rank_[w2] ) ++p2;
        else
        {
            const uint32_t d = LEExtractD(cL_[u][p1]) + LEExtractD(cL_[v][p2]);
            if(d == spd ) count += LEExtractC(cL_[u][p1]) * LEExtractC(cL_[v][p2]);
            ++p1, ++p2;
        }

    }
    while(p1 < cL_[u].size())
    {
        std::cout << "enter the count wiht p1" << std::endl;
        const uint32_t w1 = LEExtractV(cL_[u][p1]);
        const uint32_t d1 = LEExtractD(cL_[u][p1]);
        const uint32_t c1 = LEExtractC(cL_[u][p1]);
        if(w1 == v && d1 == spd) count += c1;
        ++p1;
    }
    while(p2 < cL_[v].size())
    {
        std::cout << "enter the count with p2" << std::endl;
        const uint32_t w1 = LEExtractV(cL_[v][p2]);
        const uint32_t d1 = LEExtractD(cL_[v][p2]);
        const uint32_t c1 = LEExtractC(cL_[v][p2]);
        if(w1 == u && d1 == spd) count += c1;
        ++p2;
    }
  }*/
  if(debug)std::cout << "the distance for" << u << " : " << v << " is " << spd << std::endl;
  return count;
}

/*
void USPC::IndexStat() const {
  std::vector<bool> reduced(n_, false);
  if (opt_shell_) {
    for (uint32_t v = 0; v < n_; ++v) {
      if (UINT32_MAX == shr_[v] || v != shr_[v]) {
        reduced[v] = true;
      }
    }
  }
  if (opt_equiv_) {
    for (uint32_t v = 0; v < n_; ++v) {
      if (UINT32_MAX == eqr_[v] || v != eqr_[v]) {
        reduced[v] = true;
      }
    }
  }
  uint32_t n = 0;
  for (uint32_t i = 0; i < n_; ++i) {
    if (reduced[i]) continue;
    ++n;
  }
  //
  std::vector<uint32_t> ticks;
  ticks.push_back(1);
  for (uint32_t i = 1; i <= 20; ++i) {
    ticks.push_back(ticks.back() * 2);
  }
  ASSERT(ticks.back() == (1 << 20));
  std::vector<uint32_t> cnt(ticks.size(), 0);
  // size of labels
  uint32_t max_size = 0;
  uint64_t size_L = 0, size_cL = 0, size_dL = 0;
  for (uint32_t i = 0; i < n_; ++i) {
    if (reduced[i]) continue;
    ASSERT(!dL_[i].empty());
    // update the max size
    uint32_t size = dL_[i].size() + cL_[i].size();
    if (size > max_size) max_size = size;
    // size distribution
    for (size_t i = 0; i < ticks.size(); ++i) {
      if (size <= ticks[i]) ++cnt[i];
    }
    // update total size
    size_dL += dL_[i].size();
    size_cL += cL_[i].size();
  }
  size_L = size_dL + size_cL;
  printf("L: %" PRIu64 ", dL: %" PRIu64 ", cL: %" PRIu64 "\n",
         size_L, size_dL, size_cL);
  printf("the max label size: %" PRIu32 " \n", max_size);
  for (size_t i = 0; i < ticks.size(); ++i) {
    printf("(%zu, %.4f)\n", i, static_cast<double>(cnt[i]) / n);
  }
}
*/

void USPCQuery::IndexRead(const std::string& filename) {
  ASSERT(dL_.empty() && cL_.empty() && G_.empty());
  FILE* file = fopen(filename.c_str(), "rb");
  ASSERT(fread(&n_, sizeof(n_), 1, file) == 1);
  G_.resize(n_);
  // read graph
  for (uint32_t u = 0; u < n_; ++u) {
    uint32_t s = 0;
    ASSERT(fread(&s, sizeof(s), 1, file) == 1);
    G_[u].resize(s);
    ASSERT(fread(G_[u].data(), sizeof(G_[u].back()), s, file) == s);
  }
  // initialization
  dL_.resize(n_);
  cL_.resize(n_);
  // read labels
  uint64_t num_labels = 0;
  uint64_t num_dlabels = 0;
  uint64_t num_clabels = 0;
  uint32_t d_s, c_s;
  for (uint32_t i = 0; i < n_; ++i) {
    // read distance labels
    ASSERT(fread(&d_s, sizeof(d_s), 1, file) == 1);
    dL_[i].resize(d_s);
    ASSERT(fread(dL_[i].data(), sizeof(dL_[i].back()), d_s, file) == d_s);
    num_labels += d_s;
    num_dlabels += d_s;
    // read counting labels
    ASSERT(fread(&c_s, sizeof(c_s), 1, file) == 1);
    cL_[i].resize(c_s);
    ASSERT(fread(cL_[i].data(), sizeof(cL_[i].back()), c_s, file) == c_s);
    num_labels += c_s;
    num_clabels += c_s;
  }
  printf("total # of label entries: %" PRIu64 "\n", num_labels);
  printf("total # of dlabels: %" PRIu64 "\n", num_dlabels);
  printf("total # of clabels: %" PRIu64 "\n", num_clabels);
  // order information
  order_.resize(n_);
  ASSERT(fread(order_.data(), sizeof(order_.back()), n_, file) == n_);
  rank_.resize(n_);
  OrderRank();
  // check
  for (uint32_t i = 0; i < n_; ++i) {
    for (size_t j = 1; j < dL_[i].size(); ++j) {
      ASSERT(rank_[LEExtractV(dL_[i][j])] >
             rank_[LEExtractV(dL_[i][j - 1])]);
    }
    for (size_t j = 1; j < cL_[i].size(); ++j) {
      ASSERT(rank_[LEExtractV(cL_[i][j])] >
             rank_[LEExtractV(cL_[i][j - 1])]);
    }
  }
  // multiplicity
  eqm_.resize(n_);
  ASSERT(fread(eqm_.data(), sizeof(eqm_.back()), n_, file) == n_);
  // optimization information
  ASSERT(fread(&opt_shell_, sizeof(opt_shell_), 1, file) == 1);
  if (opt_shell_) {
    std::cout << "shell reduction triggered" << std::endl;
    shr_.resize(n_);
    ASSERT(fread(shr_.data(), sizeof(shr_.back()), n_, file) == n_);
  }
  ASSERT(fread(&opt_equiv_, sizeof(opt_equiv_), 1, file) == 1);
  if (opt_equiv_) {
    std::cout << "equiv reduction triggered" << std::endl;
    eqr_.resize(n_); eqc_.resize(n_);
    ASSERT(fread(eqr_.data(), sizeof(eqr_.back()), n_, file) == n_);
    ASSERT(fread(eqc_.data(), sizeof(eqc_.back()), n_, file) == n_);
  }
  ASSERT(fread(&opt_local_, sizeof(opt_local_), 1, file) == 1);
  if (opt_local_) {
    std::cout << "ind-set reduction triggered" << std::endl;
    std::unique_ptr<bool[]> copy(new bool[n_]);
    ASSERT(fread(copy.get(), sizeof(bool), n_, file) == n_);
    local_.resize(n_);
    for (uint32_t i = 0; i < n_; ++i) local_[i] = copy[i];
    for (uint32_t u = 0; u < n_; ++u) {
      if (local_[u]) {
        ASSERT(dL_[u].empty());
        ASSERT(cL_[u].empty());
      }
    }
  }
  // reorganize the labels
  if (!opt_local_) {
    for (uint32_t i = 0; i < n_; ++i) {
      // merge
      std::vector<LabelEntry> mL;
      mL.reserve(dL_[i].size() + cL_[i].size());
      size_t di = 0, ci = 0;
      while (di < dL_[i].size() && ci < cL_[i].size()) {
        if (rank_[LEExtractV(dL_[i][di])] < rank_[LEExtractV(cL_[i][ci])]) {
          mL.push_back(dL_[i][di++]);
        } else {
          mL.push_back(cL_[i][ci++]);
        }
      }
      while (di < dL_[i].size()) mL.push_back(dL_[i][di++]);
      while (ci < cL_[i].size()) mL.push_back(cL_[i][ci++]);
      dL_[i].clear();
      cL_[i] = mL;
    }
    decltype(dL_)().swap(dL_);
    printf("labels merged.\n");
  } else {
    printf("labels not merged.\n");
  }

  // check
  bool check = false;
  ASSERT(fread(&check, sizeof(check), 1, file) == 0);
  fclose(file);

  dH_.resize(n_, UINT32_MAX);
  cH_.resize(n_, 0);

  std::cout << "index read" << std::endl;
}

void USPCIndex::IndexWrite(const std::string& filename) const {
  ASSERT(0 != n_);
  FILE* file = fopen(filename.c_str(), "wb");
  fwrite(&n_, sizeof(n_), 1, file);
  for (uint32_t u = 0; u < n_; ++u) {
    const uint32_t s = G_[u].size();
    fwrite(&s, sizeof(s), 1, file);
    fwrite(G_[u].data(), sizeof(G_[u].back()), s, file);
  }
  uint64_t num_labels = 0;
  for (uint32_t i = 0; i < n_; ++i) {
    // write distance labels
    const uint32_t d_s = dL_[i].size();
    fwrite(&d_s, sizeof(d_s), 1, file);
    fwrite(dL_[i].data(), sizeof(dL_[i].back()), d_s, file);
    num_labels += d_s;
    // write counting labels
    const uint32_t c_s = cL_[i].size();
    fwrite(&c_s, sizeof(c_s), 1, file);
    fwrite(cL_[i].data(), sizeof(cL_[i].back()), c_s, file);
    num_labels += c_s;
  }
  printf("total # of label entries: %" PRIu64 "\n", num_labels);
  // order information
  fwrite(order_.data(), sizeof(order_.back()), n_, file);
  // write multiplicity
  fwrite(eqm_.data(), sizeof(eqm_.back()), n_, file);
  // write optimization information
  fwrite(&opt_shell_, sizeof(opt_shell_), 1, file);
  if (opt_shell_) {
    fwrite(shr_.data(), sizeof(shr_.back()), n_, file);
  }
  fwrite(&opt_equiv_, sizeof(opt_equiv_), 1, file);
  if (opt_equiv_) {
    fwrite(eqr_.data(), sizeof(eqr_.back()), n_, file);
    fwrite(eqc_.data(), sizeof(eqc_.back()), n_, file);
  }
  fwrite(&opt_local_, sizeof(opt_local_), 1, file);
  if (opt_local_) {
    std::unique_ptr<bool[]> copy(new bool[n_]);
    for (uint32_t i = 0; i < n_; ++i) copy[i] = local_[i];
    fwrite(copy.get(), sizeof(bool), n_, file);
  }
  fclose(file);
}

void USPCIndex::ShellReduction() {
  std::cout << "shell reduction triggered" << std::endl;
  shr_ = std::vector<uint32_t>(n_);
  for (uint32_t u = 0; u < n_; ++u) {
    shr_[u] = (reduced_[u] ? UINT32_MAX : u);
  }
  // find the 1-shell
  std::vector<uint32_t> deg(n_, 0);
  std::vector<uint32_t> S, R;
  for (uint32_t u = 0; u < n_; ++u) {
    deg[u] = G_[u].size();
    if (1 == deg[u]) S.push_back(u);
  }
  uint32_t cnt_reduced = 0;
  while (!S.empty()) {
    const uint32_t u = S.back(); S.pop_back();
    if (0 == deg[u]) continue;
    reduced_[u] = true;
    R.push_back(u);
    ++cnt_reduced;
    for (const uint32_t v : G_[u]) {
      if (reduced_[v]) continue;
      shr_[u] = v;
      break;
    }
    if (1 == --deg[shr_[u]]) {
      S.push_back(shr_[u]);
    }
  }
  // set the representatives
  while (!R.empty()) {
    const uint32_t u = R.back(); R.pop_back();
    shr_[u] = shr_[shr_[u]];
  }
  // graph reduction
  for (uint32_t u = 0; u < n_; ++u) {
    if (reduced_[u]) G_[u].clear();
    else {
      size_t i = 0;
      for (size_t j = 0; j < G_[u].size(); ++j) {
        if (reduced_[G_[u][j]]) continue;
        G_[u][i++] = G_[u][j];
      }
      G_[u].resize(i);
      ASSERT(1 != G_[u].size());
    }
  }
  std::cout << "shell reduction: " << cnt_reduced << " "
            << "(" << static_cast<double>(cnt_reduced) / n_ << ")"
            << std::endl;
}

void USPCIndex::EquivReduction() {
  std::cout << "equiv reduction triggered" << std::endl;
  eqr_ = std::vector<uint32_t>(n_);
  eqc_ = std::vector<uint32_t>(n_, 0);
  for (uint32_t u = 0; u < n_; ++u) {
    eqr_[u] = (reduced_[u] ? UINT32_MAX : u);
  }
  // equiv-1
  std::vector<uint32_t> remain;
  for (uint32_t u = 0; u < n_; ++u) {
    if (reduced_[u]) continue;
    remain.push_back(u);
  }
  // sort the adjacency arrays
  std::sort(remain.begin(), remain.end(),
            [this](const uint32_t v1, const uint32_t v2) {
              if (G_[v1].size() < G_[v2].size()) return true;
              else if (G_[v1].size() > G_[v2].size()) return false;
              else {
                for (size_t i = 0; i < G_[v1].size(); ++i) {
                  if (G_[v1][i] < G_[v2][i]) return true;
                  else if (G_[v1][i] > G_[v2][i]) return false;
                }
                return v1 < v2;
              }
            });
  // identify the equivalence classes
  uint32_t cnt_reduced1 = 0;
  for (size_t i = 1; i < remain.size(); ++i) {
    if (G_[remain[i]] == G_[remain[i - 1]]) {
      eqr_[remain[i]] = eqr_[remain[i - 1]];
      eqc_[eqr_[remain[i]]] = G_[remain[i]].size();
      ++eqm_[eqr_[remain[i]]];
      ++cnt_reduced1;
    }
  }
  // equiv-2
  uint32_t cnt_reduced2 = 0;
  std::vector<bool> found(n_, false);
  for (uint32_t u = 0; u < n_; ++u) {
    if (reduced_[u] || found[u]) continue;
    for (const uint32_t v : G_[u]) {
      if (v < u || G_[u].size() != G_[v].size() || found[v]) continue;
      size_t pu = 0, pv = 0, cnt = 0;
      while (pu < G_[u].size() && pv < G_[v].size()) {
        if (G_[u][pu] == v) ++pu;
        else if (G_[v][pv] == u) ++pv;
        else if (G_[u][pu] == G_[v][pv]) ++pu, ++pv, ++cnt;
        else break;
      }
      if (cnt == G_[u].size() - 1) {
        found[v] = true;
        eqr_[v] = u;
        eqc_[u] = 1;
        ++eqm_[u];
        ++cnt_reduced2;
      }
    }
  }
  // update reduced_
  uint32_t cnt_reduced = 0;
  for (uint32_t u = 0; u < n_; ++u) {
    if (reduced_[u]) continue;
    if (eqr_[u] != u) {
      reduced_[u] = true;
      ++cnt_reduced;
    }
  }
  ASSERT(cnt_reduced1 + cnt_reduced2 == cnt_reduced);
  // graph reduction
  for (uint32_t u = 0; u < n_; ++u) {
    if (reduced_[u]) G_[u].clear();
    else {
      size_t j = 0;
      for (size_t i = 0; i < G_[u].size(); ++i) {
        if (reduced_[G_[u][i]]) continue;
        G_[u][j++] = G_[u][i];
      }
      G_[u].resize(j);
    }
  }
  std::cout << "equiv reduction #1: " << cnt_reduced1 << " "
            << "(" << static_cast<double>(cnt_reduced1) / n_ << ")\n"
            << "equiv reduction #2: " << cnt_reduced2 << " "
            << "(" << static_cast<double>(cnt_reduced2) / n_ << ")\n"
            << "equiv reduction total: " << cnt_reduced << " "
            << "(" << static_cast<double>(cnt_reduced) / n_ << ")"
            << std::endl;
}

uint32_t USPCIndex::Distance(const std::vector<uint32_t>& dLu,
                             const std::vector<LabelEntry>& dLv) const {
  uint32_t d = UINT32_MAX;
  for (const auto e : dLv) {
    const uint32_t v = LEExtractV(e);
    if (UINT32_MAX == dLu[v]) continue;
    const uint32_t dd = dLu[v] + LEExtractD(e);
    if (dd < d) d = dd;
  }
  return d;
}

int USPCIndex::PrunedDistance(const std::vector<uint32_t>& dLu,
                             const std::vector<LabelEntry>& dLv, uint32_t currentD, uint32_t s, uint32_t t) const {
  int result = -1;//-1 means distance in index > currentD, 0 means equals, and 1 means <
  uint32_t d = UINT32_MAX;
  for (const auto e : dLv)//auto e = dLv.rbegin(); e != dLv.rend(); e++)//const auto e : dLv) {
  { 
    const uint32_t v = LEExtractV(e);
    //if(rank_[v] > rank_[s] ) break;
    if (UINT32_MAX == dLu[v]) continue;
    const uint32_t dd = dLu[v] + LEExtractD(e);
    if (dd < currentD) return 1;
    else if(dd == currentD) result = 0;
    if(v == s)
    {
        if(LEExtractD(e) < currentD) return 1;
        else if(LEExtractD(e) == currentD) result = 0;
    }
  }
  if(dLu[t] < currentD) return 1; 
  if(dLu[t] == currentD) result = 0; 
  return result;
}
uint32_t USPCIndex::Distance(const std::vector<LabelEntry>& dLu,
                             const std::vector<LabelEntry>& dLv, uint32_t u, uint32_t v) const {
  uint32_t spd = UINT32_MAX;
  if(u == v) return 0;
  size_t p1 = 0, p2 = 0;
  if( (dLu.size() == 0 &&  dLv.size() > 0) ||  (dLu.size() > 0 &&  dLv.size() == 0) )
  {
      p1 = 0;
      p2 = 0;
      while(p1 < dLu.size())
      {
        const uint32_t w1 = LEExtractV(dLu[p1]);
        const uint32_t d1 = LEExtractD(dLu[p1]);
        if(( w1 == v || w1 == u) && d1 < spd)
        {
            spd = d1;
        }
        ++p1;
      }
      while(p2 < dLv.size())
      {
        const uint32_t w2 = LEExtractV(dLv[p2]);
        const uint32_t d2 = LEExtractD(dLv[p2]);
        if((w2 == u || w2 == v) && d2 < spd)
        {
            spd = d2;
        }
        ++p2;
      }
  }
  else
  {
  while(p1 < dLu.size() && p2 < dLv.size())
  {
    const uint32_t w1 = LEExtractV(dLu[p1]);
    const uint32_t w2 = LEExtractV(dLv[p2]);
    const uint32_t d1 = LEExtractD(dLu[p1]);
    const uint32_t d2 = LEExtractD(dLv[p2]);
    if((w1 == v || w1 == u) && d1 < spd)
    {
        spd = d1;
    }
    if( (w2 == u|| w2 == v) && d2 < spd)
    {
        spd = d2;
    }
    if(rank_[w1] < rank_[w2]) ++p1;
    else if(rank_[w1] > rank_[w2]) ++p2;
    else
    {
        const uint32_t d = LEExtractD(dLu[p1]) + LEExtractD(dLv[p2]);
        if(d < spd) spd = d;
        ++p1; ++p2;
    }
  }
  while(p1 < dLu.size())
  {
    //std::cout << "enter dLu" << dLu.size() << " : " << p1 << std::endl;
    const uint32_t w1 = LEExtractV(dLu[p1]);
    const uint32_t d1 = LEExtractD(dLu[p1]);
    if(( w1 == v || w1 == u) && d1 < spd)
    {
        spd = d1;
    }
    ++p1;
  }
  while(p2 < dLv.size())
  {
    //std::cout << "enter dLv" << dLv.size() << " : " << p2 << std::endl;
    const uint32_t w2 = LEExtractV(dLv[p2]);
    const uint32_t d2 = LEExtractD(dLv[p2]);
    if((w2 == u || w2 == v) && d2 < spd)
    {
        spd = d2;
    }
    ++p2;
  }
  }
  return spd;
}
void USPCIndex::DegreeOrder(const Graph& graph) {
  std::cout << "degree-based ordering" << std::endl;
  std::vector<uint32_t> deg(n_);
  for (uint32_t i = 0; i < n_; ++i) {
    deg[i] = graph[i].size();
    order_[i] = i;
  }
  std::sort(order_.begin(), order_.end(),
            [&deg](const uint32_t v1, const uint32_t v2) {
              return deg[v1] > deg[v2];
            });
}

/*
void USPCIndex::DegenOrder(const Graph& graph) {
  ASSERT(false);
  // initialization
  uint32_t max_deg = 0;
  std::vector<uint32_t> deg(n_);
  for (uint32_t i = 0; i < n_; ++i) {
    deg[i] = graph[i].size();
    if (deg[i] > max_deg) max_deg = deg[i];
  }
  // sort the vertices
  std::vector<uint32_t> p(max_deg + 1, 0);
  for (uint32_t i = 0; i < n_; ++i) ++p[deg[i]];
  for (uint32_t s = 0, d = 0; d <= max_deg; ++d) {
    const uint32_t c = p[d];
    p[d] = s;
    s += c;
  }
  std::vector<uint32_t> pos(n_);
  for (uint32_t i = 0; i < n_; ++i) {
    order_[p[deg[i]]] = i;
    pos[i] = p[deg[i]]++;
  }
  for (uint32_t d = max_deg; d > 0; --d) p[d] = p[d-1];
  p[0] = 0;
  // core decomposition
  uint32_t max_k = 0;
  std::vector<bool> rm(n_, false);
  for (uint32_t i = 0; i < n_; ++i) {
    const uint32_t v = order_[i];
    if (deg[v] > max_k) max_k = deg[v];
    rm[v] = true;
    ++p[deg[v]];
    for (const uint32_t u : graph[v]) {
      if (rm[u]) continue;
      const uint32_t pu = pos[u];
      const uint32_t pw = p[deg[u]];
      if (pu != pw) {
        const uint32_t w = order_[pw];
        std::swap(order_[pu], order_[pw]);
        pos[u] = pw;  pos[w] = pu;
      }
      ++p[deg[u]];
      --deg[u];
      if (i + 1 == pos[u] || deg[order_[pos[u] - 1]] != deg[u]) {
        p[deg[u]] = pos[u];
      }
    }
  }
  // reverse
  std::reverse(order_.begin(), order_.end());
}

void USPCIndex::RevDegenOrder(const Graph& graph) {
  ASSERT(false);
  // initialization
  uint32_t max_deg = 0;
  std::vector<uint32_t> deg(n_);
  for (uint32_t i = 0; i < n_; ++i) {
    deg[i] = graph[i].size();
    if (deg[i] > max_deg) max_deg = deg[i];
  }
  // sort the vertices
  std::vector<uint32_t> p(max_deg + 1, 0);
  for (uint32_t i = 0; i < n_; ++i) ++p[deg[i]];
  uint32_t s = 0;
  for (uint32_t d = max_deg; d > 0; --d) {
    const uint32_t c = p[d];
    p[d] = s;
    s += c;
  }
  p[0] = s;
  std::vector<uint32_t> pos(n_);
  for (uint32_t i = 0; i < n_; ++i) {
    order_[p[deg[i]]] = i;
    pos[i] = p[deg[i]]++;
  }
  for (uint32_t d = 0; d <= max_deg; ++d) --p[d];
  // ordering
  std::vector<bool> rm(n_, false);
  for (uint32_t i = 0; i < n_; ++i) {
    const uint32_t v = order_[i];
    rm[v] = true;
    for (const uint32_t u : graph[v]) {
      if (rm[u]) continue;
      const uint32_t pu = pos[u];
      const uint32_t pw = p[deg[u]];
      if (pu != pw) {
        const uint32_t w = order_[pw];
        std::swap(order_[pu], order_[pw]);
        pos[u] = pw;  pos[w] = pu;
      }
      --p[deg[u]];
      --deg[u];
      if (n_ - 1 == pos[u] || deg[order_[pos[u] + 1]] != deg[u]) {
        p[deg[u]] = pos[u];
      }
    }
  }
}

void USPCIndex::GreedyOrder(const Graph& graph) {
  ASSERT(false);
  // degree
  std::vector<uint32_t> d(n_, 0);
  for (uint32_t u = 0; u < n_; ++u) d[u] = graph[u].size();
  // priority queue and the current # of open wedges covered
  std::priority_queue<std::pair<uint64_t, uint32_t>> Q;
  std::vector<uint64_t> ow(n_, 0);
  // triangle listing
  Graph g = graph;
  std::vector<uint64_t> t(n_, 0);
  { // forward algorithm
    for (uint32_t u = 0; u < n_; ++u) {
      std::sort(g[u].begin(), g[u].end());
      order_[u] = u;
    }
    std::sort(order_.begin(), order_.end(),
              [&d](const uint32_t v1, const uint32_t v2) {
                return d[v1] > d[v2];
              });
    std::vector<uint32_t> eta(n_, 0);
    for (uint32_t i = 0; i < n_; ++i) eta[order_[i]] = i;
    std::vector<std::vector<uint32_t>> A(n_);
    for (const uint32_t v : order_) {
      for (const uint32_t u : g[v]) {
        if (eta[u] < eta[v]) continue;
        uint32_t pu = 0, pv = 0;
        while (pu < A[u].size() && pv < A[v].size()) {
          if (eta[A[u][pu]] < eta[A[v][pv]]) ++pu;
          else if (eta[A[u][pu]] > eta[A[v][pv]]) ++pv;
          else {
            ++t[A[u][pu]]; ++t[u]; ++t[v];
            ++pu; ++pv;
          }
        }
        A[u].push_back(v);
      }
    }
  }
  // count the # of open wedges
  for (uint32_t u = 0; u < n_; ++u) {
    ow[u] = static_cast<uint64_t>(d[u]) * (d[u] - 1) / 2;
    for (const uint32_t v : g[u]) {
      ow[u] += d[v] - 1;
    }
    ow[u] -= 3 * t[u];
    Q.push({ow[u], u});
  }
  // has been removed?
  std::vector<bool> rmv(n_, false);
  // greedy algorithm
  for (uint32_t i = 0; i < n_; ++i) {
    ASSERT(!Q.empty());
    const auto t = Q.top(); Q.pop();
    const uint64_t u_ow = t.first;
    const uint32_t u = t.second;
    if (u_ow > ow[u]) {
      --i;
      continue;
    }
    for (const uint32_t v : g[u]) {
      if (rmv[v]) continue;
      --d[v];
      // compute the # of open wedges containing both u and v
      std::unordered_set<uint32_t> triangle;
      uint64_t sub = d[v] + d[u] - 1;
      // find triangles containing both u and v
      uint32_t pu = 0, pv = 0;
      while (pu < g[u].size() && pv < g[v].size()) {
        if (g[u][pu] > g[v][pv]) ++pv;
        else if (g[u][pu] < g[v][pv]) ++pu;
        else {
          if (!rmv[g[u][pu]]) {
            sub -= 2;
            triangle.insert(g[u][pu]);
          }
          ++pu; ++pv;
        }
      }
      if (sub > 0) {
        ow[v] -= sub;
        Q.push({ow[v], v});
      }
      // update ow for vertices of distance 2
      for (const uint32_t w : g[v]) {
        if (w != u && !rmv[w] && triangle.count(w) == 0) {
          --ow[w];
          Q.push({ow[w], w});
        }
      }
    }
    order_[i] = u;
    rmv[u] = true;
  }
  ASSERT(Q.empty());
}
*/
}  // namespace spc
