#include "u_io.h"
#include <iostream>
#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <set>

#include "macros.h"

void GraphRead(const std::string& filename,
               spc::Graph& graph,
               uint32_t& n, uint32_t& m) {
  ASSERT(graph.empty());
  FILE* gfile = fopen(filename.c_str(), "r");
  fscanf(gfile, "%" SCNu32 " %" SCNu32, &n, &m);
  // check the # of vertices
  spc::NormalV(n);
  // construct graph
  graph.resize(n);
  std::set<std::pair<uint32_t, uint32_t>> edges;
  for (uint32_t e = 0; e < m && (!feof(gfile)); ++e) {
    uint32_t v1, v2;
    fscanf(gfile, "%" SCNu32 " %" SCNu32, &v1, &v2);
    if(feof(gfile)) break;
    ASSERT(v1 < n && v2 < n);
    graph[v1].push_back(v2);
    graph[v2].push_back(v1);
    //std::cout << v1 << " : " << v2 << std::endl; 
    // check
    //ASSERT(v1 != v2);
    if(v1 == v2) 
    {
        e--;
        continue;
    }
    if (v1 > v2) std::swap(v1, v2);
    if(0 != edges.count({v1, v2})) continue;
    edges.insert({v1, v2});
  }
  for (uint32_t v = 0; v < n; ++v) {
    std::sort(graph[v].begin(), graph[v].end());
  }
  fclose(gfile);
}
