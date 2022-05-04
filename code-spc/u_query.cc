#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include <fstream>
//#include <omp.h>
#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

int main(int argc, char** argv) {
  std::string lfilename;
  std::string qfilename;
  std::string afilename;
  std::string gfilename;
  int option = -1;
  while (-1 != (option = getopt(argc, argv, "l:a:q:g:"))) {
    switch (option) {
      case 'l':
        lfilename = optarg; break;
      case 'q':
        qfilename = optarg; break;
      case 'a':
        afilename = optarg; break;
      case 'g':
        gfilename = optarg; break;
    }
  }
  printf("label file name: %s\n", lfilename.c_str());
  printf("query file anme: %s\n", qfilename.c_str());
  printf("answer file name: %s\n", afilename.c_str());
  // read index
  spc::USPCQuery uspc;
  uspc.IndexRead(lfilename);
  //valid with bfs
  std::vector<std::vector<int>> adj;
  adj.resize(uspc.get_n()+1);
  std::ifstream in;
  in.open("graph/"+gfilename);
  int u,v;
  int m,n;
  int lu = -1, lv = -1;
  in >> m;
  in >> n;
  while(! in.eof())
  {
     in >> u;
     in >> v;
     if(lu == u && lv == v) continue;
     ///std::cout << u << " : " << v << std::endl;
     adj[u].push_back(v);
     adj[v].push_back(u);
     lu = u;
     lv = v;
  }
  in.close();
  uspc.ValidWithBfsQueries(10000, adj, false);


  // read queries
  FILE* file = fopen(qfilename.c_str(), "r");
  uint32_t num_queries = 0;
  fscanf(file, "%" SCNu32, &num_queries);
  std::vector<std::pair<uint32_t, uint32_t>> queries;
  for (uint32_t q = 0; q < num_queries; ++q) {
    uint32_t v1, v2;
    fscanf(file, "%" SCNu32 " %" SCNu32, &v1, &v2);
    queries.push_back({v1, v2});
  }
  fclose(file);
  // compute the results
  std::vector<uint64_t> results;
  const auto beg = std::chrono::steady_clock::now();
  for (const auto query : queries) {
    const uint32_t v1 = query.first;
    const uint32_t v2 = query.second;
    uint64_t result;
    result = uspc.Count(v1, v2);
    results.push_back(result);
  }
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  ASSERT(results.size() == num_queries);
  printf("query costs %f micro seconds in average\n",
         std::chrono::duration<double, std::micro>(dif).count() / num_queries);
  // read answers
  FILE* afile = fopen(afilename.c_str(), "r");
  uint32_t qn = 0;
  fscanf(afile, "%" SCNu32, &qn);
  ASSERT(10000 == qn);
  std::vector<uint64_t> answers;
  for (uint32_t q = 0; q < qn; ++q) {
    uint64_t cnt = 0;
    fscanf(afile, "%" SCNu64, &cnt);
    answers.push_back(cnt);
  }
  fclose(afile);
  results.resize(qn);
  for (uint32_t q = 0; q < qn; ++q) {
    if (results[q] != answers[q]) {
      std::cout << results[q] << " != " << answers[q] << std::endl;
    }
  }
  ASSERT(results == answers);
  printf("passed.\n");
}
