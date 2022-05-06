#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <string>

#include "macros.h"
#include "u_io.h"
#include "u_label.h"
#include "u_spc.h"

bool ReadBool(const char tmp) {
  ASSERT('y' == tmp || 'n' == tmp);
  return 'y' == tmp;
}

int main(int argc, char** argv) {
  // initialize the options
  std::string gfilename;
  std::string lfilename;
  std::string osname;
  spc::USPCIndex::OrderScheme os = spc::USPCIndex::OrderScheme::kInvalid;
  bool opt_shell = false;
  bool opt_equiv = false;
  bool opt_local = false;
  int numThreads = 0;
  printf("test \n");
  {
    int option = -1;
    while (-1 != (option = getopt(argc, argv, "g:l:o:s:e:i:t:"))) {
      switch (option) {
        case 'g':
          gfilename = optarg; break;
        case 'l':
          lfilename = optarg; break;
        case 'o':
          osname = optarg;
          if ("degree" == osname) {
            os = spc::USPCIndex::OrderScheme::kDegree;
          } else if ("sigpath" == osname) {
            os = spc::USPCIndex::OrderScheme::kSigPath;
          } else {
            os = spc::USPCIndex::OrderScheme::kInvalid;
          }
          break;
        case 's':
          opt_shell = ReadBool(optarg[0]);
          break;
        case 'e':
          opt_equiv = ReadBool(optarg[0]);
          break;
        case 'i':
          opt_local = ReadBool(optarg[0]);
          break;
        case 't':
          numThreads = atoi(optarg);
          //omp_set_num_threads(1);
          omp_set_num_threads(numThreads);
          break;
      }
    }
    printf("graph file: %s\n", gfilename.c_str());
    printf("label file: %s\n", lfilename.c_str());
    printf("ordering: %s\n", osname.c_str());
    printf("numer of threads is %d\n", numThreads);
    printf("opt_shell: %c\n", opt_shell ? 'y' : 'n');
    printf("opt_equiv: %c\n", opt_equiv ? 'y' : 'n');
    printf("opt_local: %c\n", opt_local ? 'y' : 'n');
  }
  // read the graph
  uint32_t n, m;
  spc::Graph graph;
  GraphRead(gfilename, graph, n, m);
  // timer starts
  const auto beg = std::chrono::steady_clock::now();
  // build index
  spc::USPCIndex spc;
  spc.set_os(os);
  spc.set_opt_shell(opt_shell);
  spc.set_opt_equiv(opt_equiv);
  spc.set_opt_local(opt_local);
  if(numThreads <= 0) numThreads = 1;
  printf("%d threads in total \n");
  //spc.BuildIndexParallel(graph,numThreads);
  spc.BuildIndexParallel_vector(graph,numThreads);
  //spc.BuildIndex(graph);
  // timer stops
  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  printf("index construction costs %f ms\n",
         std::chrono::duration<double, std::milli>(dif).count());
  // write the results
  spc.IndexWrite(lfilename);
}
