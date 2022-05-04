#ifndef SPC_U_SPC_H_
#define SPC_U_SPC_H_

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <omp.h>
#include "macros.h"
#include "u_label.h"
#include <queue>

namespace spc {
class USPC {
 public:
  // ctors
  USPC() = default;
  USPC(const USPC&) = delete;
  USPC& operator=(const USPC&) = delete;

 protected:
  void OrderRank() {
    for (uint32_t i = 0; i < n_; ++i) {
      rank_[order_[i]] = i;
    }
  }
  // core data
  uint32_t n_;
  Graph G_;
  Label dL_, cL_;
  std::vector<uint32_t> order_;
  std::vector<uint32_t> rank_;
  // index reduction
  std::vector<uint32_t> shr_;
  std::vector<uint32_t> eqr_;
  std::vector<uint32_t> eqc_;
  std::vector<uint32_t> eqm_;
  std::vector<bool> reduced_;
  std::vector<bool> local_;
  // control parameters
  bool opt_shell_ = false;
  bool opt_equiv_ = false;
  bool opt_local_ = false;
};

class USPCIndex final: private USPC {
 public:
  enum class OrderScheme {
    kDegree,   // degree-based ordering
    kSigPath,  // significant-path-based ordering
    kInvalid
  };
  // ctors
  USPCIndex() = default;
  USPCIndex(const USPCIndex&) = delete;
  USPCIndex& operator=(const USPCIndex&) = delete;
  // construct an index
  void BuildIndex(const Graph& const_graph);
  void BuildIndexParallel(const Graph& const_graph, int num_threads);
  void BuildIndexParallel_vector(const Graph& const_graph, int num_threads);
  // write the index to disk
  void IndexWrite(const std::string& filename) const;
  // mutators
  void set_opt_shell(const bool opt_shell) { opt_shell_ = opt_shell; }
  void set_opt_equiv(const bool opt_equiv) { opt_equiv_ = opt_equiv; }
  void set_opt_local(const bool opt_local) { opt_local_ = opt_local; }
  void set_os(const OrderScheme os) { os_ = os; }

 private:
  void ShellReduction();
  void EquivReduction();
  uint32_t Distance(const std::vector<uint32_t>& dLu,
                    const std::vector<LabelEntry>& dLv) const;
  uint32_t Distance(const std::vector<LabelEntry>& dLu,
                           const std::vector<LabelEntry>& dLv, uint32_t u, uint32_t v) const;
  int PrunedDistance(const std::vector<uint32_t>& dLu, 
          const std::vector<LabelEntry>& dLv, uint32_t currentD, uint32_t s, uint32_t t) const;
 
  void SigPathIndex();
  // ordering functions
  void DegreeOrder(const Graph& graph);
  // void DegenOrder(const Graph& graph);
  // void RevDegenOrder(const Graph& graph);
  // void GreedyOrder(const Graph& graph);
  void InvalidOrder(const Graph&) {
    ASSERT_INFO(false, "invalid ordering");
  }
  std::map<OrderScheme, void (USPCIndex::*)(const Graph&)> of_ = {
    {OrderScheme::kDegree,  &USPCIndex::DegreeOrder},
    {OrderScheme::kInvalid, &USPCIndex::InvalidOrder}
  };
  OrderScheme os_ = OrderScheme::kInvalid;
};

class USPCQuery final: private USPC {
 public:
  USPCQuery() = default;
  USPCQuery(const USPCQuery&) = delete;
  USPCQuery& operator=(const USPCQuery&) = delete;
  int get_n()
  { return n_;}

// query
  uint64_t Count(uint32_t v1, uint32_t v2) const;
  uint32_t Count_New(uint32_t u, uint32_t v) const;
  // read the index from disk
  void ValidWithBfsQueries(int nq, std::vector<std::vector<int>>& adj, bool newCount);//nq is the number of queries
  std::pair<int, uint32_t> SpcBfs(std::vector<std::vector<int>>& adj, uint32_t s, uint32_t t);
  void IndexRead(const std::string& filename);

 private:
  // used for IS-Counting
  mutable std::vector<uint32_t> dH_;
  mutable std::vector<uint64_t> cH_;
};
}

#endif
