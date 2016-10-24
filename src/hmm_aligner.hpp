// 
// hmm_aligner.hpp
//
// Methods to interface with aligning to a graph with a pair-HMM
//

#ifndef HMM_ALIGNER_HPP
#define HMM_ALIGNER_HPP

#include <unordered_map>
#include <vector>

#include "vg.pb.h"
#include "utility.hpp"             // nonACGTNtoN
#include "Parcours/hmm_graph.h"
#include "Parcours/common.h"       // st_uglyf 

namespace vg {
   class HmmAligner {
    public:
        HmmAligner(Graph& G);
        
        ~HmmAligner();
        
        HmmGraph hmm_graph;
    private:
        std::unordered_map<int64_t, int64_t> id_to_index;
        
    };
};

#endif
