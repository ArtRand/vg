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

namespace vg {
   class HmmAligner {
    public:
        HmmAligner(Graph& G);
        
        ~HmmAligner();
        
        int64_t Graph_K();
        
    private:
        std::unordered_map<int64_t, int64_t> id_to_index;
        // maps the nodeId from VG::Graph to VertexId in HmmGraph
        std::unordered_map<int64_t, int64_t> nodeId_to_vertexId;
        HmmGraph hmm_graph;
    };
};

#endif
