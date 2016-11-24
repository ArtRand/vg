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
        // COMMENT ME!
        HmmAligner(Graph& G, bool reversed=false);
        
        ~HmmAligner();
        
        int64_t Graph_K();

        void Align(Alignment& aln, std::vector<Alignment>* path_alignments, 
                   AlignmentParameters& p, bool get_pairs=true);

    private:
        // maps the nodeId from VG::Graph to VertexId in HmmGraph
        std::unordered_map<int64_t, int64_t> nodeId_to_vertexId;
        // TODO might need to keep track of nodes from the VG graph..?
        std::unordered_map<int64_t, Node*> node_map;

        HmmGraph hmm_graph;
        bool _reversed;

        /*
         * Internal Methods
         */
        void reverse_hmm_graph();
        void graph_aligned_pairs_to_alignment(Alignment&);
    };
};

#endif
