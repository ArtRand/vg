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

        // Aligns the sequence in aln to the graph, if path_alignments != nullptr
        // then it adds alignments from all paths, otherwise it just adds the 
        // highest scoring path (TODO add get_pairs into AlignmentParameters?)
        void Align(Alignment& aln, std::vector<Alignment>* path_alignments, 
                   AlignmentParameters& p, bool get_pairs=true);  // TODO not implemented

        void Align(Alignment&, AlignmentParameters& p, bool get_pairs=true);

        void PrintAlignedPairs();

    private:
        // Containers
        //
        // maps the nodeId from VG::Graph to VertexId in HmmGraph and vice-versa
        std::unordered_map<int64_t, int64_t> nodeId_to_vertexId;  // (Vg_node_id, hmm_graph_node_id)
        std::unordered_map<int64_t, int64_t> vertexId_to_nodeId;  // (hmm_graph_node_id, Vg_node_id)
        // interface point to Parcours
        HmmGraph hmm_graph;

        // counters and flags
        //
        
        // Internal Methods
        //
        // converts the alignement in GraphAlignedPairs to a VG alignment object, pId is the 
        // path to use (pathId)
        void makeAlignmentFromAlignedPairs(Alignment&, int64_t pId);
        std::unordered_map<int64_t, AlignedPairs> mapAlignedPairsToVgNodes(int64_t pId);
        std::vector<double> pathScoresVector();
    };
};

#endif
