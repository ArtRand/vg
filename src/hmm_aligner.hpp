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
#include "path.hpp"

#include "Parcours/hmm_graph.h"

#define ALIGNED_PAIR_LENGTH 1

namespace vg {
   class HmmAligner {
    public:
        // COMMENT ME!
        HmmAligner(Graph& G);
        
        ~HmmAligner();
        
        int64_t Graph_K();

        // Aligns the sequence in aln to the graph, if path_alignments != nullptr
        // then it adds alignments from all paths, otherwise it just adds the 
        // highest scoring path 
        void Align(Alignment& aln, std::vector<Alignment>* path_alignments, AlignmentParameters& p, 
                   bool get_pairs=true, bool ragged_end=false);

        //void Align(Alignment&, AlignmentParameters& p, bool get_pairs=true); // TODO not implemented

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
        // translates GraphAlignedPairs to an alignment in the VG aignment object, 
        // pId is the path to use (pathId)
        void makeAlignmentFromPathAlignedPairs(Alignment&, int64_t pId, bool ragged_end);
        std::unordered_map<int64_t, AlignedPairs> mapAlignedPairsToVgNodes(int64_t pId);
        std::vector<double> pathScoresVector();
    };
};

#endif
