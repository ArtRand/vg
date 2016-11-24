#include <iostream>
#include "hmm_aligner.hpp"

using namespace vg;

HmmAligner::HmmAligner(Graph& G, bool reversed) {
    st_uglyf("[HmmAligner::HmmAligner]: Building HmmGraph from Graph object, graph has %lld nodes "
             "and %lld edges\n", G.node_size(), G.edge_size());
    st_uglyf("[HmmAligner::HmmAligner]: internal graph has size %lld\n", hmm_graph.K());
    
    if (reversed) {
        // TODO
    }

    for (int64_t i = 0; i < G.node_size(); i++) {
        Node *n = G.mutable_node(i);

        st_uglyf("[HmmAligner::HmmAligner]: node %lld has sequence %s\n", n->id(), n->sequence().c_str());;
        
        // add the vertex to the HmmGraph and add the node id (from the orig. graph) tothe hash 
        int64_t vId = hmm_graph.AddVertex(&n->sequence());
        nodeId_to_vertexId[n->id()] = vId;
    }
    
    if (hmm_graph.K() != G.node_size()) throw ParcoursException(
            "[HmmAligner::HmmAligner]: VG graph and internal graph have "
             "different K() (number of vertices/nodes)");

    st_uglyf("[HmmAligner::HmmAligner]: internal graph has %lld vertices\n", hmm_graph.K());
    
    // now go through the edges
    for (int64_t i = 0; i < G.edge_size(); i++) {
        Edge *e = G.mutable_edge(i);
        
        // Parcours uses the term "vertex" so when a variable is named "*vertex*" it should
        // reflect that it belongs to a Parcours object
        int64_t fromNodeId = e->from();
        int64_t toNodeId   = e->to();
        int64_t fromVertex = nodeId_to_vertexId[fromNodeId];
        int64_t toVertex   = nodeId_to_vertexId[toNodeId];

        if (!e->from_start() && !e->to_end()) {  // normal edge
            st_uglyf("[HmmAligner::HmmAligner]: Adding end->start edge from %lld to %lld ", fromNodeId, toNodeId);
            st_uglyf("Adding arc from %lld to %lld \n", fromVertex, toVertex);
            hmm_graph.AddArc(fromVertex, toVertex);

        } else if (e->from_start() && e->to_end()) {
            // TODO needs testing
            // this is a 'normal' edge but it's orientation is reversed
            st_uglyf("[HmmAligner::HmmAligner]: Adding start->end edge from %lld to %lld ", fromNodeId, toNodeId);
            st_uglyf("Adding arc from %lld to %lld\n", fromVertex, toVertex);
            hmm_graph.AddArc(toVertex, fromVertex);
        } else {  // TODO test
            throw ParcoursException("[HmmAligner::HmmAligner] Reversing edges not supported, yet\n");
        }
    }
    _reversed = reversed;
}

HmmAligner::~HmmAligner() {};

int64_t HmmAligner::Graph_K() { return hmm_graph.K(); }

void HmmAligner::Align(Alignment& aln, std::vector<Alignment>* path_alignments, 
                       AlignmentParameters& p, bool get_pairs) {
    // step 1: align the sequence to the graph with the HMM
    std::string* sx = aln.mutable_sequence();
    if (sx->length() <= 0) return;
    hmm_graph.AlignWithFiveStateSymbolHmm((*sx), p, get_pairs);
    for (auto kv : hmm_graph.PathAlignedPairs()) {
        st_uglyf("path %lld score: %f\n", kv.first, hmm_graph.PathScores()[kv.first]);
        for (auto p : kv.second) {
            st_uglyf("x: %lld vertex: %lld offset %lld p: %f\n", 
                      std::get<1>(p), std::get<2>(p).first, std::get<2>(p).second, std::get<0>(p));
        }
    }
    // step 2: if not doing multi-path alignment, just add the one with the max score
    if (!path_alignments) {
        graph_aligned_pairs_to_alignment(aln);
    } else {
        throw ParcoursException("HmmAligner::Align multipath not implemented");  
    }
}

void HmmAligner::graph_aligned_pairs_to_alignment(Alignment& aln) {
    // setup alignment (clear, establish score, etc)
    //aln.clear_path();
    //aln.set_score()
    // Go thought 
}
