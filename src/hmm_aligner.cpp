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
        vertexId_to_nodeId[vId] = n->id();
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
}

HmmAligner::~HmmAligner() {};

int64_t HmmAligner::Graph_K() { return hmm_graph.K(); }

void HmmAligner::Align(Alignment& aln, std::vector<Alignment>* path_alignments, 
                       AlignmentParameters& p, bool get_pairs) {
    // step 1: align the sequence to the graph with the HMM
    std::string* sx = aln.mutable_sequence();
    if (sx->length() <= 0) return;
    // TODO need try/catch here:
    hmm_graph.AlignWithFiveStateSymbolHmm((*sx), p, get_pairs);
    PrintAlignedPairs();
    // step 2: if not doing multi-path alignment, just add the one with the max score
    if (!path_alignments) {
        auto max_pathId = hmm_graph.MaxScorePath();
        makeAlignmentFromAlignedPairs(aln, max_pathId);
    } else {
        throw ParcoursException("HmmAligner::Align multipath not implemented");  
    }
}

void HmmAligner::PrintAlignedPairs() {
    auto& graph_aligned_pairs = hmm_graph.PathAlignedPairs();
    if (graph_aligned_pairs.size() == 0) std::cerr << "[HmmAligner::PrintAlignedPairs] No GraphAlignedPairs\n";
    auto graph_path_scores   = hmm_graph.PathScores();
    int64_t most_prob_path   = hmm_graph.MaxScorePath();

    for (auto kv : graph_aligned_pairs) {
        std::cerr << "path: " << kv.first << " score: " << graph_path_scores[kv.first];
        if (kv.first == most_prob_path) std::cerr << " <-- Highest Score Path\n";
        else std::cerr << "\n";
        for (auto p : kv.second) {
            int64_t x_coord = std::get<1>(p);
            int64_t vid     = std::get<2>(p).first;
            int64_t offset  = std::get<2>(p).second;
            double prob     = std::get<0>(p)/PAIR_ALIGNMENT_PROB_1;
            std::cerr << "x: " << x_coord << " vertexID: " << vid << 
                " offset: " << offset << " prob: " << prob << std::endl;
        }
    }
}

void HmmAligner::makeAlignmentFromAlignedPairs(Alignment& aln, int64_t pId) {
    // setup alignment (clear, establish score, etc)
    aln.clear_path();
    aln.set_score(hmm_graph.PathScores()[pId] * PAIR_ALIGNMENT_PROB_1);
    // an Alignment has one Path, so add that
    Path *aln_path = aln.mutable_path();
    auto mapped_pairs = mapAlignedPairsToVgNodes(pId);
    auto path_deque = hmm_graph.PathMap()[pId];  // contains the vertex IDs (in order) for this path
    st_uglyf("SENTINAL vertex path: ");
    for (int64_t vid : path_deque) {
        st_uglyf(" %lld ", vid);

    }
    st_uglyf("\n");
    // TODO TODO left off here, 
    // TODO go through the vertexPath (from the HmmGraph) and make the alignment from the 
    // TODO graph aligned pairs
}

std::unordered_map<int64_t, AlignedPairs> HmmAligner::mapAlignedPairsToVgNodes(int64_t pId) {
    std::unordered_map<int64_t, AlignedPairs> mapped_pairs;  // (vg_node_id, AlignedPairs)
    // Get the GraphAligedPairs for this path
    GraphAlignedPairs& graph_aligned_pairs = hmm_graph.PathAlignedPairs()[pId];
    st_uglyf("SENTINAL mapping path %lld, got %lld pairs\n", pId, graph_aligned_pairs.size());
    if (graph_aligned_pairs.size() == 0) return mapped_pairs;  // there are no aligned pairs to this path
    // map the vertex to to the VG node, and the offset becomes the `y` in the aligned pair
    for (GraphAlignedPair& gpair : graph_aligned_pairs) {
        double posterior_prob = std::get<0>(gpair);
        int64_t x             = std::get<1>(gpair);
        int64_t vid           = std::get<2>(gpair).first;
        int64_t offset        = std::get<2>(gpair).second;
        int64_t vg_node_id    = vertexId_to_nodeId[vid];
        st_uglyf("x: %lld vid: %lld offset: %lld vg_node_id: %lld, prob: %f\n", 
                 x, vid, offset, vg_node_id, posterior_prob);
        AlignedPair aligned_pair = std::make_tuple(posterior_prob, x, offset);
        mapped_pairs[vg_node_id].push_back(aligned_pair);
    }
    return mapped_pairs;
}

std::vector<double> HmmAligner::pathScoresVector() {
    auto path_scores_map = hmm_graph.PathScores();
    std::vector<double> scores;
    scores.reserve(path_scores_map.size());
    
    std::transform(begin(path_scores_map), end(path_scores_map), begin(scores), 
                   [] (std::pair<int64_t, double> p) { return p.second; });
    
    std::sort(begin(scores), end(scores), [](double i, double j) { return i > j; });

    for (auto i : scores) st_uglyf("%lld \n", i);

    return scores;

}
