#include <iostream>
#include "hmm_aligner.hpp"

using namespace vg;

HmmAligner::HmmAligner(Graph& G) {
    st_uglyf("[HmmAligner::HmmAligner]: Building HmmGraph from Graph object, graph has %lld nodes "
             "and %lld edges\n", G.node_size(), G.edge_size());
    st_uglyf("[HmmAligner::HmmAligner]: internal graph has size %lld\n", hmm_graph.K());
    
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

void HmmAligner::Align(Alignment& aln, std::vector<Alignment>* path_alignments, AlignmentParameters& p, 
                       bool get_pairs, bool ragged_end) {
    // step 1: align the sequence to the graph with the HMM
    std::string* sx = aln.mutable_sequence();
    if (sx->length() <= 0) return;
    // TODO need try/catch here:
    hmm_graph.AlignWithFiveStateSymbolHmm((*sx), p, get_pairs, ragged_end);
    // debug
    PrintAlignedPairs();
    // step 2: if not doing multi-path alignment, just add the one with the max score
    if (!path_alignments) {
        auto max_pathId = hmm_graph.MaxScorePath();
        makeAlignmentFromPathAlignedPairs(aln, max_pathId, ragged_end);
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
        // 
        // TODO Print vertex path 
        //
        std::cerr << "path: " << kv.first << " score: " << graph_path_scores[kv.first];
        if (kv.first == most_prob_path) std::cerr << " <-- Highest Score Path\n";
        else std::cerr << "\n";
        for (auto p : kv.second) {
            int64_t x_coord = std::get<1>(p);
            int64_t vid     = std::get<2>(p).first;
            int64_t offset  = std::get<2>(p).second;
            double prob     = std::get<0>(p)/PAIR_ALIGNMENT_PROB_1;
            std::cerr << "x: " << x_coord << " vertexID: " << vid << " NodeId: " << 
                vertexId_to_nodeId[vid] << " offset: " << offset << " prob: " << prob << std::endl;
        }
    }
}

void HmmAligner::makeAlignmentFromPathAlignedPairs(Alignment& aln, int64_t pId, bool ragged_end) {
    // helper functions for modifying the alignment 
    auto do_insert = [] (Mapping *mapping, std::string& read_sequence, int64_t x, int64_t pX) {
        Edit *ed = mapping->add_edit();
        // for an insert the length in the node sequence is zero (not there)
        ed->set_from_length(0);
        // say you have an insert in the read.. 
        // 012 34  [y]
        // ACG TA  [graph sequence]
        // ||| ||
        // ACGGTA  [read sequence]
        // 012345  [x]
        // pX = 2 (G <> G), x = 4 (T <> T), length of insert = 4 - 2 - 1 = 1;
        int64_t insert_length   = x - pX - 1;
        int64_t start_of_insert = pX + 1;
        // the `to_length` is the length of the inserted sequence
        ed->set_to_length(insert_length);
        ed->set_sequence(read_sequence.substr(start_of_insert, insert_length));
        st_uglyf("SENTINAL there's an INSERT, x: %lld, pX, %lld, length: %lld\n", x, pX, insert_length);
    };
    
    auto do_delete = [] (Mapping *mapping, int64_t delete_length) {  // TODO change to 'delete_length'
        Edit *ed = mapping->add_edit();
        // the `to_length `is zero (not present in the read)
        ed->set_to_length(0);
        // say you have a delete in the read..
        // 01234  [y]
        // ACGTA  [graph sequence]
        // || ||
        // AC TA  [read sequence]
        // 01 23  [x]
        // pY = 1 (C <> C pair), y = 3, length of delete = 3 - 1 - 1 = 1
        //int64_t delete_length = y - pY - 1;
        // for a delete the `from_length` is the length of the deleted sequence (not in the read)
        ed->set_from_length(delete_length);
        st_uglyf("SENTINAL there's an DELETE of length %lld\n", delete_length);
    };

    auto do_pair = [] (Mapping *mapping, const std::string& node_sequence, std::string& read_sequence, 
                       int64_t x, int64_t y) {
        Edit *ed = mapping->add_edit();
        ed->set_from_length(ALIGNED_PAIR_LENGTH);
        ed->set_to_length(ALIGNED_PAIR_LENGTH);
        char graph_base = node_sequence.at(y);
        char read_base  = read_sequence.at(x);
        st_uglyf("SENTINAL graph seq: %c read seq: %c", graph_base, read_base);
        if (graph_base != read_base) {
            st_uglyf(" there's a SNP\n");
            ed->set_sequence(read_sequence.substr(x, ALIGNED_PAIR_LENGTH));
        } else {
            st_uglyf(" just a MATCH\n");
        }
    };
    // Main part of function starts here!
    aln.clear_path();
    aln.set_score(hmm_graph.PathScores()[pId] * PAIR_ALIGNMENT_PROB_1);
    //
    // NOTE gssw_aligner sets the alignment `query_position` to 0 here, I'm not sure what that is 
    // supposed to do, so I don't set that member, but should remember that it's NOT set
    //

    Path *aln_path    = aln.mutable_path();               // an Alignment has one Path, so get that
    string& read_seq  = *aln.mutable_sequence();
    auto mapped_pairs = mapAlignedPairsToVgNodes(pId);    // the aligned pairs mapped to the vg nodes
    auto vertex_path  = hmm_graph.PathMap()[pId];         // the vertex IDs (in order) for this path

    // 
    // first need to find fist vertex with an aligned pair in this path and set that to the 'start'
    // of the alignment slice the vertex path to the first (and last?) vertex that has aligned pairs
    // TODO TODO
    //
 
    int64_t pX = -1;                                      // previous matching X (read coordinate)
    for (int64_t vid : vertex_path) {
        st_uglyf("SENTINAL: looking at vertex %lld (vg node: %lld)\n", vid, vertexId_to_nodeId[vid]);
        int64_t vg_node_id              = vertexId_to_nodeId[vid];
        const std::string& node_seq     = *hmm_graph.VertexSequence(vid);
        AlignedPairs node_aligned_pairs = mapped_pairs[vg_node_id];
        Mapping *mapp = aln_path->add_mapping();
        mapp->mutable_position()->set_node_id(vg_node_id);
        mapp->set_rank(aln_path->mapping_size());         // set this mapping's rank (order in mappings)
        if (node_aligned_pairs.empty()) {                 // there are no aligned pairs to this node's sequence
            st_uglyf("SENTINAL DELETE THROUGH VERTEX %lld doing delete\n", vid);
            mapp->mutable_position()->set_offset(0);      // offset for an entire node delete is 0
            // y = len(node_seq) + 1 as a hack to use the same function, basically we want the delete 
            // to be the length of node_seq
            do_delete(mapp, static_cast<int>(node_seq.size()));
            st_uglyf("after adding edits mapping:\n%s\n", mapp->DebugString().c_str());
            continue;
        }
        int64_t fY = std::get<2>(node_aligned_pairs.at(0));  // first node sequence aligned-to position
        st_uglyf("fY: %lld\n", fY);
        assert(fY >= 0 && fY < static_cast<int>(node_seq.size()));
        mapp->mutable_position()->set_offset(0);  // we default to global alignment, see note below
        if (fY != 0) {                            // there is a leading delete
            //
            // NOTE need some logic here to decide about global/local alignment
            // if we're at the first vertex, might want to start the alignment
            // at the first aligned pair and not at the 0th position
            //
            do_delete(mapp, fY);                      
        } 
        // first aligned-to position in the read, for this vertex
        //int64_t fX = std::get<1>(node_aligned_pairs.at(0)); 
        //if (``)
        st_uglyf("before adding edits mapping:\n%s\n", mapp->DebugString().c_str());
        // add the aligned pairs (as edits) 
        //int64_t pY = 0;                                   // previous offset that matched
        int64_t pY = std::get<2>(node_aligned_pairs.at(0));
        for (AlignedPair& pr : node_aligned_pairs) {
            int64_t x = std::get<1>(pr);                  // read coordinate
            int64_t y = std::get<2>(pr);                  // offset in node
            // check the pairs
            assert(x >= 0 && x < static_cast<int>(read_seq.size()));
            assert(y >= 0 && y < static_cast<int>(node_seq.size()));

            if (x - pX > 1) do_insert(mapp, read_seq, x, pX);
            if (y - pY > 1) do_delete(mapp, y - pY - 1);
            do_pair(mapp, node_seq, read_seq, x, y);
            pX = x;
            pY = y;
        }
        st_uglyf("SENTINAL - final pY %lld, node seq length %d\n", pY, node_seq.size() - 1);
        
        // deal with deletes at the end of non-terminal nodes 
        // note: node_seq.size() - 1 = the number of aligned pairs you'd get if you had 
        // perfect matching
        if ((pY != static_cast<int>(node_seq.size()) - 1)) {
            int64_t delete_length = (node_seq.size() - 1) - pY;
            do_delete(mapp, delete_length);
        }
        st_uglyf("after adding edits mapping:\n%s\n", mapp->DebugString().c_str());
        
    }
    double percent_identity = identity(aln.path());
    st_uglyf("percent identity:%f\n", percent_identity);
    aln.set_identity(percent_identity);
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
        AlignedPair aligned_pair = std::make_tuple(posterior_prob, x, offset);  // 
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
