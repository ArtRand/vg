#include <iostream>
#include "hmm_aligner.hpp"

using namespace vg;

Vertex::Vertex(int64_t i, char b) {
    id = i;
    base = b;
}

int64_t Vertex::Id() const {
    return id;
}

int64_t HmmGraph::AddVertex(char base) {
    // get the next vertex id, a monotonically increasing int
    int64_t vId = next_vertex_id;

    Vertex *v = new Vertex(vId, base);

    // check, todo remove after tested
    assert(v->Id() == vId);

    // updates to state of graph
    vertex_map[vId] = v;
    vertex_list.push_back(vId);
    nVertices += 1;
    next_vertex_id += 1;
    sorted = false;

    return vId;
}

void HmmGraph::AddArc(int64_t fromId, int64_t toId) {
    adjacentcy_list[fromId].insert(toId);
}

int64_t HmmGraph::K() { return nVertices; }

// this function splits up the nodes from the graph into single-base nodes and makes a HmmGraph object
HmmAligner::HmmAligner(Graph& G) {
    st_uglyf("[HmmAligner]: Building HmmGraph from Graph object, graph has %lld nodes "
             "and %lld edges\n", G.node_size(), G.edge_size());
    st_uglyf("[HmmAligner]: internal graph has size %lld\n", hmm_graph.K());
    
    // we need to keep track of which nodes in the original the vertices in the new graph came from
    std::unordered_map<int64_t, std::vector<int64_t>> parents;  // key: parent value: chidren
    
    // go through the nodes and assemble the HmmGraph, split nodes with len(seq) > 1 into multiple nodes
    for (int64_t i = 0; i < G.node_size(); i++) {
        Node *n = G.mutable_node(i);

        // clean the sequence, got this idea from gssw_aligner.cpp
        auto cleaned_seq = nonATGCNtoN(n->sequence()); 
        
        st_uglyf("node %lld has sequence %s\n", n->id(), cleaned_seq.c_str());

        // if a single nucleotide node, just add it
        if (cleaned_seq.length() == 1) {  
            int64_t vId = hmm_graph.AddVertex(cleaned_seq.at(0));
            st_uglyf("Adding vertex %lld as a single base sequence %s\n", vId, cleaned_seq.c_str());
            parents[n->id()].push_back(vId);
            continue;
        }
        
        // if the parent node has a "sequence" (ie more than one base) split it up
        int64_t firstId = -1;
        int64_t lastId = -1;
        for (auto b : cleaned_seq)  {
            int64_t vId = hmm_graph.AddVertex(b);
            parents[n->id()].push_back(vId);
            st_uglyf("Adding vertex %lld from multi-base sequence %c\n", vId, b);
            if (firstId < 0) {  // if we're at the start of the sequence
                firstId = vId;
            }
            
            // connect up this vertex to the previous one, 
            // (the Graph object will not contain these edges)
            if (lastId >= 0) {  
                hmm_graph.AddArc(lastId, vId);
            }
            
            lastId = vId;
        }
    }
    st_uglyf("[HmmAligner]: internal graph has size %lld\n", hmm_graph.K());
    
    // in the above routine vertices spawned from "sequence" nodes are connected already and we 
    // only want to address the first/last one, these two functions quicly get the first or
    // last vertex in a block
    auto get_firstLast_vertex = [&](int64_t id, bool get_first) -> int64_t {
        std::vector<int64_t> chld = parents[id];
        if (chld.size() == 1) {
            return chld.at(0);
        } else if (get_first) {
            return *(std::min_element(std::begin(chld), std::end(chld)));
        } else {
            return *(std::max_element(std::begin(chld), std::end(chld)));
        }
    };

    // now go through the edges
    for (int64_t i = 0; i < G.edge_size(); i++) {
        Edge *e = G.mutable_edge(i);

        // these are the nodeIDs from the graph object
        int64_t fromId = e->from();
        int64_t toId = e->to();
        
        if (!e->from_start() && !e->to_end()) {  // normal edge
            st_uglyf("[HmmAligner]: Adding end->start edge ");

            int64_t fromVertex = get_firstLast_vertex(fromId, false);
            int64_t toVertex = get_firstLast_vertex(toId, true);
            
            st_uglyf("Adding arc from %lld to %lld\n", fromVertex, toVertex);
            hmm_graph.AddArc(fromVertex, toVertex);

        } else if (e->from_start() && e->to_end()) {
            // TODO needs testing
            // this is a 'normal' edge but it's orientation is reversed
            st_uglyf("[HmmAligner]: Adding start->end edge ");
            
            int64_t fromVertex = get_firstLast_vertex(toId, false);
            int64_t toVertex = get_firstLast_vertex(fromId, true);

            st_uglyf("Adding arc from %lld to %lld\n", fromVertex, toVertex);
            hmm_graph.AddArc(fromVertex, toVertex);
        }
    }

}

HmmAligner::~HmmAligner() {};
