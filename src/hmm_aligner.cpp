#include <iostream>
#include "hmm_aligner.hpp"

using namespace vg;

HmmAligner::HmmAligner(Graph& G) {
    st_uglyf("[HmmAligner]: Building HmmGraph from Graph object, graph has %lld nodes "
             "and %lld edges\n", G.node_size(), G.edge_size());
    st_uglyf("[HmmAligner]: internal graph has size %lld\n", hmm_graph.K());
    
    // mapps the nodeId from VG::Graph to VertexId in HmmGraph
    std::unordered_map<int64_t, int64_t> nodeId_to_vertexId;

    for (int64_t i = 0; i < G.node_size(); i++) {
        Node *n = G.mutable_node(i);

        st_uglyf("node %lld has sequence %s\n", n->id(), n->sequence().c_str());;
        
        // add the vertex to the HmmGraph and add the node id (from the orig. graph) to
        // the hash so we can look it up later when we're adding the edges
        int64_t vId = hmm_graph.AddVertex(&n->sequence());
        nodeId_to_vertexId[n->id()] = vId;
    }

    st_uglyf("[HmmAligner]: internal graph has %lld vertices\n", hmm_graph.K());
    
    // now go through the edges
    for (int64_t i = 0; i < G.edge_size(); i++) {
        Edge *e = G.mutable_edge(i);
        
        int64_t fromNodeId = e->from();
        int64_t toNodeId = e->to();
        int64_t fromVertex = nodeId_to_vertexId[fromNodeId];
        int64_t toVertex = nodeId_to_vertexId[toNodeId];

        if (!e->from_start() && !e->to_end()) {  // normal edge
            st_uglyf("[HmmAligner]: Adding end->start edge from %lld to %lld ", fromNodeId, toNodeId);
            st_uglyf("Adding arc from %lld to %lld \n", fromVertex, toVertex);
            hmm_graph.AddArc(fromVertex, toVertex);

        } else if (e->from_start() && e->to_end()) {
            // TODO needs testing
            // this is a 'normal' edge but it's orientation is reversed
            st_uglyf("[HmmAligner]: Adding start->end edge from %lld to %lld ", fromNodeId, toNodeId);
            st_uglyf("Adding arc from %lld to %lld\n", fromVertex, toVertex);
            hmm_graph.AddArc(toVertex, fromVertex);
        }
    }
}

HmmAligner::~HmmAligner() {};
