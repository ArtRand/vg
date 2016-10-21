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
#include "sonLib.h"        // st_uglyf 
#include "utility.hpp"     // nonACGTNtoN

/* Vertex class.
 * Basic class for vertices (aka nodes) in a DAG. 
 */
class Vertex {
public:
    // constructors
    Vertex(): id(-1), base('\0') {};  // default construct with ID = -1
    Vertex(int64_t i);
    Vertex(int64_t i, char b);
    
    // destructor TODO test for leak 
    ~Vertex() {};

    void SetId(int64_t i);  // shouldn't ever need this

    int64_t Id() const;

    char Base() const;

    bool operator < (Vertex other) const;

private:
    int64_t id;
    char base;
};

// intermediate graph structure
class HmmGraph {
public:
    HmmGraph(): nVertices(0), nb_arcs(0), next_vertex_id(0) {st_uglyf("made HmmGrph\n");};

    ~HmmGraph() {};

    int64_t AddVertex(char base);

    bool ContainsVertex(int64_t i);

    bool ContainsVertex(Vertex& v);

    void AddArc(int64_t, int64_t);

    int64_t K();

    Vertex *VertexGetter(int64_t i);

    unsigned long VertexOutDegree(int64_t i);

    unsigned long VertexInDegree(int64_t i);

    std::vector<int64_t>& Vertices();

    void TopologicalSort();

    bool TestSort();

    bool IsSorted();

private:
    // containers
    std::unordered_map<int64_t, std::set<int64_t>> adjacentcy_list;
    std::unordered_map<int64_t, Vertex*> vertex_map;  // for looking up vertices
    std::vector<int64_t> vertex_list;  // for ordering
    std::set<std::string> labels;
    std::set<std::string> seqs;
    std::vector<int64_t> starts;

    // counters and flags
    int64_t nVertices;
    int64_t nb_arcs;
    int64_t next_vertex_id;
    bool sorted = false;
};


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
