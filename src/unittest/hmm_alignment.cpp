//
// hmm_alignment.cpp
//  
// Unit tests for HMM-based aligner interacting with VG
//

#include <stdio.h>
#include <unordered_set>

#include "sonLib.h"  // for st_uglyf
#include "hmm_aligner.hpp"
#include "alignment.hpp"
#include "vg.hpp"
#include "path.hpp"
#include "catch.hpp"

using namespace google::protobuf;


namespace vg {
namespace unittest {
// Checks runs of matching alignmed pairs
static inline void CheckMappingPairs(const Path& path, int64_t mapping_idx, int64_t expected) {
    REQUIRE(path.mapping(mapping_idx).edit_size() == expected);
    // each has correct length (1) and is a match
    for (int64_t i = 0; i < path.mapping(mapping_idx).edit_size(); ++i) {
        REQUIRE(path.mapping(mapping_idx).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(mapping_idx).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(mapping_idx).edit(i).sequence().empty());
    }
}

static inline void CheckForGlobalAlignment(const Path& path, Node *final_node) {
    // alignment must start at first position in first node
    REQUIRE(path.mapping(0).position().offset() == 0);
    auto final_mapping = [&path] () -> const Mapping * {
        int64_t idx = -1;
        const Mapping *m;
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            REQUIRE(path.mapping(i).rank() != idx);
            if (path.mapping(i).rank() > idx) {
                idx = path.mapping(i).rank();
                m = &path.mapping(i);
            }
        }
        return m;
    }();
    REQUIRE(mapping_from_length(*final_mapping) == final_node->sequence().size());
}

TEST_CASE( "HmmAligner setup tests", "[hmm]" ) {
    SECTION( "HMM aligner makes correct HmmGraph with all normal edges" ) {
        VG graph;
        
        Node* n0 = graph.create_node("AGTG");
        Node* n1 = graph.create_node("C");
        Node* n2 = graph.create_node("A");
        Node* n3 = graph.create_node("TGAAGT");
        
        graph.create_edge(n0, n1);
        graph.create_edge(n0, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n3);
        
        HmmAligner hmm(graph.graph);
        
        // test that the internal graph has a node for each nucleoide in the Graph
        REQUIRE(hmm.Graph_K() == graph.node_count());
    }

    SECTION( "HMM aligner makes correct HmmGraph with backwards edges" ) {
        VG graph;

        Node *n0 = graph.create_node("A");
        Node *n1 = graph.create_node("C");
        Node *n2 = graph.create_node("G");
        Node *n3 = graph.create_node("T");

        graph.create_edge(n0, n2, false, false);
        graph.create_edge(n1, n0, true, true);
        graph.create_edge(n1, n2, false, false);
        graph.create_edge(n2, n3, false, false);

        HmmAligner hmm(graph.graph);
        REQUIRE(hmm.Graph_K() == graph.node_count());
    }

    SECTION( "HMM aligner complains with reversing edges" ) {
        VG graph;

        Node *n0 = graph.create_node("A");
        Node *n1 = graph.create_node("C");
        Node *n2 = graph.create_node("G");
        Node *n3 = graph.create_node("T");

        graph.create_edge(n0, n2, false, false);
        graph.create_edge(n1, n0, true, false);  // reversing edge
        graph.create_edge(n1, n2, false, false);
        graph.create_edge(n2, n3, false, false);

        REQUIRE_THROWS_AS(HmmAligner hmm(graph.graph), ParcoursException);
    }
}

TEST_CASE( "Hmm aligner produces correct alignment when read matches exactly (max score path)",
          "[hmm]" ) {
    VG graph;
    
    Node* vid0 = graph.create_node("AGC");
    Node* vid1 = graph.create_node("TT");
    Node* vid2 = graph.create_node("A");
    Node* vid3 = graph.create_node("CGT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
    
    string read = string("AGCACGT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    
    // has correct number of mappings (nodes mapped to)
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    
    // has corrects edits
    CheckMappingPairs(path, 0, 3);
    CheckMappingPairs(path, 1, 1);
    CheckMappingPairs(path, 2, 3);

    REQUIRE(aln.identity() == 1.0);
}

TEST_CASE("Hmm aligner produces correct alignment when there is a mismatch", "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AGCA");
    Node* vid1 = graph.create_node("TT");
    Node* vid2 = graph.create_node("A");
    Node* vid3 = graph.create_node("CGT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // AGCAACGT
                       // ||:|||||
    string read = string("AGTAACGT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();

    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // check the correct number of edits
    REQUIRE(path.mapping(0).edit_size() == 4);
    REQUIRE(path.mapping(1).edit_size() == 1);
    REQUIRE(path.mapping(2).edit_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    // check the matching region
    REQUIRE(path.mapping(0).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    REQUIRE(path.mapping(0).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).sequence().empty());
    // check the SNP
    REQUIRE(path.mapping(0).edit(2).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).sequence() == std::string("T"));
    // check the remaining matching section
    REQUIRE(path.mapping(0).edit(3).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(3).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(3).sequence().empty());
    // check the rest of the alignment
    CheckMappingPairs(path, 1, 1);
    CheckMappingPairs(path, 2, 3);
}

TEST_CASE("Hmm Aligner produces correct alignment when there is a single base deletion", "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AAGTAGCA");
    Node* vid1 = graph.create_node("TTTT");
    Node* vid2 = graph.create_node("AA");
    Node* vid3 = graph.create_node("CGTAACAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                        
                       // AAGTAGCAAACGTAACAT  [graph sequence]
                       // AAGT-GCAAACGTAACAT  [read sequence]
    string read = string("AAGTGCAAACGTAACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    // check the matches preceding the deletion
    REQUIRE(path.mapping(0).edit_size() == vid0->sequence().size());
    REQUIRE(path.mapping(0).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    REQUIRE(path.mapping(0).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).sequence().empty());
    REQUIRE(path.mapping(0).edit(2).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).sequence().empty());
    REQUIRE(path.mapping(0).edit(3).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(3).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(3).sequence().empty());
    // check the deletion 
    REQUIRE(path.mapping(0).edit(4).from_length() == 1);
    REQUIRE(path.mapping(0).edit(4).to_length() == 0);
    REQUIRE(path.mapping(0).edit(4).sequence().empty());
    // check the remaining matching section 
    REQUIRE(path.mapping(0).edit(5).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(5).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(5).sequence().empty());
    REQUIRE(path.mapping(0).edit(6).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(6).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(6).sequence().empty());
    REQUIRE(path.mapping(0).edit(7).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(7).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(7).sequence().empty());
    CheckMappingPairs(path, 1, vid2->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when there is a multi-base deletion", "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AAGTAGCA");
    Node* vid1 = graph.create_node("TTTT");
    Node* vid2 = graph.create_node("AA");
    Node* vid3 = graph.create_node("CGTAACAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // |------|--|------|  [node bounds]
                       // AAGTAGCAAACGTAACAT  [graph sequence]
                       // AAGT--CAAACGTAACAT  [read sequence]
    string read = string("AAGTCAAACGTAACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    // check the matches preceding the deletion
    REQUIRE(path.mapping(0).edit_size() == vid0->sequence().size() - 1);  // deletion becomes one edit
    REQUIRE(path.mapping(0).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    REQUIRE(path.mapping(0).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).sequence().empty());
    REQUIRE(path.mapping(0).edit(2).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).sequence().empty());
    REQUIRE(path.mapping(0).edit(3).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(3).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(3).sequence().empty());
    // check the deletion 
    REQUIRE(path.mapping(0).edit(4).from_length() == 2);
    REQUIRE(path.mapping(0).edit(4).to_length() == 0);
    REQUIRE(path.mapping(0).edit(4).sequence().empty());
    // check the remaining matching section 
    REQUIRE(path.mapping(0).edit(5).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(5).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(5).sequence().empty());
    REQUIRE(path.mapping(0).edit(6).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(6).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(6).sequence().empty());
    CheckMappingPairs(path, 1, vid2->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}
TEST_CASE("Hmm Aligner produces correct alignment when there is a single base insertion", "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AAGTAGCA");
    Node* vid1 = graph.create_node("TTTT");
    Node* vid2 = graph.create_node("AA");
    Node* vid3 = graph.create_node("CGAACAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                        
                       // AAGTAGCAAACG-AACAT  [graph sequence]
                       // AAGTAGCAAACGTAACAT  [read sequence]
    string read = string("AAGTAGCAAACGTAACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    // check the first two perfectly matching mappings
    CheckMappingPairs(path, 0, vid0->sequence().size());
    CheckMappingPairs(path, 1, vid2->sequence().size());
    // check the mapping with the insert
    REQUIRE(path.mapping(2).edit_size() == vid3->sequence().size() + 1);  // plus 1 because of the insert edit
    // matching section
    REQUIRE(path.mapping(2).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(0).sequence().empty());
    REQUIRE(path.mapping(2).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(1).sequence().empty());
    // insert
    REQUIRE(path.mapping(2).edit(2).from_length() == 0);
    REQUIRE(path.mapping(2).edit(2).to_length() == 1);
    REQUIRE(path.mapping(2).edit(2).sequence() == string("T"));
    // the rest of the matching section
    REQUIRE(path.mapping(2).edit(3).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(3).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(3).sequence().empty());
    REQUIRE(path.mapping(2).edit(4).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(4).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(4).sequence().empty());
    REQUIRE(path.mapping(2).edit(5).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(5).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(5).sequence().empty());
    REQUIRE(path.mapping(2).edit(6).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(6).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(6).sequence().empty());
    REQUIRE(path.mapping(2).edit(7).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(7).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(7).sequence().empty());
}

TEST_CASE("Hmm Aligner produces correct alignment when there is a multi-base insertion", "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AAGTAGCA");
    Node* vid1 = graph.create_node("TCATTT");
    Node* vid2 = graph.create_node("AA");
    Node* vid3 = graph.create_node("CGACAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                        
                       // AAGTAGCAAACG--ACAT  [graph sequence]
                       // AAGTAGCAAACGTAACAT  [read sequence]
                       // 012345678901234567
    string read = string("AAGTAGCAAACGTAACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    // check the first two perfectly matching mappings
    CheckMappingPairs(path, 0, vid0->sequence().size());
    CheckMappingPairs(path, 1, vid2->sequence().size());
    // check the mapping with the insert
    REQUIRE(path.mapping(2).edit_size() == vid3->sequence().size() + 1);  // plus 1 because of the insert edit
    // matching section
    REQUIRE(path.mapping(2).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(0).sequence().empty());
    REQUIRE(path.mapping(2).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(1).sequence().empty());
    // insert
    REQUIRE(path.mapping(2).edit(2).from_length() == 0);
    REQUIRE(path.mapping(2).edit(2).to_length() == 2);
    REQUIRE(path.mapping(2).edit(2).sequence() == string("TA"));
    // the rest of the matching section
    REQUIRE(path.mapping(2).edit(3).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(3).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(3).sequence().empty());
    REQUIRE(path.mapping(2).edit(4).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(4).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(4).sequence().empty());
    REQUIRE(path.mapping(2).edit(5).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(5).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(5).sequence().empty());
    REQUIRE(path.mapping(2).edit(6).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(6).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(6).sequence().empty());
}

TEST_CASE("Hmm Aligner produces correct alignment when there is a whole node deletion", "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("CATG");
    Node* vid1 = graph.create_node("CCTT");
    Node* vid2 = graph.create_node("AAA");
    Node* vid3 = graph.create_node("CTAG");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // CATGAAACTAG
                       // CATG---CTAG
    string read = string("CATGCTAG");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckMappingPairs(path, 0, vid0->sequence().size());
    REQUIRE(path.mapping(1).edit_size() == 1);                                  // only one edit, all delete
    REQUIRE(path.mapping(1).edit(0).to_length() == 0);                          // is a delete
    // delete is the length of the enire node's sequence 
    REQUIRE(path.mapping(1).edit(0).from_length() == vid2->sequence().size());
                                                                                
    CheckMappingPairs(path, 0, vid0->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when there is a deletion that crosses a node boundary", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("CATG");
    Node* vid1 = graph.create_node("CCTT");
    Node* vid2 = graph.create_node("AAA");
    Node* vid3 = graph.create_node("CTAG");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // CATGAAACTAG
                       // CAT--AACTAG
    string read = string("CATAACTAG");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    CheckForGlobalAlignment(path, vid3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    // first mapping
    // should have a mapping for each base in the sequence
    REQUIRE(path.mapping(0).edit_size() == vid0->sequence().size()); 
    REQUIRE(path.mapping(0).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    REQUIRE(path.mapping(0).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(1).sequence().empty());
    REQUIRE(path.mapping(0).edit(2).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(0).edit(2).sequence().empty());
    REQUIRE(path.mapping(0).edit(3).from_length() == 1);  // the delete at the end of the node
    REQUIRE(path.mapping(0).edit(3).to_length() == 0);    // the delete at the end of the node
    REQUIRE(path.mapping(0).edit(3).sequence().empty()); 

    // second mapping
    // should have a mapping for each base in the sequence
    REQUIRE(path.mapping(1).edit_size() == vid2->sequence().size()); 
    REQUIRE(path.mapping(1).edit(0).from_length() == 1);
    REQUIRE(path.mapping(1).edit(0).to_length() == 0);
    REQUIRE(path.mapping(1).edit(1).from_length() == 1);
    REQUIRE(path.mapping(1).edit(1).to_length() == 1);
    REQUIRE(path.mapping(1).edit(1).sequence().empty());
    REQUIRE(path.mapping(1).edit(2).from_length() == 1);
    REQUIRE(path.mapping(1).edit(2).to_length() == 1);
    REQUIRE(path.mapping(1).edit(2).sequence().empty());
    
    // third mapping
    CheckMappingPairs(path, 2, vid3->sequence().size());
    
}

TEST_CASE("Hmm Aligner produces correct alignment when the read doesn't span the entire graph", 
        "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AAGTAGCA");
    Node* vid1 = graph.create_node("TTT");
    Node* vid2 = graph.create_node("AA");
    Node* vid3 = graph.create_node("CGACAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // AAGTAGCAAAGCACAT
                       // AAGTAGCAAAGCA---
    string read = string("AAGTAGCAAACGA");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    // first two node mappings are all matches
    CheckMappingPairs(path, 0, vid0->sequence().size());
    CheckMappingPairs(path, 1, vid2->sequence().size());
    // the last mapping should only contain 3 edits, to the first three bases
    REQUIRE(path.mapping(2).edit_size() == 4);  // 3 matches and 1 delete
    REQUIRE(path.mapping(2).edit(0).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(0).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(0).sequence().empty());
    REQUIRE(path.mapping(2).edit(1).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(1).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(1).sequence().empty());
    REQUIRE(path.mapping(2).edit(2).from_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(2).to_length() == ALIGNED_PAIR_LENGTH);
    REQUIRE(path.mapping(2).edit(2).sequence().empty());
    REQUIRE(path.mapping(2).edit(3).from_length() == 3);  // final delete
    REQUIRE(path.mapping(2).edit(3).to_length() == 0);
}

TEST_CASE("Hmm Aligner produces correct alignment when it starts with a single base insertion", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AACCCAGG");
    Node* vid1 = graph.create_node("CA");
    Node* vid2 = graph.create_node("ATA");
    Node* vid3 = graph.create_node("TGAAGT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // -AACCCAGGATTAGAAGT
    string read = string("TAACCCAGGATATGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    REQUIRE(path.mapping(0).edit(0).from_length() == 0);
    REQUIRE(path.mapping(0).edit(0).to_length() == 1);
    REQUIRE(path.mapping(0).edit(0).sequence() == string("T"));
    for (int64_t i = 1; i <= 8; ++i) {
        REQUIRE(path.mapping(0).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).sequence().empty());
    }
    CheckMappingPairs(path, 1, vid2->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when it starts with a multi-base insertion", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AACCCAGG");
    Node* vid1 = graph.create_node("CA");
    Node* vid2 = graph.create_node("ATA");
    Node* vid3 = graph.create_node("TGAAGT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                       // --AACCCAGGATTAGAAGT
    string read = string("TGAACCCAGGATATGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    REQUIRE(path.mapping(0).edit(0).from_length() == 0);
    REQUIRE(path.mapping(0).edit(0).to_length() == 2);
    REQUIRE(path.mapping(0).edit(0).sequence() == string("TG"));
    for (int64_t i = 1; i <= 8; ++i) {
        REQUIRE(path.mapping(0).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).sequence().empty());
    }
    CheckMappingPairs(path, 1, vid2->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when it starts with a single base deletion", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("GAT");
    Node* vid1 = graph.create_node("A");
    Node* vid2 = graph.create_node("ATA");
    Node* vid3 = graph.create_node("CAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                      // GARATACAT
    string read = string("ATATACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    REQUIRE(path.mapping(0).edit(0).from_length() == 1);
    for (int64_t i = 1; i <= 2; ++i) {
        REQUIRE(path.mapping(0).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).sequence().empty());
    }
    CheckMappingPairs(path, 1, vid2->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when it starts with a multi-base deletion", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("GGAT");
    Node* vid1 = graph.create_node("ATTTT");
    Node* vid2 = graph.create_node("ATA");
    Node* vid3 = graph.create_node("CAT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                      // 
    string read = string("ATATACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid2->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    REQUIRE(path.mapping(0).edit(0).from_length() == 2);
    for (int64_t i = 1; i <= 2; ++i) {
        REQUIRE(path.mapping(0).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(0).edit(i).sequence().empty());
    }
    CheckMappingPairs(path, 1, vid2->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when it starts not at the first node", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AA");
    Node* vid1 = graph.create_node("CCCAGGCA");
    Node* vid2 = graph.create_node("GTGCTATA");
    Node* vid3 = graph.create_node("TGAAGT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                      // 
    string read = string("CCCAGGCATGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid1->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    REQUIRE(path.mapping(0).edit_size() == 1);
    REQUIRE(path.mapping(0).edit(0).from_length() == 2);
    REQUIRE(path.mapping(0).edit(0).to_length() == 0);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    // test the rest of the edits
    CheckMappingPairs(path, 1, vid1->sequence().size());
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when it begins with a deletion that crosses a node boundary", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("AA");
    Node* vid1 = graph.create_node("CCCAGGCA");
    Node* vid2 = graph.create_node("GTGCTATA");
    Node* vid3 = graph.create_node("TGAAGT");
    
    graph.create_edge(vid0, vid1);
    graph.create_edge(vid0, vid2);
    graph.create_edge(vid1, vid3);
    graph.create_edge(vid2, vid3);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                      // 
    string read = string("CCAGGCATGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
    // follows correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid1->id());
    REQUIRE(path.mapping(2).position().node_id() == vid3->id());
    CheckForGlobalAlignment(path, vid3);
    // has the correct edits, should skip fitst node and first position in second node
    REQUIRE(path.mapping(0).edit_size() == 1);
    REQUIRE(path.mapping(0).edit(0).from_length() == 2);
    REQUIRE(path.mapping(0).edit(0).to_length() == 0);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    REQUIRE(path.mapping(1).edit_size() == vid1->sequence().size());
    REQUIRE(path.mapping(1).edit(0).from_length() == 1);
    REQUIRE(path.mapping(1).edit(0).to_length() == 0);
    REQUIRE(path.mapping(1).edit(0).sequence().empty());
    for (int64_t i = 1; i < vid1->sequence().size(); ++i) {
        REQUIRE(path.mapping(1).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(1).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
        REQUIRE(path.mapping(1).edit(i).sequence().empty());
    }
    CheckMappingPairs(path, 2, vid3->sequence().size());
}

TEST_CASE("Hmm Aligner produces correct alignment when it begins with a deletion that ends at a node boundary", 
          "[hmm]") {
    VG graph;
    
    Node* vid0 = graph.create_node("A");
    Node* vid1 = graph.create_node("C");
    
    graph.create_edge(vid0, vid1);

    HmmAligner hmm(graph.graph);

    AlignmentParameters p;
    p.expansion   = 2;
    p.threshold   = 0.6;
    p.ignore_gaps = false;
                      // 
    string read = string("C");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true, false);
    
    const Path& path = aln.path();
    REQUIRE(path.mapping_size() == 2);
    CheckForGlobalAlignment(path, vid1);
    // correct path
    REQUIRE(path.mapping(0).position().node_id() == vid0->id());
    REQUIRE(path.mapping(1).position().node_id() == vid1->id());
    // correct edits
    REQUIRE(path.mapping(0).edit_size() == 1);
    REQUIRE(path.mapping(0).edit(0).from_length() == 1);
    REQUIRE(path.mapping(0).edit(0).to_length() == 0);
    REQUIRE(path.mapping(0).edit(0).sequence().empty());
    REQUIRE(path.mapping(1).edit_size() == 1);
    REQUIRE(path.mapping(1).edit(0).from_length() == 1);
    REQUIRE(path.mapping(1).edit(0).to_length() == 1);
    REQUIRE(path.mapping(1).edit(0).sequence().empty());

}

// 
// TODO test when read ends with big mismatch
// TODO test when the read is soft clipped (ends with insertion)
// TODO N-match
//

}}
