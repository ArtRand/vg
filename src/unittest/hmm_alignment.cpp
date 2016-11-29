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
#include "gssw_aligner.hpp"
#include "vg.hpp"
#include "path.hpp"
#include "catch.hpp"
#include "json2pb.h"

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
    hmm.Align(aln, nullptr, p, true);
    
    const Path& path = aln.path();
    
    // has correct number of mappings (nodes mapped to)
    REQUIRE(path.mapping_size() == 3);

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
    hmm.Align(aln, nullptr, p, true);
    
    const Path& path = aln.path();

    REQUIRE(path.mapping_size() == 3);
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
    hmm.Align(aln, nullptr, p, true);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
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
                        
                       // AAGTAGCAAACGTAACAT  [graph sequence]
                       // AAGT--CAAACGTAACAT  [read sequence]
    string read = string("AAGTCAAACGTAACAT");
    Alignment aln;
    aln.set_sequence(read);
    hmm.Align(aln, nullptr, p, true);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
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
TEST_CASE("Hmm Aligner produces correct alignment when there is an insertion", "[hmm][current]") {
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
    hmm.Align(aln, nullptr, p, true);
    
    const Path& path = aln.path();
    // maps to 3 nodes
    REQUIRE(path.mapping_size() == 3);
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
}}
