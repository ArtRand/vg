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
            auto check_mapping_pairs = [&path] (int64_t mapping_idx, int64_t expected) {
                // has correct number of aligned pairs
                REQUIRE(path.mapping(mapping_idx).edit_size() == expected);
                // each has correct length (1) and is a match
                for (int64_t i = 0; i < path.mapping(mapping_idx).edit_size(); ++i) {
                    REQUIRE(path.mapping(mapping_idx).edit(i).from_length() == ALIGNED_PAIR_LENGTH);
                    REQUIRE(path.mapping(mapping_idx).edit(i).to_length() == ALIGNED_PAIR_LENGTH);
                    REQUIRE(path.mapping(mapping_idx).edit(i).sequence().empty());
                }
            };
            check_mapping_pairs(0, 3);
            check_mapping_pairs(1, 1);
            check_mapping_pairs(2, 3);
            REQUIRE(aln.identity() == 1.0);
        }

        TEST_CASE("Hmm aligner produces correct alignment when there is a mismatch", "[hmm][current]") {
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
            // check the matching region
            
        }
    }
}
