
/***************************************************************************
 *   Copyright (c) 2025 Théo Veilleux-Trinh <theo.veilleux.trinh@proton.me>*
 *                                                                         *
 *   This file is part of the FreeCAD CAx development system.              *
 *                                                                         *
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Library General Public           *
 *   License as published by the Free Software Foundation; either          *
 *   version 2 of the License, or (at your option) any later version.      *
 *                                                                         *
 *   This library  is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Library General Public License for more details.                  *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this library; see the file COPYING.LIB. If not,    *
 *   write to the Free Software Foundation, Inc., 59 Temple Place,         *
 *   Suite 330, Boston, MA  02111-1307, USA                                *
 *                                                                         *
 ***************************************************************************/

#include "Decomposition.h"
#include "Mod/Sketcher/App/planegcs/Substitution.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

#include <unordered_set>
#include <iostream>

#define USE_STRONGLY_CONNECTED_COMPONENTS

namespace GCS
{
using namespace boost;

#ifdef USE_STRONGLY_CONNECTED_COMPONENTS

using graph_t = adjacency_list<vecS, vecS, bidirectionalS>;
using vertex_t = graph_traits<graph_t>::vertex_descriptor;
using edge_t = graph_traits<graph_t>::edge_descriptor;

// Helper class to find the maximum matching of
// the constraint-unknowns system
class BipartiteGraph
{
public:
    struct BPMOutput
    {
        std::vector<std::pair<int, int>> matching;
        std::vector<int> unsaturatedA;
        std::vector<int> unsaturatedB;
    };

private:
    size_t A {0};
    size_t B {0};
    std::vector<std::unordered_set<int>> B_to_A {};

public:
    // Edges go from A->B
    BipartiteGraph(int szA, int szB, std::vector<std::pair<int, int>> edges)
        : A(szA)
        , B(szB)
        , B_to_A(szB)
    {
        for (auto edge : edges) {
            std::cerr << "Edge " << edge.first << " :: " << edge.second << "\n";
            B_to_A[edge.second].insert(edge.first);
        }
    }

    // Returns matched pairs from A->B
    BPMOutput maxBPM()
    {
        std::vector<int> matched(B, -1);
        std::unordered_set<int> unsaturatedA;

        for (size_t a = 0; a < A; ++a) {
            std::vector<int> seen(B, 0);
            bpm(a, matched, seen);
            unsaturatedA.insert(a);
        }

        BPMOutput out;
        for (size_t b = 0; b < B; ++b) {
            if (matched[b] != -1) {
                out.matching.emplace_back(matched[b], b);
                unsaturatedA.erase(matched[b]);  // This A vertex is saturated
            }
            else {
                out.unsaturatedB.push_back(b);
            }
        }
        out.unsaturatedA.assign(unsaturatedA.begin(), unsaturatedA.end());
        return out;
    }

private:
    bool bpm(int A_index, std::vector<int>& matched, std::vector<int>& seen)
    {

        for (size_t b = 0; b < B; ++b) {
            if (B_to_A[b].count(A_index) == 1 && !seen[b]) {
                seen[b] = true;

                if (matched[b] < 0 || bpm(matched[b], matched, seen)) {
                    matched[b] = A_index;
                    return true;
                }
            }
        }
        return false;
    }
};


// Visitor class for the breadth first search to build the
// under and over constrained subsystems
class AggregationVisitor: public default_bfs_visitor
{
public:
    std::set<vertex_t>* vertices {nullptr};
    AggregationVisitor(std::set<vertex_t>* vertices_)
        : vertices(vertices_)
    {}
    template<typename Vertex, typename Graph>
    void discover_vertex(Vertex vert, const Graph& /*g*/) const
    {
        // std::cout<<"vert: "<<vert<<"\n";
        vertices->insert(vert);
        // vertices.insert(vert);
    }
};


SystemDecomposition::SystemDecomposition(size_t numParams,
                                         const std::vector<Constraint*>& constraints,
                                         const std::map<Constraint*, VEC_pD>& c2p,
                                         const UMAP_pD_I& pIndex)
    : parameterComponent(numParams)
{
    // This constructor implements algorithms detailed in <source>
    // in the source paper, the underconstrained graph is called G3
    // and the overconstrained graph is called G2
    //
    // Constraints are also in the group Y (y1, y2, y3..) and parameters
    // are in the group X (x1, x2, x3..)

    // The graph G contains all the the constraints and all the parameters
    graph_t G(numParams + constraints.size());

    // std::cerr << "System size: " << numParams + constraints.size() << "\n";

    std::vector<std::pair<int, int>> bipartiteEdges;
    for (size_t c = 0; c < constraints.size(); ++c) {
        Constraint* constr = constraints[c];

        for (const auto param : constr->paramsIndex()) {
            if (param == -1) {
                continue;
            }
            // The constraints are placed after the parameters in the constraints list
            add_edge(c + numParams, param, G);
            bipartiteEdges.emplace_back(c, param);
        }
    }
    BipartiteGraph bpGraph(constraints.size(), numParams, bipartiteEdges);
    BipartiteGraph::BPMOutput matchingOut = bpGraph.maxBPM();

    // Edges that are part of the maximum matching of the bipartite graph get a back edge
    for (auto match : matchingOut.matching) {
        // std::cerr << "Match: " << match.first << " :: " << match.second << "\n";
        // a match goes from equation -> unknown but both sets start at 0 so we have to offset the
        // equation index
        add_edge(match.second, match.first + numParams, G);
    }

    // Build the over constrained subsystem description (G2 in the paper)
    std::set<vertex_t> OCVertices;
    AggregationVisitor OCVisitor(&OCVertices);

    for (auto unsaturated : matchingOut.unsaturatedA) {
        breadth_first_search(G, unsaturated + numParams, visitor(OCVisitor));
    }
    // std::cerr << "OC vertices: \n";
    for (auto vert : OCVertices) {
        // Remove all edges from or to this vertex so that it is not tangled in well constrained
        // strongly connected components
        clear_vertex(vert, G);

        if (vert >= numParams) {
            // std::cerr << "[eq]" << vert - numParams << "\n";
            overConstrained.equations.push_back(vert - numParams);
        }
        else {
            // std::cerr << "[va]" << vert << "\n";
            overConstrained.unknowns.push_back(vert);
        }
    }

    // Build the under constrained subsystem description (G3 in the paper)
    std::set<vertex_t> UCVertices;
    AggregationVisitor UCVisitor(&UCVertices);

    for (auto unsaturated : matchingOut.unsaturatedB) {
        reverse_graph<graph_t> rg(G);
        breadth_first_search(rg, unsaturated, visitor(UCVisitor));
    }

    // std::cerr << "UC vertices: \n";
    for (auto vert : UCVertices) {
        // Remove all edges from or to this vertex so that it is not tangled in well constrained
        // strongly connected components
        clear_vertex(vert, G);
        if (vert >= numParams) {
            // std::cerr << "[eq]" << vert - numParams << "\n";
            underConstrained.equations.push_back(vert - numParams);
        }
        else {
            // std::cerr << "[va]" << vert << "\n";
            underConstrained.unknowns.push_back(vert);
        }
    }

    // Define the vertex index property map using boost::property_map and
    // boost::associative_property_map
    std::map<vertex_t, size_t> index_map;
    int id = 0;
    for (auto v : make_iterator_range(vertices(G))) {
        index_map.emplace(v, id++);
    }

    std::vector<size_t> components(num_vertices(G));
    auto vi_map = make_assoc_property_map(index_map);
    int num_scc = strong_components(G, make_iterator_property_map(components.begin(), vi_map));

    // std::cerr << "WC: \n";
    std::vector<std::vector<int>> vertexPerComponent;
    for (auto v : make_iterator_range(vertices(G))) {
        size_t ind = vi_map[v];

        if (UCVertices.count(ind) == 1 || OCVertices.count(ind) == 1) {
            continue;
        }

        size_t comp = components[vi_map[v]];
        if (comp >= vertexPerComponent.size()) {
            vertexPerComponent.resize(comp + 1);
        }
        vertexPerComponent[comp].push_back(ind);
    }
    for (auto vertices : vertexPerComponent) {
        if (vertices.empty()) {
            continue;
        }
        wellConstrained.push_back({});
        for (auto vertex : vertices) {
            if (vertex >= numParams) {
                // std::cerr << "[eq]" << ind - numParams << " :: " << comp << "\n";
                wellConstrained.back().equations.push_back(vertex - numParams);
            }
            else {
                // std::cerr << "[va]" << ind << " :: " << comp << "\n";
                wellConstrained.back().unknowns.push_back(vertex);
                parameterComponent[vertex] = wellConstrained.size() - 1;
            }
        }
    }

    for (auto pind : underConstrained.unknowns) {
        parameterComponent[pind] = wellConstrained.size();
    }
    for (auto pind : overConstrained.unknowns) {
        parameterComponent[pind] = wellConstrained.size() + 1;
    }
}
#else

SystemDecomposition::SystemDecomposition(size_t numParams,
                                         const std::vector<Constraint*>& constraints,
                                         const std::map<Constraint*, VEC_pD>& c2p,
                                         const UMAP_pD_I& pIndex)
{
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

    // partitioning into decoupled components
    Graph g;
    for (int i = 0; i < int(numParams + constraints.size()); i++) {
        boost::add_vertex(g);
    }

    int cvtid = int(numParams);
    for (const auto constr : constraints) {
        auto cparams = c2p.find(constr);

        if (cparams == c2p.end()) {
            continue;
        }

        for (const auto param : cparams->second) {
            auto it = pIndex.find(param);
            if (it != pIndex.end()) {
                boost::add_edge(cvtid, it->second, g);
            }
        }
        ++cvtid;
    }

    VEC_I components(boost::num_vertices(g));
    int componentsSize = 0;
    if (!components.empty()) {
        componentsSize = boost::connected_components(g, &components[0]);
    }

    wellConstrained.resize(componentsSize);

    for (size_t i = 0; i < numParams; ++i) {
        wellConstrained[components[i]].unknowns.push_back(i);
    }
    for (size_t i = 0; i < constraints.size(); ++i) {
        wellConstrained[components[i + numParams]].equations.push_back(i);
    }
}

#endif

size_t SystemDecomposition::size() const
{
    return wellConstrained.size() + (!overConstrained.empty()) + (!underConstrained.empty());
}

std::vector<SubsystemPrecursor>
SystemDecomposition::makeSubsystemPrecursors(const Substitution& substitution) const
{
    std::vector<SubsystemPrecursor> out;
    out.reserve(size());

    for (const auto& subsyst : wellConstrained) {
        out.push_back(subsyst.makeSubsystemPrecursor(substitution));
    }
    if (!underConstrained.empty()) {
        out.push_back(underConstrained.makeSubsystemPrecursor(substitution));
    }
    if (!overConstrained.empty()) {
        out.push_back(overConstrained.makeSubsystemPrecursor(substitution));
    }
    return out;
}

bool SubsystemDescription::empty() const
{
    return equations.empty() || unknowns.empty();
}
SubsystemPrecursor
SubsystemDescription::makeSubsystemPrecursor(const Substitution& substitution) const
{
    SubsystemPrecursor out;
    out.constraints.reserve(equations.size());
    out.parameters.reserve(unknowns.size());

    for (auto unknown : unknowns) {
        out.parameters.push_back(substitution.parameters[unknown]);
    }

    for (auto eq : equations) {
        Constraint* constr = substitution.constraints[eq];
        out.constraints.push_back(constr);

        auto origParams = constr->origParams();
        auto paramSubstIndices = constr->paramsIndex();
        for (size_t i = 0; i < origParams.size(); ++i) {
            if (paramSubstIndices[i] != -1) {
                out.reductionMap[origParams[i]] = substitution.parameters[paramSubstIndices[i]];
            }
        }
    }
    return out;
}

}  // namespace GCS
