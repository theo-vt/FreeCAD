
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

// The goal of this file is to break down the constraint system in smaller constraints systems
// that can be solved in order. This is desirable because the numerical solvers have a complexity
// of the order of O(n^3), so breaking systems down can have a massive impact on performance

#ifndef PLANEGCS_DECOMPOSITION_H
#define PLANEGCS_DECOMPOSITION_H

#include "Util.h"
#include "Constraints.h"
#include <vector>

namespace GCS
{

struct SubsystemDescription
{
    std::vector<int> unknowns;   // Indices of unknowns in the GCS
    std::vector<int> equations;  // Indices of equations in the GCS

    bool empty() const;

    std::vector<Constraint*> makeClist(const std::vector<Constraint*>& constraints,
                                       const std::set<Constraint*>& reducedConstraints) const;
    std::vector<double*> makePlist(const std::vector<double*>& parameters) const;
};

struct SystemDecomposition
{
public:
    std::vector<SubsystemDescription> wellConstrained;
    SubsystemDescription overConstrained;
    SubsystemDescription underConstrained;
    std::vector<int> parameterComponent;

public:
    SystemDecomposition(const VEC_pD& parameters,
                        const std::vector<Constraint*>& constraints,
                        const std::map<Constraint*, VEC_pD>& c2p,
                        const MAP_pD_I& pIndex);

    size_t size() const;

    std::vector<std::vector<Constraint*>>
    makeClists(const std::vector<Constraint*>& constraints,
               const std::set<Constraint*>& reducedConstraints) const;
    std::vector<std::vector<double*>> makePlists(const std::vector<double*>& parameters) const;
};

}  // namespace GCS

#endif
