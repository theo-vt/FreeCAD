
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

// The goal of this file is to handle substitutions, that is constraints which can be implicitly
// solved Currently the only substitution mecanism is reduction, so if two parameters are set to be
// equal (e.g. if a line is horizonal, y1=y2), then only one parameter y1 will be forwarded to the
// solver and all instances of y2 will be replaced by y1 in every constraints


#ifndef PLANEGCS_SUBSTITUTION_H
#define PLANEGCS_SUBSTITUTION_H

#include "Util.h"
#include "Constraints.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace GCS
{

struct Substitution
{
    std::vector<Constraint*> constraints;  // Constraints that could not be reduced

    std::vector<double*> parameters;  // Parameters that could not be reduced
    std::vector<std::vector<double*>>
        reducedParameter;  // For every parameter in the parameters vector there is a vector of
                           // equal valued reduced parameters

    std::unordered_map<double*, int> parameterToIndex;  // For every parameter in GCS, an index to
                                                        // the parameter in the substitution

    Substitution(const std::vector<double*>& initialParameters,
                 const std::vector<Constraint*>& constraints,
                 const UMAP_pD_I& paramToIndex);
    Substitution() = default;
};

}  // namespace GCS

#endif
