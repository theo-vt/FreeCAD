
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

#include "Substitution.h"

#include "Constraints.h"
#include "Util.h"

#include <boost/parameter/parameters.hpp>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <optional>

namespace GCS
{

// Each parameter is represented as it's index in the gcs
// plist, a substitution blob is a set of parameters which have
// the same value
struct SubstitutionBlob
{
    std::vector<int> parameters;
    std::optional<double> constValue;

    void add(int param)
    {
        parameters.push_back(param);
    }
    void merge(const SubstitutionBlob& other)
    {
        parameters.reserve((parameters.size() + other.parameters.size()));
        for (auto param : other.parameters) {
            parameters.push_back(param);
        }
    }
};

struct SubstitutionFactory
{
    std::vector<SubstitutionBlob> substitutionBlobs;
    std::unordered_map<int, int>
        paramToBlob;  // Map between parameters (index) and substitution blob (index)

    bool addEqual(int a, int b)
    {
        auto foundA = paramToBlob.find(a);
        auto foundB = paramToBlob.find(b);

        // Neither parameter are in a blob, so create a new one
        if (foundA == paramToBlob.end() && foundB == paramToBlob.end()) {
            substitutionBlobs.push_back({});
            substitutionBlobs.back().add(a);
            substitutionBlobs.back().add(b);

            paramToBlob[a] = substitutionBlobs.size() - 1;
            paramToBlob[b] = substitutionBlobs.size() - 1;

            return true;
        }
        // Both parameters are in a blob so merge them
        if (foundA != paramToBlob.end() && foundB != paramToBlob.end()) {
            int blobA = foundA->second;
            int blobB = foundB->second;

            if (blobA == blobB) {
                return false;  // This substitution already exist, maybe the constraint is redundant
            }

            if (blobA > blobB) {
                std::swap(blobA, blobB);
            }

            substitutionBlobs[blobA].merge(substitutionBlobs[blobB]);

            reassignBlob(blobB, blobA);
            for (size_t i = blobB + 1; i < substitutionBlobs.size(); ++i) {
                reassignBlob(i, i - 1);
            }
            substitutionBlobs.erase(substitutionBlobs.begin() + blobB);

            return true;
        }

        // One of the two parameters is in a blob, so add the other to it
        if (foundB != paramToBlob.end()) {
            std::swap(foundA, foundB);
            std::swap(a, b);
        }
        substitutionBlobs[foundA->second].add(foundB->second);
        paramToBlob[b] = foundA->second;

        return true;
    }
    bool areEqual(int a, int b)
    {
        auto foundA = paramToBlob.find(a);
        auto foundB = paramToBlob.find(b);

        return foundA != paramToBlob.end() && foundB != paramToBlob.end()
            && foundA->second == foundB->second;
    }
    bool addConst(int param, double value)
    {
        auto foundParam = paramToBlob.find(param);

        // Create a new const substitution blob with a single parameter
        if (foundParam == paramToBlob.end()) {
            substitutionBlobs.push_back({});
            substitutionBlobs.back().add(param);
            substitutionBlobs.back().constValue = value;

            return true;
        }

        // make the substitution blob const

        // The blob is already const, not allowed
        if (substitutionBlobs[foundParam->second].constValue.has_value()) {
            return false;
        }

        substitutionBlobs[foundParam->second].constValue = value;
        return true;
    }
    bool isConst(int param)
    {
        auto foundParam = paramToBlob.find(param);

        if (foundParam != paramToBlob.end()) {
            return substitutionBlobs[foundParam->second].constValue.has_value();
        }
        return false;
    }
    void reassignBlob(size_t src, size_t dst)
    {
        for (auto param : substitutionBlobs[src].parameters) {
            paramToBlob[param] = dst;
        }
    }
};

enum class ConstraintSubstitutionAttempt
{
    Yes,     // Substitution is possible
    No,      // Substitution is not possible
    Maybe,   // Substitution could be possible if some parameters are right
    Unknown  // Not checked yet
};

ConstraintSubstitutionAttempt substitutionForEqual(Constraint* constr,
                                                   const UMAP_pD_I& paramToIndex,
                                                   SubstitutionFactory& factory)
{
    const auto foundP1 = paramToIndex.find(constr->params()[0]);
    const auto foundP2 = paramToIndex.find(constr->params()[1]);

    if (foundP1 == paramToIndex.end() && foundP2 != paramToIndex.end()) {
        factory.addConst(foundP2->second, *constr->params()[0]);
        return ConstraintSubstitutionAttempt::Yes;
    }
    if (foundP1 == paramToIndex.end() && foundP2 != paramToIndex.end()) {
        factory.addConst(foundP2->second, *constr->params()[0]);
        return ConstraintSubstitutionAttempt::Yes;
    }

    if (foundP1 == paramToIndex.end() && foundP2 == paramToIndex.end()) {
        return ConstraintSubstitutionAttempt::No;
    }

    if (factory.addEqual(foundP1->second, foundP2->second)) {
        return ConstraintSubstitutionAttempt::Yes;
    }
    return ConstraintSubstitutionAttempt::No;
}
ConstraintSubstitutionAttempt substitutionForPerpendicular(Constraint* constr,
                                                           const UMAP_pD_I& paramToIndex,
                                                           SubstitutionFactory& factory)
{
    const auto l1x1 = paramToIndex.find(constr->params()[0]);
    const auto l1y1 = paramToIndex.find(constr->params()[1]);
    const auto l1x2 = paramToIndex.find(constr->params()[2]);
    const auto l1y2 = paramToIndex.find(constr->params()[3]);
    const auto l2x1 = paramToIndex.find(constr->params()[4]);
    const auto l2y1 = paramToIndex.find(constr->params()[5]);
    const auto l2x2 = paramToIndex.find(constr->params()[6]);
    const auto l2y2 = paramToIndex.find(constr->params()[7]);

    if (l1x1 == paramToIndex.end() || l1y1 == paramToIndex.end() || l1x2 == paramToIndex.end()
        || l1y2 == paramToIndex.end() || l2x1 == paramToIndex.end() || l2y1 == paramToIndex.end()
        || l2x2 == paramToIndex.end() || l2y2 == paramToIndex.end()) {
        return ConstraintSubstitutionAttempt::No;
    }

    // First line is vertical, second line should be horizontal
    if (factory.areEqual(l1x1->second, l1x2->second)) {
        factory.addEqual(l2y1->second, l2y2->second);
        return ConstraintSubstitutionAttempt::Yes;
    }

    // First line is horizontal, second line should be vertical
    if (factory.areEqual(l1y1->second, l1y2->second)) {
        factory.addEqual(l2x1->second, l2x2->second);
        return ConstraintSubstitutionAttempt::Yes;
    }

    // Second line is vertical, first line should be horizontal
    if (factory.areEqual(l2x1->second, l2x2->second)) {
        factory.addEqual(l1y1->second, l1y2->second);
        return ConstraintSubstitutionAttempt::Yes;
    }

    // Second line is horizontal, fist line should be vertical
    if (factory.areEqual(l2y1->second, l2y2->second)) {
        factory.addEqual(l1x1->second, l1x2->second);
        return ConstraintSubstitutionAttempt::Yes;
    }
    return ConstraintSubstitutionAttempt::Maybe;
}


Substitution::Substitution(const std::vector<double*>& initialParameters,
                           const std::vector<Constraint*>& initialConstraints,
                           const UMAP_pD_I& paramToIndex)
{
    SubstitutionFactory factory;

    std::vector<ConstraintSubstitutionAttempt> attempts(initialConstraints.size(),
                                                        ConstraintSubstitutionAttempt::Unknown);
    // bool done = false;

    for (size_t i = 0; i < initialConstraints.size(); ++i) {
        auto constr = initialConstraints[i];

        if (constr->getTypeId() != Equal) {
            continue;
        }

        attempts[i] = substitutionForEqual(constr, paramToIndex, factory);
    }

    // while (!done) {
    //     done = true; // Done until proven false by having a successful substituttion
    //     for (size_t i = 0; i < initialConstraints.size(); ++i) {
    //         auto constr = initialConstraints[i];

    //         // No use trying again
    //         if (attempts[i] == ConstraintSubstitutionAttempt::No || attempts[i] ==
    //         ConstraintSubstitutionAttempt::Yes) {
    //             continue;
    //         }

    //         // No substitution for temporary constraints
    //         if (constr->getTag() < 0) {
    //             continue;
    //         }

    //         ConstraintSubstitutionAttempt attempt = ConstraintSubstitutionAttempt::No;
    //         switch (constr->getTypeId()) {
    //         case Perpendicular:
    //             attempt = substitutionForEqual(constr, paramToIndex, factory);
    //         }
    //         attempts[i] = attempt;
    //         if (attempt == ConstraintSubstitutionAttempt::Yes) {
    //             done = false;
    //         }
    //     }
    // }

    parameters.reserve(factory.substitutionBlobs.size());
    for (const auto& blob : factory.substitutionBlobs) {
        if (!blob.constValue.has_value()) {
            parameters.push_back(initialParameters[blob.parameters[0]]);
            reducedParameter.push_back({});

            parameterToIndex[parameters.back()] = parameters.size() - 1;
            for (size_t i = 1; i < blob.parameters.size(); ++i) {
                double* param = initialParameters[blob.parameters[i]];
                reducedParameter.back().push_back(param);
                parameterToIndex[param] = parameters.size() - 1;
            }
        }
        else {
            constParams.emplace_back(VEC_pD {}, *blob.constValue);
            for (auto blobParam : blob.parameters) {
                double* param = initialParameters[blobParam];
                constParams.back().first.push_back(param);
                parameterToIndex[param] = -1;
            }
        }
    }


    for (size_t i = 0; i < initialParameters.size(); ++i) {
        auto foundBlob = factory.paramToBlob.find(i);

        if (foundBlob == factory.paramToBlob.end()) {
            parameters.push_back(initialParameters[i]);
            parameterToIndex[initialParameters[i]] = parameters.size() - 1;
        }
    }
    parameterVals.reserve(parameters.size());
    for (auto parameter : parameters) {
        parameterVals.push_back(*parameter);
    }
    for (auto paramAndIndex : parameterToIndex) {
        if (paramAndIndex.second > 0) {
            parameterToVal[paramAndIndex.first] = &parameterVals[paramAndIndex.second];
        }
    }

    for (size_t i = 0; i < attempts.size(); ++i) {
        if (attempts[i] == ConstraintSubstitutionAttempt::Yes) {
            continue;  // This constraint was substituted, no need to solve for it
        }
        constraints.push_back(initialConstraints[i]);
        initialConstraints[i]->fillParamIndices(parameterToIndex);
    }
    std::cerr << "#Constraints " << constraints.size() << " (from " << initialConstraints.size()
              << ")\n";
}

void Substitution::apply()
{
    for (size_t i = 0; i < parameters.size(); ++i) {
        *parameters[i] = parameterVals[i];
    }

    for (size_t r = 0; r < reducedParameter.size(); ++r) {
        for (size_t p = 0; p < reducedParameter[r].size(); ++p) {
            *reducedParameter[r][p] = *parameters[r];
        }
    }
}
void Substitution::applyConst()
{
    for (auto constSet : constParams) {
        for (auto param : constSet.first) {
            *param = constSet.second;
        }
    }
}

}  // namespace GCS
