/***************************************************************************
 *   Copyright (c) 2011 Konstantinos Poulios <logari81@gmail.com>          *
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

#ifdef _MSC_VER
#pragma warning(disable : 4251)
#endif

#include <iostream>
#include <iterator>

#include "SubSystem.h"


namespace GCS
{

// SubSystem
SubSystem::SubSystem(const std::vector<Constraint*>& clist_, const VEC_pD& params)
    : clist(clist_)
    , plist(params)
{
    UMAP_pD_pD dummymap;
    initialize(dummymap);
}

SubSystem::SubSystem(const std::vector<Constraint*>& clist_,
                     const VEC_pD& params,
                     const UMAP_pD_pD& reductionmap)
    : clist(clist_)
    , plist(params)
{
    initialize(reductionmap);
}

SubSystem::~SubSystem()
{}

void SubSystem::initialize(const UMAP_pD_pD& reductionmap)
{
    csize = static_cast<int>(clist.size());
    psize = static_cast<int>(plist.size());


    pvals.resize(psize);

    for (int j = 0; j < psize; j++) {
        pvals[j] = *plist[j];
        pmap[plist[j]] = &pvals[j];
    }

    p2c.clear();
    for (const auto constr : clist) {
        for (const auto p : constr->origParams()) {
            p2c[p].push_back(constr);
        }
    }
}

void SubSystem::redirectParams()
{
    // copying values to pvals
    for (auto p : pmap) {
        *(p.second) = *(p.first);
    }

    // redirect constraints to point to pvals
    for (auto constr : clist) {
        constr->redirectParams(pmap);
    }
}

void SubSystem::revertParams()
{
    for (auto constr : clist) {
        constr->revertParams();
    }
}

void SubSystem::getParamMap(UMAP_pD_pD& pmapOut)
{
    pmapOut = pmap;
}

void SubSystem::getParamList(VEC_pD& plistOut)
{
    plistOut = plist;
}

void SubSystem::getParams(VEC_pD& params, Eigen::VectorXd& xOut)
{
    if (xOut.size() != int(params.size())) {
        xOut.setZero(params.size());
    }

    for (int j = 0; j < int(params.size()); j++) {
        auto pmapfind = pmap.find(params[j]);
        if (pmapfind != pmap.end()) {
            xOut[j] = *(pmapfind->second);
        }
    }
}

void SubSystem::getParams(Eigen::VectorXd& xOut)
{
    if (xOut.size() != psize) {
        xOut.setZero(psize);
    }

    for (int i = 0; i < psize; i++) {
        xOut[i] = pvals[i];
    }
}

void SubSystem::setParams(VEC_pD& params, Eigen::VectorXd& xIn)
{
    assert(xIn.size() == int(params.size()));
    for (int j = 0; j < int(params.size()); j++) {
        auto pmapfind = pmap.find(params[j]);
        if (pmapfind != pmap.end()) {
            *(pmapfind->second) = xIn[j];
        }
    }
}

void SubSystem::setParams(Eigen::VectorXd& xIn)
{
    assert(xIn.size() == psize);
    for (int i = 0; i < psize; i++) {
        pvals[i] = xIn[i];
    }
}

void SubSystem::getConstraintList(std::vector<Constraint*>& clist_)
{
    clist_ = clist;
}

double SubSystem::error()
{
    double err = 0.;
    for (std::vector<Constraint*>::const_iterator constr = clist.begin(); constr != clist.end();
         ++constr) {
        double tmp = (*constr)->error();
        err += tmp * tmp;
    }
    err *= 0.5;
    return err;
}

void SubSystem::calcResidual(Eigen::VectorXd& r)
{
    assert(r.size() == csize);

    int i = 0;
    for (std::vector<Constraint*>::const_iterator constr = clist.begin(); constr != clist.end();
         ++constr, i++) {
        r[i] = (*constr)->error();
    }
}

void SubSystem::calcResidual(Eigen::VectorXd& r, double& err)
{
    assert(r.size() == csize);

    int i = 0;
    err = 0.;
    for (std::vector<Constraint*>::const_iterator constr = clist.begin(); constr != clist.end();
         ++constr, i++) {
        r[i] = (*constr)->error();
        err += r[i] * r[i];
    }
    err *= 0.5;
}

void SubSystem::calcJacobi(VEC_pD& params, Eigen::MatrixXd& jacobi)
{
    jacobi.setZero(csize, params.size());
    for (int j = 0; j < int(params.size()); j++) {
        const auto pmapfind = pmap.find(params[j]);
        if (pmapfind != pmap.end()) {
            for (int i = 0; i < csize; i++) {
                jacobi(i, j) = clist[i]->grad(pmapfind->second);
            }
        }
    }
}

void SubSystem::calcJacobi(Eigen::MatrixXd& jacobi)
{
    calcJacobi(plist, jacobi);
}

void SubSystem::calcGrad(VEC_pD& params, Eigen::VectorXd& grad)
{
    assert(grad.size() == int(params.size()));

    grad.setZero();
    for (int j = 0; j < int(params.size()); j++) {
        const auto pmapfind = pmap.find(params[j]);
        if (pmapfind != pmap.end()) {
            std::vector<Constraint*> constrs = p2c[pmapfind->second];
            for (const auto constr : constrs) {
                grad[j] += constr->error() * constr->grad(pmapfind->second);
            }
        }
    }
}

void SubSystem::calcGrad(Eigen::VectorXd& grad)
{
    calcGrad(plist, grad);
}

double SubSystem::maxStep(VEC_pD& params, Eigen::VectorXd& xdir)
{
    assert(xdir.size() == int(params.size()));

    MAP_pD_D dir;
    for (int j = 0; j < int(params.size()); j++) {
        const auto pmapfind = pmap.find(params[j]);
        if (pmapfind != pmap.end()) {
            dir[pmapfind->second] = xdir[j];
        }
    }

    double alpha = 1e10;
    for (auto constr : clist) {
        alpha = constr->maxStep(dir, alpha);
    }

    return alpha;
}

double SubSystem::maxStep(Eigen::VectorXd& xdir)
{
    return maxStep(plist, xdir);
}

void SubSystem::applySolution()
{
    for (const auto p2p : pmap) {
        *(p2p.first) = *(p2p.second);
    }
}

void SubSystem::analyse(Eigen::MatrixXd& /*J*/, Eigen::MatrixXd& /*ker*/, Eigen::MatrixXd& /*img*/)
{}

void SubSystem::report()
{}

void SubSystem::printResidual()
{
    Eigen::VectorXd r(csize);
    double err = 0.;
    for (size_t i = 0; i < clist.size(); ++i) {
        Constraint* constr = clist[i];
        r[i] = constr->error();
        err += r[i] * r[i];
    }
    err *= 0.5;
    std::cout << "Residual r = " << r.transpose() << std::endl;
    std::cout << "Residual err = " << err << std::endl;
}


}  // namespace GCS
