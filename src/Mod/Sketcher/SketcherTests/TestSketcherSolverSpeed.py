# **************************************************************************
#   Copyright (c) 2025 Théo Veilleux-Trinh <theo.veilleux.trinh@proton.me>*
#                                                                         *
#   This file is part of the FreeCAD CAx development system.              *
#                                                                         *
#   This program is free software; you can redistribute it and/or modify  *
#   it under the terms of the GNU Lesser General Public License (LGPL)    *
#   as published by the Free Software Foundation; either version 2 of     *
#   the License, or (at your option) any later version.                   *
#   for detail see the LICENCE text file.                                 *
#                                                                         *
#   FreeCAD is distributed in the hope that it will be useful,            *
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#   GNU Library General Public License for more details.                  *
#                                                                         *
#   You should have received a copy of the GNU Library General Public     *
#   License along with FreeCAD; if not, write to the Free Software        *
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
#   USA                                                                   *
# **************************************************************************


import os, tempfile, unittest
import FreeCAD, Part, Sketcher
from Part import Precision

App = FreeCAD

xy_normal = FreeCAD.Vector(0, 0, 1)


class TestSketcherSolverSpeed(unittest.TestCase):

    def speedTest(sketchBuilder, numIt):
        timer = Timer()
        for _ in range(num_it):
            sketch = sketch_builder()
            timer.start()
            solve()
            timer.stop()

        return timer.elapsed()

    def assertSuccessfulSolve(self, sketch, msg=None):
        status = sketch.solve()
        # TODO: can we get the solver's messages somehow to improve the message?
        self.assertTrue(status == 0, msg=msg or "solver didn't converge")
