# SPDX-License-Identifier: LGPL-2.1-or-later

# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2011 Yorik van Havre <yorik@uncreated.net>              *
# *                                                                         *
# *   This file is part of FreeCAD.                                         *
# *                                                                         *
# *   FreeCAD is free software: you can redistribute it and/or modify it    *
# *   under the terms of the GNU Lesser General Public License as           *
# *   published by the Free Software Foundation, either version 2.1 of the  *
# *   License, or (at your option) any later version.                       *
# *                                                                         *
# *   FreeCAD is distributed in the hope that it will be useful, but        *
# *   WITHOUT ANY WARRANTY; without even the implied warranty of            *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
# *   Lesser General Public License for more details.                       *
# *                                                                         *
# *   You should have received a copy of the GNU Lesser General Public      *
# *   License along with FreeCAD. If not, see                               *
# *   <https://www.gnu.org/licenses/>.                                      *
# *                                                                         *
# ***************************************************************************

__title__  = "FreeCAD Panel"
__author__ = "Yorik van Havre"
__url__    = "https://www.freecad.org"

## @package ArchPanel
#  \ingroup ARCH
#  \brief The Panel object and tools
#
#  This module provides tools to build Panel objects.
#  Panels consist of a closed shape that gets extruded to
#  produce a flat object.

import math

import FreeCAD
import ArchCommands
import ArchComponent
import Draft
import DraftVecUtils
import Part

from FreeCAD import Vector
from draftutils import params

if FreeCAD.GuiUp:
    from PySide import QtCore, QtGui
    from PySide.QtCore import QT_TRANSLATE_NOOP
    import FreeCADGui
    from draftutils.translate import translate
else:
    # \cond
    def translate(ctxt,txt):
        return txt
    def QT_TRANSLATE_NOOP(ctxt,txt):
        return txt
    # \endcond


class _Panel(ArchComponent.Component):

    "The Panel object"

    def __init__(self,obj):

        ArchComponent.Component.__init__(self,obj)
        self.Type = "Panel"
        self.setProperties(obj)
        obj.IfcType = "Plate"

    def setProperties(self,obj):

        pl = obj.PropertiesList
        if not "Length" in pl:
            obj.addProperty("App::PropertyLength","Length","Panel",   QT_TRANSLATE_NOOP("App::Property","The length of this element, if not based on a profile"), locked=True)
        if not "Width" in pl:
            obj.addProperty("App::PropertyLength","Width","Panel",    QT_TRANSLATE_NOOP("App::Property","The width of this element, if not based on a profile"), locked=True)
        if not "Thickness" in pl:
            obj.addProperty("App::PropertyLength","Thickness","Panel",QT_TRANSLATE_NOOP("App::Property","The thickness or extrusion depth of this element"), locked=True)
        if not "Sheets" in pl:
            obj.addProperty("App::PropertyInteger","Sheets","Panel",  QT_TRANSLATE_NOOP("App::Property","The number of sheets to use"), locked=True)
            obj.Sheets = 1
        if not "Offset" in pl:
            obj.addProperty("App::PropertyDistance","Offset","Panel",   QT_TRANSLATE_NOOP("App::Property","The offset between this panel and its baseline"), locked=True)
        if not "WaveLength" in pl:
            obj.addProperty("App::PropertyLength","WaveLength","Panel", QT_TRANSLATE_NOOP("App::Property","The length of waves for corrugated elements"), locked=True)
        if not "WaveHeight" in pl:
            obj.addProperty("App::PropertyLength","WaveHeight","Panel", QT_TRANSLATE_NOOP("App::Property","The height of waves for corrugated elements"), locked=True)
        if not "WaveOffset" in pl:
            obj.addProperty("App::PropertyDistance","WaveOffset","Panel", QT_TRANSLATE_NOOP("App::Property","The horizontal offset of waves for corrugated elements"), locked=True)
        if not "WaveDirection" in pl:
            obj.addProperty("App::PropertyAngle","WaveDirection","Panel", QT_TRANSLATE_NOOP("App::Property","The direction of waves for corrugated elements"), locked=True)
        if not "WaveType" in pl:
            obj.addProperty("App::PropertyEnumeration","WaveType","Panel", QT_TRANSLATE_NOOP("App::Property","The type of waves for corrugated elements"), locked=True)
            obj.WaveType = ["Curved","Trapezoidal","Spikes"]
        if not "WaveBottom" in pl:
            obj.addProperty("App::PropertyBool","WaveBottom","Panel", QT_TRANSLATE_NOOP("App::Property","If the wave also affects the bottom side or not"), locked=True)
        if not "Area" in pl:
            obj.addProperty("App::PropertyArea","Area","Panel",       QT_TRANSLATE_NOOP("App::Property","The area of this panel"), locked=True)
        if not "FaceMaker" in pl:
            obj.addProperty("App::PropertyEnumeration","FaceMaker","Panel",QT_TRANSLATE_NOOP("App::Property","The facemaker type to use to build the profile of this object"), locked=True)
            obj.FaceMaker = ["None","Simple","Cheese","Bullseye"]
        if not "Normal" in pl:
            obj.addProperty("App::PropertyVector","Normal","Panel",QT_TRANSLATE_NOOP("App::Property","The normal extrusion direction of this object (keep (0,0,0) for automatic normal)"), locked=True)
        obj.setEditorMode("VerticalArea",2)
        obj.setEditorMode("HorizontalArea",2)

    def onDocumentRestored(self,obj):

        ArchComponent.Component.onDocumentRestored(self,obj)
        self.setProperties(obj)

    def loads(self,state):

        self.Type = "Panel"

    def execute(self,obj):

        "creates the panel shape"

        if self.clone(obj):
            return

        layers = []
        length = 0
        width = 0
        thickness = 0

        # base tests
        if obj.Base:
            if hasattr(obj.Base,'Shape'):
                if obj.Base.Shape.isNull():
                    return
            elif obj.Base.isDerivedFrom("Mesh::Feature"):
                if not obj.Base.Mesh.isSolid():
                    return
        else:
            if obj.Length.Value:
                length = obj.Length.Value
            else:
                return
            if obj.Width.Value:
                width = obj.Width.Value
            else:
                return
        if obj.Thickness.Value:
            thickness = obj.Thickness.Value
        else:
            if not obj.Base:
                return
            elif hasattr(obj.Base,'Shape'):
                if not obj.Base.Shape.Solids:
                    return
        if hasattr(obj,"Material"):
            if obj.Material:
                if hasattr(obj.Material,"Materials"):
                    varwidth = 0
                    thicknesses = [t for t in obj.Material.Thicknesses if t >= 0]
                    restwidth = thickness - sum(thicknesses)
                    if restwidth > 0:
                        varwidth = [t for t in thicknesses if t == 0]
                        if varwidth:
                            varwidth = restwidth/len(varwidth)
                    for t in obj.Material.Thicknesses:
                        if t:
                            layers.append(t)
                        elif varwidth:
                            layers.append(varwidth)
        # creating base shape
        pl = obj.Placement
        base = None
        normal = None
        if hasattr(obj,"Normal"):
            if obj.Normal.Length > 0:
                normal = Vector(obj.Normal)
                normal.normalize()
                normal.multiply(thickness)
        baseprofile = None
        if obj.Base:
            base = obj.Base.Shape.copy()
            if not base.Solids:
               # p = FreeCAD.Placement(obj.Base.Placement)
                if base.Faces:
                    baseprofile = base
                    if not normal:
                        normal = baseprofile.Faces[0].normalAt(0,0).multiply(thickness)
                    if layers:
                        layeroffset = 0
                        shps = []
                        for l in layers:
                            if l >= 0:
                                n = Vector(normal).normalize().multiply(abs(l))
                                b = base.extrude(n)
                                if layeroffset:
                                    o = Vector(normal).normalize().multiply(layeroffset)
                                    b.translate(o)
                                shps.append(b)
                            layeroffset += abs(l)
                        base = Part.makeCompound(shps)
                    else:
                        base = base.extrude(normal)
                elif base.Wires:
                    fm = False
                    if hasattr(obj,"FaceMaker"):
                        if obj.FaceMaker != "None":
                            try:
                                baseprofile = Part.makeFace(base.Wires,"Part::FaceMaker"+str(obj.FaceMaker))
                                fm = True
                            except Exception:
                                FreeCAD.Console.PrintError(translate("Arch","Facemaker returned an error")+"\n")
                                return
                    if not fm:
                        closed = True
                        for w in base.Wires:
                            if not w.isClosed():
                                closed = False
                        if closed:
                            baseprofile = ArchCommands.makeFace(base.Wires)
                    if not normal:
                        normal = baseprofile.normalAt(0,0).multiply(thickness)
                    if layers:
                        layeroffset = 0
                        shps = []
                        for l in layers:
                            if l >= 0:
                                n = Vector(normal).normalize().multiply(abs(l))
                                b = baseprofile.extrude(n)
                                if layeroffset:
                                    o = Vector(normal).normalize().multiply(layeroffset)
                                    b.translate(o)
                                shps.append(b)
                            layeroffset += abs(l)
                        base = Part.makeCompound(shps)
                    else:
                        base = baseprofile.extrude(normal)
                elif obj.Base.isDerivedFrom("Mesh::Feature"):
                    if obj.Base.Mesh.isSolid():
                        if obj.Base.Mesh.countComponents() == 1:
                            sh = ArchCommands.getShapeFromMesh(obj.Base.Mesh)
                            if sh.isClosed() and sh.isValid() and sh.Solids:
                                base = sh
        else:
            if layers:
                shps = []
                layeroffset = 0
                for l in layers:
                    if l >= 0:
                        if normal:
                            n = Vector(normal).normalize().multiply(l)
                        else:
                            n = Vector(0,0,1).multiply(abs(l))
                        l2 = length/2 or 0.5
                        w2 = width/2 or 0.5
                        v1 = Vector(-l2,-w2,layeroffset)
                        v2 = Vector(l2,-w2,layeroffset)
                        v3 = Vector(l2,w2,layeroffset)
                        v4 = Vector(-l2,w2,layeroffset)
                        base = Part.makePolygon([v1,v2,v3,v4,v1])
                        baseprofile = Part.Face(base)
                        base = baseprofile.extrude(n)
                        shps.append(base)
                    layeroffset += abs(l)
                base = Part.makeCompound(shps)
            else:
                if not normal:
                    normal = Vector(0,0,1).multiply(thickness)
                l2 = length/2 or 0.5
                w2 = width/2 or 0.5
                v1 = Vector(-l2,-w2,0)
                v2 = Vector(l2,-w2,0)
                v3 = Vector(l2,w2,0)
                v4 = Vector(-l2,w2,0)
                base = Part.makePolygon([v1,v2,v3,v4,v1])
                baseprofile = Part.Face(base)
                base = baseprofile.extrude(normal)

        if hasattr(obj,"Area"):
            if baseprofile:
                obj.Area = baseprofile.Area

        if hasattr(obj,"WaveLength"):
            if baseprofile and obj.WaveLength.Value and obj.WaveHeight.Value:
                # corrugated element
                bb = baseprofile.BoundBox
                bb.enlarge(bb.DiagonalLength)
                downsegment = None
                if hasattr(obj,"WaveBottom"):
                    if not obj.WaveBottom:
                        if obj.WaveType == "Curved":
                            if obj.Thickness.Value > obj.WaveHeight.Value:
                                downsegment = obj.Thickness.Value
                            else:
                                downsegment = obj.WaveHeight.Value + obj.Thickness.Value
                        else:
                            downsegment = obj.Thickness.Value
                p1 = Vector(0,0,0)
                p5 = Vector(obj.WaveLength.Value*2,0,0)
                if obj.WaveType == "Curved":
                    p2 = Vector(obj.WaveLength.Value/2,0,obj.WaveHeight.Value)
                    p3 = Vector(obj.WaveLength.Value,0,0)
                    e1 = Part.Arc(p1,p2,p3).toShape()
                    p4 = Vector(obj.WaveLength.Value*1.5,0,-obj.WaveHeight.Value)
                    e2 = Part.Arc(p3,p4,p5).toShape()
                    upsegment = Part.Wire([e1,e2])
                    if not downsegment:
                        if obj.Thickness.Value < e1.Curve.Radius:
                            c3 = e1.Curve.copy()
                            c3.Radius = e1.Curve.Radius-obj.Thickness.Value
                            e3 = Part.Arc(c3,e1.FirstParameter,e1.LastParameter).toShape()
                            c4 = e2.Curve.copy()
                            c4.Radius = e2.Curve.Radius+obj.Thickness.Value
                            e4 = Part.Arc(c4,e2.FirstParameter,e2.LastParameter).toShape()
                            downsegment = Part.Wire([e3,e4])
                        else:
                            r = e2.Curve.Radius+obj.Thickness.Value
                            z = math.sqrt(r**2 - obj.WaveLength.Value**2)
                            p6 = e2.Curve.Center.add(Vector(-obj.WaveLength,0,-z))
                            p7 = e2.Curve.Center.add(Vector(0,0,-r))
                            p8 = e2.Curve.Center.add(Vector(obj.WaveLength,0,-z))
                            downsegment = Part.Arc(p6,p7,p8).toShape()

                elif obj.WaveType == "Trapezoidal":
                    p2 = Vector(obj.WaveLength.Value/4,0,obj.WaveHeight.Value)
                    p3 = Vector(obj.WaveLength.Value,0,obj.WaveHeight.Value)
                    p4 = Vector(obj.WaveLength.Value*1.25,0,0)
                    upsegment = Part.makePolygon([p1,p2,p3,p4,p5])
                    if not downsegment:
                        a = ((p1.sub(p2)).getAngle(p3.sub(p2)))/2
                        tx = obj.Thickness.Value/math.tan(a)
                        d1 = Vector(tx,0,-obj.Thickness.Value)
                        d2 = Vector(-tx,0,-obj.Thickness.Value)
                        p6 = p1.add(d1)
                        if tx >= p3.sub(p2).Length/2:
                            d3 = p2.sub(p1)
                            d3.normalize()
                            d3.multiply((0.625*obj.WaveLength.Value)/d3.x)
                            d4 = Vector(d3.x,0,-d3.z)
                            p7 = p6.add(d3)
                            p8 = p7.add(d4)
                            p9 = p5.add(d1)
                            downsegment = Part.makePolygon([p6,p7,p8,p9])
                        elif tx <= 0.625*obj.WaveLength.Value:
                            p7 = p2.add(d1)
                            p8 = p3.add(d2)
                            p9 = p4.add(d2)
                            p10 = p5.add(d1)
                            downsegment = Part.makePolygon([p6,p7,p8,p9,p10])
                        else:
                            downsegment = obj.Thickness.Value

                else: # spike
                    p2 = Vector(obj.WaveHeight.Value,0,obj.WaveHeight.Value)
                    p3 = Vector(obj.WaveHeight.Value*2,0,0)
                    upsegment = Part.makePolygon([p1,p2,p3,p5])
                    if not downsegment:
                        downsegment = obj.Thickness.Value

                upsegment.translate(Vector(bb.getPoint(0).x,bb.getPoint(0).y,bb.Center.z))
                if isinstance(downsegment,Part.Shape):
                    downsegment.translate(Vector(bb.getPoint(0).x,bb.getPoint(0).y,bb.Center.z))
                if hasattr(obj,"WaveOffset"):
                    if obj.WaveOffset.Value:
                        upsegment.translate(Vector(obj.WaveOffset.Value,0,0))
                        if isinstance(downsegment,Part.Shape):
                            downsegment.translate(Vector(obj.WaveOffset.Value,0,0))

                upedges = []
                downedges = []
                for i in range(int(bb.XLength/(obj.WaveLength.Value*2))):
                    w1 = upsegment.copy()
                    w1.translate(Vector(obj.WaveLength.Value*2*i,0,0))
                    upedges.extend(w1.Edges)
                    if isinstance(downsegment,Part.Shape):
                        w2 = downsegment.copy()
                        w2.translate(Vector(obj.WaveLength.Value*2*i,0,0))
                        downedges.extend(w2.Edges)
                upwire = Part.Wire(upedges)
                FreeCAD.upwire = upwire # REMOVE
                if isinstance(downsegment,Part.Shape):
                    downwire = Part.Wire(downedges)
                    FreeCAD.downwire = downwire # REMOVE
                    e1 = Part.LineSegment(upwire.Vertexes[0].Point,downwire.Vertexes[0].Point).toShape()
                    e2 = Part.LineSegment(upwire.Vertexes[-1].Point,downwire.Vertexes[-1].Point).toShape()
                    basewire = Part.Wire(upwire.Edges+[e1,e2]+downwire.Edges)
                else:
                    z = obj.Thickness.Value
                    if obj.WaveType == "Curved":
                        z += obj.WaveHeight.Value
                    p1 = upwire.Vertexes[0].Point
                    p2 = p1.add(Vector(0,0,-z))
                    p3 = Vector(upwire.Vertexes[-1].Point.x,upwire.Vertexes[-1].Point.y,p2.z)
                    p4 = upwire.Vertexes[-1].Point
                    w = Part.makePolygon([p1,p2,p3,p4])
                    basewire = Part.Wire(upwire.Edges+w.Edges)

                FreeCAD.basewire = basewire
                if not basewire.isClosed():
                    print("Error closing base wire - check FreeCAD.basewire")
                    return

                baseface = Part.Face(basewire)
                base = baseface.extrude(Vector(0,bb.YLength,0))
                rot = FreeCAD.Rotation(FreeCAD.Vector(0,0,1),normal)
                base.rotate(bb.Center,rot.Axis,math.degrees(rot.Angle))
                if obj.WaveDirection.Value:
                    base.rotate(bb.Center,normal,obj.WaveDirection.Value)
                n1 = normal.negative().normalize().multiply(obj.WaveHeight.Value*2)
                self.vol = baseprofile.copy()
                self.vol.translate(n1)
                self.vol = self.vol.extrude(n1.negative().multiply(2))
                base = self.vol.common(base)
                base = base.removeSplitter()
                if not base:
                    FreeCAD.Console.PrintError(translate("Arch","Error computing shape of")+" "+obj.Label+"\n")
                    return False

        if base and (obj.Sheets > 1) and normal and thickness:
            bases = [base]
            for i in range(1,obj.Sheets):
                n = FreeCAD.Vector(normal).normalize().multiply(i*thickness)
                b = base.copy()
                b.translate(n)
                bases.append(b)
            base = Part.makeCompound(bases)

        if base and normal and hasattr(obj,"Offset"):
            if obj.Offset.Value:
                v = DraftVecUtils.scaleTo(normal,obj.Offset.Value)
                base.translate(v)

        # process subshapes
        base = self.processSubShapes(obj,base,pl)

        # applying
        if base:
            if not base.isNull():
                if base.isValid() and base.Solids:
                    if len(base.Solids) == 1:
                        if base.Volume < 0:
                            base.reverse()
                        if base.Volume < 0:
                            FreeCAD.Console.PrintError(translate("Arch","Couldn't compute a shape"))
                            return
                        base = base.removeSplitter()
                    obj.Shape = base
                    if not pl.isNull():
                        obj.Placement = pl


class _ViewProviderPanel(ArchComponent.ViewProviderComponent):

    "A View Provider for the Panel object"

    def __init__(self,vobj):

        ArchComponent.ViewProviderComponent.__init__(self,vobj)
        vobj.ShapeColor = ArchCommands.getDefaultColor("Panel")

    def getIcon(self):

        #import Arch_rc
        if hasattr(self,"Object"):
            if hasattr(self.Object,"CloneOf"):
                if self.Object.CloneOf:
                    return ":/icons/Arch_Panel_Clone.svg"
        return ":/icons/Arch_Panel_Tree.svg"

    def updateData(self,obj,prop):

        if prop in ["Placement","Shape","Material"]:
            if hasattr(obj,"Material"):
                if obj.Material:
                    if hasattr(obj.Material,"Materials"):
                        activematerials = [obj.Material.Materials[i] for i in range(len(obj.Material.Materials)) if obj.Material.Thicknesses[i] >= 0]
                        if len(activematerials) == len(obj.Shape.Solids):
                            cols = []
                            for i,mat in enumerate(activematerials):
                                c = obj.ViewObject.ShapeColor
                                c = (c[0],c[1],c[2],1.0-obj.ViewObject.Transparency/100.0)
                                if 'DiffuseColor' in mat.Material:
                                    if "(" in mat.Material['DiffuseColor']:
                                        c = tuple([float(f) for f in mat.Material['DiffuseColor'].strip("()").split(",")])
                                if 'Transparency' in mat.Material:
                                    c = (c[0],c[1],c[2],1.0-float(mat.Material['Transparency']))
                                cols.extend([c for j in range(len(obj.Shape.Solids[i].Faces))])
                            if obj.ViewObject.DiffuseColor != cols:
                                obj.ViewObject.DiffuseColor = cols
        ArchComponent.ViewProviderComponent.updateData(self,obj,prop)


class PanelCut(Draft.DraftObject):

    "A flat, 2D view of an Arch Panel"

    def __init__(self, obj):
        Draft.DraftObject.__init__(self,obj)
        obj.Proxy = self

        # setProperties of ArchComponent will be overwritten
        # thus setProperties from ArchComponent will be explicit called to get the properties
        ArchComponent.ViewProviderComponent.setProperties(self, obj)

        self.setProperties(obj)

    def setProperties(self,obj):

        pl = obj.PropertiesList
        if not "Source" in pl:
            obj.addProperty("App::PropertyLink","Source","PanelCut",QT_TRANSLATE_NOOP("App::Property","The linked object"), locked=True)
        if not "TagText" in pl:
            obj.addProperty("App::PropertyString","TagText","PanelCut",QT_TRANSLATE_NOOP("App::Property","The text to display. Can be %tag%, %label% or %description% to display the panel tag or label"), locked=True)
            obj.TagText = "%tag%"
        if not "TagSize" in pl:
            obj.addProperty("App::PropertyLength","TagSize","PanelCut",QT_TRANSLATE_NOOP("App::Property","The size of the tag text"), locked=True)
            obj.TagSize = 10
        if not "TagPosition" in pl:
            obj.addProperty("App::PropertyVector","TagPosition","PanelCut",QT_TRANSLATE_NOOP("App::Property","The position of the tag text. Keep (0,0,0) for center position"), locked=True)
        if not "TagRotation" in pl:
            obj.addProperty("App::PropertyAngle","TagRotation","PanelCut",QT_TRANSLATE_NOOP("App::Property","The rotation of the tag text"), locked=True)
        if not "FontFile" in pl:
            obj.addProperty("App::PropertyFile","FontFile","PanelCut",QT_TRANSLATE_NOOP("App::Property","The font of the tag text"), locked=True)
            obj.FontFile = params.get_param("FontFile")
        if not "MakeFace" in pl:
            obj.addProperty("App::PropertyBool","MakeFace","PanelCut",QT_TRANSLATE_NOOP("App::Property","If True, the object is rendered as a face, if possible."), locked=True)
        if not "AllowedAngles" in pl:
            obj.addProperty("App::PropertyFloatList","AllowedAngles","PanelCut",QT_TRANSLATE_NOOP("App::Property","The allowed angles this object can be rotated to when placed on sheets"), locked=True)
        self.Type = "PanelCut"
        if not "CutOffset" in pl:
            obj.addProperty("App::PropertyDistance","CutOffset","PanelCut",QT_TRANSLATE_NOOP("App::Property","An offset value to move the cut plane from the center point"), locked=True)

    def onDocumentRestored(self,obj):

        self.setProperties(obj)

    def execute(self, obj):

        pl = obj.Placement
        if obj.Source:
            base = None
            n = None
            if Draft.getType(obj.Source) == "Panel":
                import DraftGeomUtils
                import Part
                baseobj = None
                if obj.Source.CloneOf:
                    baseobj = obj.Source.CloneOf.Base
                if obj.Source.Base:
                    baseobj = obj.Source.Base
                if baseobj:
                    if hasattr(baseobj,'Shape'):
                        if baseobj.Shape.Solids:
                            center = baseobj.Shape.BoundBox.Center
                            diag = baseobj.Shape.BoundBox.DiagonalLength
                            if obj.Source.Normal.Length:
                                n = obj.Source.Normal
                            elif baseobj.isDerivedFrom("Part::Extrusion"):
                                n = baseobj.Dir
                            if not n:
                                n = Vector(0,0,1)
                            if hasattr(obj,"CutOffset") and obj.CutOffset.Value:
                                l = obj.CutOffset.Value
                                d = Vector(n)
                                d.multiply(l)
                                center = center.add(d)
                            plane = Part.makePlane(diag,diag,center,n)
                            plane.translate(center.sub(plane.BoundBox.Center))
                            wires = []
                            for sol in baseobj.Shape.Solids:
                                s = sol.section(plane)
                                wires.extend(DraftGeomUtils.findWires(s.Edges))
                            if wires:
                                base = self.buildCut(obj,wires)
                        else:
                            base = self.buildCut(obj,baseobj.Shape.Wires)
                            for w in base.Wires:
                                n = DraftGeomUtils.getNormal(w)
                                if n:
                                    break
                            if not n:
                                n = Vector(0,0,1)
                        if base and n:
                            base.translate(base.BoundBox.Center.negative())
                            r = FreeCAD.Rotation(n,Vector(0,0,1))
                            base.rotate(Vector(0,0,0),r.Axis,math.degrees(r.Angle))
                    elif baseobj.isDerivedFrom("Mesh::Feature"):
                        return
                else:
                    l2 = obj.Source.Length/2
                    w2 = obj.Source.Width/2
                    v1 = Vector(-l2,-w2,0)
                    v2 = Vector(l2,-w2,0)
                    v3 = Vector(l2,w2,0)
                    v4 = Vector(-l2,w2,0)
                    base = Part.makePolygon([v1,v2,v3,v4,v1])
                if base:
                    self.outline = base
                    if obj.FontFile and obj.TagText and obj.TagSize.Value:
                        if obj.TagPosition.Length == 0:
                            pos = base.BoundBox.Center
                        else:
                            pos = obj.TagPosition
                        if obj.TagText == "%tag%":
                            string = obj.Source.Tag
                        elif obj.TagText == "%label%":
                            string = obj.Source.Label
                        elif obj.TagText == "%description%":
                            string = obj.Source.Description
                        else:
                            string = obj.TagText
                        chars = []
                        for char in Part.makeWireString(string,obj.FontFile,obj.TagSize.Value,0):
                            chars.extend(char)
                        textshape = Part.Compound(chars)
                        textshape.translate(pos.sub(textshape.BoundBox.Center))
                        textshape.rotate(textshape.BoundBox.Center,Vector(0,0,1),obj.TagRotation.Value)
                        self.tag = textshape
                        base = Part.Compound([base,textshape])
                    else:
                        base = Part.Compound([base])
                    obj.Shape = base
                    obj.Placement = pl

    def buildCut(self,obj,wires):

        """buildCut(obj,wires): builds the object shape"""

        import Part
        if hasattr(obj,"MakeFace"):
            if obj.MakeFace:
                face = None
                if len(wires) > 1:
                    d = 0
                    ow = None
                    for w in wires:
                        if w.BoundBox.DiagonalLength > d:
                            d = w.BoundBox.DiagonalLength
                            ow = w
                    if ow:
                        face = Part.Face(ow)
                        for w in wires:
                            if w.hashCode() != ow.hashCode():
                                wface = Part.Face(w)
                                face = face.cut(wface)
                else:
                    face = Part.Face(wires[0])
                if face:
                    return face
        return Part.makeCompound(wires)

    def getWires(self, obj):

        """getWires(obj): returns a tuple containing 3 shapes
        that define the panel outline, the panel holes, and
        tags (engravings): (outline,holes,tags). Any of these can
        be None if nonexistent"""

        tag = None
        outl = None
        inl = None
        if not hasattr(self,"outline"):
            self.execute(obj)
        if not hasattr(self,"outline"):
            return None
        outl = self.outline.copy()
        if hasattr(self,"tag"):
            tag = self.tag.copy()
        if tag:
            tag.Placement = obj.Placement.multiply(tag.Placement)

        outl = self.outline.copy()
        outl.Placement = obj.Placement.multiply(outl.Placement)
        if len(outl.Wires) > 1:
            # separate outline
            d = 0
            ow = None
            for w in outl.Wires:
                if w.BoundBox.DiagonalLength > d:
                    d = w.BoundBox.DiagonalLength
                    ow = w
            if ow:
                inl = Part.Compound([w for w in outl.Wires if w.hashCode() != ow.hashCode()])
                outl = Part.Compound([ow])
        else:
            inl = None
            outl = Part.Compound([outl.Wires[0]])
        return (outl, inl, tag)


class ViewProviderPanelCut(Draft.ViewProviderDraft):

    "a view provider for the panel cut object"

    def __init__(self,vobj):

        Draft.ViewProviderDraft.__init__(self,vobj)
        self.setProperties(vobj)

    def setProperties(self,vobj):

        pl = vobj.PropertiesList
        if not "Margin" in pl:
            vobj.addProperty("App::PropertyLength","Margin","Arch",QT_TRANSLATE_NOOP("App::Property","A margin inside the boundary"), locked=True)
        if not "ShowMargin" in pl:
            vobj.addProperty("App::PropertyBool","ShowMargin","Arch",QT_TRANSLATE_NOOP("App::Property","Turns the display of the margin on/off"), locked=True)

    def onDocumentRestored(self,vobj):

        self.setProperties(vobj)

    def attach(self,vobj):

        Draft.ViewProviderDraft.attach(self,vobj)
        from pivy import coin
        self.coords = coin.SoCoordinate3()
        self.lineset = coin.SoLineSet()
        self.lineset.numVertices.setValue(-1)
        lineStyle = coin.SoDrawStyle()
        lineStyle.linePattern = 0x0f0f
        self.color = coin.SoBaseColor()
        self.switch = coin.SoSwitch()
        sep = coin.SoSeparator()
        self.switch.whichChild = -1
        sep.addChild(self.color)
        sep.addChild(lineStyle)
        sep.addChild(self.coords)
        sep.addChild(self.lineset)
        self.switch.addChild(sep)
        vobj.Annotation.addChild(self.switch)
        self.onChanged(vobj,"ShowMargin")
        self.onChanged(vobj,"LineColor")

    def onChanged(self,vobj,prop):

        if prop in ["Margin","ShowMargin"]:
            if hasattr(vobj,"Margin") and hasattr(vobj,"ShowMargin"):
                if (vobj.Margin.Value > 0) and vobj.Object.Shape and vobj.ShowMargin:
                    self.lineset.numVertices.setValue(-1)
                    if vobj.Object.Shape.Wires:
                        d = 0
                        dw = None
                        for w in vobj.Object.Shape.Wires:
                            if w.BoundBox.DiagonalLength > d:
                                d = w.BoundBox.DiagonalLength
                                dw = w
                        if dw:
                            ow = dw.makeOffset2D(vobj.Margin.Value)
                            verts = []
                            for v in ow.OrderedVertexes:
                                v = vobj.Object.Placement.inverse().multVec(v.Point)
                                verts.append((v.x,v.y,v.z))
                            if dw.isClosed():
                                verts.append(verts[0])
                        self.coords.point.setValues(verts)
                        self.lineset.numVertices.setValue(len(verts))
                        self.switch.whichChild = 0
                else:
                    self.switch.whichChild = -1
        elif prop == "LineColor":
            if hasattr(vobj,"LineColor"):
                c = vobj.LineColor
                self.color.rgb.setValue(c[0],c[1],c[2])
        Draft.ViewProviderDraft.onChanged(self,vobj,prop)

    def updateData(self,obj,prop):

        if prop in ["Shape"]:
            self.onChanged(obj.ViewObject,"Margin")
        Draft.ViewProviderDraft.updateData(self,obj,prop)

    def doubleClicked(self,vobj):

        # See setEdit in ViewProviderDraft.
        FreeCADGui.runCommand("Std_TransformManip")
        return True

class PanelSheet(Draft.DraftObject):

    "A collection of Panel cuts under a sheet"

    def __init__(self, obj):

        Draft.DraftObject.__init__(self,obj)
        obj.Proxy = self
        self.setProperties(obj)

    def setProperties(self,obj):

        pl = obj.PropertiesList
        if not "Group" in pl:
            obj.addProperty("App::PropertyLinkList","Group","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The linked Panel cuts"), locked=True)
        if not "TagText" in pl:
            obj.addProperty("App::PropertyString","TagText","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The tag text to display"), locked=True)
        if not "TagSize" in pl:
            obj.addProperty("App::PropertyLength","TagSize","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The size of the tag text"), locked=True)
            obj.TagSize = 10
        if not "TagPosition" in pl:
            obj.addProperty("App::PropertyVector","TagPosition","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The position of the tag text. Keep (0,0,0) for center position"), locked=True)
        if not "TagRotation" in pl:
            obj.addProperty("App::PropertyAngle","TagRotation","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The rotation of the tag text"), locked=True)
        if not "FontFile" in pl:
            obj.addProperty("App::PropertyFile","FontFile","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The font of the tag text"), locked=True)
            obj.FontFile = params.get_param("FontFile")
        if not "Width" in pl:
            obj.addProperty("App::PropertyLength","Width","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The width of the sheet"), locked=True)
            obj.Width = params.get_param_arch("PanelLength")
        if not "Height" in pl:
            obj.addProperty("App::PropertyLength","Height","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The height of the sheet"), locked=True)
            obj.Height = params.get_param_arch("PanelWidth")
        if not "FillRatio" in pl:
            obj.addProperty("App::PropertyPercent","FillRatio","PanelSheet",QT_TRANSLATE_NOOP("App::Property","The fill ratio of this sheet"), locked=True)
            obj.setEditorMode("FillRatio",2)
        if not "MakeFace" in pl:
            obj.addProperty("App::PropertyBool","MakeFace","PanelSheet",QT_TRANSLATE_NOOP("App::Property","If True, the object is rendered as a face, if possible."), locked=True)
        if not "GrainDirection" in pl:
            obj.addProperty("App::PropertyAngle","GrainDirection","PanelSheet",QT_TRANSLATE_NOOP("App::Property","Specifies an angle for the wood grain (Clockwise, 0 is North)"), locked=True)
        if not "Scale" in pl:
            obj.addProperty("App::PropertyFloat","Scale","PanelSheet", QT_TRANSLATE_NOOP("App::Property","Specifies the scale applied to each panel view."), locked=True)
            obj.Scale = 1.0
        if not "Rotations" in pl:
            obj.addProperty("App::PropertyFloatList","Rotations","PanelSheet", QT_TRANSLATE_NOOP("App::Property","A list of possible rotations for the nester"), locked=True)
        self.Type = "PanelSheet"

    def onDocumentRestored(self, obj):

        self.setProperties(obj)

    def execute(self, obj):

        import Part
        self.sheettag = None
        self.sheetborder = None
        pl = obj.Placement
        if obj.Width.Value and obj.Height.Value:
            l2 = obj.Width.Value/2
            w2 = obj.Height.Value/2
            v1 = Vector(-l2,-w2,0)
            v2 = Vector(l2,-w2,0)
            v3 = Vector(l2,w2,0)
            v4 = Vector(-l2,w2,0)
            base = Part.makePolygon([v1,v2,v3,v4,v1])
            if hasattr(obj,"MakeFace"):
                if obj.MakeFace:
                    base = Part.Face(base)
            self.sheetborder = base
            wires = []
            area = obj.Width.Value * obj.Height.Value
            subarea = 0
            for v in obj.Group:
                if hasattr(v,'Shape'):
                    wires.extend(v.Shape.Wires)
                    if Draft.getType(v) == "PanelCut":
                        if v.Source:
                            subarea += v.Source.Area.Value
                    else:
                        for w in v.Shape.Wires:
                            if w.isClosed():
                                f = Part.Face(w)
                                subarea += f.Area
            if wires:
                base = Part.Compound([base]+wires)
            if obj.FontFile and obj.TagText and obj.TagSize.Value:
                chars = []
                for char in Part.makeWireString(obj.TagText,obj.FontFile,obj.TagSize.Value,0):
                    chars.extend(char)
                textshape = Part.Compound(chars)
                textshape.translate(obj.TagPosition)
                textshape.rotate(textshape.BoundBox.Center,Vector(0,0,1),obj.TagRotation.Value)
                self.sheettag = textshape
                base = Part.Compound([base,textshape])
            base.scale(obj.Scale, FreeCAD.Vector())
            obj.Shape = base
            obj.Placement = pl
            obj.FillRatio = int((subarea/area)*100)

    def getOutlines(self,obj,transform=False):

        """getOutlines(obj,transform=False): returns a list of compounds whose wires define the
        outlines of the panels in this sheet. If transform is True, the placement of
        the sheet will be added to each wire"""

        outp = []
        for p in obj.Group:
            ispanel = False
            if hasattr(p,"Proxy"):
                if hasattr(p.Proxy,"getWires"):
                    ispanel = True
                    w = p.Proxy.getWires(p)
                    if w[0]:
                        w = w[0]
                        w.scale(obj.Scale, FreeCAD.Vector())
                        if transform:
                            w.Placement = obj.Placement.multiply(w.Placement)
                        outp.append(w)
            if not ispanel:
                if hasattr(p,'Shape'):
                    for w in p.Shape.Wires:
                        w.scale(obj.Scale, FreeCAD.Vector())
                        if transform:
                            w.Placement = obj.Placement.multiply(w.Placement)
                        outp.append(w)
        return outp

    def getHoles(self,obj,transform=False):

        """getHoles(obj,transform=False): returns a list of compound whose wires define the
        holes contained in the panels in this sheet. If transform is True, the placement of
        the sheet will be added to each wire"""

        outp = []
        for p in obj.Group:
            if hasattr(p,"Proxy"):
                if hasattr(p.Proxy,"getWires"):
                    w = p.Proxy.getWires(p)
                    if w[1]:
                        w = w[1]
                        w.scale(obj.Scale, FreeCAD.Vector())
                        if transform:
                            w.Placement = obj.Placement.multiply(w.Placement)
                        outp.append(w)
        return outp

    def getTags(self,obj,transform=False):

        """getTags(obj,transform=False): returns a list of compounds whose wires define the
        tags (engravings) contained in the panels in this sheet and the sheet intself.
        If transform is True, the placement of the sheet will be added to each wire.
        Warning, the wires returned by this function may not be closed,
        depending on the font"""

        outp = []
        for p in obj.Group:
            if hasattr(p,"Proxy"):
                if hasattr(p.Proxy,"getWires"):
                    w = p.Proxy.getWires(p)
                    if w[2]:
                        w = w[2]
                        w.scale(obj.Scale, FreeCAD.Vector())
                        if transform:
                            w.Placement = obj.Placement.multiply(w.Placement)
                        outp.append(w)
        if self.sheettag is not None:
            w = self.sheettag.copy()
            w.scale(obj.Scale, FreeCAD.Vector())
            if transform:
                w.Placement = obj.Placement.multiply(w.Placement)
            outp.append(w)

        return outp


class ViewProviderPanelSheet(Draft.ViewProviderDraft):

    "a view provider for the panel sheet object"

    def __init__(self,vobj):

        Draft.ViewProviderDraft.__init__(self,vobj)
        self.setProperties(vobj)
        vobj.PatternSize = 0.0035

    def setProperties(self,vobj):

        pl = vobj.PropertiesList
        if not "Margin" in pl:
            vobj.addProperty("App::PropertyLength","Margin","PanelSheet",QT_TRANSLATE_NOOP("App::Property","A margin inside the boundary"), locked=True)
        if not "ShowMargin" in pl:
            vobj.addProperty("App::PropertyBool","ShowMargin","PanelSheet",QT_TRANSLATE_NOOP("App::Property","Turns the display of the margin on/off"), locked=True)
        if not "ShowGrain" in pl:
            vobj.addProperty("App::PropertyBool","ShowGrain","PanelSheet",QT_TRANSLATE_NOOP("App::Property","Turns the display of the wood grain texture on/off"), locked=True)

    def onDocumentRestored(self,vobj):

        self.setProperties(vobj)

    def getIcon(self):

        return ":/icons/Arch_Panel_Sheet_Tree.svg"

    def setEdit(self, vobj, mode):
        if mode == 1 or mode == 2:
            return None

        taskd = SheetTaskPanel(vobj.Object)
        taskd.update()
        FreeCADGui.Control.showDialog(taskd)
        return True

    def unsetEdit(self, vobj, mode):
        if mode == 1 or mode == 2:
            return None

        FreeCADGui.Control.closeDialog()
        return True

    def attach(self,vobj):

        Draft.ViewProviderDraft.attach(self,vobj)
        from pivy import coin
        self.coords = coin.SoCoordinate3()
        self.lineset = coin.SoLineSet()
        self.lineset.numVertices.setValue(-1)
        lineStyle = coin.SoDrawStyle()
        lineStyle.linePattern = 0x0f0f
        self.color = coin.SoBaseColor()
        self.switch = coin.SoSwitch()
        sep = coin.SoSeparator()
        self.switch.whichChild = -1
        sep.addChild(self.color)
        sep.addChild(lineStyle)
        sep.addChild(self.coords)
        sep.addChild(self.lineset)
        self.switch.addChild(sep)
        vobj.Annotation.addChild(self.switch)
        self.onChanged(vobj,"ShowMargin")
        self.onChanged(vobj,"LineColor")

    def onChanged(self,vobj,prop):

        if prop in ["Margin","ShowMargin"]:
            if hasattr(vobj,"Margin") and hasattr(vobj,"ShowMargin"):
                if (vobj.Margin.Value > 0) and (vobj.Margin.Value < vobj.Object.Width.Value/2) and (vobj.Margin.Value < vobj.Object.Height.Value/2):
                    l2 = vobj.Object.Width.Value/2
                    w2 = vobj.Object.Height.Value/2
                    v = vobj.Margin.Value
                    v1 = (-l2+v,-w2+v,0)
                    v2 = (l2-v,-w2+v,0)
                    v3 = (l2-v,w2-v,0)
                    v4 = (-l2+v,w2-v,0)
                    self.coords.point.setValues([v1,v2,v3,v4,v1])
                    self.lineset.numVertices.setValue(5)
                if vobj.ShowMargin:
                    self.switch.whichChild = 0
                else:
                    self.switch.whichChild = -1
        elif prop == "LineColor":
            if hasattr(vobj,"LineColor"):
                c = vobj.LineColor
                self.color.rgb.setValue(c[0],c[1],c[2])
        elif prop == "ShowGrain":
            if hasattr(vobj,"ShowGrain"):
                if vobj.ShowGrain:
                    vobj.Pattern = "woodgrain"
                else:
                    vobj.Pattern = "None"
        Draft.ViewProviderDraft.onChanged(self,vobj,prop)


    def updateData(self,obj,prop):

        if prop in ["Width","Height"]:
            self.onChanged(obj.ViewObject,"Margin")
        elif prop == "GrainDirection":
            if hasattr(self,"texcoords"):
                if self.texcoords:
                    s = FreeCAD.Vector(self.texcoords.directionS.getValue().getValue()).Length
                    vS  = DraftVecUtils.rotate(FreeCAD.Vector(s,0,0),-math.radians(obj.GrainDirection.Value))
                    vT  = DraftVecUtils.rotate(FreeCAD.Vector(0,s,0),-math.radians(obj.GrainDirection.Value))
                    self.texcoords.directionS.setValue(vS.x,vS.y,vS.z)
                    self.texcoords.directionT.setValue(vT.x,vT.y,vT.z)
        Draft.ViewProviderDraft.updateData(self,obj,prop)


class SheetTaskPanel(ArchComponent.ComponentTaskPanel):

    def __init__(self,obj):

        ArchComponent.ComponentTaskPanel.__init__(self)
        self.obj = obj
        self.optwid = QtGui.QWidget()
        self.optwid.setWindowTitle(QtGui.QApplication.translate("Arch", "Tools", None))
        lay = QtGui.QVBoxLayout(self.optwid)
        self.editButton = QtGui.QPushButton(self.optwid)
        self.editButton.setIcon(QtGui.QIcon(":/icons/Draft_Edit.svg"))
        self.editButton.setText(QtGui.QApplication.translate("Arch", "Edit views positions", None))
        lay.addWidget(self.editButton)
        QtCore.QObject.connect(self.editButton, QtCore.SIGNAL("clicked()"), self.editNodes)
        self.form = [self.form,self.optwid]

    def editNodes(self):

        FreeCADGui.Control.closeDialog()
        FreeCADGui.runCommand("Draft_Edit")
