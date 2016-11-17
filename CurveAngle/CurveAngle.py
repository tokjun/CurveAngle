import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
import math

#
# 3D Slicer module to check if the models under the model hierarchy is intersecting
# the trajoectory given by the fiducial node.
#

class CurveAngle(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "CurveAngle" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Informatics"]
    self.parent.dependencies = []
    self.parent.contributors = ["Sarah Elbenna (ENSEEIHT), Junichi Tokuda (BWH)"]
    self.parent.helpText = """
    This module traces a curve and lists structures intersecting with it.
    """
    self.parent.acknowledgementText = """
    This module was created using a template developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
    """ # replace with organization, grant and thanks.

#
# CurveAngleWidget
#

class CurveAngleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    ####################
    # For debugging
    #
    # Reload and Test area
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload && Test"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "CurveAngle Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)
    #
    ####################
    #
    # Entry Angle area
    #
    angleCollapsibleButton = ctk.ctkCollapsibleButton()
    angleCollapsibleButton.text = "Entry Angle"
    angleCollapsibleButton.collapsed = True
    self.layout.addWidget(angleCollapsibleButton)
    angleFormLayout = qt.QFormLayout(angleCollapsibleButton)

    # - input fiducials (trajectory) selector
    distanceLayout = qt.QVBoxLayout()

    self.inputFiducialSelector = slicer.qMRMLNodeComboBox()
    self.inputFiducialSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.inputFiducialSelector.selectNodeUponCreation = True
    self.inputFiducialSelector.addEnabled = True
    self.inputFiducialSelector.removeEnabled = True
    self.inputFiducialSelector.noneEnabled = False
    self.inputFiducialSelector.showHidden = False
    self.inputFiducialSelector.showChildNodeTypes = False
    self.inputFiducialSelector.setMRMLScene( slicer.mrmlScene )
    self.inputFiducialSelector.setToolTip( "Pick the trajectory." )
    angleFormLayout.addRow("Trajectory: ", self.inputFiducialSelector)

    #  - Model selector
    self.inputModelSelector = slicer.qMRMLNodeComboBox()
    self.inputModelSelector.nodeTypes = ["vtkMRMLModelHierarchyNode"]
    self.inputModelSelector.selectNodeUponCreation = True
    self.inputModelSelector.addEnabled = True
    self.inputModelSelector.removeEnabled = True
    self.inputModelSelector.noneEnabled = False
    self.inputModelSelector.showHidden = False
    self.inputModelSelector.showChildNodeTypes = False
    self.inputModelSelector.setMRMLScene( slicer.mrmlScene )
    self.inputModelSelector.setToolTip( "Select a 3D Model" )
    angleFormLayout.addRow("Model:", self.inputModelSelector)

    self.modelNode = None
    self.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onModelSelected)

    self.angleTable = qt.QTableWidget(1, 2)
    self.angleTable.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
    self.angleTable.setSelectionMode(qt.QAbstractItemView.SingleSelection)
    self.angleTableHeaders = ["Model", "Entry Angle"]
    self.angleTable.setHorizontalHeaderLabels(self.angleTableHeaders)
    self.angleTable.horizontalHeader().setStretchLastSection(True)
    angleFormLayout.addWidget(self.angleTable)

    self.extrapolateCheckBox = qt.QCheckBox()
    self.extrapolateCheckBox.checked = 0
    self.extrapolateCheckBox.setToolTip("Extrapolate the first and last segment to calculate the distance")
    self.extrapolateCheckBox.connect('toggled(bool)', self.updateAngleTable)
    self.extrapolateCheckBox.text = 'Extrapolate curves to measure the distances'

    self.showErrorVectorCheckBox = qt.QCheckBox()
    self.showErrorVectorCheckBox.checked = 0
    self.showErrorVectorCheckBox.setToolTip("Show error vectors, which is defined by the target point and the closest point on the curve. The vector is perpendicular to the curve, unless the closest point is one end of the curve.")
    self.showErrorVectorCheckBox.connect('toggled(bool)', self.updateAngleTable)
    self.showErrorVectorCheckBox.text = 'Show error vectors'

    distanceLayout.addWidget(self.extrapolateCheckBox)
    distanceLayout.addWidget(self.showErrorVectorCheckBox)
    angleFormLayout.addRow(distanceLayout)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    angleFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)


    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputModelSelector.currentNode() and self.inputFiducialSelector.currentNode()

  def onApplyButton(self):
    logic = CurveAngleLogic()
    #logic.EntryAngle(self.inputModelSelector.currentNode(), self.inputFiducialSelector.currentNode())
    logic.CheckIntersections(self.inputModelSelector.currentNode(), self.inputFiducialSelector.currentNode())

  def onReload(self,moduleName="CurveAngle"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)



  def onModelSelected(self):

    # Remove observer if previous node exists
    if self.modelNode and self.tag:
      self.modelNode.RemoveObserver(self.tag)

    # Update selected node, add observer, and update control points
    if self.inputModelSelector.currentNode():
      self.modelNode = self.inputModelSelector.currentNode()
      self.tag = self.modelNode.AddObserver('ModifiedEvent', self.onAngleUpdated)
    else:
      self.modelNode = None
      self.tag = None
    self.updateAngleTable()

  def onAngleUpdated(self,caller,event):
    if caller.IsA('vtkMRMLModelHierarchyNode') and event == 'ModifiedEvent':
      #self.updateAngleTable()
      pass


  def updateAngleTable(self):

    ##    logic = CurveAngleLogic()

    if not self.modelNode:
      self.angleTable.clear()
      self.angleTable.setHorizontalHeaderLabels(self.angleTableHeaders)

    else:

      self.angleTableData = []
      nOfControlPoints = 0
      if self.inputModelSelector:
        nOfControlPoints = self.inputModelSelector.GetNumberOfChildrenNodes()

      if self.angleTable.rowCount != nOfControlPoints:
        self.angleTable.setRowCount(nOfControlPoints)
        row = ['aa', 'bb']
        self.angleTable.setItem(i, 0, row[0])
        self.angleTable.setItem(i, 1, row[1])
        self.angleTableData.append(row)
    self.angleTable.show()

#
# CurveAngleLogic
#

class CurveAngleLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def CheckIntersections(self, inputModelNode, inputFiducialNode):
    """
    Run the actual algorithm
    """

    logging.info('Processing started')
    print ("EntryAngle() is called.")

    if inputModelNode == None:
      return None
    if inputFiducialNode == None:
      return None

    nOfModels = inputModelNode.GetNumberOfChildrenNodes()
    objectIDs = []
    objectNames = []
    entryAngles = []
    totalLengthInObject = []

    for i in range(nOfModels):
        chnode = inputModelNode.GetNthChildNode(i)
        if chnode == None:
            continue
        mnode = chnode.GetAssociatedNode()
        if mnode == None:
            continue

        name = mnode.GetName()
        objectPoly = mnode.GetPolyData()
        if objectPoly == None:
            continue

        print "Processing object: %s" % name

        intersectingPoints = vtk.vtkPoints()
        trajectoryPoints = vtk.vtkPoints()

        idList = vtk.vtkIdList()
        pos0 = [0.0] * 3
        pos1 = [0.0] * 3

        nFiducials = inputFiducialNode.GetNumberOfFiducials()
        posStart = None
        posEnd = None

        # Look for points inside the object
        for j in range(nFiducials):
            inputFiducialNode.GetNthFiducialPosition(j, pos0)
            trajectoryPoints.InsertNextPoint(pos0)

        trajectoryPoly = vtk.vtkPolyData()
        trajectoryPoly.SetPoints(trajectoryPoints)
        enclosed = vtk.vtkSelectEnclosedPoints()
        enclosed.SetInputData(trajectoryPoly)
        enclosed.SetSurfaceData(objectPoly)
        enclosed.Update()

        lengthInObject = 0.0

        isInside = False

        angles = []
        for j in range(nFiducials-1):

            inputFiducialNode.GetNthFiducialPosition(j, pos0)
            inputFiducialNode.GetNthFiducialPosition(j+1, pos1)

            inout0 = enclosed.IsInside(j)
            inout1 = enclosed.IsInside(j+1)

            print "Point %d: from (%f, %f, %f) (%d) to (%f, %f, %f) (%d)" % (j, pos0[0], pos0[1], pos0[2], inout0, pos1[0], pos1[1], pos1[2], inout1)

            traj = np.array(pos1) - np.array(pos0)

            if inout0 and inout1:
                ## Both in the object
                lSegment = np.linalg.norm(traj)
                lengthInObject = lengthInObject + lSegment

            if inout0 != inout1:
                ## Potential intersection

                bspTree = vtk.vtkModifiedBSPTree()
                bspTree.SetDataSet(objectPoly)
                bspTree.BuildLocator()
                tolerance = 0.001
                bspTree.IntersectWithLine(pos0, pos1, tolerance, intersectingPoints, idList)

                if intersectingPoints.GetNumberOfPoints < 1:
                    continue

                if idList.GetNumberOfIds() < 1:
                    continue

                cell = objectPoly.GetCell(idList.GetId(0))
                if cell == None:
                    continue

                isInside = True

                # Get intersecting point and measure the length inside the boject
                p = [0.0]*3
                intersectingPoints.GetPoint(0, p)
                if inout0:
                    segmentInObject = np.array(p) - np.array(pos0)
                    lengthInObject = lengthInObject + np.linalg.norm(segmentInObject)
                elif inout1:
                    segmentInObject = np.array(p) - np.array(pos1)
                    lengthInObject = lengthInObject + np.linalg.norm(segmentInObject)

                points2 = cell.GetPoints()
                if points2 == None:
                    continue
                p0 = [0.0] * 3
                p1 = [0.0] * 3
                p2 = [0.0] * 3
                points2.GetPoint(0, p0)
                points2.GetPoint(1, p1)
                points2.GetPoint(2, p2)
                npap0 = np.array(p0)
                npap1 = np.array(p1)
                npap2 = np.array(p2)
                v0 = npap1-npap0
                v1 = npap2-npap0
                v = np.cross(v0, v1)
                norm = np.linalg.norm(v,ord=1)
                normVec = v/norm
                angle = vtk.vtkMath.AngleBetweenVectors(normVec, traj)
                angles.append(angle)

                print "  -- Intersecting at (%f, %f, %f) with angle %f" % (p[0], p[1], p[2], angle)

        print "Length in object = %f" % lengthInObject

        if isInside:
            objectIDs.append(i)
            objectNames.append(name)
            entryAngles.append(angles)
            totalLengthInObject.append(lengthInObject)

    return (objectIDs, entryAngles, totalLengthInObject)


class CurveAngleTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_CurveAngle1()

  def test_CurveAngle1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = CurveAngleLogic()
    self.delayDisplay('Test passed!')
