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

class PathCollisionAnalysis(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "PathCollisionAnalysis" # TODO make this more human readable by adding spaces
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
# PathCollisionAnalysisWidget
#

class PathCollisionAnalysisWidget(ScriptedLoadableModuleWidget):
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
    reloadCollapsibleButton.collapsed = True
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "PathCollisionAnalysis Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)
    #
    ####################
    #
    # Main area
    #
    mainCollapsibleButton = ctk.ctkCollapsibleButton()
    mainCollapsibleButton.text = "Collision Analysis"
    mainCollapsibleButton.collapsed = False
    self.layout.addWidget(mainCollapsibleButton)
    mainFormLayout = qt.QFormLayout(mainCollapsibleButton)

    # - input fiducials (trajectory) selector
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
    mainFormLayout.addRow("Trajectory: ", self.inputFiducialSelector)

    #  - Model selector
    self.inputModelHierarchySelector = slicer.qMRMLNodeComboBox()
    self.inputModelHierarchySelector.nodeTypes = ["vtkMRMLModelHierarchyNode"]
    self.inputModelHierarchySelector.selectNodeUponCreation = True
    self.inputModelHierarchySelector.addEnabled = True
    self.inputModelHierarchySelector.removeEnabled = True
    self.inputModelHierarchySelector.noneEnabled = False
    self.inputModelHierarchySelector.showHidden = False
    self.inputModelHierarchySelector.showChildNodeTypes = False
    self.inputModelHierarchySelector.setMRMLScene( slicer.mrmlScene )
    self.inputModelHierarchySelector.setToolTip( "Select a 3D Model" )
    mainFormLayout.addRow("Model:", self.inputModelHierarchySelector)

    self.modelNode = None
    #self.inputModelHierarchySelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onModelSelected)

    self.intersectionTable = qt.QTableWidget(1, 4)
    self.intersectionTable.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
    self.intersectionTable.setSelectionMode(qt.QAbstractItemView.SingleSelection)
    self.intersectionTableHeader = ["Model", "Entry 1", "Entry 2", "Length"]
    self.intersectionTable.setHorizontalHeaderLabels(self.intersectionTableHeader)
    self.intersectionTable.horizontalHeader().setStretchLastSection(True)
    mainFormLayout.addWidget(self.intersectionTable)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    mainFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputModelHierarchySelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

    # variables
    self.objectIDs = None
    self.objectNames = None
    self.normalVectors = None
    self.entryAngles = None
    self.totalLengthInObject = None

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputModelHierarchySelector.currentNode() and self.inputFiducialSelector.currentNode()

  def onApplyButton(self):
    logic = PathCollisionAnalysisLogic()
    #logic.EntryAngle(self.inputModelSelector.currentNode(), self.inputFiducialSelector.currentNode())
    [self.objectIDs, self.objectNames, self.normalVectors, self.entryAngles, self.totalLengthInObject] = logic.CheckIntersections(self.inputModelHierarchySelector.currentNode(), self.inputFiducialSelector.currentNode())
    self.updateIntersectionTable()

  def updateIntersectionTable(self):
    ##    logic = PathCollisionAnalysisLogic()
    modelHierarchyNode = self.inputModelHierarchySelector.currentNode()
    fiducialNode = self.inputFiducialSelector.currentNode()

    if (not modelHierarchyNode) or (not fiducialNode):
        self.intersectionTable.clear()
        self.intersectionTable.setHorizontalHeaderLabels(self.intersectionTableHeader)

    else:
        nObjects = len(self.objectIDs)
        self.intersectionTable.setRowCount(nObjects)
        for i in range(nObjects):
            # "Model", "Entry 1", "Entry 2", "Length"
            self.intersectionTable.setItem(i, 0, qt.QTableWidgetItem(self.objectNames[i]))
            angles = self.entryAngles[i]
            normals = self.normalVectors[i]
            nEntry = 2
            if len(angles) < 2:
                nEntry = 1
            for j in range(2):
                if j < nEntry:
                    lb = "%f (%f, %f, %f)" % (angles[j], normals[j][0], normals[j][1], normals[j][2])
                    self.intersectionTable.setItem(i, j+1, qt.QTableWidgetItem(lb))
                else:
                    self.intersectionTable.setItem(i, j+1, qt.QTableWidgetItem("--"))

            self.intersectionTable.setItem(i, 3, qt.QTableWidgetItem("%f" % self.totalLengthInObject[i]))

    self.intersectionTable.show()


  def onReload(self,moduleName="PathCollisionAnalysis"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)

#
# PathCollisionAnalysisLogic
#

class PathCollisionAnalysisLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def CheckIntersections(self, inputModelHierarchyNode, inputFiducialNode):
    """
    Look for intersections between the path (inputFiducialNode) and the models under the hierarchy node.
    Return (objectIDs, objectNames, normalVectors, entryAngles, curvatures, radiusVectors, totalLengthInObject),
    where the elements are the arrays of object IDs, names, normal vectors at the entry points, angles between the path
    and the normal vectors, curvatures at the entry points, normal vectors to the center, and the total length of the path
    in the object.
    """
    if inputModelHierarchyNode == None:
      return None
    if inputFiducialNode == None:
      return None

    nOfModels = inputModelHierarchyNode.GetNumberOfChildrenNodes()
    objectIDs = []
    objectNames = []
    entryAngles = []
    normalVectors = []
    totalLengthInObject = []

    for i in range(nOfModels):
        chnode = inputModelHierarchyNode.GetNthChildNode(i)
        if chnode == None:
            continue
        mnode = chnode.GetAssociatedNode()
        if mnode == None:
            continue

        name = mnode.GetName()
        objectPoly = mnode.GetPolyData()

        if objectPoly == None:
            continue

        triangle = vtk.vtkTriangleFilter()
        triangle.SetInputData(objectPoly)
        triangle.Update()
        objectTrianglePoly = triangle.GetOutput()

        # print "Processing object: %s" % name

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
        enclosed.SetSurfaceData(objectTrianglePoly)
        enclosed.SetTolerance(0.0001) # Very important to get consistent result.
        enclosed.Update()

        lengthInObject = 0.0

        isInside = False

        angles = []
        normals = []
        curvatures = []
        radiusNormals = []

        surfaceNormals = vtk.vtkPolyDataNormals()
        surfaceNormals.SetInputData(objectTrianglePoly)
        surfaceNormals.ComputeCellNormalsOn()
        surfaceNormals.Update()
        surfaceNormalsOutput = surfaceNormals.GetOutput()

        # extract the cell data
        surfaceNormalsCellData = surfaceNormalsOutput.GetCellData();
        sNormals = surfaceNormalsCellData.GetNormals()

        for j in range(nFiducials-1):

            inputFiducialNode.GetNthFiducialPosition(j, pos0)
            inputFiducialNode.GetNthFiducialPosition(j+1, pos1)

            isInside0 = enclosed.IsInside(j)
            isInside1 = enclosed.IsInside(j+1)

            ## For debug
            #print "Point %d: from (%f, %f, %f) (%d) to (%f, %f, %f) (%d)" % (j, pos0[0], pos0[1], pos0[2], isInside0, pos1[0], pos1[1], pos1[2], isInside1)

            # A vector that represents the trajectory between pos0 and pos1
            # The orientation will be adjusted later to direct to the outside of the object.
            trajVec = np.array(pos1) - np.array(pos0)

            if isInside0 and isInside1:
                ## Both in the object
                lSegment = np.linalg.norm(trajVec)
                lengthInObject = lengthInObject + lSegment
                isInside = True

            intersectingPoint = [0.0]*3

            if isInside0 != isInside1:
                ## Potential intersection
                bspTree = vtk.vtkModifiedBSPTree()
                bspTree.SetDataSet(objectTrianglePoly)
                bspTree.BuildLocator()
                tolerance = 0.0001
                pCoord = [0.0]*3
                t = vtk.mutable(0)
                subID = vtk.mutable(0)
                cellID = vtk.mutable(0)
                fIntersect = bspTree.IntersectWithLine(pos0, pos1, tolerance, t, intersectingPoint, pCoord, subID, cellID)
                # idList = vtk.vtkIdList()
                # intersectingPoints = vtk.vtkPoints()
                # fIntersect = bspTree.IntersectWithLine(pos0, pos1, tolerance, intersectingPoints, idList)

                if fIntersect == 0:
                    ## If this happens, consider smaller tolerance
                    continue

                isInside = True

                # Get intersecting point and measure the length inside the boject
                # intersectingPoints.GetPoint(0, intersectingPoint)
                if isInside0:
                    segmentInObject = np.array(intersectingPoint) - np.array(pos0)
                    lengthInObject = lengthInObject + np.linalg.norm(segmentInObject)
                elif isInside1:
                    segmentInObject = np.array(intersectingPoint) - np.array(pos1)
                    lengthInObject = lengthInObject + np.linalg.norm(segmentInObject)
                    trajVec = - trajVec

                # cellID = idList.GetId(0)
                cell = objectTrianglePoly.GetCell(cellID)
                if cell == None:
                    continue

                cellNormal = [0.0]*3
                sNormals.GetTuple(cellID, cellNormal)

                # Check cell type?
                # if cell.GetCellType() == vtk.VTK_TRIANGLE:
                #     subID = 0
                # elif cell.GetCellType() == vtk.VTK_TRIANGLE_STRIP:
                #     print "Triangle Strip"

                # # Get subID -- no need since the cells have already converted to triangles
                # cell.IntersectWithLine(pos0, pos1, tolerance, t, intersectingPoint, pCoord, subID)
                points = cell.GetPoints()
                if points == None:
                    print "continue 4"
                    continue
                p0 = [0.0] * 3
                p1 = [0.0] * 3
                p2 = [0.0] * 3

                # Get point (when subID is used)
                # points.GetPoint(subID + 0, p0)
                # points.GetPoint(subID + 1, p1)
                # points.GetPoint(subID + 2, p2)

                points.GetPoint(0, p0)
                points.GetPoint(1, p1)
                points.GetPoint(2, p2)

                # print (intersectingPoint, p0, p1, p2)
                # npap0 = np.array(p0)
                # npap1 = np.array(p1)
                # npap2 = np.array(p2)
                # v0 = npap1-npap0
                # v1 = npap2-npap0
                # v = np.cross(v0, v1)
                # norm = np.linalg.norm(v,ord=1)
                # normVec = v / norm
                # print "Normal = (%f, %f, %f) / (%f, %f, %f)" %  (normVec[0], normVec[1], normVec[2], cellNormal[0], cellNormal[1], cellNormal[2])

                # Compute average normal
                #clippedModel = clip.GetOutput()
                # cellsNormal = clippedModel.GetCell(cellID).GetPointData().GetNormals()
                #
                # averageNormal = [0.0, 0.0, 0.0]
                # nOfNormals = 0;
                #
                # for cellIndex in range(0, cellsNormal.GetNumberOfTuples()):
                #     cellNormal = [0.0, 0.0, 0.0]
                #     cellsNormal.GetTuple(cellIndex, cellNormal)
                #
                # if not(math.isnan(cellNormal[0]) or math.isnan(cellNormal[1]) or math.isnan(cellNormal[2])):
                #     averageNormal[0] = averageNormal[0] + cellNormal[0]
                #     averageNormal[1] = averageNormal[1] + cellNormal[1]
                #     averageNormal[2] = averageNormal[2] + cellNormal[2]
                #     nOfNormals = nOfNormals + 1


                # Calculate the entry angle. Entry angle is zero, when the trajectory is perpendicular
                # to the surface.
                # angle = vtk.vtkMath.AngleBetweenVectors(normVec, trajVec) * 180.0 / math.pi
                angle = vtk.vtkMath.AngleBetweenVectors(cellNormal, trajVec) * 180.0 / math.pi
                angles.append(angle)
                normals.append(cellNormal)

                trajNorm = np.linalg.norm(trajVec, ord=1)
                nTrajVec = trajVec
                if trajNorm > 0.0:
                    nTrajVec = trajVec / trajNorm

                # # Caluclate the entry vector defined by: <v_e> = <n_t> - <n>
                # # where <n_t> and <n> are the trajectory vector and the normal.
                # entryVec = nTrajVec - normVec
                #print "  -- Intersecting at (%f, %f, %f) with angle %f and normal vector (%f, %f, %f)" % (p[0], p[1], p[2], angle, normVec[0], normVec[1], normVec[2])

        if isInside:
            objectIDs.append(i)
            objectNames.append(name)
            normalVectors.append(normals)
            entryAngles.append(angles)
            totalLengthInObject.append(lengthInObject)

    return (objectIDs, objectNames, normalVectors, entryAngles, totalLengthInObject)

  def ComputeCurvature(self, p1, p2, p3):
    # Given a curve that runs trhough points p1, p2, and p3, the function calculates
    # the curvature at p2, and the normal vector to the curve center.
    # p1, p2, and p3 must be Numpy arrays

    # Curvature
    v12 = p2 - p1
    v23 = p3 - p2
    pT  = v12 / numpy.linalg.norm(v12)
    ds = numpy.linalg.norm(v23)
    T  = v23 / ds
    kappa = numpy.linalg.norm(T-pT) / ds

    # Normal vector
    v13 = p3 - p1
    n13 = v13 / numpy.linalg.norm(v13)
    pr12 = numpy.inner(v12, n13) * n13 # Projection of v12 on n13
    vc = pr12 - v12
    nvc = vc / numpy.linalg.norm(vc)

    return (kappa, vc)


  def ComputeNormal(self, poly, center, radius):
    pass


class PathCollisionAnalysisTest(ScriptedLoadableModuleTest):
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
    self.test_PathCollisionAnalysis1()

  def test_PathCollisionAnalysis1(self):
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
    logic = PathCollisionAnalysisLogic()
    self.delayDisplay('Test passed!')
