import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
import math

#
# 
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
    self.parent.contributors = ["Sarah Elbenna (ENSEEIHT), Junichi Tokuda (BWH)"] # replace with "Firstname Lastname (Organization)"
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
    angleFormLayout.addRow("Trajectory: ", self.inputFiducialSelector)


  

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
    logic.EntryAngle(self.inputModelSelector.currentNode(), self.inputFiducialSelector.currentNode())
    

  def onReload(self,moduleName="CurveAngle"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)


 
    
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



  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputLabelMapNode, inputFiducialNode):
    """Validates if the output is not the same as input
    """
    if not inputLabelMapNode:
      logging.debug('isValidInputOutputData failed: no input label map node defined')
      return False
    if not inputFiducialNode:
      logging.debug('isValidInputOutputData failed: no input fiducial  node defined')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def EntryAngle(self, inputModelNode, inputFiducialNode):
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

    for i in range(nOfModels):
      
      chnode = inputModelNode.GetNthChildNode(i)
      if chnode == None:
        continue
      
      mnode = chnode.GetAssociatedNode()
      if mnode == None:
        continue
      
      poly = mnode.GetPolyData()
      if poly == None:
        continue
      
      points = vtk.vtkPoints()
      idList = vtk.vtkIdList()
      pos0 = [0.0, 0.0, 0.0]
      posN = [0.0, 0.0, 0.0]
      
      inputFiducialNode.GetNthFiducialPosition(0, pos0)
      inputFiducialNode.GetNthFiducialPosition(12, posN)
      
      traj = [posN[0]- pos0[0], posN[1]-pos0[1], posN[2]-pos0[2]]
      inter1= [0.0, 0.0, 0.0]
      inter2= [0.0, 0.0, 0.0]
      bspTree = vtk.vtkModifiedBSPTree()
      bspTree.SetDataSet(poly)
      bspTree.BuildLocator()
      inputFiducialNode.GetNthFiducialPosition(0, inter1)
      inputFiducialNode.GetNthFiducialPosition(14, inter2)
      tolerance = 0.001
      bspTree.IntersectWithLine(inter1, inter2, tolerance, points, idList)

      if points.GetNumberOfPoints < 1:
        continue
      
      #points.GetPoint(0)

      if idList.GetNumberOfIds() < 1:
        continue
      
      cell0 = poly.GetCell(idList.GetId(0))
      if cell0 == None:
        continue
      
      p0 = cell0.GetPoints()
      if p0 == None:
        continue
      
      x0 = p0.GetPoint(1)[0]- p0.GetPoint(0)[0]
      y0 = p0.GetPoint(1)[1]- p0.GetPoint(0)[1]
      z0 = p0.GetPoint(1)[2]- p0.GetPoint(0)[2]
      v0 = [x0, y0, z0]
      x1 = p0.GetPoint(2)[0]- p0.GetPoint(0)[0]
      y1 = p0.GetPoint(2)[1]- p0.GetPoint(0)[1]
      z1 = p0.GetPoint(2)[2]- p0.GetPoint(0)[2]
      v1 = [x1, y1, z1]
      v0xv1 = np.cross(v0, v1)
      norm = math.sqrt(v0xv1[0]*v0xv1[0] + v0xv1[1]*v0xv1[1] + v0xv1[2]*v0xv1[2] )
      
      normal = (1/norm)*v0xv1
      angle = vtk.vtkMath.AngleBetweenVectors(normal, traj)
      
      print angle 

    

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
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
