#!/usr/bin/env -S pvpython --enable-bt

import argparse
from paraview.simple import *

parser = argparse.ArgumentParser(description='Make animation from solidificationTest')
parser.add_argument('foamfile', type=str, help='file *.foam')
parser.add_argument('-g', '--gap', type=float, default=.2e-6, help='gap between pictures')
parser.add_argument('-z', '--zoom', type=float, default=1.4, help='camera zoom')
parser.add_argument('-f', '--font', type=int, default=18, help='font size')
parser.add_argument('-r', '--rate', type=int, default=5, help='frame rate')
parser.add_argument('-x', type=int, default=1200, help='resolution along the x axis')
parser.add_argument('-y', type=int, default=1200, help='resolution along the y axis')
parser.add_argument('-s', '--screenshot', action='store_true', help='make a screenshot instead')
args = parser.parse_args()

source = OpenFOAMReader(FileName=args.foamfile)
view = GetActiveViewOrCreate('RenderView')

# show source in view
sourceDisplay = Show(source, view, 'UnstructuredGridRepresentation')
ColorBy(sourceDisplay, ('POINTS', 'phase'))
sourceDisplay.RescaleTransferFunctionToDataRange(True, False)
sourceDisplay.SetScalarBarVisibility(view, True)

xmin, xmax, ymin, ymax, zmin, zmax = GetActiveSource().GetDataInformation().GetBounds()
print("Bounds:", xmin, xmax, ymin, ymax, zmin, zmax)

# create a new 'Transform'
transform = Transform(Input=source)
transform.Transform = 'Transform'
transform.Transform.Translate = [xmax + args.gap, 0.0, 0.0]
transformDisplay = Show(transform, view, 'UnstructuredGridRepresentation')
ColorBy(transformDisplay, ('POINTS', 'concentrationCr'))
transformDisplay.SetScalarBarVisibility(view, True)

# current camera placement for view
view.ViewSize = [args.x, args.y]
view.InteractionMode = '2D'
view.CameraViewUp = [1, 0, 0] # rotate clockwise 90 degrees
camFar = 1.0
bounds_cx = (xmin + xmax)/2
bounds_cy = (ymin + ymax)/2
bounds_cz = 0
pos = ymax - ymin
view.CameraPosition = [bounds_cx, bounds_cy, -pos*camFar]
view.CameraFocalPoint = [bounds_cx, bounds_cy, bounds_cz]
view.CenterAxesVisibility = 0
view.OrientationAxesVisibility = 0
view.CameraParallelScale = 0
view.ResetCamera()

LUT = GetColorTransferFunction('phase')
LUT.RescaleTransferFunction(0, 1.0)
LUTColorBar = GetScalarBar(LUT, view)
LUTColorBar.Position = [0.01, 0.05]
LUTColorBar.ScalarBarLength = 0.4
LUTColorBar.WindowLocation = 'AnyLocation'
LUTColorBar.LabelFontSize = LUTColorBar.TitleFontSize = args.font
LUTColorBar.LabelFormat = LUTColorBar.RangeLabelFormat = '%-#6.1g'

LUT = GetColorTransferFunction('concentrationCr')
LUTColorBar = GetScalarBar(LUT, view)
LUTColorBar.Position = [0.01, 0.55]
LUTColorBar.ScalarBarLength = 0.4
LUTColorBar.WindowLocation = 'AnyLocation'
LUTColorBar.LabelFontSize = LUTColorBar.TitleFontSize = args.font
LUTColorBar.LabelFormat = LUTColorBar.RangeLabelFormat = '%-#6.4g'
LUTColorBar.AutomaticLabelFormat = 0

cam = GetActiveCamera()
cam.Zoom(args.zoom)

animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()
animationScene.GoToLast()

# rescale color and/or opacity maps used to exactly fit the current source range
transformDisplay.RescaleTransferFunctionToDataRange(False, True)

print("We are ready!")

if args.screenshot:
    SaveScreenshot('screenshot.png', view)
else:
    SaveAnimation('video.avi', view, FrameRate=args.rate, Compression=0)
