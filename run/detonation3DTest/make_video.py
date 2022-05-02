#!/usr/bin/env -S pvpython --enable-bt

import argparse
from paraview.simple import *

str2pair = lambda s: [float(item) for item in s.split(':')]
str2pair2 = lambda s: [int(item) for item in s.split(':')]

parser = argparse.ArgumentParser(description='Make animation for detonation3DTest')
parser.add_argument('foamfile', type=str, help='file *.foam')
parser.add_argument('-f', '--font', type=int, default=18, help='font size')
parser.add_argument('-a', '--rate', type=int, default=5, help='frame rate')
parser.add_argument('-d', '--display-field', type=str, default='normalisedGradRho', help='field to be display')
parser.add_argument('-n', '--display-field-name', type=str, default='normalised gradient of density', help='name of the field to be display')
parser.add_argument('-r', '--display-field-range', type=str2pair, default='0.05:0.5', help='range of the field to be display')
parser.add_argument('-i', '--iso-field', type=str, default='lambda', help='field for isosurface')
parser.add_argument('-v', '--value', type=float, default=0.01, help='value for isosurface')
parser.add_argument('--dpi', type=str2pair2, default='600:600', help='resolution along the x:y axes')
parser.add_argument('-s', '--screenshot', action='store_true', help='make a screenshot instead')
args = parser.parse_args()

# 1. Read the case
reader = OpenFOAMReader(FileName=args.foamfile)
view = GetActiveViewOrCreate('RenderView')

# 2. Show the geometry outline to measure the dimensions
source = FindSource(args.foamfile)
sourceDisplay = Show(source, view, 'UnstructuredGridRepresentation')
sourceDisplay.SetRepresentationType('Outline')
GetActiveSource().CaseType = 'Decomposed Case'
xmin, xmax, ymin, ymax, zmin, zmax = GetActiveSource().GetDataInformation().GetBounds()
y0, z0 = (ymin + ymax)/2, (zmin + zmax)/2
print("Bounds:", xmin, xmax, ymin, ymax, zmin, zmax)
Hide(source, view)

# 3. Show a contour isosurface
contour = Contour(Input=source)
contour.ContourBy = ['POINTS', args.iso_field]
contour.Isosurfaces = [args.value]
contourDisplay = Show(contour, view, 'GeometryRepresentation')
contourDisplay.Representation = 'Surface'
ColorBy(contourDisplay, ('POINTS', args.display_field))
contourDisplay.SetScalarBarVisibility(view, True)

# 4. Configure the color bar
LUT = GetColorTransferFunction(args.display_field)
LUT.RescaleTransferFunction(*args.display_field_range)
LUTColorBar = GetScalarBar(LUT, view)
LUTColorBar.Orientation = 'Horizontal'
LUTColorBar.WindowLocation = 'Any Location'
LUTColorBar.Position = [0.15, 0.05]
LUTColorBar.ScalarBarLength = 0.7
LUTColorBar.Title = args.display_field_name
LUTColorBar.LabelFontSize = LUTColorBar.TitleFontSize = args.font
LUTColorBar.LabelFormat = LUTColorBar.RangeLabelFormat = '%-#6.1g'

# 5. Set the canvas size
layout = GetLayout()
layout.SetSize(*args.dpi)

# 6. Configure the View
view.InteractionMode = '2D'
view.CameraViewUp = [0, 0, 1]
view.CameraPosition = [1e10, y0, z0]
view.CameraFocalPoint = [0, y0, z0]
view.CameraParallelScale = max(y0, z0)
view.OrientationAxesVisibility = 0
view.Update()

# 7. Update the AnimationScene
animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()
animationScene.GoToLast()

print("We are ready!")

if args.screenshot:
    SaveScreenshot('screenshot.png', view)
else:
    SaveAnimation('video.avi', view, FrameRate=args.rate)
    print("Video is saved!")
