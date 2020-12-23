# Reload the original source

from paraview.simple import *

# Find the original source and reload it
for key, value in GetSources().items():
    if key[0].endswith('foam'):
        print("Reload " + key[0]);
        ReloadFiles(value[0])

# Go to the last time step
GetAnimationScene().GoToLast()

