Paraview macros
========

Python scripts for automatization of Paraview chores.

Shortcuts for macros are especially useful, e.g., use `Command+R` for `R.py`.

Installation
-------

For MacOS, create symlinks to these files:
```shell
mkdir -p ~/.config/ParaView/Macros
for f in $(ls *.py); do ln -s `pwd`/$f ~/.config/ParaView/Macros; done
``
