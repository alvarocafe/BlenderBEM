# BlenderBEM
This is a Blender addon to solve boundary element method equations using Python and Julia solvers.

## Install
First, you'll need Blender: https://www.blender.org/download/

Then you can already use BlenderBEM by downloading this repository and running the script BlenderBEML.py. This script uses a Cython implementation of the BEM for Laplace problems.

If you wish to run BlenderBEM for Helmholtz problems, you'll need to install Julia: https://julialang.org/

Then, add Julia to the PATH. In GNU/Linux, this may be achieved by making a symbolic link to a folder that is in your PATH.

`sudo ln -s /path/to/your/julia-1.5.3/bin/ /usr/bin/julia`

You'll need to install and build PyCall in Julia:

`]add PyCall`

`]build PyCall`

Install pip in the Blender Python binary by downloading the get_pip.py script from: https://github.com/pypa/get-pip

And run it with:

`/path/to/your/blender-2.91.0-linux64/2.91/python/bin/python3.7m get_pip.py`

Now, using the Python binary from Blender, you'll install PyJulia:

`/path/to/your/blender-2.91.0-linux64/2.91/python/bin/python3.7m -m pip install julia`

Open the iterative console of the Blender Python:

`/path/to/your/blender-2.91.0-linux64/2.91/python/bin/python3.7m`

Import Julia and install it:

`import julia`

`julia.install()`

You can now run the BlenderBEM.py script in Blender to use the Julia solver. 


## How to use

First, add a mesh object in Blender and save the project in the same folder as the BlenderBEM script. Select the object and click "Prepare Mesh" to start the method and let BlenderBEM know which object will be used for the analysis.

Now, change to Object Mode by pressing the Tab key and select some polygons with the mouse. Use the Shift key together with the left mouse button or press b to box select. Press the "Submit polys to potential" to submit those polygons to the boundary condition potential specified.

When the boundary conditions are all set, go back to Object Mode by pressing the Tab key and press the "Run Laplace" button to solve the problem. 

A demo of the tool can be seen here:

[![Elementos de contorno no Blender!! (BlenderBEM)](https://img.youtube.com/vi/WVS4Ix-wXA8/0.jpg)](https://www.youtube.com/watch?v=WVS4Ix-wXA8)
