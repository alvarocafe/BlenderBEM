# Boundary element method base (BEM_base)
# Author: Álvaro Campos Ferreira
# Copyright (C) 2018 Álvaro Campos Ferreira
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
# BEM_base follows the following scheme:
#Start----Problem--------------Method-----------------post-processing
#--^You are here!----------------------------------------------------
# This is the start of the program, S. From here, you'll need to define the
#problem P which will be solved and the method M for building and solving
#the matricial system. The results then must be post-processed.
# This is the main program of the BEM_base. Add here your new implementation!
# Returns the potential and flux on boundary points and then on domain points.
# 2-dimensional elements.
include("src/const2D/const2D.jl")
include("src/wavenurbs2D/wavenurbs2D.jl")
#include("src/nurbs2D/nurbs2D.jl")
# 3-dimensional elements.
include("src/const3D_tri/const3D_tri.jl")
#include("src/waveconst3d/waveconst3d.jl")
include("src/const3D_tri/const3D_tri_POT.jl")
include("src/transconst2D/transconst2D.jl")
# Tests
include("tests/wave_tests.jl")
include("tests/pot_tests.jl")
using .const2D, .const3D_tri, .potconst3d, SpecialFunctions, KrylovMethods, Distributed
# Main function
function BEM_base(file,PONTOS_int=[],BCFace = [],k=1, equation = "wave")
    println("Importing mesh...")
    @time mshinfo = const3D_tri.lermsh(file,3) #Read the mesh generated 
    NOS_GEO,ELEM,elemint,CDC = mshinfo
    # Choose the right kernel for the element type of the mesh. 
    if equation == "wave"
	if size(ELEM,2) == 5
	    u,q,uint,qint = const3D_tri.solve(mshinfo,PONTOS_int,BCFace,k)
	else
	    u,q,uint,qint = const3D_quad.solve(mshinfo,PONTOS_int,BCFace,k)
	end
    elseif equation == "heat"
	if size(ELEM,2) == 5
	    u,q,uint,qint = const3D_tri_POT(mshinfo,PONTOS_int,BCFace,k)
	else
	    u,q,uint,qint = const3D_quad_POT(mshinfo,PONTOS_int,BCFace,k)
	end
    elseif typeof(equation) != String
	println("Error: the equation which will be solved must be specified.")
    end
    return u,q,uint,qint
end
