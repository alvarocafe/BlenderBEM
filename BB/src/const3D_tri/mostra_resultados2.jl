function mostra_resultados2(XYZ,tri,T)
nelem=size(tri,1)
x=zeros(3)
y=zeros(3)
z=zeros(3)
zc=zeros(nelem,1)
pc=[zeros(3,3)];
triang=zeros(3,3)
for elem =1:nelem
    no1=tri[elem,2]
    no2=tri[elem,3]
    no3=tri[elem,4]
    x[1]=XYZ[no1,2]
    y[1]=XYZ[no1,3]
    z[1]=XYZ[no1,4]
    x[2]=XYZ[no2,2]
    y[2]=XYZ[no2,3]
    z[2]=XYZ[no2,4]
    x[3]=XYZ[no3,2]
    y[3]=XYZ[no3,3]
    z[3]=XYZ[no3,4]
    triang=[[x[1] y[1] z[1]
    x[2] y[2] z[2]
    x[3] y[3] z[3]]]
    append!(pc,triang)
end
fig = plt.figure()
ax = mp.Axes3D(fig)
q = ar.Poly3DCollection(pc[2:end], linewidths=1)
ax[:add_collection3d](q)
m = cm.ScalarMappable(cmap=cm.jet)
b=m[:to_rgba](T[1:nelem])
q[:set_facecolor](b[:,1:3])
m[:set_array]([minimum(T),maximum(T)])
m[:set_clim](vmin=minimum(T),vmax=maximum(T))
plt.colorbar(m, orientation="vertical",shrink=0.9)
ax[:set_xlabel]("x")
ax[:set_ylabel]("y")
ax[:set_zlabel]("z")
ax[:set_xlim](minimum(XYZ[:,2]),maximum(XYZ[:,2]))
ax[:set_ylim](minimum(XYZ[:,3]),maximum(XYZ[:,3]))
ax[:set_zlim](minimum(XYZ[:,4]),maximum(XYZ[:,4]))
ax[:set_aspect]("equal")
ax[:view_init](elev=18., azim=43.)
return
end
