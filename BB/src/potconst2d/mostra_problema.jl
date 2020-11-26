using PyPlot
plt=PyPlot

function mostra_problema(ELEM,NOS_GEO,NOS)
nelem=size(ELEM,1);
figure();
ind=collect(1:nelem);
xmax=maximum(NOS_GEO[:,2])
ymax=maximum(NOS_GEO[:,3])
xmin=minimum(NOS_GEO[:,2])
ymin=minimum(NOS_GEO[:,3])
deltax=xmax-xmin
deltay=ymax-ymin
dmax=sqrt(deltax^2+deltay^2)
ax = plt.gca() # get current axes
plt.grid("on")
plt.plot(NOS[:,2],NOS[:,3],"ro",markersize=4);	# Plot the node of the elements
end

