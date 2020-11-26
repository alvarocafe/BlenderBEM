function mostra_geo(NOS_GEO,ELEM)
# geometria da estrutura

nelem=size(ELEM,1);
xmin = mini(NOS_GEO[:,2]);
xmax = maxi(NOS_GEO[:,2]);
ymin = mini(NOS_GEO[:,3]);
ymax = maxi(NOS_GEO[:,3]);
zmin = mini(NOS_GEO[:,4]);
zmax = maxi(NOS_GEO[:,4]);
NOS=zeros(nelem,4);

for i = 1:nelem
    no1=ELEM[i,2];
    no2=ELEM[i,3];
    no3=ELEM[i,4];
    no4=ELEM[i,5];

    xno1=NOS_GEO[no1,2];
    yno1=NOS_GEO[no1,3];
    zno1=NOS_GEO[no1,4];

    xno2=NOS_GEO[no2,2];
    yno2=NOS_GEO[no2,3];
    zno2=NOS_GEO[no2,4];

    xno3=NOS_GEO[no3,2];
    yno3=NOS_GEO[no3,3];
    zno3=NOS_GEO[no3,4];

    xno4=NOS_GEO[no4,2];
    yno4=NOS_GEO[no4,3];
    zno4=NOS_GEO[no4,4];

    xnos=[xno1,xno2,xno3,xno4];
    ynos=[yno1,yno2,yno3,yno4];
    znos=[zno1,zno2,zno3,zno4];
    #xlabel ('{\it x}')
    #ylabel ('{\it y}')
    #zlabel ('{\it z}')
    plot3D(xnos,ynos,znos,linestyle="none",marker="o");	#Parte do programa que mostra geometria
    #axis equal
    #axis([xmin xmax ymin ymax zmin zmax])
    #hold on;
    xno=sum(xnos)/4;
    yno=sum(ynos)/4;
    zno=sum(znos)/4;
    #text(xno,yno,zno,num2str(i))
    NOS[i,:]=[i,xno,yno,zno];
end;
#view(3)
#hold off;

return NOS
end

