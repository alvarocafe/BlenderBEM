function ACAF(Tree,block,fHeG,arg,erro=1e-5)
#arg = [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC]
  n=size(block,1)
  Aaca=Array{Any}(n,2)
  #Baca=Array{Any}(n,2)
  b=zeros(size(arg[1],1))
  for i=1:n
    b1=Tree[block[i,1]]
    b2=Tree[block[i,2]]
    if block[i,3]==0
      # if 0==0
      Aaca[i,1],B=fHeG(b1,b2,arg)
      b[b1]+=B*arg[7][b2,3]
    else
      # println(i)

      INDB1=[]
      INDB2=[]
      B1=zeros(0,length(b2))
      B2=zeros(length(b1),0)

      ind1=trues(length(b1))
      ind2=trues(length(b2))

      #linha 1
      indaref=1
      aref=0*ind1
      for indaref =1 :length(ind1)
        aref,btemp=fHeG(b1[indaref],b2,arg)
        aref=aref'
        push!(INDB1,indaref)
        B1=[B1;btemp]
        if norm(aref)>1e-10
          break
        end
        ind1[indaref]=0
      end
      if ind1==falses(ind1)
        Aaca[i,1]=ind1*0
        Aaca[i,2]=ind2'*0
      else
        #coluna 1
        indbref=1
        bref=0*ind2

        for ii =1 :length(ind2)
          indbref=indmin(abs.(aref[ind2]))
          indbref+=cumsum(ind2.==0)[ind2][indbref]
          bref,btemp=fHeG(b1,b2[indbref],arg)
          push!(INDB2,indbref)
          B2=[B2 btemp]
          if norm(bref)>1e-10
            break
          end
          ind2[indbref]=0
        end

        arefmax=indmax(abs.(aref))
        brefmax=indmax(abs.(bref))
        Umax=0
        Vmax=0
        nmin=min(length(ind1),length(ind2))
        U=zeros(length(ind1),nmin)
        V=zeros(nmin,length(ind2))
        Ap=zeros(length(ind1),length(ind2))
        norma0=0.0
        cont=0
        for cont=1:nmin-1
          if abs.(aref[arefmax])>abs.(bref[brefmax])
            U[:,cont],btemp=fHeG(b1,b2[arefmax],arg)
            U[:,cont]=U[:,cont]-Ap[:,arefmax]
            push!(INDB2,arefmax)
            B2=[B2 btemp]
            Umax=indmax(abs.(U[:,cont]))
            V[cont,:],btemp=fHeG(b1[Umax],b2,arg)
            V[cont,:]=(V[cont,:][:]-Ap[Umax,:])'/U[Umax,cont]
            push!(INDB1,Umax)
            B1=[B1;btemp]
            # V[cont,:]=(fHeG(b1[Umax],b2,arg)[1,:]-Ap[Umax,:])'/U[Umax,cont]
            Vmax=arefmax
          else
            V[cont,:],btemp=fHeG(b1[brefmax],b2,arg)
            V[cont,:]=(V[cont,:][:]-Ap[brefmax,:])
            push!(INDB1,brefmax)
            B1=[B1;btemp]
            Vmax=indmax(abs.(V[cont,:]))
            U[:,cont],btemp=fHeG(b1,b2[Vmax],arg)
            U[:,cont]=(U[:,cont]-Ap[:,Vmax])/V[cont,Vmax]
            push!(INDB2,Vmax)
            B2=[B2 btemp]
            Umax=brefmax
          end
          Ap=U*V
          norma1=vecnorm(Ap)
          if abs.((norma1-norma0)/norma1) < erro
            break
          else
            norma0=norma1
          end
          ind1[Umax]=0
          ind2[Vmax]=0
          if indaref==Umax && indbref==Vmax

            for indaref =1 :sum(ind1)
              indaref=indmax(ind1)
              aref,btemp=fHeG(b1[indaref],b2,arg)
              aref=aref-Ap[indaref,:]'
              push!(INDB1,indaref)
              B1=[B1;btemp]
              if norm(aref)>1e-10
                break
              end
              ind1[indaref]=0
            end



            for indbref =1 :sum(ind2)
              indbref=indmin(abs.(aref[ind2]))
              indbref=indbref+cumsum(ind2.==0)[ind2][indbref]
              bref,btemp=fHeG(b1,b2[indbref],arg)
              bref=bref-Ap[:,indbref]

              push!(INDB2,indbref)
              B2=[B2 btemp]
              if norm(bref)>1e-10
                break
              end
              ind2[indbref]=0
            end

          elseif indaref==Umax
            bref=bref-U[:,cont]*V[cont,indbref]'
            #  indaref=indmin(abs.(bref[ind1]))
            #  indaref=indaref+cumsum(ind1.==0)[ind1][indaref]
            #  aref=fHeG(b1[indaref],b2)'-Ap[indaref,:]'

            for indaref =1 :sum(ind1)
              indaref=indmin(abs.(bref[ind1]))
              indaref=indaref+cumsum(ind1.==0)[ind1][indaref]
              aref,btemp=fHeG(b1[indaref],b2,arg)
              aref=aref[:]-Ap[indaref,:][:]
              push!(INDB1,indaref)
              B1=[B1;btemp]
              if norm(aref)>1e-10
                break
              end
              ind1[indaref]=0
            end


          elseif indbref==Vmax
            aref=aref[:]-U[indaref,cont]*V[cont,:]
            # indbref=indmin(abs.(aref[ind2]))
            # indbref=indbref+cumsum(ind2.==0)[ind2][indbref]
            # bref=fHeG(b1,b2[indbref])-Ap[:,indbref]

            for indbref =1 :sum(ind2)
              indbref=indmin(abs.(aref[ind2]))
              indbref=indbref+cumsum(ind2.==0)[ind2][indbref]
              bref,btemp=fHeG(b1,b2[indbref],arg)
              bref=bref-Ap[:,indbref]

              push!(INDB2,indbref)
              B2=[B2 btemp]
              if norm(bref)>1e-10
                break
              end
              ind2[indbref]=0
            end
          else
            aref=aref[:]-U[indaref,cont]*V[cont,:]
            bref=bref-U[:,cont]*V[cont,indbref]'
          end
          arefmax=indmax(abs.(aref))
          brefmax=indmax(abs.(bref))
        end
        Aaca[i,1]=U[:,1:cont]
        Aaca[i,2]=V[1:cont,:]
        # Aaca[i,1]=U
        # Aaca[i,2]=V
      end
      max1=ind2sub(size(B1),indmax(abs.(B1[:,INDB2])))
      maxv=B1[max1[1],INDB2[max1[2]]]
      Vb=B1[max1[1],:]'
      Ub=B2[:,max1[2]]/maxv
      Bp=Ub*Vb
      B1-=Bp[INDB1,:]
      B2-=Bp[:,INDB2]
      for i=1:size(B1,1)-1
        max1=ind2sub(size(B1),indmax(abs.(B1[:,INDB2])))
        maxv=B1[max1[1],INDB2[max1[2]]]
        if abs.(maxv)<1e-12
          break
        end
        Vb=[Vb; B1[max1[1],:]']
        Ub=[Ub B2[:,max1[2]]/maxv]
        lastBp=Ub[:,end]*Vb[end,:]'
        B1-=lastBp[INDB1,:]
        B2-=lastBp[:,INDB2]

      end
      b[b1]+=Ub*Vb*arg[7][b2,3]

      # println(Bp)

    end

    # println(INDB2)

  end

  return Aaca,b
end
function erroblocos(hmat,A,block,Tree)
  for i =1:length(block[:,3])
    if block[i,3]==1
      n=vecnorm(A[Tree[block[i,1]],Tree[block[i,2]]]-hmat[i,1]*hmat[i,2])
      println("erro no bloco $i = $n")
    end
  end
end

function matvec(hmat,b,block,Tree)
  v=b*0
  for i =1:length(block[:,3])
    if block[i,3]==1
      v[Tree[block[i,1]]]+=hmat[i,1]*(hmat[i,2]*b[Tree[block[i,2]]])
    else
      v[Tree[block[i,1]]]+=hmat[i,1]*b[Tree[block[i,2]]]
    end
  end
return  v
end
function montacheia(hmat,block,Tree,n)
A=zeros(n,n)
  for i =1:length(block[:,3])
    if block[i,3]==1
      A[Tree[block[i,1]],Tree[block[i,2]]]=hmat[i,1]*hmat[i,2]
    else
      A[Tree[block[i,1]],Tree[block[i,2]]]=hmat[i,1]
    end
  end
  A
end
function tamanho(hmat,block,Tree)
A=0
  for i =1:length(block[:,3])
    if block[i,3]==1
      A+=length(hmat[i,1])+length(hmat[i,2])
    else
      A+=length(hmat[i,1])
    end
  end
return  A
end
