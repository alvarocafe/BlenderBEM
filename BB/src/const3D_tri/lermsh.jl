function lermsh(nome,dim=2,CCSeg=false)
    dados=DelimitedFiles.readdlm(nome)
    index=zeros(Int,2);idenx=findall(dados.=="\$Nodes");index[1]=idenx[1][1];index[2]=idenx[1][2]
    if(isempty(index))
        index=findall(dados.=="\$ParametricNodes")
    end
    npontos=dados[index[1][1]+1,1]
    pontos=Array{Float64}(dados[index[1][1]+2:index[1][1]+npontos+1,2:4])
    index[1,1]=index[1,1]+4+npontos
    nelem=dados[index[1][1],1]
    elemtipo=dados[index[1][1]+1:index[1][1]+nelem,2]
    elemint=Array{Int}[]
    if dim==2
        if sum(elemtipo.==8)!= 0 # 3-node second order line
            elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==8),[6:8;5]])
            # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==8),4])
        elseif sum(elemtipo.==1)!= 0 # 2-node line.
            elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==1),[6:7;5]])
            # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==1),4])
        end

    elseif dim==3
        if sum(elemtipo.==2)!= 0 #3-node triangle.
            elemcont=Array{Int}(dados[index[1][1].+findall(elemtipo.==2),[6:8;5]])
            elemnum=Array{Int}(dados[index[1][1].+findall(elemtipo.==2),1])
            # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==2),4])
        elseif sum(elemtipo.==9)!= 0 #6-node second order triangle
            elemcont=Array{Int}(dados[index[1][1].+findall(elemtipo.==9),[6:11;5]])
            # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==9),4])
        end
        if sum(elemtipo.==4)!= 0 #4-node tetrahedron.
            elemcont=Array{Int}(dados[index[1][1].+findall(elemtipo.==4),[6:9;5]])
            elemint=Array{Int}(dados[index[1][1].+findall(elemtipo.==4),[6:9;5]])
        elseif sum(elemtipo.==11)!= 0 #10-node second order tetrahedron
            elemcont=Array{Int}(dados[index[1][1].+findall(elemtipo.==11),[6:15;5]])
            elemint=Array{Int}(dados[index[1][1].+findall(elemtipo.==11),[6:15;5]])
        elseif sum(elemtipo.==3)!= 0 #4-node quadrangle
            elemcont=Array{Int}(dados[index[1][1].+findall(elemtipo.==3),[6:9;5]])    	
        end

    end
    if CCSeg==false
        # return pontos,elemcont,elemint
        return [1:size(pontos,1) pontos],[1:size(elemcont,1) elemcont],[1:size(elemint,1) elemint],elemtipo


    else
        CDC=Array{Float64}(elemcont)
        for i=1:size(CCSeg,1)
            index=find(seg.==i)
            CDC[index,:]=CCSeg[i,:]
        end
        return pontos,elemcont,elemint,CDC
    end
end

function salva(nome,T,elemtipo,tipo=2)
    dados=readdlm(nome)
    index=findfirst(dados[:,1].=="\$ElementData")
    if isempty(index)
        open(nome, "a") do f
            writedlm(f,["\$ElementData"
                        1
                        "Temperatura"
                        1
                        0.0
                        3
                        0
                        1
                        size(T,1)])
            writedlm(f,[(1:size(T,1))+findfirst(elemtipo,tipo)-1 T])
            write(f,"\$EndElementData")
        end
    else
        open(nome, "w") do f
            writedlm(f,dados[1:index-1,:])
            writedlm(f,["\$ElementData"
                        1
                        "Temperatura"
                        1
                        0.0
                        3
                        0
                        1
                        size(T,1)])
            writedlm(f,[string.(collect((1:size(T,1))+findfirst(elemtipo,tipo)-1)) T])
            write(f,"\$EndElementData")
        end
    end
end
