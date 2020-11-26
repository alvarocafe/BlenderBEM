function divnode(X,t)
    # Realiza divisao binária dos nós da malha;
    # X = matrix que contêm as coordenadas dos nós;{t,2}
    # t = vetor com os números dos nós; [t]
    n = length(t)   # Quantidade de nós
    x = X[t,:] # Matrix com as coordenadas dos nós que pertencem ao bloco a ser dividido; {t,2}
    c = mean(x,1)   # Vetor com as coordenadas do centro geometrico do conjunto de nós;{1,2}
    #mean(x,1) = média ao longo da 1 dimensão da matrix (1dim = linhas).
    covx = cov(x)                           # Calcula matrix de covariancia de x
    eig_valx,eig_vecx  = eig(covx)          # Calucla autovalores e autovetores da matriz de covariancia de x
    ref = eig_vecx[:,indmax(eig_valx)]      # Define como referencia o autovetor relacionado ao maior autovalor
    # a direcao desse autovetor e a direcao e maior cvariabilidade dos dados
    attcond = zeros(n)
        for i=1:n
            attcond[i] = (x.-c)[i,:]'*ref  # Condicao que divide os nos em dois blococ diferentes.
        end
    x1 = t[attcond.>=0]         # Bloco tal que a condicao e >= 0
    x2 = t[attcond.<0]          # Bloco tal que a condicao e < 0
    diam = 2*maximum(sqrt.(((x.-c).*(x.-c))[:,1]+((x.-c).*(x.-c))[:,2]))
    # Calcula diametro do conjunto de dados, centralizando eles; dia = 2*norma do ponto mais distante do centro
    return x1,x2,diam,c
end

function cluster(X, max_elem,η = 1.0)
    # X = Coordenadas (x,y) dos nós
    # max_elem = Define máximo de nós em cada folha, tal que: max_elem/2 =< nós em cada folha < max_elem

    m,n = size(X)                     # Tamanho da matriz contendo as coordernadas de cada nós {m,2}
    max_clt = ceil(Int,2*m/max_elem)  # Define limite superior para tamanho (nº de linhas) das matrizes e vetores utilizados
    child1 = zeros(1,2*max_clt)
    child2 = zeros(1,2*max_clt)
    t = collect(1:m)                # Nós
    inode = 1                       # Começa a contagem de nós da árvore
    ileaf = 1                       # Começa a contagem de folhas da árvore
    nodes = Array{Any}(2*max_clt)   # Aloca um vetor de vetores para guardar os nós da malha pertencentes a cada nó da árvore
    # Nodes[i] = vetor com os nós da malha pertencentes ao nó i da árvore
    leaves = Array{Any}(2*max_clt)  # Aloca um vetor para guardar as folhas
    child = zeros(Int,2*max_clt,2)  # Aloca uma matriz para guardar os filhos de cada nó.
    # Child[i,:] = filhos do nó i
    nodes[1] = t                     # O 1º nó da árvore (raiz) contem todos os nós da malha.
    center_row = zeros(2*max_clt,2)  # Aloca um vetor para guardar o centro geometrico de cada bloco de nós da malha.
    diam = zeros(2*max_clt)          # Aloca um vetor ppara guardar o diametro de cada bloco de nós da malha.
    i = 1                            # Começa contagem de nós
    while inode >= ileaf             # Enquanto o quantidade de nós for maior que a de folhas.
    # Observe que a condição só não vai ser satisfeita quando a árvore estiver completa.
        t1,t2,d,c = divnode(X,nodes[i])      # Executa a rotina que divide os nós da malha.
        center_row[i,:] = c;                 # Salva centro geometrico do nó i da árvore
        diam[i] = d;                         # Salva diametro do nó i da árvore
            if length(t1)> max_elem          # Se a quantidade de nós em t1 for maior que max_elem, definir como nó
                inode = inode + 1            # Chama proximo nó
                nodes[inode] = t1            # Define t1 como um nó
                child[i,1] = inode           # Define t1 como filho do nó i
            else                             # Se a quantidade de nós for menor que max_elem, definir como folha
                leaves[ileaf] = t1           # Define t1 como uma folha
                ileaf = ileaf + 1            # Chama proxima folha
                child1[i] = ileaf            # Define t1 como folha do nó i
            end
            # Realiza o mesmo para t2---------------------------------------
            if length(t2) > max_elem
                inode = inode + 1
                nodes[inode] = t2
                child[i,2] = inode
            else
                leaves[ileaf] = t2
                ileaf = ileaf + 1
                child2[i] = ileaf
            end
            # --------------------------------------------------------------
        i = i + 1
    end

    Tree = Array{Any}(inode+ileaf-1)    # Define tamanho da árvore
    for i=1:inode                       # Para todos os nós
        Tree[i] = nodes[i]              # Insere os nós na árvore
        if child1[i] > 0                # Se aquele nó tem folhas
            child[i,1] = child1[i] + inode - 1   # Adiciona as folhas pares na matriz child
        end
        if child2[i] > 0                         # Se aquele nó tem folhas
            child[i,2] = child2[i] + inode - 1   # Adiciona as folhas impares na matriz child
        end
    end

    for i=1:ileaf-1     # Para todos as folhas
        Tree[inode+i] = leaves[i]   # Insere as folhas na árvore
        t1,t2,d,c = divnode(X,leaves[i]) # Calcula o diam e c das folhas
        # Havia sido calculado somente os dos bloco que foram divididos, ou seja, os nós.
        center_row[inode+i,:] = c   # Adicona o c das folhas na matrixc center_row
        diam[inode+i] = d           # Adicona diam das folhas na matrix diam
        child[inode+i,:] = [0 0]    # Completa a matrix child com pares [0,0], pois folhas nao tem filhos
    end

    admiss = zeros(inode+ileaf-1,inode+ileaf-1)
    # Cria matriz para alocar o resultado da aplicadao da condicao de admissiblidade entre os blocos
    for i=1:inode+ileaf-1     # Para todos os nós da malha
        for j=1:inode+ileaf-1
            admiss[i,j] = η*norm(center_row[i]-center_row[j],2) - max(diam[i],diam[j])
            # Condicao de adimissiblidade, para satisfazer deve ser > 0
        end
    end
    allow = admiss.>=0; # Salva blocos onde a condicao e satisfeita
    block = blocks(Tree,child,allow) # Funcao que retorna os blocos admissiveis
    return  Tree,block
end

function blocks(Tree,child,allow)
    fc1 = [2; 2; 3; 3]   # Primeiros Blocos a serem avaliados
    fc2 = [2; 3; 2; 3]   # Primeiros Blocos a serem avaliados
    # fc1(1) e fc(2) formam blocos a serem analisados -> (22, 23, 32, 33)
    block = zeros(Int,0,3)
    # Matrix que aloca os blocos admissiveis [:,1:2]
    # e se atende a condicao de admissiblidade [:,3]
    c1 = 0;     # Contador
    while c1 < length(fc1)/2
        for i=1:2
            if allow[fc1[c1*2+i],fc2[c1*2+i]]==1    # Se blocos são admissiveis
                block = vcat(block,[fc1[c1*2+i] fc2[c1*2+i] 1])
                # Adicionar blocos e identificador 1 (admissivel) a proxima linha matrix block
            else # Se blocos não são admissiveis
                if child[fc1[c1*2+i],1]==0 && child[fc2[c1*2+i],1]==0
                    # Se ambos os blocos não tem filhos, ou seja, se ambos sao folhas
                    block = vcat(block,[fc1[c1*2+i] fc2[c1*2+i] 0])
                    # Adicionar blocos e identificador 0 (não admissivel) a proxima linha matrix block
                else
                    if length(Tree[fc1[c1*2+i]])>=length(Tree[fc2[c1*2+i]])
                        # Se a quantidade de elementos no bloco Tree[fc1[...]]] e >=  Tree[fc2[...]]
                        fc1 = [fc1; child[fc1[c1*2+i],:]]        # Adiciona filhos a fc1[...]
                        fc2 = [fc2; fc2[c1*2+i]; fc2[c1*2+i]]    # Repete elemento de fc2[...]
                    else
                        fc1 = [fc1; fc1[c1*2+i]; fc1[c1*2+i]]    # Repete elemento de fc1[...]
                        fc2 = [fc2; child[fc2[c1*2+i],:]]        # Adiciona filhos a fc2[...]
                    end
                end
            end
        end
        c1 = c1 + 1  # Atualiza contador
    end
    return block  #Matriz que contem os blocos analisados e sua condicao de admiiblidade
end
