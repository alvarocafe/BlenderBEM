# Matriz de transferência para cavidades acústicas concatenadas
# Autor: Alvaro Campos Ferreira - e-mail: alvaro.campos.ferreira@gmail.com
# Linguagem: Julia
function transfer_matrix(c,ω,l,s,bc)
    w=(ω/10):(ω/100):ω; #faixa de frequencia analisada rad/s
    n=size(l,2);    #numero de cavidades acopladas

    ga=zeros(n,size(w,1)); #gamma é o número de onda para cada valor da frequência

    for i=1:n
        ga[i,:]=w.*(l[1,i] ./c);
    end

    m=2;   #Numero de divisoes por elemento
    global b11, b12, b21,b22


    b11 = cos.(ga[n,:]);
    b12 = (c ./((complex(0,1).*s[1,n])).*sin.(ga[n,:]));
    b21 = (-(1 ./c).*(complex(0,1).*s[1,n]).*sin.(ga[n,:]));
    b22 = cos.(ga[n,:]);

    for i = n-1:-1:1
        a11 = cos.(ga[i,:]);
        a12 = (c/(complex(0,1)*s[1,i])).*sin.(ga[i,:]);
        a21 = -((complex(0,1).*s[1,i]) ./c).*sin.(ga[i,:]);
        a22 = copy(a11);
        a = b11 .*a11 + b12.*a21;
        b = b11 .*a12 + b12.*a22;
        c = b21 .*a11 + b22.*a21;
        d = b21 .*a12 + b22.*a22;
        b11 = a; b12=b; b21=c; b22=d;
    end

    bXX = b21;
    function findlocalmaxima(signal::Vector)
        inds = Int[]
        if length(signal)>1
            if signal[1]>signal[2]
                push!(inds,1)
            end
            for i=2:length(signal)-1
                if signal[i-1]<signal[i]>signal[i+1]
                    push!(inds,i)
                end
            end
            if signal[end]>signal[end-1]
                push!(inds,length(signal))
            end
        end
        inds
    end


    loc = findlocalmaxima(1.0 ./abs.(bXX));

    # Plota a solucao das matrizes de influencia
    #semilogy(w ./(2*pi),abs(1  ./imag(bXX)),w ./(2*pi),abs(1  ./real(bXX)));
    #grid on;
    #title('Matrix solution');
    #xlabel('frequency (Hz)');
    #ylabel('frequency response');

    # gap=complex(zeros(n,m*n));
    # p=complex(zeros(size(loc,1),n*m));
    # q=complex(zeros(size(loc,1),n*m));
    # lin = zeros(1,n*m);

    # if loc != 0
    #     for k=1:length(loc)
    #         p[k,1]=1;
    #         q[k,1]=0;
    #         vez=1;
    #         lin[vez]=0;
    #         for i=1:n
    #             for j=1:m-1
    #                 vez=vez+1;
    #                 gap[i,j]=(w[loc[k]] ./c).*(l[1,i] ./m);
    #                 lin[vez]=lin[vez-1]+l[1,i] ./m;
    #                 b11=cos.(gap[i,j]);
    #                 b12=(c ./((complex(0,1).*s[1,i])).*sin.(gap[i,j]));
    #                 b21=(-(1 ./c).*(complex(0,1).*s[1,i]).*sin.(gap[i,j]));
    #                 b22=cos.(gap[i,j]);
    #                 p[k,vez]= b11*p[k,vez-1] + b12*q[k,vez-1];
    #                 q[k,vez]= b21*p[k,vez-1] + b22*q[k,vez-1];
    #             end
    #         end
    #         #         Plota a forma modal
    #         #         figure(k+1)
    #         #         plot(lin,(p[k,:])/max(abs(p[k,:])));
    #         #         grid on;
    #         #         title('Modal shape');
    #         #         xlabel('distance [m]');
    #         #         ylabel('modal shape [-1,1]');
    #         #         pause
    #     end
    # end
    #        figure(10)
    # plot(lin,(p[k,:])); # ./max(abs(p[k,:])
    # grid on;
    # hold on;
    # title('Modal shape');
    # xlabel('distance [m]');
    # ylabel('modal shape [-1,1]');
    return loc, loc
end
