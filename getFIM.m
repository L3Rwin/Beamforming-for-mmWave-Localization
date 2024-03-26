function J_ = getFIM(X, T, h, Ts, H_path, AOD, AOA, G, N, Nr, Nt, power_clk,power_noise)
% calculate FIM J and transform it to J~ via getTmat function

% calculate derivative of H with respect to eta.
D=zeros(Nr,Nt,5*G,N);
for n=1:N
    %derivative with respect to AOD 
    for g=1:G
        D(:,:,g,n)=H_path(:,:,g,n)*cos(AOD(g))*(-diag((0:Nt-1))*1j*pi);
    end

    %derivative with respect to AOA
    for g=1:G
        D(:,:,G+g,n)=(-diag((0:Nr-1))*1j*pi)*cos(AOA(g))*H_path(:,:,g,n);
    end

    %derivative with respect to h
    for g=1:G
        D(:,:,2*G+g,n)=H_path(:,:,g,n)/h(g);
        D(:,:,3*G+g,n)=1j*H_path(:,:,g,n)/h(g);
    end
    
    %derivative with respect to tao
    for g=1:G
        D(:,:,4*G+g,n)=H_path(:,:,g,n)*(-1j*2*pi*(n-1)/(N*Ts));
    end
end
J=zeros(5*G,5*G);
for i=1:5*G
    for j=1:5*G
            sum=0;
        for n=1:N
            sum=sum+real(trace(X*conj(D(:,:,i,n)).'*D(:,:,j,n)));
        end
    end
    J(i,j)=sum*2/power_noise;
end

Jprior=zeros(4*G+2,4*G+2);
Jprior(4*G+2,4*G+2)=(1/power_clk);
J_=T*J*T'+Jprior;
end