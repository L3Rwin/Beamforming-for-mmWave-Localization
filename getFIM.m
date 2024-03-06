function J_ = getFIM(X, T, h, Ts, H_path, G, N, Nr, Nt, power_clk)
% calculate FIM J and transform it to J~ via getTmat function

% calculate derivative of H with respect to eta.
D=zeros(Nr,Nt,5*G,N);
for n=1:N
    %derivative with respect to AOD 
    for g=1:G
        D(:,:,g,n)=H_path(:,:,g,n)*(-diag((1:Nt))*1j*pi);
    end

    %derivative with respect to AOA
    for g=1:G
        D(:,:,G+g,n)=(-diag((1:Nr))*1j*pi)*H_path(:,:,g,n);
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
        for k=1:N
            J(i,j)=J(i,j)+real(trace(X*conj(D(:,:,i,n)).'*D(:,:,j,n)));
        end
    end
end

J_=T*J*T';
J_(4*G+2,4*G+2)=J_(4*G+2,4*G+2)+(1/power_clk);
end