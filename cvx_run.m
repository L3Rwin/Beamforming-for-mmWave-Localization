function [PEB, status] = cvx_run(D,Jprior,Tmats,idMat,power_noise, N, G, Nt, Nsc, Ptot)
%calculate via cvx
    cvx_begin sdp
    variable t
    variable u(2,N)
    variable X(Nt,Nt)
    minimize t
    expression J(5*G,5*G,N)
    for nN=1:N
        for i=1:5*G
            for j=1:5*G
                sum=0;
                for n=1:Nsc
                    sum=sum+real(trace(X*conj(D(:,:,i,n,nN)).'*D(:,:,j,n,nN)));
                end
                J(i,j,nN)=sum*2/power_noise;
            end
        end
        expression Jj(4*G+2,4*G+2,nN)
        Jj(:,:,nN)=Tmats(:,:,nN)*J*Tmats(:,:,nN).'+Jprior;
    end
    subject to
        for nN=1:N
            [Jj(:,:,nN),idMat(:,1); idMat(1,:),u(1,nN)] == semidefinite(4*G+3);
            [Jj(:,:,nN),idMat(:,2); idMat(2,:),u(2,nN)] == semidefinite(4*G+3);
            u(1,nN)+u(2,nN)<=t
        end
        trace(X) == Ptot/Nsc;
        X == semidefinite(Nt);       
cvx_end 

PEB = getPEB(Jj,'1')
status = cvx_status;


%saveData2MatFile('clkBiasResults.mat',[sigma_in_meter;PEB],3)
end