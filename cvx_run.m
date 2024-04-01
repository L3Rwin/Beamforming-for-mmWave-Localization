function [PEB, status] = cvx_run(D,Jprior,T,idMat,power_noise, G, Nt, Nsc, Ptot)
%calculate via cvx
    cvx_begin sdp
    %cvx_save_prefs
    variable u(2,1)
    variable X(Nt,Nt)
    minimize(ones(1,2)*u)
    expression J(5*G,5*G)
    for i=1:5*G
        for j=1:5*G
            sum=0;
            for n=1:Nsc
                sum=sum+real(trace(X*conj(D(:,:,i,n)).'*D(:,:,j,n)));
            end
            J(i,j)=sum*2/power_noise;
        end
    end
    expression Jj(4*G+2,4*G+2)
    Jj=T*J*T'+Jprior;
    subject to
        [Jj,idMat(:,1); idMat(1,:),u(1,1)] == semidefinite(4*G+3);
        [Jj,idMat(:,2); idMat(2,:),u(2,1)] == semidefinite(4*G+3);
        trace(X) == Ptot/Nsc;
        X == semidefinite(Nt);       
cvx_end 

PEB = getPEB(Jj,'1')
status = cvx_status;


%saveData2MatFile('clkBiasResults.mat',[sigma_in_meter;PEB],3)
end