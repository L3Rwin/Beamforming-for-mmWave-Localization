function PEB = getPEB(Jj, status)
% get positiong error bound (PEB)
    PEB=-1;
    if length(status) ~= length('Failed')
        invMat=inv(Jj);
        PEB=sqrt(trace(invMat(1:2,1:2)));
    end
end