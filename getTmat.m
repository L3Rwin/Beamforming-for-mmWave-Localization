function T = getTmat(G,p,q,AOA,AOD,SP,h,c)
T=zeros(4*G+2,5*G);

% solve the partial derivative of AOD with respect to p
T(1:2,1)=[-sin(AOD(1)); cos(AOD(1))]/norm(p-q);
% g>0, the derivative is 0.

% solve the partial derivative of AOA with respect to p
T(1:2,G+1)=[-sin(AOD(1)); cos(AOD(1))]/norm(p-q);
for j=1:G-1
    T(1:2,G+1+j)=[sin(AOA(1+j)); -cos(AOA(1+j))]/norm(SP(j,:)'-p);
end

% solve the partial derivative of tao with respect to p
T(1:2,4*G+1)=[cos(AOD(1)); sin(AOD(1))]/c;
for j=1:G-1
    T(1:2,4*G+1+j)=-[cos(AOA(1+j)); sin(AOA(1+j))]/c;
end

% solve the partial derivative with respect to alpha
T(3,G+1:2*G)=-ones(1,G);
% other terms are zeros

% solve the partial derivative of AOD with respect to r
for j=1:G-1
    T(2+2*j:3+2*j,1+j)=[-sin(AOD(1+j)); cos(AOD(1+j))]/norm(SP(j,:)'-q);
end

% solve the partial derivative of AOA with respect to r
for j=1:G-1
    T(2+2*j:3+2*j,G+1+j)=[-sin(AOA(1+j)); cos(AOA(1+j))]/norm(SP(j,:)'-p);
end

% solve the partial derivative of tao with respect to r
for j=1:G-1
    T(2+2*j:3+2*j,4*G+1+j)=[cos(AOD(1+j))+cos(AOA(1+j));sin(AOD(1+j))+sin(AOA(1+j))]/c;
end

% solve the partial derivative of h with respect to h
for j=1:G
    T(2*G+1+j,2*G+j)=1;
    T(3*G+1+j,3*G+j)=1;
end

% solve the partial derivative of tao with respect to delta t
T(4*G+2,4*G+1:5*G)=ones(1,G);
