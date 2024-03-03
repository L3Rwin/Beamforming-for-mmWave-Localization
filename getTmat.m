function T = getTmat(G,p,q,AOA,AOD,SP,h,c)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
T=zeros(5*G,4*G+2);

% solve the partial derivative of AOD with respect to p
T(1:2,0)=[-sin(AOD(1)); cos(AOD(1))]/norm(p-q);
% g>0, the derivative is 0.

% solve the partial derivative of AOA with respect to p
T(1:2,G+1)=[-sin(AOD(1)); cos(AOD(1))]/norm(p-q);
for j=1:G-1
    T(1:2,G+1+j)=[sin(AOA(1+j)); -cos(AOA(1+j))]/norm(SP(j,:)-q);
end

% solve the partial derivative of h with respect to pJ
J= h(1)*(p-q)/power(norm(p-q),2);
T(1:2,2*G+1)=-real(J);
T(1:2,3*G+1)=-imag(J);
for j=1:G-1
    r=SP(j,:);
    J=-h(1+j)*(p-r)/((norm(p-r)+norm(q-r))*norm(p-r));
    T(1:2,2*G+1+j)=real(J);
    T(1:2,3*G+1+j)=imag(J);
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
    T(4+j:5+j,1+j)=[-sin(AOD(1+j)); cos(AOD(1+j))]/norm(SP(j,:)-q);
end

% solve the partial derivative of AOA with respect to r
for j=1:G-1
    T(4+j:5+j,G+1+j)=[-sin(AOA(1+j)); cos(AOA(1+j))]/norm(SP(j,:)-p);
end