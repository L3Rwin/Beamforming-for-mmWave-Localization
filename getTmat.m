function T = getTmat(G,p,q,AOA,AOD,SP,h,c)
%Transformation from J and J~
T=zeros(4*G+2,5*G);

% solve the partial derivative of AOD with respect to p
T(1:2,1)=[-sin(AOD(1)); cos(AOD(1))]/norm(p-q);
% g>0, the derivative is 0.

% solve the partial derivative of AOA with respect to p
T(1:2,G+1)=[-sin(AOD(1)); cos(AOD(1))]/norm(p-q);
for j=1:G-1
    T(1:2,G+1+j)=[sin(AOA(1+j)); -cos(AOA(1+j))]/norm(SP(j,:)-q);
end

% solve the partial derivative of h with respect to p
t= -h(1)*(p-q)/power(norm(p-q),2);
T(1:2,2*G+1)=real(t);
T(1:2,3*G+1)=imag(t);
for j=1:G-1
    r=SP(j,:);
    t=-h(1+j)*(p-r)/((norm(p-r)+norm(q-r))*norm(p-r));
    T(1:2,2*G+1+j)=real(t);
    T(1:2,3*G+1+j)=imag(t);
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
    T(2+2*j:3+2*j,1+j)=[-sin(AOD(1+j)); cos(AOD(1+j))]/norm(SP(j,:)-q);
end

% solve the partial derivative of AOA with respect to r
for j=1:G-1
    T(2+2*j:3+2*j,G+1+j)=[-sin(AOA(1+j)); cos(AOA(1+j))]/norm(SP(j,:)-p);
end

% solve the partial derivative of h with respect to r
for j=1:G-1
    r=SP(j,:);
    t=-h(1+j)*(((r-q)/norm(q-r))+((r-p)/norm(p-r)))/(norm(p-r)+norm(q-r));
    T(2+2*j:3+2*j,2*G+1+j)=real(t);
    T(2+2*j:3+2*j,3*G+1+j)=imag(t);
end

% solve the partial derivative of tao with respect to r
for j=1:G-1
    r=SP(j,:);
    T(2+2*j:3+2*j,4*G+1+j)=(((r-q)/norm(q-r))+((r-p)/norm(p-r)))/c;
end

% solve the partial derivative of h with respect to h
for j=1:G
    T(2*G+1+j,2*G+j)=1;
    T(3*G+1+j,3*G+j)=1;
    T(2*G+1+j,3*G+j)=real(h(j))/imag(h(j));
    T(2*G+1+j,3*G+j)=imag(h(j))/real(h(j));
end

% solve the partial derivative of tao with respect to delta t
T(4*G+2,5*G)=1;
