close all
clear all

%% System parameters
G=2;                % number of paths (including LOS)
N=10;             % number of subcarriers 
Nt=32;              % number of TX antennas
Nr=16;              % number of RX antennas
c=300;              % speed of light meter / us
Rs=N*120;           % total BW in MHz
Ts=1/Rs;            % sampling period in usL
fc=28e3;            % frequence of carrier in MHz
posRx=[25 10]';     % RX (user) position, TX is assumed to be in [0, 0]
alpha=0.2;          % user orientation
gamma=0.1;          % NLOS path reflection coefficient
power_clk=0.01;     % Clock bias uncertainty
Ptot=0.1;

%% generate scatter points
SP=rand(G-1,2)*20-10;      % random points uniformly placed in a 20 m x 20 m area 


%% Compute Channel Parameters for G paths
TOA(1)=norm(posRx)/c;                                   % LOS TOA
AOD(1)=atan2(posRx(2),posRx(1));                        % LOS AOD
AOA(1)=atan2(posRx(2),posRx(1))-pi-alpha;               % LOS AOA
for g=1:G-1
    AOD(g+1)=atan2(SP(g,2),SP(g,1));
    AOA(g+1)=atan2(SP(g,2)-posRx(2),SP(g,1)-posRx(1))-alpha;
    TOA(g+1)=(randn(1)/sqrt(power_clk))+(norm(SP(g,:))+norm(posRx-SP(g,:)))/c;     % note: max distance should be below (N*Ts*c)
end

h=zeros(1,G);% some high channel gains
h(1)=exp(1j*2*pi*rand(1,1))*c/(4*pi*fc*norm(posRx));
for g=1:G-1
    h(g+1)=exp(1j*2*pi*rand(1,1))*gamma*c/(4*pi*fc*(norm(posRx-SP(g,:))+norm(SP(g,:))));
end

%% Generate channel
H=zeros(Nr,Nt,N);
% Contributions of each path to H
H_path=zeros(Nr,Nt,G,N);

for n=1:N
    for g=1:G 
        H_path(:,:,g,n)=h(g)*exp(-1j*2*pi*TOA(g)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA(g)))*sqrt(Nt)*getResponse(Nt,sin(AOD(g)))';
        H(:,:,n)=H(:,:,n)+H_path(:,:,g,n);
    end   
end

T = getTmat(G,posRx,[0,0]',AOA,AOD,SP,h,c);
Ut=[getResponse(Nt,AOD),(-diag((1:Nt))*1j*pi)*getResponse(Nt,AOD)];
idMat=eye(4*G+2,4*G+2);
Jprior(4*G+2,4*G+2)=(1/power_clk);

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

cvx_begin sdp
    variable lambda(2*G,2*G)
    variable u(1,2)
    minimize(u*ones(2,1))
    expression Xopt(Nt,Nt)
    Xopt=Ut*lambda*conj(Ut).';
    expression J(5*G,5*G)
    for i=1:5*G
        for j=1:5*G
            sum=0;
            for k=1:N
                sum=sum+real(trace(Xopt*conj(D(:,:,i,n)).'*D(:,:,j,n)));
            end
            J(i,j)=sum;
        end
    end
    expression Jj(4*G+2,4*G+2)
    Jj=T*J*T'+Jprior;
    subject to
        [Jj,idMat(:,1); idMat(1,:),u(1,1)] == hermitian_semidefinite(4*G+3);
        [Jj,idMat(:,2); idMat(2,:),u(1,2)] == hermitian_semidefinite(4*G+3);
        trace(Xopt)==Ptot/N;
        lambda==hermitian_semidefinite(2*G);
cvx_end
