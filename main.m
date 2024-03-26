close all
clear all

%% System parameters
G=2;                % number of paths (including LOS)
Nsc=100;             % number of subcarriers 
Nt=32;              % number of TX antennas
Nr=16;              % number of RX antennas
c=3e8;              % speed of light meter / s
Rs=Nsc*12e7;        % total BW in Hz
Ts=1/Rs;            % sampling period in usL
fc=28e9;            % frequence of carrier in Hz
posRx=[25 10]';     % RX (user) position, TX is assumed to be in [0, 0]
gamma=1e-2;        % NLOS path reflection coefficient
power_clk=1e-10;     % Clock bias uncertainty
M=16;               % number of pilot frames 
Ptot=0.1*M;         % power of beamforming vector
power_noise=power(10,0.1*(8-30-174))*Nsc*12e7;

%% generate scatter points
SP=[15,25];      


%% Compute Channel Parameters for G paths
alpha=rand(1,1)*2*pi-pi;                                % User Orientation
TOA(1)=norm(posRx)/c+(randn(1)*sqrt(power_clk));        % LOS TOA
AOD(1)=atan2(posRx(2),posRx(1));                        % LOS AOD
AOA(1)=atan2(posRx(2),posRx(1))-pi-alpha;               % LOS AOA
for g=1:G-1
    AOD(g+1)=atan2(SP(g,2),SP(g,1));
    AOA(g+1)=atan2(SP(g,2)-posRx(2),SP(g,1)-posRx(1))-alpha;
    TOA(g+1)=(norm(SP(g,:))+norm(posRx-SP(g,:)))/c+(randn(1)*sqrt(power_clk));     % note: max distance should be below (N*Ts*c)
end

h=zeros(1,G);% some high channel gains
h(1)=exp(1j*2*pi*rand(1,1))*c/(4*pi*fc*norm(posRx));
for g=1:G-1
    h(g+1)=exp(1j*2*pi*rand(1,1))*gamma*c/(4*pi*fc*(norm(posRx-SP(g,:))+norm(SP(g,:))));
end

%h=10*ones(1,Nsc);

%% Generate channel
H=zeros(Nr,Nt,Nsc);

% Contributions of each path to H
H_path=zeros(Nr,Nt,G,Nsc);

for n=1:Nsc
    for g=1:G 
        H_path(:,:,g,n)=h(g)*exp(-1j*2*pi*TOA(g)*(n-1)/(Nsc*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA(g)))*sqrt(Nt)*getResponse(Nt,sin(AOD(g)))';
        H(:,:,n)=H(:,:,n)+H_path(:,:,g,n);
    end   
end

%% Preparation before CVX

T = getTmat(G,posRx,[0,0]',AOA,AOD,SP,h,c);
Ut=[getResponse(Nt,sin(AOD)),(-diag((0:Nt-1))*1j*pi)*getResponse(Nt,sin(AOD))*diag(cos(AOD))]*sqrt(Nt);
idMat=eye(4*G+2,4*G+2);
Jprior=zeros(4*G+2,4*G+2);
Jprior(4*G+2,4*G+2)=(1/power_clk);

%% Calculate derivative of H with respect to eta

D=zeros(Nr,Nt,5*G,Nsc);
for n=1:Nsc
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
        D(:,:,4*G+g,n)=H_path(:,:,g,n)*(-1j*2*pi*(n-1)/(Nsc*Ts));
    end
end
 
%% CVX solve
%{
cvx_begin sdp
    %cvx_solver mosek
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

%}
cvx_begin sdp
    %cvx_solver mosek
    %cvx_save_prefs
    variable lambda(2*G,2*G)
    variable u(2,1) 
    minimize(ones(1,2)*u)
    expression Xopt(Nt,Nt)
    Xopt=Ut*lambda*conj(Ut).';
    expression J(5*G,5*G)
    for i=1:5*G
        for j=1:5*G
            sum=0;
            for n=1:Nsc
                sum=sum+real(trace(Xopt*conj(D(:,:,i,n)).'*D(:,:,j,n)));
            end
            J(i,j)=sum*2/power_noise;
        end
    end
    expression Jj(4*G+2,4*G+2)
    Jj=T*J*T'+Jprior;
    subject to
        [Jj,idMat(:,1); idMat(1,:),u(1,1)] == semidefinite(4*G+3);
        [Jj,idMat(:,2); idMat(2,:),u(2,1)] == semidefinite(4*G+3);
        trace(Xopt) == Ptot/Nsc;
        lambda == semidefinite(2*G);
cvx_end
