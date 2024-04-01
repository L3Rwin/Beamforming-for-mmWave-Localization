close all
clear all

%% System parameters
G=2;                % number of paths (including LOS)
Nsc=10;            % number of subcarriers 
Nt=32;              % number of TX antennas
Nr=16;              % number of RX antennas
c=3e8;              % speed of light meter / s
Rs=Nsc*12e7;        % total BW in Hz
Ts=1/Rs;            % sampling period in usL
fc=28e9;            % frequence of carrier in Hz
cenRx=[25 10]';     % RX (user) position, TX is assumed to be in [0, 0]
gamma=0.2;         % NLOS path reflection coefficientzzr
sigma_in_meter=10   ;
power_clk=power(sigma_in_meter/c,2);% Clock bias uncertainty
M=16;               % number of pilot frames 
Ptot=0.1*M;         % power of beamforming vector
power_noise=power(10,0.1*(8-30-174))*Nsc*12e7;
shrinkFactor=1e-1;   % shrink the optimization to avoid large absolute size
N=16;               % number of grid point 
spUncertain=5;       % Uncertain of SP
ueUncertain=0.3;       % Uncertain of UE

%% generate scatter points

SPc=[15,25];  
idMat=eye(4*G+2,4*G+2);
Jprior=zeros(4*G+2,4*G+2);
Jprior(4*G+2,4*G+2)=1/power_clk;

%% space declaration
TOA=zeros(N,G);
AOD=zeros(N,G);
AOA=zeros(N,G);
h=zeros(N,G);% some high channel gains
H=zeros(Nr,Nt,Nsc,N);
H_path=zeros(Nr,Nt,G,Nsc,N); % Contributions of each path to H
Tmats=zeros(4*G+2,5*G,N);
for nN=1:N
    biasWord=dec2bin(nN-1,4);
    SP=SPc+[str2num(biasWord(1))-0.5,str2num(biasWord(2))-0.5]*spUncertain;
    posRx=cenRx+[str2num(biasWord(3))-0.5,str2num(biasWord(4))-0.5]*ueUncertain;

    %% Compute Channel Parameters for G paths
    alpha=0.2;                                % User Orientation
    TOA(nN,1)=norm(posRx)/c+(randn(1)*sqrt(power_clk));        % LOS TOA
    AOD(nN,1)=atan2(posRx(2),posRx(1));                        % LOS AOD
    AOA(nN,1)=atan2(posRx(2),posRx(1))-alpha;                  % LOS AOA
    for g=1:G-1
        AOD(nN,g+1)=atan2(SP(g,2),SP(g,1));
        AOA(nN,g+1)=atan2(SP(g,2)-posRx(2),SP(g,1)-posRx(1))-alpha;
        TOA(nN,g+1)=(norm(SP(g,:))+norm(posRx-SP(g,:)))/c+(randn(1)*sqrt(power_clk));     % note: max distance should be below (N*Ts*c)
    end

    h(nN,1)=exp(1j*2*pi*(rand(1,1)-0.5))*c/(4*pi*fc*norm(posRx));
    for g=1:G-1
        h(nN,g+1)=exp(1j*2*pi*(rand(1,1)-0.5))*gamma*c/(4*pi*fc*(norm(posRx-SP(g,:))+norm(SP(g,:))));
    end

    %% Generate channel
    for n=1:Nsc
        for g=1:G 
            H_path(:,:,g,n,nN)=h(nN,g)*exp(-1j*2*pi*TOA(nN,g)*(n-1)/(Nsc*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA(nN,g)))*sqrt(Nt)*getResponse(Nt,sin(AOD(nN,g)))';
            H(:,:,n,nN)=H(:,:,n,nN)+H_path(:,:,g,n,nN);
        end   
    end

    %% Preparation before CVX
    Tmats(:,:,nN) = getTmat(G,posRx,[0,0]',AOA,AOD,SP,h,c);

    %% Calculate derivative of H with respect to eta
    D=zeros(Nr,Nt,5*G,Nsc,N);
    for n=1:Nsc
        %derivative with respect to AOD 
        for g=1:G
            D(:,:,g,n,nN)=H_path(:,:,g,n,nN)*cos(AOD(nN,g))*(-diag((0:Nt-1))*1j*pi);
        end

        %derivative with respect to AOA
        for g=1:G
            D(:,:,G+g,n,nN)=(-diag((0:Nr-1))*1j*pi)*cos(AOA(nN,g))*H_path(:,:,g,n,nN);
        end

        %derivative with respect to h
        for g=1:G
            D(:,:,2*G+g,n,nN)=H_path(:,:,g,n,nN)/h(nN,g);
            D(:,:,3*G+g,n,nN)=1j*H_path(:,:,g,n,nN)/h(nN,g);
        end
    
        %derivative with respect to tao
        for g=1:G
            D(:,:,4*G+g,n,nN)=H_path(:,:,g,n,nN)*(-1j*2*pi*(n-1)/(Nsc*Ts));
        end
    end
end


 
%% CVX solve

[PEB, status]=cvx_run(D,Jprior,Tmats,idMat,power_noise, N, G, Nt, Nsc, Ptot);
%[PEB, status]=cvx_run(D,Jprior,T,idMat,power_noise, G, Nt, Nsc, Ptot);
