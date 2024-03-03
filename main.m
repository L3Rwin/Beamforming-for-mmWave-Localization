close all
clear all

%% System parameters
G=4;                % number of paths (including LOS)
Rs=100;             % total BW in MHz
N=10;               % number of subcarriers 
Nt=32;              % number of TX antennas
Nr=Nt;              % number of RX antennas
Nb=Nt*2;            % number of beams in dictionary; 
Ns=20;              % number of beams sent
c=300;              % speed of light meter / us
Ts=1/Rs;            % sampling period in usL
fc=28e3;            % frequence of carrier in MHz
posRx=[5 1]';       % RX (user) position, TX is assumed to be in [0, 0]
alpha=0.2;          % user orientation
sigma=0.1;          % noise standard deviation
gamma=0.1;           % NLOS path reflection coefficient

%% generate scatter points
SP=rand(G-1,2)*20-10;      % random points uniformly placed in a 20 m x 20 m area 


%% Compute Channel Parameters for G paths
TOA(1)=norm(posRx)/c;                                   % LOS TOA
AOD(1)=atan2(posRx(2),posRx(1));                        % LOS AOD
AOA(1)=atan2(posRx(2),posRx(1))-pi-alpha;               % LOS AOA
for g=1:G-1
    AOD(g+1)=atan2(SP(g,2),SP(g,1));
    AOA(g+1)=atan2(SP(g,2)-posRx(2),SP(g,1)-posRx(1))-alpha;
    TOA(g+1)=(norm(SP(g,:))+norm(posRx-SP(g,:)))/c;     % note: max distance should be below (N*Ts*c)
end

h=zeros(1,G);% some high channel gains
h(1)=exp(1j*2*pi*rand(1,1))*c/(4*pi*fc*norm(posRx));
for g=1:G-1
    h(g+1)=exp(1j*2*pi*rand(1,1))*gamma*c/(4*pi*fc*(norm(posRx-SP(g,:))+norm(SP(g,:))));
end


%% Create dictionary 
Ut=zeros(Nt,Nb);
Ur=zeros(Nr,Nb);
aa=-Nb/2:Nb/2-1;
aa=2*aa/Nb;                                 % dictionary of spatial frequencies
for m=1:Nb
    Ut(:,m)=getResponse(Nt,aa(m))*sqrt(Nt);  
    Ur(:,m)=getResponse(Nr,aa(m))*sqrt(Nr);
end


%% Generate channel: eq. (1)-(5) from the paper
H=zeros(Nr,Nt,N);
for n=1:N
    for g=1:G       
        H(:,:,n)=H(:,:,n)+h(g)*exp(-1j*2*pi*TOA(g)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA(g)))*sqrt(Nt)*getResponse(Nt,sin(AOD(g)))';
    end   
end


%% Visualize the beamspace channel for 1 subcarrier in AOA/AOD space
Hb=Ur'*H(:,:,1)*Ut;
mesh(asin(aa),asin(aa),abs(Hb));
xlabel('AOD'); ylabel('AOA');



