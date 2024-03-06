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
gamma=0.1;          % NLOS path reflection coefficient
power_clk=0.1;      % Clock bias uncertainty

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

%% Generate channel: eq. (1)-(5) from the paper
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

