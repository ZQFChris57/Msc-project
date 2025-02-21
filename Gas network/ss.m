clear all;
close all;
clc


%% computational Parameters definition
lambda1 = 0.5;
lambda2 = 1-lambda1;
alfa_p = 0.7;
alfa_g = 0.3;
toll_g = 0.01; %[kg/s]
toll_p = 50; %[Pa]
%% max number of loops that are acceptable
max1 = 1e4;
max2 = 1e4;

%% physical Parameters definition
T = 273.15+5; %K
UGC = 8.3144598*1e3; %J/kmol/k
MM = 16; %molar mass
RR = UGC./MM;

%% Upload the topology and features in the network
nodei = xlsread('EN.xlsx', 'BRANCHES', 'B2:B32');
nodeo = xlsread('EN.xlsx', 'BRANCHES', 'C2:C32');
Diametri = xlsread('EN.xlsx', 'BRANCHES', 'D2:D32').*1e-3; %m
L = xlsread('EN.xlsx', 'BRANCHES', 'E2:E32'); %m
epsi = xlsread('EN.xlsx', 'BRANCHES', 'F2:F32').*1e-3; %m
%% Calculatie the cross section f each branch
S = pi*(Diametri.^2)/4;

%% the boundary condtition
Gexe = xlsread('EN.xlsx', 'NODES', 'B2:B32'); %kg/s
Gguess = xlsread('EN.xlsx', 'BRANCHES', 'G2:G32'); %kg/s
Pguess = xlsread('EN.xlsx', 'NODES', 'C2:C32').*1e5; %Pa

%% Find number of nodes considering one of the vectors associated to the nodes
nNodes = size(Gexe,1); % give the number of rows for the incidence matrix
nBranches = size(nodei,1); %give the number of branches for the incidence matrix
%%
GR = digraph(nodei,nodeo);
A = -incidence(GR);
figure;
plot(GR)

%%
%incidence matrix initialization
A = zeros(nNodes,nBranches);

%fill the matrix A moving by colomns, each colomun is a branch and inlet node is 1, while outlet is -1.
for i = 1:nBranches
	in = nodei(i,1); %colume is the one indicated by number of Branch
	out = nodeo(i,1); %outlet ndoe is in nodeo, nodeo and nodei are vectors 
	A(in,i) = 1;
	A(out,i) = -1;
end
Aplus = (A>0).*1;
Aminus = (A<0).*1;

%Considering the inital guesses to enter the algorithm (START), these numbers are alerady from t excel
P2guess = Pguess.^2;% square of the pressure
G0 = Gguess;
G_1 = Gguess;

AT=A.'; %transpose of the incidence matrix转置

%Now entering in SIMPLE algorithm
abserrorp = 100;
i = 0;
%START of the first loop, big cycle
while abserrorp > toll_p && i < max1 % SIMPLE algorithm for cycle
	%% START of the second loop, fixed point
	abserror = 1;
	k = 0;
	while abserror > toll_g && k < max2
		Gr = lambda1*G0+lambda2*G_1; %calculate the related point of mass flow rate first
		%Now, we need momentum equation to find new values of mass flow rates
		%Thus to calcluare Fluid-dynamic resistances vector (Rf)
		f = Frictionfactoraverage(epsi,Gr,Diametri);
		pm = 2/3*(Aplus'*Pguess.^2+Aminus'*Pguess.^2+(Aplus'*Pguess).*(Aminus'*Pguess))./((Aplus'*Pguess)+(Aminus'*Pguess));
		ZZ = Papay(pm,T);
		c2 = ZZ.*RR.*T;
		Rf = (f.*L./Diametri.*c2)./(S.^2).*abs(Gr); %阻力向量
		Y = diag(Rf.^(-1));
		%write momentum equation, second step of fixed point algorithm
		G = Y*AT*P2guess; %second step of FP, new guess value
		%Now, have the mass flow rates, and to calculate mass flow rate from previous iternation
		abserror = max(abs(G-G0));
        G_1 = G0;
		G0 = G;
		k = k+1;
	end %end of FPA
		if (k == max2 && abserror > toll_g) %converage check
			disp('Maximum number of iteration on FP reached');
		end
		if abserror < toll_g
			disp('abs less than tollerance in fpa');
		end
	%%
	%Now to proceed with calculations of P' and G', at first to calculate P'
	H = A*Y*AT;
	b = -A*G0-Gexe; %before I had put G now I try to put G0
	%Boundary conditions setting, we set pressure, so we need to modify H and b
	P2guess(1) = (6*10^5).^2; % impose pressure at node 1 in thermal plant as minimun pressure, return pipe
	H(1,:) = (zeros);
	H(1,1) = 1; %(diagonal term is setted to 1)对角线项设置为 1
	b(1) = 0;
	%P1 = P' pressure corrections
	P1 = H\b; % results for the corrections	
	%Then to calculate new corrections on mass flow rates G'
	G1 = Y*AT*P1;%G'=Y*AT*P'
	%Now, can obtain new values for mass flow rates and pressures
	Pnew = P2guess+P1*alfa_p; %alfa are underrelaxtion factors for simple cycle
	Gnew = G+G1*alfa_g;

	abserrorp = max(abs(Pnew-P2guess));
	
	P2guess = Pnew;
	Pguess = sqrt(P2guess);
	G0 = Gnew;
	G_1 = Gnew;
	i = i+1;
	disp("Hello");
	if (i == max1)
		disp('Reached max value of cycles for SIMPLE cycle part');
	end
end %end of SIMPLE
