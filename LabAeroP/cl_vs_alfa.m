%CL vs ALFA considerando un profilo naca aperto
clear all
clc

%dati iniziali:codice profilo ed estremi e passo dell'intervallo in cui
%varia l'incidenza
k='0012'; 
r=5; 
b=0.5; %i valori di r,b,l possono essere modificati per variare l'intervallo dell'angolo di incidenza
l=0; 
alfa=(l:b:r);
u=length(alfa);
alfar=alfa*(pi/180);%passaggio da gradi a radianti

%inizializzazione del cl e cd
cl=zeros(u,1);
cd=zeros(u,1);

%estrazione dati relativi al profilo
filename       = strcat('data/NACA',k,'.dat');
format long;
airfoil_struct = importdata(filename,' ',1);
coord          = airfoil_struct.('data');
x0             = coord(:,1);
y0             = coord(:,2);

for j=1:u
%estrazione file sul cp 
filename       = strcat('data/naca',k,'_cp_a',num2str(alfa(j)*10),'.txt');
analysis_struct = importdata(filename,' ',6); 
analysis        = analysis_struct.('data');
cp_alfa          = analysis(:,2);

%rotazione cordinate profilo
x=(x0-0.25)*cos(-alfar(j))-y0*sin(-alfar(j))+0.25;
y=(x0-0.25)*sin(-alfar(j))+y0*cos(-alfar(j));

%inizializzazione e calcolo delle quantità relative ai singoli pannelli
np=length(x);
nP=np-1;
deltax=zeros(nP,1);
deltay=zeros(nP,1);
cp_alfaP=zeros(nP,1);
for i=1:nP
deltax(i) = x(i+1)-x(i);
deltay(i)=y(i+1)-y(i);
cp_alfaP(i)=(cp_alfa(i)+cp_alfa(i+1))/2;
end
deltal=sqrt(deltax.^2+deltay.^2);

%prodotti scalari n*ey, n*ex
costheta=deltax./deltal; 
sentheta=deltay./deltal;

%calcolo cl
clp=cp_alfaP.*costheta.*deltal;
cl(j)=sum(clp);

%calcolo cd
cdp=cp_alfaP.*sentheta.*deltal;
cd(j)=sum(cdp);
end
disp('        alpha                cl                    cd')
G=[alfa', cl, cd]
figure(1);
plot(alfa,cl,'r.','MarkerSize',20);