function[COEFF]=coeff_aero(k,alfa) 
%k=numero serie del profilo e deve esssere dato some stringa es:'0012'
%alpha deve essere dato in gradi
%restitutisce la matrice le cui due colonne sono cl e cd

%estrazione dati relativi al profilo
filename       = strcat('data/NACA',k,'.dat');
format long;
airfoil_struct = importdata(filename,' ',1);
coord          = airfoil_struct.('data');
x0             = coord(:,1);
y0             = coord(:,2);
clear filename;

%estrazione file sul cp 
filename       = strcat('data/naca',k,'_cp_a',num2str(alfa*10),'.txt');
analysis_struct = importdata(filename,' ',6); 
analysis        = analysis_struct.('data');
cp_alfa          = analysis(:,2);

%rotazione cordinate profilo
alfar=alfa*(pi/180);
x=(x0-0.25)*cos(-alfar)-y0*sin(-alfar)+0.25;
y=(x0-0.25)*sin(-alfar)+y0*cos(-alfar);

%quantità relative ai singoli pannelli
np=length(x0);
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
clP=cp_alfaP.*costheta.*deltal;
cl=sum(clP);
%calcolo cd
cdP=cp_alfaP.*sentheta.*deltal;
cd=sum(cdP);
COEFF=[cl,cd];
end