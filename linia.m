clear all
clc
%PARAMETRY MODELU
n = 400;
Rc = 1 ; %rezystancja ca³ej linii
Gc = 10e3;
Lc = 0.04e-5; %indukcyjnoœæ ca³kowita
Cc = 5e-9;
R1 = 20; %rezystancja na pocz¹tku linii
R2 = 50; %rezystancja na koñcu linii

%parrametry jednego elementu
R = Rc/n;
G = n*Gc; %mówi jaka jest rezystancja miêdzy "+", a "-", to siê wymna¿a
C = Cc/n;
L = Lc/n;

%PARAMETRY MODELU
Tk = 0.25e-6;
dT = 0.5e-11;
T = 0:dT:Tk;
U = zeros(size(T));
Tim = 1e-8;
U(1:round(Tim/dT)) = 1; %wpisanie 1(impulsu) do wektora pobudzeñ, 

%generowanie macierzy A
N = 2*n;
A = sparse([],[],[],N,N,3*N);
A(1,1) = -1/(G*C); A(1,2) = 1/C; A(1,4) = -1/C;
A(2,1) = -1/L; A(2,2) = -(R+R1)/L;
K = 3:2:N-3;
%wprowadzenie zmiennej pomocniczej k, bo elementy zmieniaj¹ siê co 2
for i=K,
    A(i,i) = -1/(G*C); A(i,i+1) = 1/C;
    A(i,i+3) = -1/C; %bo trzeba przeskoczyæ napiêcie
    A(i+1,i-2) = 1/L; A(i+1,i) = -1/L; %i-2, bo U_n, trzeba przeskoczyæ napiêcie
    A(i+1,i+1) = -R/L;
end;
%koniec linii
A(N-1,N-1) = -1/(G*C)-1/(R2*C); %napiêcie na koñcu linii
A(N-1,N) = 1/C;
A(N,N-3) = 1/L; A(N,N-1) = -1/L;
A(N,N-1) = -1/L;
A(N,N) = -R/L;

%macierz B
B = sparse([],[],[],N,1,1);
B(2,1) = 1/L;

C = sparse([],[],[],1,N,1);
C(1,1) = 1;

D = 0;

%wektor pocz¹tkowy
X0 = zeros(N,1);

%SYMULACJA
%zwraca wartoœæ odpowiedzi modelu stanowego na pobudzenie z wektora
%pobudzeñ
Nt = size(T,2);
X = zeros(N,Nt);
Y = zeros(size(T));
As = dT*A;
Bs = dT*B;
Cs = C;
for i=1:Nt-1;
    %metoda ca³kowania trapezów
    Xp = X(:,i) + As*X(:,i) + Bs*U(i);
    X(:,i+1) = X(:,i) + 0.5*(As*X(:,i) + Bs*U(i) + As*Xp + Bs*U(i));
    Y(i+1) = Cs*X(:,i+1);
end

plot(T,Y); grid on; pause;

%ANIMACJA
L = (1/n)*[1:n];
P = mod(1:N,2)==1;

Nt = size(T,2);
Tindx = 1:Nt;
Tindx = Tindx(mod(1:Nt,10)==0);
for i = Tindx
    plot(L,X(P,i));
    axis([0 1 -0.1 0.6]);
    grid on;
    pause(0.005);
end;