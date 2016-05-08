clear

// Loi normale
function [x]=normale(y,m,s2)
  x=%e^(-(y-m).^2/2/s2)/sqrt(2*%pi*s2)
endfunction;

// Test du chi 2

// Cette fonction retourne la p valeur P(chi2>zeta_n) 
// du test du chi2 d'adequation de loi

// N est un vecteur ligne des occurences observees
// p0 est un vecteur ligne correspondant a la loi sous H0

function[proba]=test_chi2(N,p0)
  n=sum(N);// taille de l'echantillon observe
  // calcul de zeta_ n
  zeta_n=n*sum(((N/n-p0).^2)./p0);
  // nombre de degres de liberte  (= nombre de classes dans N-1)
  d= length(N)-1;
  // on calcule la proba pour un chi 2 à d-1 degres d'etre superieur a zeta
  [p,q]=cdfchi("PQ",zeta_n,d);
  proba=q;
endfunction;

// Ouvrir le fichier de données (nombre de crabes par intervalle)
x=fscanfMat('crabe.txt');
x=x;

// intervalles
y=.580+.002+.004*[0:28];
yM=y+.002;
ym=y-.002;
Max=25;


// Dessiner la loi normale correspondante

nbCrabes = sum(x);
esp = (y * x) / nbCrabes;
var = ((y + 0.002 - esp).*(y + 0.002 - esp))*x / (nbCrabes);
//plot(y, normale(y, esp, var))

// Tracer l'histogramme
//bar(y,x/4)

// P-valeur du test du chi 2
disp("P-valeur du test du chi 2:")
disp(test_chi2(x'/4, normale(y ,esp , var)))

// Données
pi0=[1; 3; 4]/2/2/2;
pi=pi0;
mu=[.57; .67; .60];
s2=[1; 1; 1]/10000;

rho=ones(2,1000);

// Algorithme EM pour les crabes
//------------------------------

N=1000;
R=zeros(8,N+1);
R(:,1)=[mu(1);mu(2);mu(3);pi(1);pi(2);s2(1);s2(2);s2(3)];

nbPopulation = 3
Y = zeros(1, 1000)
crabId = 0
yLength = 29
firstSegmentCrabId = 0
for id=1:yLength
  yi = y(id)
  xi = x(id)
  for idInSegment=1:xi
    crabId = idInSegment + firstSegmentCrabId
    Y(crabId) = yi
  end
  firstSegmentCrabId = xi + firstSegmentCrabId
end

function [vecteurF]= f(Y, mu, s2)
  vecteurF = zeros(1000, nbPopulation)
  for crabId=1:1000
    for i=1:nbPopulation
      vecteurF(crabId, i) = normale(Y(crabId), mu(i), s2(i))
    end
  end
endfunction

function [vecteurFTheta]= fTheta(Y, mu, s2, pi)
  vecteurFTheta = f(Y, mu, s2) * pi
endfunction

for k=1:N;
  if modulo(k, 10) == 0 then
    disp("pi :")
    disp(pi)
    disp("mu :")
    disp(mu)
    disp("sigma :")
    disp(sqrt(s2))
  end
  rho = diag(1 ./ fTheta(Y)) * f(Y, mu, s2) * diag(pi)
  pi = (sum(rho, 1) / N)'
  mu = Y * rho ./ sum(rho, 1)
  for j=1:nbPopulation
    matS2 = ((Y - mu(j)).^2 * rho) ./ sum(rho, 1)
    s2(j) = matS2(j)
  end
  s2 = s2'
end;

// Affichages
bar(y,x/4)
plot(y, pi(1)*normale(y, mu(1), s2(1)),"red")
plot(y, pi(2)*normale(y, mu(2), s2(2)),"red")
plot(y, pi(3)*normale(y, mu(3), s2(3)),"red")
plot(y, pi(1)*normale(y, mu(1), s2(1)) + ...
        pi(2)*normale(y, mu(2), s2(2)) + ...
        pi(3)*normale(y, mu(3), s2(3)),"red")
