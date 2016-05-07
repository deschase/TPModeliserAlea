clear

// Loi normale
function [x]=normale(y,m,s2)
  x=%e^(-(y-m).^2/2/s2)/sqrt(2*%pi*s2)
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

nbCrabes = 1000;
esp = y * x / nbCrabes;
var = ((y - esp).^2) * x / nbCrabes;
//plot(normale(y ,esp , var))

// Tracer l'histogramme
//plot(x)

// Données
pi0=[1; 3]/2/2;
pi=pi0;
mu=[.57; .67];
s2=[1 ;1]/10000;

rho=ones(2,1000);

// Algorithme EM pour les crabes
//------------------------------

N=1000;
R=zeros(5,N+1);
R(:,1)=[mu(1);mu(2);pi(1);s2(1);s2(2)];

nbPopulation = 2
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
    disp("s2 :")
    disp(s2)
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
// A FAIRE
