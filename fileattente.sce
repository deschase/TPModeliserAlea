function [x]=expo(lambda)
    rand("uniform");
    nb = rand()
    x = -(1/lambda)*log(nb);
endfunction;

function [X_t]=queue(lambda, mu, tpsmax)
    X_t = zeros(2,1)
    tps = 0;
    i = 2;
    while tps < tpsmax
        if (X_t(2,i-1) == 0.) then
            temp = expo(lambda);
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            X_t(2,i) = X_t(2,i-1)+1;
        else
            temp = expo(lambda + mu);
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            rand("uniform");
            r = rand();
            if r < lambda/(lambda + mu) then
                X_t(2,i) = X_t(2,i-1)+1;
            else
                X_t(2,i) = X_t(2,i-1)-1;
            end
        end
        i = i+1;
    end
endfunction

function [dist] = distribution(n)
    dist = 0;
    if (X_t(2,i) == n) then
        dist = dist + X_t(1,1)
    end
    for i =2:size(X_t(1,:),"c")
        if (X_t(2,i) == n) then
            dist = dist + X_t(1,i) - X_t(1,i-1);
        end
    dist = 1/100 * dist
    end
endfunction

function [esp] = espTheorique(lambda, mu)
  rho = lambda / mu
  esp = rho / (1 - rho)
endfunction

function [esp] = espPratique(X, borneinf)
  esp = integChemin(X, Id, borneinf)
endfunction

function [var] = varTheorique(lambda, mu)
  rho = lambda / mu
  var = rho / (1 - (rho^2))
endfunction

function [var] = varPratique(X,borneinf)
  esp = integChemin(X, Id, borneinf)
  espXCarre = integChemin(X, carre, borneinf)
  var = espXCarre - esp^2
endfunction

function [proba] = probaTheorique(lambda, mu, X)
  tailleMax = max(X(2, :))
  proba = zeros(1, tailleMax + 1)
  rho = lambda / mu
  for i=1:(tailleMax + 1)
    proba(i) = (rho ** (i-1)) * (1 - rho)
  end                      
endfunction                
                           
function [proba] = probaPratique(X,borneinf)
  tailleMax = max(X(2, :))
  proba = zeros(1, tailleMax + 1)
  for i=1:(tailleMax + 1)
    proba(i) = integChemin(X, indicatrice(i),borneinf)
  end
endfunction

function [integ] = integChemin(X, f, borneinf)
  integ = 0
  [nbCols, nbChangements] = size(X)
  tpsmax = X(1, nbChangements)
  for i= borneinf:nbChangements
    dt = X(1, i) - X(1, i - 1)
    integ = integ + (f(X(2, i)) * dt)
  end
  integ = integ / (tpsmax - borneinf)
endfunction

function [res] = Id(x)
  res = x
endfunction

function [res] = carre(x)
  res = x^2
endfunction

function [ind] = indicatrice(x)
  execstr("function [res] = ind(y); if " + string(x) + " == y then; res = 1; else; res = 0; end; endfunction")
endfunction



lambda = 1
mu = 3

tpsmin = 1000
tpsmax = 2000
X_t = queue(lambda, mu, tpsmax);
disp("espérance théorique")
disp(espTheorique(lambda, mu))
disp("espérance pratique")
disp(espPratique(X_t, tpsmin))
disp("variance théorique")
disp(varTheorique(lambda, mu))
disp("variance pratique")
disp(varPratique(X_t, tpsmin))
disp("Lois")
absc = 0:max(X_t(2, :))
plot(absc, probaTheorique(lambda, mu, X_t), '--r+')
plot(absc, probaPratique(X_t, tpsmin), '--mo')

// Calcul de la probabilité du nombre de personnes dans la file à l'arrivée du client numéro client.

numtronc=100
P = zeros(numtronc,numtronc);
for i = 1:(numtronc-1)
    for j = 1:(i+1)
        if j == 1 then
            P(i,j) = (mu/(lambda+mu))^(i-j+1)
        else
            P(i,j) = (lambda/(lambda+mu))*(mu/(lambda+mu))^(i-j+1)
        end
    end
end

client = 50;
P_n = P^(client);
disp("Distribution de probabilité de la taille de la fileà l'arrrivée du client nu")
disp(P_n(1,:)')

