function [x]=unif(lambda)
    rand("uniform");
    nb = lambda*rand()
    x = nb;
endfunction;

function [X_t]=queue(lambda, mu, tpsmax)
    X_t = zeros(2,1)
    tps = 0;
    i = 2;
    while tps < tpsmax
        if (X_t(2,i-1) == 0.) then
            temp = unif(2/lambda);
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            X_t(2,i) = X_t(2,i-1)+1;
        else
            temp = unif(2/(lambda +mu));
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            rand("uniform");
            r = (lambda + mu)*rand();
            if r < lambda then
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

function [esp] = espPratique(X)
  esp = integChemin(X, Id)
endfunction


function [var] = varPratique(X)
  esp = integChemin(X, Id)
  espXCarre = integChemin(X, carre)
  var = espXCarre - esp^2
endfunction
             
                           
function [proba] = probaPratique(X)
  tailleMax = max(X(2, :))
  proba = zeros(1, tailleMax + 1)
  for i=1:(tailleMax + 1)
    proba(i) = integChemin(X, indicatrice(i-1))
  end
endfunction

function [integ] = integChemin(X, f)
  integ = 0
  [nbCols, nbChangements] = size(X)
  tpsmax = X(1, nbChangements)
  for i= 2:nbChangements
    dt = X(1, i) - X(1, i - 1)
    integ = integ + (f(X(2, i-1)) * dt)
  end
  integ = integ / (tpsmax)
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
mu = 2

tpsmax = 200
X_t = queue(lambda, mu, tpsmax);
disp("espérance pratique")
disp(espPratique(X_t))
disp("variance pratique")
disp(varPratique(X_t))
disp("Lois")
absc = 0:max(X_t(2, :))
//plot2d2(X_t(1,:),X_t(2,:))
plot(absc, probaPratique(X_t), '--mo')

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
disp("Distribution de probabilité de la taille de la fileà l arrivée du client numéro 50")
disp(P_n(1,:)')

