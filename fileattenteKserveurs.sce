
function [x]=expo(lambda)
    rand("uniform");
    nb = rand()
    x = -(1/lambda)*log(nb);
endfunction;

function [X_t]=queue(lambda, mu, K, tpsmax)
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
            temp = expo(lambda + min(K,X_t(2,i-1))*mu);
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            rand("uniform");
            r = rand();
            if r < lambda/(lambda + min(K, X_t(2,i-1))*mu) then
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

function [esp] = espTheorique(lambda, mu, K)
  rho = lambda / (K*mu)
  esp = K*rho + rho**(K+1)/((1 - rho)**2)*K**K/factorial(K)
endfunction

function [esp] = espPratique(X)
  esp = integChemin(X, Id)
endfunction

function [proba] = probaTheorique(lambda, mu, K, X)
  tailleMax = max(X(2, :))
  proba = zeros(1, tailleMax + 1)
  rho = lambda / (K*mu)
  somme = 0
  for i = 0:(K-1)
      somme = somme + (K*rho)**i/(factorial(i))
  end
  proba(1) = 1/(somme + 1/(1-rho)*(K*rho)**K/(factorial(K)))
  for i=1:min(K-1,tailleMax)
      proba(i+1) = proba(1)*((rho*K)**i)/(factorial(i))
  end  
  for j = K:(tailleMax)
      proba(j+1) = proba(1)*(rho**j)*(K**K)/(factorial(K))
  end
  disp(proba)
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



lambda = 20
mu = 1
K = 10

tpsmax = 100
X_t = queue(lambda, mu, K, tpsmax);
disp("espérance théorique")
disp(espTheorique(lambda, mu, K))
disp("espérance pratique")
disp(espPratique(X_t))
disp("Lois")
absc = 0:max(X_t(2, :))
plot2d2(X_t(1,:),X_t(2,:))
//plot(absc, probaTheorique(lambda, mu, K, X_t), '--r+')
//plot(absc, probaPratique(X_t), '--mo')

// Calcul de la probabilité du nombre de personnes dans la file à l'arrivée du client numéro client.

numtronc=100
P = zeros(numtronc,numtronc);
for i = 1:(numtronc-1)
    for j = 1:(i+1)
        if j == 1 then
            P(i,j) = (min(K, i)*mu/(lambda+min(K, i)*mu))^(i-j+1)
        else
            P(i,j) = (lambda/(lambda+min(K, i)*mu))*(min(K, i)*mu/(lambda+min(K, i)*mu))^(i-j+1)
        end
    end
end

client = 50;
P_n = P^(client);
disp("Distribution de probabilité de la taille de la fileà l arrivée du client numéro 50")
disp(P_n(1,:)')

