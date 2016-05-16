function [x]=expo(lambda)
    rand("uniform");
    nb = rand()
    x=-(1/lambda)*log(nb);
endfunction;

function [A_ij]=geninf(lambda,mu,K,i,j)
    if (j == i+1) then
        A_ij = lambda;
    elseif (i == j+1) then
        A_ij = min(i,K);
    elseif (i == j) then
        A_ij = -lambda -min(i,K);
    else
        A_ij = 0;
    end
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

//disp(size(X_t(1,:),"c"))
//X_t = queue(2,1,100);
//k=[1:size(X_t(1,:),"c")]
//for i = 1:size(X_t(1,:),"c")
//    distr(1,i) = distribution(size(X_t(1,:),"c"))
//end
//disp(distribution(0))
//plot2d2(k,distr)

function [esp] = espTheorique(lambda, mu)
  rho = lambda / mu
  esp = rho / (1 - rho)
endfunction

function [esp] = espPratique(X)
  esp = integChemin(X, Id)
endfunction

function [var] = varTheorique(lambda, mu)
  rho = lambda / mu
  var = rho / (1 - (rho^2))
endfunction

function [var] = varPratique(X)
  esp = integChemin(X, Id)
  espXCarre = integChemin(X, carre)
  var = espXCarre - esp^2
endfunction

function [proba] = probaTheorique(lambda, mu, X)
  tailleMax = max(X(2, :))
  proba = zeros(1, tailleMax + 1)
  rho = lambda / mu
  for i=1:(tailleMax + 1)
    proba(i) = (rho ** i) * (1 - rho)
  end
endfunction

function [proba] = probaPratique(X)
  tailleMax = max(X(2, :))
  proba = zeros(1, tailleMax + 1)
  for i=1:(tailleMax + 1)
    proba(i) = integChemin(X, indicatrice(i))
  end
endfunction

function [integ] = integChemin(X, f)
  integ = 0
  [nbCols, nbChangements] = size(X)
  tpsmax = X(1, nbChangements)
  for i=2:nbChangements
    dt = X(1, i) - X(1, i - 1)
    integ = integ + (f(X(2, i)) * dt)
  end
  integ = integ / tpsmax
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
tpsmax = 1000
X_t = queue(lambda, mu, tpsmax);
disp("espérance théorique")
disp(espTheorique(lambda, mu))
disp("espérance pratique")
disp(espPratique(X_t))
disp("variance théorique")
disp(varTheorique(lambda, mu))
disp("variance pratique")
disp(varPratique(X_t))
disp("Lois")
absc = 0:max(X_t(2, :))
plot(absc, probaTheorique(lambda, mu, X_t), '--r+')
plot(absc, probaPratique(X_t), '--mo')
