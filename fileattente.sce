function [x]=expo(lamb)
    rand("uniform");
    nb = rand()
    x=-(1/lamb)*log(nb);
endfunction;

function [A_ij]=geninf(lamb,mu,K,i,j)
    if (j == i+1) then
        A_ij = lamb;
    elseif (i == j+1) then
        A_ij = min(i,K);
    elseif (i == j) then
        A_ij = -lamb -min(i,K);
    else
        A_ij = 0;
    end
endfunction;

function [X_t]=queue(lamb,mu,tpsmax)
    X_t = zeros(2,1)
    tps = 0;
    i = 2;
    while tps < tpsmax 
        if (X_t(2,i-1) == 0.) then
            temp = expo(lamb);
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            X_t(2,i) = X_t(2,i-1)+1;
        else
            temp = expo(lamb + mu);
            tps = tps + temp
            X_t(1,i) = X_t(1,i-1) + temp;
            rand("uniform");
            r = rand();
            if r < lamb/(lamb + mu) then
                X_t(2,i) = X_t(2,i-1)+1;
            else
                X_t(2,i) = X_t(2,i-1)-1;
            end
        end
        i = i+1;
    end
endfunction

function [dist] = distribution(t)
    X_t = queue(1,2,t)
    dist = 1/t * intg(0,t,X_t);
endfunction

X_t = queue(2,1,30);
n=[1:100]
plot(distribution(n))
