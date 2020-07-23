%% Transition Density Function
function TDF = Density(muX, sigmaX, K, J)
    syms a  b  c
    syms xs ys zs
    syms x  y  z
    syms h t s 
    %Change X to Y
    fX2Y=int(1/sigmaX,x);
    fY2X=subs((finverse(fX2Y)), x,y);
    %Y's Drift and Diffusion
    muY_temp=muX/sigmaX-sym('1')/sym('2')*diff(sigmaX,x);
    muY=simplify(subs(muY_temp, x, fY2X));
    sigmaY=sym('1');
    
    %Change Y to Z
    fY2Z=h^(-1/2)*(y-ys);
    
    syms Htemp Expectation 
    sym Beta
    clear Beta Htemp Expectation 
    %slides 25 (7)
    for n=1:K
         HTemp=subs(Hermite(n), z, fY2Z);
         Expectation=HTemp;
         for k=1:J 
           HTemp=muY*diff(HTemp,y,1)+sym('1')/sym('2')*sigmaY*diff(HTemp, y, 2);
           %h = t - s
           Expectation=Expectation + h^k/factorial(k)*HTemp;
         end
         Beta{n}= sym('1')/factorial(n-1) * subs(Expectation, y, ys);
    end
    
    pZ=sym('0');
    for m=1:K
      pZ=pZ+Beta{m}*Hermite(m);
    end
    
    pZ=exp(-z^2/2)/sqrt(2*pi)*pZ;
    pY=(h^(-1/2))*subs(pZ, z, fY2Z);
    pX=(sigmaX^(-1))*subs(pY, y, fX2Y);
    pX=subs(pX, ys, subs(fX2Y, x, xs)) ;
    TDF=simplify(pX);
end
