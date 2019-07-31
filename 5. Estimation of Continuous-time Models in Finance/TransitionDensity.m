function pX = TransitionDensity(muX,sigmaX, K, J)
%%%%% Transformation X(t) to Y(t)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a  b  c
syms xs ys zs
syms x  y  z
syms h t s 

fX2Y=int(1/sigmaX,x);
fY2X=subs((finverse(fX2Y)), x,y);

%%%%%  Drift and Diffusion for Y(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

muY_temp=muX/sigmaX-sym('1')/sym('2')*diff(sigmaX,x,1);
muY=subs(muY_temp, x, fY2X);
muY=simplify(muY);

sigmaY=sym('1');


%%%%%%  Transformation Y(t) to Z(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fY2Z=h^(-1/2)*(y-ys);
fZ2Y=h^(1/2)*z+ys;


%%%%% Generating Beta   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms Htemp Expectation Beta_t  
clear Beta_t Htemp Expectation  

for k=1:K
     HTemp=subs(Hermite(k), z, fY2Z);
     Expectation=HTemp;

     for j=1:J 
       HTemp=muY*diff(HTemp,y,1)+sym('1')/sym('2')*diff(HTemp, y, 2);
       Expectation=Expectation + h^j/factorial(j)*HTemp;
     end
     Beta_t{k}= sym('1')/factorial(k-1) * subs(Expectation, y, ys);
end

%%%%% Geberating pZ With Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pZ=sym('0');

for m=1:K
  pZ=pZ+Beta_t{m}*Hermite(m);
end
findsym(pZ)

%%%%% Generating pY pX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pZ=exp(-z^2/2)/sqrt(2*pi)*pZ;
pY=(h^(-1/2))*subs(pZ, z, fY2Z);
pX=(sigmaX^(-1))*subs(pY, y, fX2Y);
pX=subs(pX, ys, subs(fX2Y, x, xs));
pX=simplify(pX);
end

