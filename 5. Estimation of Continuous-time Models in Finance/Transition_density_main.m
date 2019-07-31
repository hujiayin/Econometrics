%global xs ys zs pX pY pZ 
syms a  b  c
syms xs ys zs
syms x  y  z
syms h t s 
 
%%%%% Drift and Diffusion for X(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vasicek
muX=a*(b-x)
sigmaX=c
%sigmaX=c*sqrt(x)
K = 3
J = 4
density_v = TransitionDensity(muX, sigmaX, K, J);

%%%%% pX   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g1=subs(density_v, {a,b,c,h,xs}, {1,1,2,1/250,1})

%%%%%% Ploting Exact Density for Vasicek %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamm=sigmaX*sqrt(1-exp(-2*a*h))
density_ve=(pi*gamm^2/a)^(-1/2)*exp( -(x-b-(xs-b)*exp(-a*h))^2 *a/(gamm^2) )
g2=subs(density_ve, {a,b,c,h,xs},{1,1,2,1/250,1})
g2=simplify(g2)
gDiff=g1-g2

fig=figure
subplot(2,1,1)
ez1=fplot(g1,[0, 2], 'r-')
hold on
ez2=fplot(g2,[0, 2], 'b:')
legend('Transformation Density','Actual Density')
%%%%%% Plot Density Difference   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
fplot(gDiff, [0,2])
legend('Difference')


% Black Scholes%%%%%%%%%%%%%%%%%%%%
muX_b=a*x
sigmaX_b=b*x
K = 5
J = 6
density_b = TransitionDensity(muX_b, sigmaX_b, K, J);
g3=subs(density_b, {a,b,h,xs},{1,1,1/250,2})


%%%%% Exact Density for Black Schoes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%muX_be=a*x
%sigmaX_be=b*x
densityBS=(sym('1')/sqrt(2*pi*b^2*h)/x)*exp(-(log(x)-log(xs)-(a-b^2)/2*h)^2/(2*b^2*h))
g4=subs(densityBS, {a,b,h,xs},{1,1,1/250,2})
gDiff_bs=g3-g4

fig_bs=figure
subplot(2,1,1)
ezbs1=fplot(g3,[0, 4], 'r-')
hold on
ezbs2=fplot(g4,[0, 4], 'b:')
legend('Transformation Density','Actual Density')
%%%%%% Plot Density Difference   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
fplot(gDiff_bs, [0, 4])
legend('Difference')

syms r k
muX_b= (r - 1/2 * b^2)*x
sigmaX_b = b*x
K = 6
J = 7
density_b = TransitionDensity(muX_b, sigmaX_b, K, J);

syms dif
difxk = piecewise(dif>0, dif, 0)
difxk = subs(difxk, dif, x-k)
expect = difxk * density_b;
expect_num = subs(expect, {k,b,h,r,xs},{34, 0.1, 0.25, 0.05, 32})
int_expect = vpa(int(expect_num, x, 0, Inf))
price_bs = eval(subs(exp(-r*h) * int_expect, {r, h}, {0.05, 0.25}))
%price = simplify(exp(-r*h) * vpa(int(expect,x,[0,1000])));
%price_numerical = eval(subs(price, {k,b,h,r,xs},{30, 0.1, 0.25, 0.05, 31}))


d1 = (log(xs/k) + (r + 1/2 * b^2) * h)/ (b * sqrt(h));
d2 = (log(xs/k) + (r - 1/2 * b^2) * h)/ (b * sqrt(h));
price_c = xs *  normcdf(d1) - k * exp(-r * h) * normcdf(d2);
price_formula = eval(subs(price_c, {k,b,h,r,xs},{34, 0.1, 0.25, 0.05, 32}))
% 1.5232

% CEV Model
syms a  b  c d
syms xs ys zs
syms x  y  z
syms h t s 
muX_cev=a*(b-x)
sigmaX_cev=c*x^d
K = 2
J = 3
density_cev = TransitionDensity(muX_cev, sigmaX_cev, K, J);
g_cev=subs(density_cev, {h},{1/250})


