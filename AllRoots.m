function RZs = AllRoots(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz)
%Compute all roots of the four-species system, sorted vaguely by the number
%of nonzero entries in increasing order.

%1 Species Alive

%Subsidy only
RZs{1} = [I/Gamma,0,0,0];

%2 Species Alive

%Subsidy-Predator
RZs{2} = [-dy/(Q*dy-1),0,(I*Q*dy+Gamma*dy-I)/(Phi*dy*(Q*dy-1)),0];

%Subsidy-Prey
RZs{3} = [I/Gamma,K,0,0];

%3 Species Alive

%Subsidy-Prey-Predator - here we don't have an analytic expression, but we
%do know that s* satisfies a polynomial A*s*^2+B*s+C=0. It has coefficients 
%given by:
A = -Phi*(Q*dy-1)/((P*dy-1)*K);B = -Gamma-Phi*(K+dy/(P*dy-1))/K; C =I;
s = roots([A,B,C]);s = sort(s,'comparisonmethod','real');

%In theory, P(1) should *always* be negative so never be feasible, but do
%check this as you sweep parameters.
xs = @(s)-(Q*dy*s+dy-s)/(P*dy-1); ys=@(s)(-1+(P-Q)*s)*((K*P+Q*s+1)*dy-K-s)/((P*dy-1)^2*K);
RZs{4} = [s(2),xs(s(2)),ys(s(2)),0]; RZs{5} = [s(1),xs(s(1)),ys(s(1)),0];


%Subsidy-Predator-Superpredator - again here we can reduce to a third
%degree polynomial, and I'm fairly certain only one of these is ever
%feasible, but please check as you go along.
A = Q*Gamma*(R*dz-Psi); B=(-I*Q*R+Gamma*R-Phi)*dz+Psi*(I*Q-Gamma);
C = -I*R*dz+I*Psi;
s = roots([A,B,C]);s = sort(s,'comparisonmethod','real');

%Again, I suspect RZs{6} will sometimes be positive and the other two
%always negative, but check as you sweep parameters. You may have to
%determine which, if any, of the following are feasible.
ys = @(s)-(Gamma*s-I)*(Q*s+1)/(Phi*s); zs=@(s)(Gamma*Q*R*s^2-I*Q*R*s+Gamma*R*s-I*R-Phi*s)*(Q*dy*s+dy-s)/(Phi*s*(Q*s+1));
RZs{6} = [s(2),0,ys(s(2)),zs(s(2))];
RZs{7} = [s(1),0,ys(s(1)),zs(s(1))];

%Subsidy-Prey-Predator-Superpredator
ys = dz/(-R*dz+Psi);
ss = @(x)-(-K*P*R*dz*x+P*R*dz*x^2+K*P*Psi*x-P*Psi*x^2-K*R*dz+R*dz*x+K*Psi-K*dz-Psi*x)/(Q*(-K*R*dz+R*dz*x+K*Psi-Psi*x));
zs = @(s,x)-(P*dy*x+Q*dy*s+dy-s-x)*Psi/((P*x+Q*s+1)*(-R*dz+Psi));
A = Phi*(R*dz-Psi)*P; B=-2*P*(Phi+(1/2)*Gamma)*(R*dz-Psi)*K+Phi*(R*dz-Psi);
C = P*(R*dz-Psi)*(Phi+Gamma)*K^2-(((2*R+1)*Phi+(I*Q+Gamma)*R)*dz-(I*Q+Gamma+2*Phi)*Psi)*K;
D = (((R+1)*Phi+(I*Q+Gamma)*R+Gamma)*dz-(I*Q+Gamma+Phi)*Psi)*K^2;

xs = roots([A,B,C,D]);xs = sort(xs,'comparisonmethod','real');

RZs{8} = [ss(xs(3)),xs(3),ys,zs(ss(xs(3)),xs(3))];
RZs{9} = [ss(xs(2)),xs(2),ys,zs(ss(xs(2)),xs(2))];
RZs{10} = [ss(xs(1)),xs(1),ys,zs(ss(xs(1)),xs(1))];
end

