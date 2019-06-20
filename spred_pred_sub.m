
myrm2(0,ss3,a1,a2,b1,b2,d1,d2)

%% Determining stability by Jacobian

npt = 500; %Size of parameter sampling range (not too large or it'll take forever!)

%Starting parameter values
%a1=3; a2=4; b1=0.5; b2=1; d1=1; d2=1; %Default
a1=5; a2=1/10; b1=3; b2=2; d1=2/5; d2=1/100; %chaos
minpr = 0.2; maxpr = 7;
minpr2 = 2; maxpr2 = 5;

prange = linspace(minpr,maxpr,npt); %Parameter range (can change)
prange2 = linspace(minpr2,maxpr2,npt);



mina1 = 0.2; maxa1=15; mind1=0.1; maxd1=6;

prangea1 = linspace(mina1,maxa1,npt); %a1 prange
pranged1 = linspace(mind1,maxd1,npt); %d1 prange

M= zeros(npt,npt,3); %Initial colour matrix
%MM = zeros(npt,npt);

%Expressions for steady states in terms of parameters (and other steady
%states)

xs2 = @(a1,b1,d1) (d1./(a1-b1.*d1));
ys2 = @(a1,b1,d1) (1-(d1./(a1-b1.*d1))).*(1+b1.*(d1./(a1-b1.*d1)))./(a1.*(d1./(a1-b1.*d1)));

xs3 = @(a1,b1,a2,b2,d2) ((b1-1-a1.*(d2./(a2-d2.*b2))+sqrt((1+a1.*(d2./(a2-d2.*b2))-b1).^2+4*b1))./(2*b1));
ys3 = @(a2,b2,d2) (d2./(a2-d2.*b2));
zs3 = @(a1,b1,d1,a2,b2,d2) (((a1.*((b1-1-a1.*(d2./(a2-d2.*b2))+sqrt((1+a1.*(d2./(a2-d2.*b2))-b1).^2+4*b1))./(2*b1)))./(1+b1.*((b1-1-a1.*(d2./(a2-d2.*b2))+sqrt((1+a1.*(d2./(a2-d2.*b2))-b1).^2+4*b1))./(2*b1)))-d1).*((1+b2.*(d2./(a2-d2.*b2)))./(a2)));

J = @(a1,b1,d1,a2,b2,d2,xs,ys,zs) [-1-(a1.*ys)./((1+b1.*xs).^2), -(a1.*xs)./(1+b1.*xs), 0;
    (a1.*ys)./((1+b1.*xs).^2), (a1.*xs)./(1+b1.*xs)-(a2.*zs)./((1+b2.*ys).^2)-d1, -(a2.*ys)./(1+b2.*ys);
    0, (a2.*zs)./((1+b2.*ys).^2), (a2.*ys)./(1+b2.*ys)-d2];

%Determine value of Jacobian eigenvalues at each fixed point, over a range
%of parameter values

stab = zeros(3,1);

for i=1:npt
    a1 = prangea1(i);
    for j=1:npt
        d1 = pranged1(j);
        
        ss1 = [1; 0; 0];
        x2 = xs2(a1,b1,d1); y2 = ys2(a1,b1,d1);
        ss2 = [x2; y2; 0];
        y3 = ys3(a2,b2,d2); x3 = xs3(a1,b1,a2,b2,d2); z3 = zs3(a1,b1,d1,a2,b2,d2);
        ss3 = [x3; y3; z3];
        
        J1 = J(a1,b1,d1,a2,b2,d2,ss1(1),ss1(2),ss1(3));
        J2 = J(a1,b1,d1,a2,b2,d2,ss2(1),ss2(2),ss2(3));
        J3 = J(a1,b1,d1,a2,b2,d2,ss3(1),ss3(2),ss3(3));
        
        if(any(isnan(ss1)) || any(isnan(ss2)) || any(isnan(ss3)) || any(isinf(ss1)) || any(isinf(ss2)) || any(isinf(ss3)))
            continue
        end
        
        
        %1 if stable (and feasible), 0 if not
        
        stab(1) = ~any(real(eig(J1))>0);
        stab(2) = ~any(real(eig(J2))>0);
        stab(3) = ~any(real(eig(J3))>0);
        
        if(any(ss2<0) || norm(imag(ss2))>1e-14)
            stab(2)=0;
        end
        if(any(ss3<0) || norm(imag(ss3))>1e-14)
            stab(3)=0;
        end
        
        
        if(stab(1)) %Subsidy
            M(npt+1-i,j,1)=116/255;
            M(npt+1-i,j,2)=173/255;
            M(npt+1-i,j,3)=209/255;
        end
        if(stab(2)) %Predator Subsidy
            M(npt+1-i,j,1)=254/255;
            M(npt+1-i,j,2)=224/255;
            M(npt+1-i,j,3)=144/255;
        end
        if(stab(3)) %Superpredator predator subsidy
            M(npt+1-i,j,1)=215/255;
            M(npt+1-i,j,2)=48/255;
            M(npt+1-i,j,3)=39/255;
        end
        
        if(sum(stab)>1)
            d1
            a1
        end
        
    end
end

%Colour key:
%Red - Subsidy survives
%Green - prey/predator coexistence
%Blue - three species coexistence
%If we have multi-stability, it will be a combination of these colours

imagesc(M);
set(gca,'XTick',linspace(0,npt,9));
set(gca,'XTickLabel', linspace(minpr,maxpr,9),'FontSize',16);
set(gca,'YTick',linspace(0,npt,9));
set(gca,'YTickLabel', linspace(maxpr,minpr,9),'FontSize',16);
xlabel('d1','FontSize',16)
ylabel('a1','FontSize',16)



%% Single run of ODE45

a1=2; a2=1/10; b1=3; b2=2; d1=0.495; d2=1/100; %chaos


options = odeset('RelTol',1e-11,'AbsTol',1e-11);

t_fin = 50000; %Final time

dt=1;
tspan = [0:dt:t_fin]; %Time span (consider reducing if it takes too long to run, though be wary of 
y0 = [0.5 0.5 0.5]; %Initial state

[t,x] = ode45(@(t,y) myrm2(t,y,a1,a2,b1,b2,d1,d2), tspan, y0, options);



plot(t,x)
legend('subsidy','predator','super-predator')


xs2 = @(a1,b1,d1) (d1./(a1-b1.*d1));
ys2 = @(a1,b1,d1) (1-(d1./(a1-b1.*d1))).*(1+b1.*(d1./(a1-b1.*d1)))./(a1.*(d1./(a1-b1.*d1)));

xs3 = @(a1,b1,a2,b2,d2) ((b1-1-a1.*(d2./(a2-d2.*b2))+sqrt((1+a1.*(d2./(a2-d2.*b2))-b1).^2+4*b1))./(2*b1));
ys3 = @(a2,b2,d2) (d2./(a2-d2.*b2));
zs3 = @(a1,b1,d1,a2,b2,d2) (((a1.*((b1-1-a1.*(d2./(a2-d2.*b2))+sqrt((1+a1.*(d2./(a2-d2.*b2))-b1).^2+4*b1))./(2*b1)))./(1+b1.*((b1-1-a1.*(d2./(a2-d2.*b2))+sqrt((1+a1.*(d2./(a2-d2.*b2))-b1).^2+4*b1))./(2*b1)))-d1).*((1+b2.*(d2./(a2-d2.*b2)))./(a2)));

J = @(a1,b1,d1,a2,b2,d2,xs,ys,zs) [-1-(a1.*ys)./((1+b1.*xs).^2), -(a1.*xs)./(1+b1.*xs), 0;
    (a1.*ys)./((1+b1.*xs).^2), (a1.*xs)./(1+b1.*xs)-(a2.*zs)./((1+b2.*ys).^2)-d1, -(a2.*ys)./(1+b2.*ys);
    0, (a2.*zs)./((1+b2.*ys).^2), (a2.*ys)./(1+b2.*ys)-d2];


        ss1 = [1; 0; 0];
        x2 = xs2(a1,b1,d1); y2 = ys2(a1,b1,d1);
        ss2 = [x2; y2; 0];
        y3 = ys3(a2,b2,d2); x3 = xs3(a1,b1,a2,b2,d2); z3 = zs3(a1,b1,d1,a2,b2,d2);
        ss3 = [x3; y3; z3];

        
        J1 = J(a1,b1,d1,a2,b2,d2,ss1(1),ss1(2),ss1(3));
        J2 = J(a1,b1,d1,a2,b2,d2,ss2(1),ss2(2),ss2(3));
        J3 = J(a1,b1,d1,a2,b2,d2,ss3(1),ss3(2),ss3(3));
        
             
        stab(1) = ~any(real(eig(J1))>0);
        stab(2) = ~any(real(eig(J2))>0);
        stab(3) = ~any(real(eig(J3))>0);
        
        if(any(ss2<0) || norm(imag(ss2))>1e-14)
            stab(2)=0;
        end
        if(any(ss3<0) || norm(imag(ss3))>1e-14)
            stab(3)=0;
        end

        





