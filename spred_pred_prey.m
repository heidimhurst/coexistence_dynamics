close()
hold on;

%% Determining stability by Jacobian

npt = 500; %Size of parameter sampling range (not too large or it'll take forever!)

%Starting parameter values

%a1=3; a2=4; b1=0.5; b2=1; d1=1.7; d2=0.1; %Weird stuff

%a1=4; a2=4; b1=10; b2=1; d1=0.1; d2=1; %All exist
a1=3; a2=4; b1=0.5; b2=1; d1=1; d2=1; %Default
%a1=5; a2=1/10; b1=3; b2=2; d1=2/5; d2=1/100; %chaos
minpr = 0.2; maxpr = 7;
%minpr2 = 2; maxpr2 = 5;

prange = linspace(minpr,maxpr,npt); %Parameter range (can change)
%prange2 = linspace(minpr2,maxpr2,npt);



mina1 = 0.1; maxa1=15; mind1=0.1; maxd1=6;

prangea1 = linspace(mina1,maxa1,npt); %a1 prange
pranged1 = linspace(mind1,maxd1,npt); %d1 prange

M= zeros(npt,npt,3); %Initial colour matrix
%MM = zeros(npt,npt);

%Expressions for steady states in terms of parameters (and other steady
%states)

xs3 = @(a1,b1,d1) d1./(a1-b1.*d1);
ys3 = @(a1,b1,d1) (1-d1./(a1-b1.*d1)).*(1+b1.*d1./(a1-b1.*d1))./a1;

xs4 = @(a1,b1,a2,b2,d2) ((b1-1+sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1));
ys4 = @(a2,b2,d2) (d2./(a2-d2.*b2));
zs4 = @(a1,b1,d1,a2,b2,d2) ((a1.*((b1-1+sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))./(1+b1.*((b1-1+sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))-d1).*(1+b2.*(d2./(a2-d2.*b2)))./a2;

xs5 = @(a1,b1,a2,b2,d2) ((b1-1-sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1));
ys5 = @(a2,b2,d2) (d2./(a2-d2.*b2));
zs5 = @(a1,b1,d1,a2,b2,d2) ((a1.*((b1-1-sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))./(1+b1.*((b1-1-sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))-d1).*(1+b2.*(d2./(a2-d2.*b2)))./a2;

J = @(a1,b1,d1,a2,b2,d2,xs,ys,zs) [1-2*xs-((1+b1.*xs).*a1.*ys-a1.*b1.*xs.*ys)./((1+b1.*xs).^2), -(a1.*xs)./(1+b1.*xs), 0;
    ((1+b1.*xs).*a1.*ys-a1.*b1.*xs.*ys)./((1+b1.*xs).^2), (a1.*xs)./(1+b1.*xs)-((1+b2.*ys).*a2.*zs-a2.*b2.*ys.*zs)./((1+b2.*ys).^2)-d1, -(a2.*ys)./(1+b2.*ys);
    0, ((1+b2.*ys).*a2.*zs-a2.*b2.*ys.*zs)./((1+b2.*ys).^2), (a2.*ys)./(1+b2.*ys)-d2];

%Determine value of Jacobian eigenvalues at each fixed point, over a range
%of parameter values

stab = zeros(5,1);

for i=1:npt
    a1 = prange(i);
    for j=1:npt
        a2 = prange(j);
        
        ss1 = [0; 0; 0]; %Don't need to check since it's always unstable
        ss2 = [1; 0; 0];
        x3 = xs3(a1,b1,d1); y3 = ys3(a1,b1,d1);
        ss3 = [x3; y3; 0];
        y4 = ys4(a2,b2,d2); x4 = xs4(a1,b1,a2,b2,d2); z4 = zs4(a1,b1,d1,a2,b2,d2);
        ss4 = [x4; y4; z4];
        y5 = ys5(a2,b2,d2); x5 = xs5(a1,b1,a2,b2,d2); z5 = zs5(a1,b1,d1,a2,b2,d2);
        ss5 = [x5; y5; z5];
        
        if(any(isnan(ss2)) || any(isnan(ss3)) || any(isnan(ss4)) || any(isnan(ss5)) || any(isinf(ss2)) || any(isinf(ss3)) || any(isinf(ss4))|| any(isinf(ss5)))
            continue
        end
        
        J1 = J(a1,b1,d1,a2,b2,d2,ss1(1),ss1(2),ss1(3));
        J2 = J(a1,b1,d1,a2,b2,d2,ss2(1),ss2(2),ss2(3));
        J3 = J(a1,b1,d1,a2,b2,d2,ss3(1),ss3(2),ss3(3));
        J4 = J(a1,b1,d1,a2,b2,d2,ss4(1),ss4(2),ss4(3));
        J5 = J(a1,b1,d1,a2,b2,d2,ss5(1),ss5(2),ss5(3));
        
        %1 if stable (and feasible), 0 if not
        
        stab(1) = ~any(real(eig(J1))>0);
        stab(2) = ~any(real(eig(J2))>0);
        stab(3) = ~any(real(eig(J3))>0);
        stab(4) = ~any(real(eig(J4))>0);
        stab(5) = ~any(real(eig(J5))>0);
        
        if(any(ss3<0) || norm(imag(ss3))>1e-14)
            stab(3)=0;
        end
        if(any(ss4<0) || norm(imag(ss4))>1e-14)
            stab(4)=0;
        end
        if(any(ss5<0) || norm(imag(ss5))>1e-14)
            stab(5)=0;
        end
        
        
        if(stab(2)) %Prey only
            M(npt+1-i,j,1)=49/255;
            M(npt+1-i,j,2)=54/255;
            M(npt+1-i,j,3)=149/255;
        end
        if(stab(3)) %Prey Predator
            M(npt+1-i,j,1)=244/255;
            M(npt+1-i,j,2)=109/255;
            M(npt+1-i,j,3)=67/255;
        end
        if(stab(4)) %3spec
            M(npt+1-i,j,1)=93/255;
            M(npt+1-i,j,2)=0/255;
            M(npt+1-i,j,3)=11/255;
        end
        if(stab(5)) %3spec
            M(npt+1-i,j,1)=93/255;
            M(npt+1-i,j,2)=0/255;
            M(npt+1-i,j,3)=11/255;
        end
        
        if(sum(stab)>1)
            d1
            a1
        end
        
    end
end

%Colour key:
%Dark red (Maroon) - prey survives
%Green - prey/predator coexistence
%Blue - three species coexistence
%Yellowish-brown - three species coexistence (negative solution)
%If we have multi-stability, it will be a combination of these colours

imagesc(M);
set(gca,'XTick',linspace(0,npt,9));
set(gca,'XTickLabel', linspace(minpr,maxpr,9),'fontsize',16);
set(gca,'YTick',linspace(0,npt,9));
set(gca,'YTickLabel', linspace(maxpr,minpr,9),'fontsize',16);
xlabel('a2','fontsize',16)
ylabel('a1','fontsize',16)

%% Single run of ODE45

a1=5; a2=1/10; b1=3; b2=2; d1=2/5; d2=1/100; %chaos

options = odeset('RelTol',1e-11,'AbsTol',1e-11);

t_fin = 50000; %Final time

dt=1;
tspan = [0:dt:t_fin]; %Time span (consider reducing if it takes too long) 
y0 =[0.01 0.01 0.01]; %Initial state

[t,x] = ode45(@(t,y) myrm(t,y,a1,a2,b1,b2,d1,d2), tspan, y0, options);



plot(t,x)
legend('prey','predator','super-predator')

%Testing steady states

xs3 = @(a1,b1,d1) d1./(a1-b1.*d1);
ys3 = @(a1,b1,d1) (1-d1./(a1-b1.*d1)).*(1+b1.*d1./(a1-b1.*d1))./a1;

xs4 = @(a1,b1,a2,b2,d2) ((b1-1+sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1));
ys4 = @(a2,b2,d2) (d2./(a2-d2.*b2));
zs4 = @(a1,b1,d1,a2,b2,d2) ((a1.*((b1-1+sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))./(1+b1.*((b1-1+sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))-d1).*(1+b2.*(d2./(a2-d2.*b2)))./a2;

xs5 = @(a1,b1,a2,b2,d2) ((b1-1-sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1));
ys5 = @(a2,b2,d2) (d2./(a2-d2.*b2));
zs5 = @(a1,b1,d1,a2,b2,d2) ((a1.*((b1-1-sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))./(1+b1.*((b1-1-sqrt((1-b1).^2+4*b1.*(1-a1.*(d2./(a2-d2.*b2)))))./(2*b1)))-d1).*(1+b2.*(d2./(a2-d2.*b2)))./a2;

J = @(a1,b1,d1,a2,b2,d2,xs,ys,zs) [1-2*xs-((1+b1.*xs).*a1.*ys-a1.*b1.*xs.*ys)./((1+b1.*xs).^2), -(a1.*xs)./(1+b1.*xs), 0;
    ((1+b1.*xs).*a1.*ys-a1.*b1.*xs.*ys)./((1+b1.*xs).^2), (a1.*xs)./(1+b1.*xs)-((1+b2.*ys).*a2.*zs-a2.*b2.*ys.*zs)./((1+b2.*ys).^2)-d1, -(a2.*ys)./(1+b2.*ys);
    0, ((1+b2.*ys).*a2.*zs-a2.*b2.*ys.*zs)./((1+b2.*ys).^2), (a2.*ys)./(1+b2.*ys)-d2];


        ss1 = [0; 0; 0]; %Don't need to check since it's always unstable
        ss2 = [1; 0; 0];
        x3 = xs3(a1,b1,d1); y3 = ys3(a1,b1,d1);
        ss3 = [x3; y3; 0];
        y4 = ys4(a2,b2,d2); x4 = xs4(a1,b1,a2,b2,d2); z4 = zs4(a1,b1,d1,a2,b2,d2);
        ss4 = [x4; y4; z4];
        y5 = ys5(a2,b2,d2); x5 = xs5(a1,b1,a2,b2,d2); z5 = zs5(a1,b1,d1,a2,b2,d2);
        ss5 = [x5; y5; z5];
        
               
        J1 = J(a1,b1,d1,a2,b2,d2,ss1(1),ss1(2),ss1(3));
        J2 = J(a1,b1,d1,a2,b2,d2,ss2(1),ss2(2),ss2(3));
        J3 = J(a1,b1,d1,a2,b2,d2,ss3(1),ss3(2),ss3(3));
        J4 = J(a1,b1,d1,a2,b2,d2,ss4(1),ss4(2),ss4(3));
        J5 = J(a1,b1,d1,a2,b2,d2,ss5(1),ss5(2),ss5(3));
        
        %1 if stable (and feasible), 0 if not
        
        stab(1) = ~any(real(eig(J1))>0);
        stab(2) = ~any(real(eig(J2))>0);
        stab(3) = ~any(real(eig(J3))>0);
        stab(4) = ~any(real(eig(J4))>0);
        stab(5) = ~any(real(eig(J5))>0);
        
        if(any(ss3<0) || norm(imag(ss3))>1e-14)
            stab(3)=0;
        end
        if(any(ss4<0) || norm(imag(ss4))>1e-14)
            stab(4)=0;
        end
        if(any(ss5<0) || norm(imag(ss5))>1e-14)
            stab(5)=0;
        end


%%
a1=13; a2=1/10; b1=3; b2=2; d1=2.5; d2=1/100; %chaos









