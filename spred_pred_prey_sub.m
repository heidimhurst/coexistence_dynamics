close()
hold on;

myrm3(0,stt,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2)

%% Determining stability by Jacobian

npt = 250; %Size of parameter sampling range (not too large or it'll take forever!)

%Starting parameter values

%a1=3; a2=4; b1=0.5; b2=1; d1=1.7; d2=0.1; %Weird stuff

%a1=4; a2=4; b1=10; b2=1; d1=0.1; d2=1; %All exist
%a1=3; a2=4; a3=4; b1=0.5; b2=1; b3=0.1; c1=2; c2=3; c3=1; d1=1; d2=0.1; %??
a1=4; a2=4; a3=1; b1=3; b2=2; b3=0.2; c1=1; c2=1; c3=1; d1=2/5; d2=1; %GoodDefaultStuff
%a1=4; a2=4; a3=1; b1=3; b2=2; b3=0.1; c1=20; c2=1; c3=1; d1=0.5; d2=0.5; %Sequence

minpr = 0.2; maxpr = 15;
minpr2 = 2; maxpr2 = 5;

prange = linspace(minpr,maxpr,npt); %Parameter range (can change)
prange2 = linspace(minpr2,maxpr2,npt);



mina1 = 0.1; maxa1=15; mind1=0.1; maxd1=6;

prangea1 = linspace(mina1,maxa1,npt); %a1 prange
pranged1 = linspace(mind1,maxd1,npt); %d1 prange

M= zeros(npt,npt,3); %Initial colour matrix
%MM = zeros(npt,npt);

%Expressions for steady states in terms of parameters (and other steady
%states)

ss3 = @(a1,b1,d1) d1/(a1-b1*d1);
ys3 = @(a1,b1,c1,ss) (b1*ss+1)*(c1-ss)/(a1*ss);

ss4c = @(a1,a2,b1,b2,c1,c2,c3,d1) [a1.^2.*c2-a1.*c2.*d1.*b1; a1.*a2.*c2.*c3+a2.^2-a1.*c2.*d1-d1.*b2.*a1.*c2.*c3-d1.*b2.*a2; -a2.^2.*c1+d1.*b2.*a2.*c1];
xs4 = @(a1,a2,c1,c2,c3,s) (a1.*c2.*c3.*s-a2.*(c1-s))./(a1.*c2.*s);
ys4 = @(a2,b1,b2,c2,c3,s,x) (c3-x).*(1+b1.*s+b2.*x).*c2./a2;

ys6 = @(a3,b3,d2) d2/(a3-b3*d2);
ss6 = @(a1,b1,c1,ys) (b1*c1-a1*ys-1+sqrt(b1^2*c1^2+a1^2*ys^2+2*b1*c1-2*(a1*b1*c1-a1)*ys+1))/(2*b1);
zs6 = @(a1,b1,a3,b3,d1,d2,ss) ((b1*d1-a1)*ss+d1)/(b3*d2+(b1*b3*d2-a3*b1)*ss-a3);

%ys7 = @(a3,b3,d2) d2/(a3-b3*d2);
%ss7 = @(a1,b1,c1,ys) (b1*c1-a1*ys-1-sqrt(b1^2*c1^2+a1^2*ys^2+2*b1*c1-2*(a1*b1*c1-a1)*ys+1))/(2*b1);
%zs7 = @(a1,b1,a3,b3,d1,d2,ss) ((b1*d1-a1)*ss+d1)/(b3*d2+(b1*b3*d2-a3*b1)*ss-a3);

ys8 = @(a3,b3,d2) d2./(a3-b3.*d2);
ss8 = @(a2,a3,b2,b3,c1,c2,c3,d2,x) (b2.*b3.*c2.*c3.*d2.*x-b2.*b3.*c2.*d2.*x.^2-a3.*b2.*c2.*c3.*x+a3.*b2.*c2.*x.^2+b3.*c2.*c3.*d2-b3.*c2.*d2.*x-a3.*c2.*c3+a3.*c2.*x+a2.*d2)./(b1.*c2.*(-b3.*c3.*d2+b3.*d2.*x+a3.*c3-a3.*x));
zs8 = @(a1,a2,a3,b1,b2,b3,d1,d2,s,x) (-b1.*d1.*s-b2.*d1.*x+a1.*s+a2.*x-d1)./((b1.*s+b2.*x+1).*(-b3.*d2+a3));
xs8c = @(a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2) [-a1.*b2.*(b3.*d2-a3).*c2.^2;(-(-2*a1.*c3.*b2+a1).*(b3.*d2-a3).*c2.^2+b2.*(b3.*d2-a3).*a2.*c2); (-(a1.*c3.^2.*b2-2*a1.*c3).*(b3.*d2-a3).*c2.^2+((b1.*b3.*c1+a1+b3).*d2+(-b1.*c1-1).*a3-b2.*(b3.*d2-a3).*c3).*a2.*c2); -a1.*c3.^2.*(b3.*d2-a3).*c2.^2-((b1.*b3.*c1+a1+b3).*d2+(-b1.*c1-1).*a3).*c3.*a2.*c2-a2.^2.*d2];

J = @(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,ss,xs,ys,zs)[-1-((1+b1.*ss+b2.*xs).*a1.*ys-a1.*b1.*ss.*ys)./(1+b1.*ss+b2.*xs).^2, (a1.*b2.*ss.*ys)./((1+b1.*ss+b2.*xs).^2),(-a1.*ss)./((1+b1.*ss+b2.*xs)),0;
    (a2.*b1.*xs.*ys)./((1+b1.*ss+b2.*xs).^2), c2.*c3-2*xs.*c2-((1+b1.*ss+b2.*xs).*a2.*ys-a2.*b2.*xs.*ys)./((1+b1.*ss+b2.*xs).^2), -(a2.*xs)./((1+b1.*ss+b2.*xs)), 0;
    ((1+b1.*ss+b2.*xs).*a1.*ys-(a1.*ss.*ys+a2.*xs.*ys).*b1)./(1+b1.*ss+b2.*xs).^2, ((1+b1.*ss+b2.*xs).*a2.*ys-(a1.*ss.*ys+a2.*xs.*ys).*b2)./(1+b1.*ss+b2.*xs).^2, (a1.*ss+a2.*xs)./(1+b1.*ss+b2.*xs) - ((1+b3.*ys).*a3.*zs-a3.*b3.*ys.*zs)./(1+b3.*ys).^2-d1, -(a3.*ys)./(1+b3.*ys);
    0, 0, ((1+b3.*ys).*a3.*zs-a3.*b3.*ys.*zs)./(1+b3.*ys).^2, (a3.*ys)./(1+b3.*ys)-d2];
            
            
stab = zeros(8,1);

for i=1:npt
    a1 = prange(i);
    for j=1:npt
        a2 = prange(j);
        
        st1 = [c1; 0; 0; 0]; 
        st2 = [c1; c3; 0; 0];
        s3 = ss3(a1,b1,d1); y3 = ys3(a1,b1,c1,s3);
        st3 = [s3; 0; y3; 0];
        
        ss4t = roots(ss4c(a1,a2,b1,b2,c1,c2,c3,d1));
        xs4t =  xs4(a1,a2,c1,c2,c3,ss4t); ys4t = ys4(a2,b1,b2,c2,c3,ss4t,xs4t);
        
        st4 = [ss4t(1); xs4t(1); ys4t(1); 0];
        if(length(ss4t)==1)
            st5 = [-1; -1; -1; -1];
        else
            st5 = [ss4t(2); xs4t(2); ys4t(2); 0];
        end
        
        y6 = ys6(a3,b3,d2); s6 = ss6(a1,b1,c1,y6); z6 = zs6(a1,b1,a3,b3,d1,d2,s6);
        st6 = [s6; 0; y6; z6];
 %       y7 = ys7(a3,b3,d2); s7 = ss7(a1,b1,c1,y6); z7 = zs7(a1,b1,a3,b3,d1,d2,s6);
 %       st7 = [s7; 0; y7; z7];
        
        xs8t = roots(xs8c(a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2));
        ys8t = ys8(a3,b3,d2); ss8t = ss8(a2,a3,b2,b3,c1,c2,c3,d2,xs8t); zs8t = zs8(a1,a2,a3,b1,b2,b3,d1,d2,ss8t,xs8t);
        st8 = [ss8t(1); xs8t(1); ys8t; zs8t(1)];
        if(length(xs8t)==1)
            st9 = [-1; -1; -1; -1];
            st10 = [-1; -1; -1; -1];
        elseif(length(xs8t)==2)
            st9 = [ss8t(2); xs8t(2); ys8t; zs8t(2)];
            st10 = [-1; -1; -1; -1];
        else
            st9 = [ss8t(2); xs8t(2); ys8t; zs8t(2)];
            st10 = [ss8t(3); xs8t(3); ys8t; zs8t(3)];
        end
        
        
        
        J1 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st1(1),st1(2),st1(3),st1(4));
        J2 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st2(1),st2(2),st2(3),st2(4));
        J3 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st3(1),st3(2),st3(3),st3(4));
        J4 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st4(1),st4(2),st4(3),st4(4));
        J5 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st5(1),st5(2),st5(3),st5(4));
        J6 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st6(1),st6(2),st6(3),st6(4));
    %    J7 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st7(1),st7(2),st7(3),st7(4));
        J8 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st8(1),st8(2),st8(3),st8(4));
        J9 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st9(1),st9(2),st9(3),st9(4));
        J10 = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st10(1),st10(2),st10(3),st10(4));
        
        if(any(isnan(st3)) || any(isnan(st4)) || any(isnan(st5)) || any(isnan(st6)) || any(isnan(st8)) || any(isnan(st9)) || any(isnan(st10)) || any(isinf(st3)) || any(isinf(st4)) || any(isinf(st5))|| any(isinf(st6))|| any(isinf(st8)) || any(isinf(st9)) || any(isinf(st10)))
            continue
        end
        
        %1 if stable (and feasible), 0 if not
        
        stab(1) = ~any(real(eig(J1))>0);
        stab(2) = ~any(real(eig(J2))>0);
        stab(3) = ~any(real(eig(J3))>0);
        stab(4) = ~any(real(eig(J4))>0);
        stab(5) = ~any(real(eig(J5))>0);
        stab(6) = ~any(real(eig(J6))>0);
    %    stab(7) = ~any(real(eig(J7))>0);
        stab(8) = ~any(real(eig(J8))>0);
        stab(9) = ~any(real(eig(J9))>0);
        stab(10) = ~any(real(eig(J10))>0);
        
        if(any(st3<0) || norm(imag(st3))>1e-14)
            stab(3)=0;
        end
        if(any(st4<0) || norm(imag(st4))>1e-14)
            stab(4)=0;
        end
        if(any(st5<0) || norm(imag(st5))>1e-14)
            stab(5)=0;
        end
        if(any(st6<0) || norm(imag(st6))>1e-14)
            stab(6)=0;
        end
  %      if(any(st7<0) || norm(imag(st7))>1e-14)
  %          stab(7)=0;
  %      end
        if(any(st8<0) || norm(imag(st8))>1e-14)
            stab(8)=0;
        end
        if(any(st9<0) || norm(imag(st9))>1e-14)
            stab(9)=0;
        end
        if(any(st10<0) || norm(imag(st10))>1e-14)
            stab(10)=0;
        end
        
        if(stab(1)) %Subsidy only
            M(npt+1-i,j,1)=116/255;
            M(npt+1-i,j,2)=173/255;
            M(npt+1-i,j,3)=209/255;
        end
        if(stab(2)) %Subsidy prey
            M(npt+1-i,j,1)=69/255;
            M(npt+1-i,j,2)=117/255;
            M(npt+1-i,j,3)=180/255;
        end
        if(stab(3)) %Subsidy predator
            M(npt+1-i,j,1)=254/255;
            M(npt+1-i,j,2)=224/255;
            M(npt+1-i,j,3)=144/255;
        end
        if(stab(4))  %Subsidy prey predator
            M(npt+1-i,j,1)=253/255;
            M(npt+1-i,j,2)=174/255;
            M(npt+1-i,j,3)=97/255;
        end
        if(stab(5)) %Subsidy prey predator
            M(npt+1-i,j,1)=253/255;
            M(npt+1-i,j,2)=174/255;
            M(npt+1-i,j,3)=97/255;
        end
        if(stab(6)) %Subsidy Predator superpredator
            M(npt+1-i,j,1)=215/255;
            M(npt+1-i,j,2)=48/255;
            M(npt+1-i,j,3)=39/255;
        end
 %       if(stab(7))
 %           M(npt+1-i,j,1)=0.6;
 %           M(npt+1-i,j,3)=0.6;
 %       end
        if(stab(8)) %Quadruple trouble
            M(npt+1-i,j,1)=145/255;
            M(npt+1-i,j,2)=0/255;
            M(npt+1-i,j,3)=88/255;
        end
        if(stab(9)) %Quadruple trouble
            M(npt+1-i,j,1)=145/255;
            M(npt+1-i,j,2)=0/255;
            M(npt+1-i,j,3)=88/255;
        end
        if(stab(10)) %Quadruple trouble
            M(npt+1-i,j,1)=145/255;
            M(npt+1-i,j,2)=0/255;
            M(npt+1-i,j,3)=88/255;
        end
        
        if(sum(stab)>1)
            M(npt+1-i,j,1)=1;
            M(npt+1-i,j,2)=1;
            M(npt+1-i,j,3)=1;
        end
 
%        if(sum(stab)>1)
%            a1
%            a2
%        end
        
    end
end

imagesc(M);
set(gca,'XTick',linspace(0,npt,9));
set(gca,'XTickLabel', linspace(minpr,maxpr,9),'FontSize',16);
set(gca,'YTick',linspace(0,npt,9));
set(gca,'YTickLabel', linspace(maxpr,minpr,9),'FontSize',16);
xlabel('a2','FontSize',16)
ylabel('a1','FontSize',16)
title('c1=20','FontSize',16)

%%

imagesc(M);
set(gca,'XTick',linspace(0,npt,9));
set(gca,'XTickLabel', linspace(minpr,maxpr,9));
set(gca,'YTick',linspace(0,npt,9));
set(gca,'YTickLabel', linspace(maxpr,minpr,9));
xlabel('c1')
ylabel('a1')
title('Bifurcation Diagram')

            
 

%% Single run of ODE45

%a1=1; a2=8; a3=1; b1=3; b2=2; b3=0.1; c1=1; c2=1; c3=1; d1=0.5; d2=0.5; %GoodDefaultStuff
a1=0.5; a2=5; a3=1/10; b1=0.5; b2=3; b3=2; d1=2/5; d2=1/100; c1=0.5; c2=1; c3=1; %Chaos
options = odeset('RelTol',1e-11,'AbsTol',1e-11);

t_fin = 50000; %Final time

dt=0.2;
tspan = [0:dt:t_fin]; %Time span (consider reducing if it takes too long to run, though be wary of 
y0 = [0.5 0.5 0.5 0.5]; %Initial state

[t,x] = ode45(@(t,y) myrm3(t,y,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2), tspan, y0, options);

%%
clf;
%x = x(10000:end,:);
plot_traj(xx)
title("Combined Model")
xlabel("Prey")
ylabel("Predator")
zlabel("Superpredator")


%%


plot(t,x(:,1)); hold on; plot(t,x(:,4))
legend('subsidy','super-predator')

%%

fsolve(@(x) myrm3(0,x,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2), [0.5 2 10 10])

%% Multistab


a3=1; b1=3; b2=2; b3=0.2; c1=1; c2=1; c3=1; d1=2/5; d2=1; %GoodDefaultStuff
a1=4;
%a2=1.5;
%a2=1.47;
a2=1.44;

npt = 50; %Size of parameter sampling range (not too large or it'll take forever!)

M= zeros(npt,npt,3); %Initial colour matrix

t_fin = 4000; %Final time

dt=1;
tspan = 1:dt:t_fin; %Time span (consider reducing if it takes too long to run, though be wary of 

preyrange = linspace(0.1,1,npt);
predrange = linspace(0.1,1,npt);


eps = 1e-4; %Closeness of convergence


for i=1:npt
    for j=1:npt 
        y0=[0.5, preyrange(npt+1-i), predrange(j), 0.5];
        [t,x] = ode45(@(t,y) myrm3(t,y,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2), tspan, y0, options);
        fin = x(end,:);
        hist = x((end-100:end),:);
        %Check steady states against (1) known value, (2) are they
        %unchanged from the second last entry to the last
        if(norm(fin(2)-0)<eps && (max(hist(:,2))-min(hist(:,2)))<eps) %No Prey
            M(npt+1-i,j,1)=215/255;
            M(npt+1-i,j,2)=48/255;
            M(npt+1-i,j,3)=39/255;
        else %Prey
            M(npt+1-i,j,1)=145/255;
            M(npt+1-i,j,2)=0/255;
            M(npt+1-i,j,3)=88/255;
        end
                
    end
end

imagesc(M);
set(gca,'XTick',linspace(0,npt,6));
set(gca,'XTickLabel', linspace(0,1,6),'FontSize',16);
set(gca,'YTick',linspace(0,npt,6));
set(gca,'YTickLabel', linspace(1,0,6),'Fontsize',16);
xlabel('y0 (Initial Predator Density)','Fontsize',16)
ylabel('x0 (Initial Prey Density)','Fontsize',16)




%%

        y0=[0.5, 0.8, 0.6, 0.5];
        [t,x] = ode45(@(t,y) myrm3(t,y,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2), tspan, y0, options);
        plot(t,x)
        
        %%
function outplot = plot_traj(x1)
    c = 1:length(x1(:,1));
    colormap(jet);
    patch([x1(:,1)' nan],[x1(:,2)' nan],[x1(:,3)' nan],[c nan],'FaceColor','none','EdgeColor','interp')
    hcb = colorbar;
    title(hcb, "time")
    view(3)
end








