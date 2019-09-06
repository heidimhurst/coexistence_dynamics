close()
hold on;

myrm3(0,stt,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2)

%% Determining stability by Jacobian

npt = 250; %Size of parameter sampling range (not too large or it'll take forever!)

%Starting parameter values

%a1=3; a2=4; b1=0.5; b2=1; d1=1.7; d2=0.1; %Weird stuff

%a1=4; a2=4; b1=10; b2=1; d1=0.1; d2=1; %All exist
%a1=3; a2=4; a3=4; b1=0.5; b2=1; b3=0.1; c1=2; c2=3; c3=1; d1=1; d2=0.1; %??
%a1=4; a2=4; a3=1; b1=3; b2=2; b3=0.2; c1=1; c2=1; c3=1; d1=2/5; d2=1; %GoodDefaultStuff
%a1=4; a2=4; a3=1; b1=3; b2=2; b3=0.1; c1=20; c2=1; c3=1; d1=0.5; d2=0.5; %Sequence

K = 1; I = 0.1; Gamma = 1; Phi = 1; Psi = 0.1; P = 0.5; Q = 2; R = 1.25; dy = 0.1; dz = 0.01;

prange = linspace(minpr,maxpr,npt); %Parameter range (can change)
prange2 = linspace(minpr2,maxpr2,npt);

M= zeros(npt,npt,1); %Initial colour matrix
%MM = zeros(npt,npt);

J = @(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,ss,xs,ys,zs) [-(Gamma-(Phi*ys)/(1+P*xs+Q*ss))-ss*(Q*Phi*ys)/(1+P*xs+Q*ss)^2,...
                                                  -(sx*Phi*ys*P)/(1+P*xs+Q*ss)^2,...
                                                  ss*Phi/(1+P*xs+Q*ss),...
                                                  0;...
                                                  xs*ys*Q/(1+P*xs+Q*ss)^2,...
                                                  1-xs/K-ys/(1+P*xs+Q*ss)+xs*(-1/K+P*ys/(1+P*xs+Q*ss)^2),...
                                                  -xs/(1+P*xs+Q*ss),...
                                                  0;...
                                                  (ys*(1+P*xs+Q*ss)-ys*(xs+ss)*Q)/(1+P*xs+Q*ss)^2,...
                                                  (ys*(1+P*xs+Q*ss)+ys*(xs+ss)*P)/(1+P*xs+Q*ss)^2,...
                                                  (xs+ss)/(1+P*xs+Q*ss)-dy-zs/(1+R*ys)+ys*(R*zs)/(1+R*ys)^2,...
                                                  -ys/(1+R*ys);...
                                                  0,...
                                                  0,...
                                                  zs*((Psi*(1+R*ys)-Psi*ys*R)/(1+R*ys)^2),...
                                                  Psi*ys/(1+R*ys)-dz];       
            
stab = zeros(10,1);

for i=1:npt
    K = prange(i);
    for j=1:npt
        I= prange(j);
        
        % First steady state
        s1 = I./Gamma;
        x1 = 0; y1 = 0; z1 = 0;
        
        st1 = [s1; x1; y1; z1];
        
        % Second steady state
        s2 = I./Gamma;
        x2 = K;
        y2 = 0; z2 = 0;
        
        st2 = [s2; x2; y2; z2];
        
        %Third steady state
        
        s3 = dy./(1-dy.*Q);
        y3 = (s3.*Gamma - I).*(1+Q.*s3)./(Phi.*s3);
        x3 = 0; z3 = 0;
        
        st3 = [s3; x3; y3; z3];
        
        %Fourth/fifth steady states
        xPPS = [(Phi./K) .* (1- P); + ((Gamma - Phi) .* (1- dy - dz.*P) - dy .* (Phi./K)); ...
        I - dy.*(Gamma - Phi) - dy.*Q.*I];

        rPPS = roots(xPPS);

        if(length(rPPS) == 1)
            x4 = rPPS;
            x5 = -1;
            s4 = I./(Gamma - Phi.*(1-x4./K));
            s5 = -1;
            y4 = (1 - x4./K).*(1+P.*x4+Q.*s4);
            y5 = -1;
        else
            x4 = rPPS(1);
            x5 = rPPS(2);
            s4 = I./(Gamma - Phi*(1-x4./K));
            s5 = I./(Gamma - Phi*(1-x5./K));
            y4 = (1 - x4./K).*(1+P.*x4+Q.*s4);
            y5 = (1 - x5./K).*(1+P.*x5+Q.*s5);
        end

        z4 = 0;
        z5 = 0;
        
        st4 = [s4; x4; y4; z4];
        st5 = [s5; x5; y5; z5];
        
        %Sixth steady state
        
        y6 = dz./(Psi - R*dz);
        s6 = (Phi.*y6 + Q.*I - Gamma + sqrt((Gamma - Q.*I - Phi.*y6).^2 + 4.*I.*Gamma.*Q))/(2.*Gamma.*Q);
        z6 = (s6./(1+Q.*s6) - dy).*(1+R.*y6);
        x6=0;
        
        st6 = [s6; x6; y6; z6];
        
        %Seventh steady state is never stable
        st7 = [-1; -1; -1; -1];
        
        %Eighth-tenth steady states
        x4S = [dz.*P.*Phi.*(R.*dz-Phi); dz.*(-(-2.*P.*Phi.^2 + 2.*P.*(R.*dz - Gamma./2).*Phi + dz.*Gamma.*P.*R).*K + Phi.*(R.*dz-Phi)); ...
            dz.*((-P.*Phi.^2 + P.*(R.*dz-Gamma).*Phi + dz.*Gamma.*P.*R).*K.^2 - (-2*Phi.^2 + ((2*R+1).*dz - I.*Q - Gamma) .* Phi + (I.*Q + Gamma).*R.*dz).*K); ...
            dz.*(-Phi.^2 + ((R+1).*dz - I.*Q - Gamma).*Phi + dz.*((I.*Q + Gamma).*R + Gamma)).*K.^2];

        r4S = roots(x4S);

        if(length(r4S)==1)
            x8 = r4S;
            x9 = -1;
            x10 = -1;
        elseif(length(r4S)==2)
            x8 = r4S(1);
            x9 = r4S(2);
            x10 = -1;
        else
            x8 = r4S(1);
            x9 = r4S(2);
            x10 = r4S(3);
        end

        y8 = dz./(Psi - dz.*R);
        y9 = y8; y10 = y8;
        s8 = (-K.*P.*R.*dz.*x8 + P.*R.*dz.*x8.^2 + K.*P.*Phi.*x8 - P.*Phi.*x8.^2 - K.*R.*dz + R.*dz.*x8 + K.*Phi - dz.*K - Phi.*x8)./...
            Q.*(-K.*R.*dz + R.*dz.*x8 + K.*Phi - Phi.*x8);
        s9 = (-K.*P.*R.*dz.*x9 + P.*R.*dz.*x9.^2 + K.*P.*Phi.*x9 - P.*Phi.*x9.^2 - K.*R.*dz + R.*dz.*x9 + K.*Phi - dz.*K - Phi.*x9)./...
            Q.*(-K.*R.*dz + R.*dz.*x9 + K.*Phi - Phi.*x9);
        s10 = (-K.*P.*R.*dz.*x10 + P.*R.*dz.*x10.^2 + K.*P.*Phi.*x10 - P.*Phi.*x10.^2 - K.*R.*dz + R.*dz.*x10 + K.*Phi - dz.*K - Phi.*x10)./...
            Q.*(-K.*R.*dz + R.*dz.*x10 + K.*Phi - Phi.*x10);
        z8 = - ((P.*dy.*x8 + Q.*dy.*s8 + dy - s8 - x8).*Phi)./((P.*x8 + Q.*s8 + 1).*(Phi - R.*dz));
        z9 = - ((P.*dy.*x9 + Q.*dy.*s9 + dy - s9 - x9).*Phi)./((P.*x9 + Q.*s9 + 1).*(Phi - R.*dz));
        z10 = - ((P.*dy.*x10 + Q.*dy.*s10 + dy - s10 - x10).*Phi)./((P.*x10 + Q.*s10 + 1).*(Phi - R.*dz));
        
        st8 = [s8; x8; y8; z8];
        st9 = [s9; x9; y9; z9];
        st10 = [s10; x10; y10; z10];
        
        
        J1  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st1(1),st1(2),st1(3),st1(4));
        J2  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st2(1),st2(2),st2(3),st2(4));
        J3  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st3(1),st3(2),st3(3),st3(4));
        J4  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st4(1),st4(2),st4(3),st4(4));
        J5  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st5(1),st5(2),st5(3),st5(4));
        J6  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st6(1),st6(2),st6(3),st6(4));
       %J7  = J(a1,b1,c1,a2,b2,c2,a3,b3,c3,d1,d2,st7(1),st7(2),st7(3),st7(4));
        J8  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st8(1),st8(2),st8(3),st8(4));
        J9  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st9(1),st9(2),st9(3),st9(4));
        J10 = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st10(1),st10(2),st10(3),st10(4));
        
        if(any(isnan(st3)) || any(isnan(st4)) || any(isnan(st5)) || any(isnan(st6)) || any(isnan(st8)) || any(isnan(st9)) || any(isnan(st10)) || any(isinf(st3)) || any(isinf(st4)) || any(isinf(st5))|| any(isinf(st6))|| any(isinf(st8)) || any(isinf(st9)) || any(isinf(st10)))
            continue
        end
        
        %1 if stable (and feasible), 0 if not
        
        stab(1)  = ~any(real(eig(J1))>0);
        stab(2)  = ~any(real(eig(J2))>0);
        stab(3)  = ~any(real(eig(J3))>0);
        stab(4)  = ~any(real(eig(J4))>0);
        stab(5)  = ~any(real(eig(J5))>0);
        stab(6)  = ~any(real(eig(J6))>0);
       %stab(7)  = ~any(real(eig(J7))>0);
        stab(8)  = ~any(real(eig(J8))>0);
        stab(9)  = ~any(real(eig(J9))>0);
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
  %     if(any(st7<0) || norm(imag(st7))>1e-14)
  %         stab(7)=0;
  %     end
        if(any(st8<0) || norm(imag(st8))>1e-14)
            stab(8)=0;
        end
        if(any(st9<0) || norm(imag(st9))>1e-14)
            stab(9)=0;
        end
        if(any(st10<0) || norm(imag(st10))>1e-14)
            stab(10)=0;
        end
        
        for stab_index=1:length(stab)
            if(stab(stab_index))
                M(npt+1-i,j,1) = M(npt+1-i,j,1) + 10^(stab_index-1);
            end
        end
        
%         if(stab(1)) %Subsidy only
% %             M(npt+1-i,j,1)=116/255;
% %             M(npt+1-i,j,2)=173/255;
% %             M(npt+1-i,j,3)=209/255;
%             M(npt+1-i,j,1) = M(npt+1-i,j,1) + 
%         end
%         if(stab(2)) %Subsidy prey
%             M(npt+1-i,j,1)=69/255;
%             M(npt+1-i,j,2)=117/255;
%             M(npt+1-i,j,3)=180/255;
%         end
%         if(stab(3)) %Subsidy predator
%             M(npt+1-i,j,1)=254/255;
%             M(npt+1-i,j,2)=224/255;
%             M(npt+1-i,j,3)=144/255;
%         end
%         if(stab(4))  %Subsidy prey predator
%             M(npt+1-i,j,1)=253/255;
%             M(npt+1-i,j,2)=174/255;
%             M(npt+1-i,j,3)=97/255;
%         end
%         if(stab(5)) %Subsidy prey predator
%             M(npt+1-i,j,1)=253/255;
%             M(npt+1-i,j,2)=174/255;
%             M(npt+1-i,j,3)=97/255;
%         end
%         if(stab(6)) %Subsidy Predator superpredator
%             M(npt+1-i,j,1)=215/255;
%             M(npt+1-i,j,2)=48/255;
%             M(npt+1-i,j,3)=39/255;
%         end
%  %       if(stab(7))
%  %           M(npt+1-i,j,1)=0.6;
%  %           M(npt+1-i,j,3)=0.6;
%  %       end
%         if(stab(8)) %Quadruple trouble
%             M(npt+1-i,j,1)=145/255;
%             M(npt+1-i,j,2)=0/255;
%             M(npt+1-i,j,3)=88/255;
%         end
%         if(stab(9)) %Quadruple trouble
%             M(npt+1-i,j,1)=145/255;
%             M(npt+1-i,j,2)=0/255;
%             M(npt+1-i,j,3)=88/255;
%         end
%         if(stab(10)) %Quadruple trouble
%             M(npt+1-i,j,1)=145/255;
%             M(npt+1-i,j,2)=0/255;
%             M(npt+1-i,j,3)=88/255;
%         end
%         
%         if(sum(stab)>1)
%             M(npt+1-i,j,1)=1;
%             M(npt+1-i,j,2)=1;
%             M(npt+1-i,j,3)=1;
%         end
 
%        if(sum(stab)>1)
%            a1
%            a2
%        end
        
    end
end

%%
% save output to file based on
outfile = "spps_" + npts + "_" + date;
save(outfile)

%% visualize output

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








