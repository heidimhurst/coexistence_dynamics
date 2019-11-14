% Using the fixed points function from Andrew


%% Determining stability by Jacobian

npt = 250; %Size of parameter sampling range (not too large or it'll take forever!)

%Starting parameter values

%a1=3; a2=4; b1=0.5; b2=1; d1=1.7; d2=0.1; %Weird stuff

%a1=4; a2=4; b1=10; b2=1; d1=0.1; d2=1; %All exist
%a1=3; a2=4; a3=4; b1=0.5; b2=1; b3=0.1; c1=2; c2=3; c3=1; d1=1; d2=0.1; %??
%a1=4; a2=4; a3=1; b1=3; b2=2; b3=0.2; c1=1; c2=1; c3=1; d1=2/5; d2=1; %GoodDefaultStuff
%a1=4; a2=4; a3=1; b1=3; b2=2; b3=0.1; c1=20; c2=1; c3=1; d1=0.5; d2=0.5; %Sequence

% K = 1.1; I = 0.2; Gamma = 0.9; Phi = 1; Psi = 0.11; P = 0.5; Q = 2.2; R = 1.25; dy = 0.1; dz = 0.1;
K = 1.1; I = 0.2; Gamma = 0.9; Phi = 1; Psi = 0.11; P = 0.7; Q = 0.8; R = 1.25; dy = 0.1; dz = 0.1;

minpr = 0.01; maxpr = 1;

prange = linspace(minpr,maxpr,npt); %Parameter range (can change)
%prange2 = linspace(minpr2,maxpr2,npt);

M= zeros(npt,npt,1); %Initial colour matrix
%MM = zeros(npt,npt);

% J = @(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,ss,xs,ys,zs) [-(Gamma+(Phi*ys)/(1+P*xs+Q*ss))+ss*(Q*Phi*ys)/(1+P*xs+Q*ss)^2,...
%                                                   +(xs*Phi*ys*P)/(1+P*xs+Q*ss)^2,...
%                                                   -ss*Phi/(1+P*xs+Q*ss),...
%                                                   0;...
%                                                   xs*ys*Q/(1+P*xs+Q*ss)^2,...
%                                                   1-xs/K-ys/(1+P*xs+Q*ss)+xs*(-1/K+P*ys/(1+P*xs+Q*ss)^2),...
%                                                   -xs/(1+P*xs+Q*ss),...
%                                                   0;...
%                                                   (ys*(1+P*xs+Q*ss)-ys*(xs+ss)*Q)/(1+P*xs+Q*ss)^2,...
%                                                   (ys*(1+P*xs+Q*ss)-ys*(xs+ss)*P)/(1+P*xs+Q*ss)^2,...
%                                                   (xs+ss)/(1+P*xs+Q*ss)-dy-zs/(1+R*ys)+ys*(R*zs)/(1+R*ys)^2,...
%                                                   -ys/(1+R*ys);...
%                                                   0,...
%                                                   0,...
%                                                   zs*((Psi*(1+R*ys)-Psi*ys*R)/(1+R*ys)^2),...
%                                                   Psi*ys/(1+R*ys)-dz]; 
                                              
J = @(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,ss,xs,ys,zs) [-Gamma - (Phi*ys*(1+P*xs))/((1+P*xs+Q*ss)^2),...
                                                  +(ss*Phi*ys*P)/(1+P*xs+Q*ss)^2,...
                                                  -ss*Phi/(1+P*xs+Q*ss),...
                                                  0;...
                                                  xs*ys*Q/(1+P*xs+Q*ss)^2,...
                                                  1-2*xs/K-(ys+Q*ss*ys)/((1+P*xs+Q*ss).^2),...
                                                  -xs/(1+P*xs+Q*ss),...
                                                  0;...
                                                  (ys+P*xs*ys-Q*ys*xs)/((1+P*xs+Q*ss)^2),...
                                                  (1+Q*ss*ys-P*ss*ys)/((1+P*xs+Q*ss)^2),...
                                                  (xs+ss)/(1+P*xs+Q*ss)-dy-zs/((1+R*ys)^2),...
                                                  -ys/(1+R*ys);...
                                                  0,...
                                                  0,...
                                                  Psi*zs/((1+R*ys)^2),...
                                                  Psi*ys/(1+R*ys)-dz];  
            
stab = zeros(10,1);

savedvars = zeros(npt,npt,2);

for i=1:npt
    K = prange(i);
    for j=1:npt
        I= prange(j);
        
        RZs = AllRoots(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz);
        
        st1 = RZs{1}';
        st2 = RZs{2}';
        st3 = RZs{3}';
        st4 = RZs{4}';
        st5 = RZs{5}';
        st6 = RZs{6}';
        st7 = RZs{7}';
        st8 = RZs{8}';
        st9 = RZs{9}';
        st10 = RZs{10}';
        
        
        J1  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st1(1),st1(2),st1(3),st1(4));
        J2  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st2(1),st2(2),st2(3),st2(4));
        J3  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st3(1),st3(2),st3(3),st3(4));
        J4  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st4(1),st4(2),st4(3),st4(4));
        J5  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st5(1),st5(2),st5(3),st5(4));
        J6  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st6(1),st6(2),st6(3),st6(4));
        J7  = J(K,I,Gamma,Phi,Psi,P,Q,R,dy,dz,st7(1),st7(2),st7(3),st7(4));
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

        savedvars(npt+1-i,j,1) = prange(i);
        savedvars(npt+1-i,j,2) = prange(j);
        
        if mod(j,30)==0
            j
        end
        
    end
end

%%
% save output to file based on
outfile = "spps_" + npt + "_" + date;
save(outfile)

%% visualize output - color based on binary interpretation
% create output matrix mapping binary matrix to indices
outmat = zeros(npt,npt);
list_of = unique(M); % ALTERNATIVELY: construct independently to cover all cases

for i=1:length(list_of)
    outmat(M==list_of(i)) = i;
end

% plot outmat
imshow(outmat, [min(min(outmat)),max(max(outmat))])
colormap(lines(length(list_of)))
c = colorbar;
c.Ticks = linspace(1.5,length(list_of)-0.5,length(list_of)) 
c.TickLabels = {"none","ss2","ss3","ss4","ss5","ss3,5 (multi)","ss3,4,5 (multi)"};