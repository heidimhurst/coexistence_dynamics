%% Determining stability by simulation

npt = 100; %Size of parameter sampling range (not too large or it'll take forever!)


M= zeros(npt,npt,3); %Initial colour matrix

t_fin = 200; %Final time

tspan = [0 t_fin]; %Time span (consider reducing if it takes too long to run, though be wary of 
y0 = [0.6 0.6 0.6]; %Initial state

eps = 1e-2; %Closeness of convergence


%Parameter values
a1=3; a2=4; b1=0.5; b2=1; d1=1; d2=1; %Default
prange = linspace(0.2,7,npt); %Parameter range (can change)

options = odeset('RelTol',1e-11,'AbsTol',1e-11);

for i=1:npt
    a1=prange(i); %<-- Change this to one of the other parameters
    for j=1:npt 
        a2=prange(j); %<-- Change this to one of hte other parameters
        [t,x] = ode45(@(t,y) myrm(t,y,a1,a2,b1,b2,d1,d2), tspan, y0, options);
        fin = x(end,:);
        hist = x(end-1,:);
        %Check steady states against (1) known value, (2) are they
        %unchanged from the second last entry to the last
        if(norm(fin(1)-0)<eps && norm(fin-hist)<eps) %Extinction: Black
            M(npt+1-i,j,1)=0;
            M(npt+1-i,j,2)=0;
            M(npt+1-i,j,3)=0;
        elseif(norm(fin(1)-1)<eps && norm(fin-hist)<eps) %Prey: Red
            M(npt+1-i,j,1)=255;
        elseif(norm(fin(1)-(d1./(a1-b1*d1)))<eps && norm(fin-hist)<eps) %Pred/Prey: Blue
            M(npt+1-i,j,3)=255;
        elseif(norm(fin(2)-(d2./(a2-b2*d2)))<eps && norm(fin-hist)<eps) %Spred/Pred/Prey: Green
            M(npt+1-i,j,2)=255;
        else %Who knows?!?: Purple
            M(npt+1-i,j,1)=255;
            M(npt+1-i,j,3)=255;
        end
                
    end
end

%Plot image
imagesc(M);
set(gca,'XTick',[0 20 40 60 80 100]);
set(gca,'XTickLabel', [0.2 1.56 2.92 4.28 5.64 7]);
set(gca,'YTick',[0 20 40 60 80 100]);
set(gca,'YTickLabel', [7 5.64 4.28 2.92 1.56 0.2]);
xlabel('d2')
ylabel('d1')
title('Bifurcation Diagram')



%plot(t,x(:,1),'r-',t,x(:,2),'g-',t,x(:,3),'b-')
%title('Tritrophic RM Model');
%xlabel('Time t');
%ylabel('Solution');
%legend('x','y','z')