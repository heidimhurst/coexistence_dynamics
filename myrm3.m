function dxdt = myrm3(t,x,K,I,Gamma,Phi,Psi,P,Q,R,dy,dz)
% x(1) = s, x(2) = x, x(3) = y, x(4) = z

    dxdt = [ I - x(1).*(Gamma - Phi.*x(3)./(1+P.*x(2)+Q.*x(1)));  % ds/dt
            x(2).*(1 - x(2)./K - x(3)./(1+P.*x(2)+Q.*x(1))); % dx/dt
            x(3).*((x(2)+x(1))./(1+P.*x(2)+Q.*x(1)) - dy - x(4)./(1+R.*x(3))); % dy/dt
            x(4).*(Psi.*x(3)./(1+R.*x(3))- dz)]; % dz/dt

end