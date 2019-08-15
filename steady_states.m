%% Parameters are I, Gamma, Phi, Psi, P, Q, R, dy, dz

%% 1/ Subsidy only

s1 = I./Gamma;
x1 = 0; y1 = 0; z1 = 0;

%% 2/ Prey-Subsidy

s2 = I./Gamma;
x2 = K;
y2 = 0; z2 = 0;


%% 3/ Predator-Subsidy

s3 = dy./(1-dy.*Q);
y3 = (s3.*Gamma - I).*(1+Q.*s3)./(Phi.*s3);
x3 = 0; z3 = 0;

%% 4-5/ Predator-Prey-Subsidy

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

%% 6/ SuperPredator-Predator-Subsidy

y6 = dz./(Psi - R*dz);
s6 = (Phi.*y6 + Q.*I - Gamma + sqrt((Gamma - Q.*I - Phi.*y6).^2 + 4.*I.*Gamma.*Q))/(2.*Gamma.*Q);
z6 = (s6./(1+Q.*s6) - dy).*(1+R.*y6);
x6=0;

%% 7-8-9/ Four Species Coexistence

x4S = [dz.*P.*Phi.*(R.*dz-Phi); dz.*(-(-2.*P.*Phi.^2 + 2.*P.*(R.*dz - Gamma./2).*Phi + dz.*Gamma.*P.*R).*K + Phi.*(R.*dz-Phi)); ...
    dz.*((-P.*Phi.^2 + P.*(R.*dz-Gamma).*Phi + dz.*Gamma.*P.*R).*K.^2 - (-2*Phi.^2 + ((2*R+1).*dz - I.*Q - Gamma) .* Phi + (I.*Q + Gamma).*R.*dz).*K); ...
    dz.*(-Phi.^2 + ((R+1).*dz - I.*Q - Gamma).*Phi + dz.*((I.*Q + Gamma).*R + Gamma)).*K.^2];

r4S = roots(x4S);

if(length(r4S)==1)
    x7 = r4S;
    x8 = -1;
    x9 = -1;
elseif(length(r4S)==2)
    x7 = r4S(1);
    x8 = r4S(2);
    x9 = -1;
else
    x7 = r4S(1);
    x8 = r4S(2);
    x9 = r4S(3);
end

y7 = dz./(Psi - dz.*R);
y8 = y7; y9 = y7;
s7 = (-K.*P.*R.*dz.*x7 + P.*R.*dz.*x7.^2 + K.*P.*Phi.*x7 - P.*Phi.*x7.^2 - K.*R.*dz + R.*dz.*x7 + K.*Phi - dz.*K - Phi.*x7)./...
    Q.*(-K.*R.*dz + R.*dz.*x7 + K.*Phi - Phi.*x7);
s8 = (-K.*P.*R.*dz.*x8 + P.*R.*dz.*x8.^2 + K.*P.*Phi.*x8 - P.*Phi.*x8.^2 - K.*R.*dz + R.*dz.*x8 + K.*Phi - dz.*K - Phi.*x8)./...
    Q.*(-K.*R.*dz + R.*dz.*x8 + K.*Phi - Phi.*x8);
s9 = (-K.*P.*R.*dz.*x9 + P.*R.*dz.*x9.^2 + K.*P.*Phi.*x9 - P.*Phi.*x9.^2 - K.*R.*dz + R.*dz.*x9 + K.*Phi - dz.*K - Phi.*x9)./...
    Q.*(-K.*R.*dz + R.*dz.*x9 + K.*Phi - Phi.*x9);
z7 = - ((P.*dy.*x7 + Q.*dy.*s7 + dy - s7 - x7).*Phi)./((P.*x7 + Q.*s7 + 1).*(Phi - R.*dz));
z8 = - ((P.*dy.*x8 + Q.*dy.*s8 + dy - s8 - x8).*Phi)./((P.*x8 + Q.*s8 + 1).*(Phi - R.*dz));
z9 = - ((P.*dy.*x9 + Q.*dy.*s9 + dy - s9 - x9).*Phi)./((P.*x9 + Q.*s9 + 1).*(Phi - R.*dz));




