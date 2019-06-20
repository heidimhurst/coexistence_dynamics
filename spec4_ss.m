%% s-x-y-z steady states

a1=1.2; a2=4; a3=4.1; b1=0.1; b2=2; b3=0.2; c1=1.1; c2=2; c3=1; d1=0.1; d2=2; %Default - two ss

ysc = @(a3,b3,d2) d2./(a3-b3.*d2);
ssc = @(a2,a3,b2,b3,c1,c2,c3,d2,x) (b2.*b3.*c2.*c3.*d2.*x-b2.*b3.*c2.*d2.*x.^2-a3.*b2.*c2.*c3.*x+a3.*b2.*c2.*x.^2+b3.*c2.*c3.*d2-b3.*c2.*d2.*x-a3.*c2.*c3+a3.*c2.*x+a2.*d2)./(b1.*c2.*(-b3.*c3.*d2+b3.*d2.*x+a3.*c3-a3.*x));
zsc = @(a1,a2,a3,b1,b2,b3,d1,d2,s,x) (-b1.*d1.*s-b2.*d1.*x+a1.*s+a2.*x-d1)./((b1.*s+b2.*x+1).*(-b3.*d2+a3));
xcoeff = [-a1.*b2.*(b3.*d2-a3).*c2.^2;(-(-2*a1.*c3.*b2+a1).*(b3.*d2-a3).*c2.^2+b2.*(b3.*d2-a3).*a2.*c2); (-(a1.*c3.^2.*b2-2*a1.*c3).*(b3.*d2-a3).*c2.^2+((b1.*b3.*c1+a1+b3).*d2+(-b1.*c1-1).*a3-b2.*(b3.*d2-a3).*c3).*a2.*c2); -a1.*c3.^2.*(b3.*d2-a3).*c2.^2-((b1.*b3.*c1+a1+b3).*d2+(-b1.*c1-1).*a3).*c3.*a2.*c2-a2.^2.*d2];
xxx = roots(xcoeff);
yyy=ysc(a3,b3,d2);
sss=ssc(a2,a3,b2,b3,c1,c2,c3,d2,xxx);
zzz=zsc(a1,a2,a3,b1,b2,b3,d1,d2,sss,xxx);


stt1 = [sss(2) xxx(2) yyy zzz(2)];
stt2 = [sss(3) xxx(3) yyy zzz(3)];
%%

a1=3; a2=5; a3=8; b1=5; b2=2; b3=1; c1=1.3; c2=5; c3=3; d1=0.5; d2=2; %Default

ysc = @(a3,b3,d2) d2./(a3-b3.*d2);
ssc = @(a2,a3,b2,b3,c1,c2,c3,d2,x) (b2.*b3.*c2.*c3.*d2.*x-b2.*b3.*c2.*d2.*x.^2-a3.*b2.*c2.*c3.*x+a3.*b2.*c2.*x.^2+b3.*c2.*c3.*d2-b3.*c2.*d2.*x-a3.*c2.*c3+a3.*c2.*x+a2.*d2)./(b1.*c2.*(-b3.*c3.*d2+b3.*d2.*x+a3.*c3-a3.*x));
zsc = @(a1,a2,a3,b1,b2,b3,d1,d2,s,x) (-b1.*d1.*s-b2.*d1.*x+a1.*s+a2.*x-d1)./((b1.*s+b2.*x+1).*(-b3.*d2+a3));
xcoeff = [-a1.*b2.*(b3.*d2-a3).*c2.^2;(-(-2*a1.*c3.*b2+a1).*(b3.*d2-a3).*c2.^2+b2.*(b3.*d2-a3).*a2.*c2); (-(a1.*c3.^2.*b2-2*a1.*c3).*(b3.*d2-a3).*c2.^2+((b1.*b3.*c1+a1+b3).*d2+(-b1.*c1-1).*a3-b2.*(b3.*d2-a3).*c3).*a2.*c2); -a1.*c3.^2.*(b3.*d2-a3).*c2.^2-((b1.*b3.*c1+a1+b3).*d2+(-b1.*c1-1).*a3).*c3.*a2.*c2-a2.^2.*d2];
xxx = roots(xcoeff);
yyy=ysc(a3,b3,d2);
sss=ssc(a2,a3,b2,b3,c1,c2,c3,d2,xxx);
zzz=zsc(a1,a2,a3,b1,b2,b3,d1,d2,sss,xxx);


stt1 = [sss(2) xxx(2) yyy zzz(2)];
stt2 = [sss(3) xxx(3) yyy zzz(3)];

%% s-x-y steady states


a1=0.4; a2=2.5; a3=1; b1=3; b2=2; b3=0.2; c1=1; c2=1; c3=1; d1=2/5; d2=1; %chaos

scoeffs = [a1.^2.*c2-a1.*c2.*d1.*b1; a1.*a2.*c2.*c3+a2.^2-a1.*c2.*d1-d1.*b2.*a1.*c2.*c3-d1.*b2.*a2; -a2.^2.*c1+d1.*b2.*a2.*c1];
x = @(a1,c1,c2,c3,s) (a1.*c2.*c3.*s-a2.*(c1-s))./(a1.*c2.*s);
y = @(a2,b1,b2,c2,c3,s,x) (c3-x).*(1+b1.*s+b2.*x).*c2./a2;

ss = roots(scoeffs);
xx = x(a1,c1,c2,c3,ss);
yy = y(a2,b1,b2,c2,c3,ss,xx);

stt = [ss(1) xx(1) yy(1) 0];








