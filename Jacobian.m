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