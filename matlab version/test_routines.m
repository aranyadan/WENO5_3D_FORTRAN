%% Align
sz = [151,151]; Re=4000;xpos=151;

X=reshape(X,sz); Y=reshape(Y,sz);dy = Y(1,2)-Y(1,1);
u=reshape(u,sz);v=reshape(v,sz);p=reshape(p,sz);rho=reshape(rho,sz);T=reshape(T,sz);

upl = u(xpos,:);ypl=Y(xpos,:)-dy/2;plot(upl,ypl);
neu = sqrt(Re)*(ypl)./sqrt(X(xpos,1));
ind=sum(neu<=5);
    
% Exact soln
neuex = [0,0.4:0.4:4.8,5.0];
uex = [0,0.13277,0.26471,0.39378,0.51676,0.62977,0.72899,0.81152,0.87609,0.92333,0.95552,0.97587,0.98779,0.99155];

clf;plot(upl(2:ind),neu(2:ind)); grid;hold on;scatter(uex,neuex,'r');