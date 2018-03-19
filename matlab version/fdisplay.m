function [] = fdisplay(X,Y,A)
% h=surf(X,Y,A);
h = contour(X,Y,A,20);
h2 = colorbar;
% set(h2,'ylim',[0 1]);
view(0,90);
% set(h,'LineStyle','none');
end