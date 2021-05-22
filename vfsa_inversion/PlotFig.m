function PlotFig(rho,dGcal,xm,zm,iter,RMS)
%PLOTFIG Summary of this function goes here
%   Detailed explanation goes here

  global x;
  global dGobs;
  global m;

  global nx;
  global nz;

  figure(3);
  subplot(2,1,1);
  hold off;
  plot(x,dGobs,'*','color','r','MarkerSize',5);
  hold on
  plot(x,dGcal,'-','color','b','LineWidth',2);
  xlabel(['\bf \fontsize{10}\fontname{Times}Distance (m)']);
  ylabel(['\bf \fontsize{10}\fontname{Times}Gravity Effect (mGal)']);
  leg = legend('Obs','Cal'); set(leg,'Location','South','fontsize',8);
  title(['\bf \fontsize{12}\fontname{Times}VFSA || Misfit : ', num2str(RMS),...
    ' || iter : ', num2str(iter)]);

  xlim([min(x)-m max(x)+m]);

  subplot(2,1,2);
  hold off;
  hold on;

  for i=1:nz
    for j= 1:nx
        xleft(i,j) = xm(i,j);
        xright(i,j) = xm(i,j+1);
        xfill = [xleft(i,j) xright(i,j) xright(i,j) xleft(i,j)];
        zfill = [zm(i,j) zm(i,j+1) zm(i+1,j+1) zm(i+1,j)];
        color = rho(i,j);
        patch(xfill,zfill,color,'EdgeColor','none');
    end
  end

  ylabel(['\bf \fontsize{10}\fontname{Times}Depth (m)']);
  set(gca,'ydir','reverse');
  xlim([min(x) max(x)]);
  axis tight;

  % hp = get(subplot(2,1,2),'Position');
  % % % bla,atas,
  % hcb = colorbar('Position', [hp(1)+hp(3)+0.022  hp(1)+0.022  0.025  hp(2)+hp(3)-0.25]);
  % title(hcb,'g/cm^3','fontweight','bold','fontsize',10);
  % caxis([0 1]);

end