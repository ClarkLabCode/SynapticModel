function ConfAxis(fontSize)
ax = gca;
ax.YLabel.FontSize = fontSize;
ax.XLabel.FontSize = fontSize;
set(gca,'FontSize',fontSize)
set(gca,'LineWidth',2);
set(gca,'box','off');
end