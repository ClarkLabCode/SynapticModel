function ErrorBarChart(x,y,el,eu, corder)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    bar(x,y, 'EdgeColor','none');
else
    bar(x,y, 'FaceColor','flat','cData', corder);
end

hold on;
errorbar(x,y,y-el,eu-y, 'k','Linewidth',2,'CapSize', 0, 'LineStyle','None');

end

