function ActualvsModel(gdata, gpredicted)
%ACTUALVSMODEL Plot actual glucose data with predicted glucose values

figure;
hold on;
%Plot scatterplot predicted glucose vs actual glucose values
scatter(gpredicted,gdata,15)
xlabel ('G predicted (mM)');
ylabel('G observed (mM)');
plot(linspace(5.5, 9, 500),linspace(5.5, 9, 500),'LineWidth',1)
%calculate correlation coefficient of predicted and actual glucose values
tmp=corrcoef(gdata, gpredicted);
str=sprintf('r= %1.2f',tmp(1,2));
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
axis square
hold off;
end

