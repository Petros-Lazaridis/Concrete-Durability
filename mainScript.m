
%φλορινα  Τ = [ .5, 2.7, 6.7, 11.6, 16.8, 21, 23.1, 22.5, 18.4, 12.6, 7, 2.2 ]
% θε/νικη Τ = [ 5.2, 6.7, 9.7, 14.2, 19.6, 24.4, 26.6, 26, 21.8, 16.2, 11, 6.9 ] 
%ηρακλειο Τ = [ 0.5, 2.7, 6.7, 11.6, 16.8, 21, 23.1, 22.5, 18.4, 12.6, 7, 2.2 ]
%alaska
  Temp = [-4.9, -3, 1.1, 6.9, 13.3, 17.1, 18.6, 17.5, 12.8, 4.7, -2.3, -4]

%U = clorideExposure (.04,3*36000000, 0.004, 360000, 3.5*10^(-12), 0.3, 20, 4);

c = 1 % m
te = 50 %years
dx = 0.005 % in m
dt = .05 % in months
D = 1.2*10^(-12) %m^2/sec
m = 0.3 
maxSurfConc = .8 %%wt conc
maxSurfConcFirstYear = 0 %years to build max syrface concentration
ct = 0.05
cmin = [3 4 5]%cm
%country = 'Alaska'
country = '\Theta\epsilon\sigma\sigma\alpha\lambda\o\nu\iota\kappa\eta'


tic
[U, D, r] = clorideExposure (c, te, dx, dt, D, m, Temp, maxSurfConc,...
 maxSurfConcFirstYear );
toc

size(U)

figure
%set(gcf, 'Units','centimeters', 'Position',[0 0 13 11])
plot( [0:dx:c], U(:,1:200:end) )
hold 
plot([0 c], [ct ct], 'linestyle', '--', 'linewidth', 2, 'color', [0 0 0]) 
set( gca, 'TickLength', [0 0], 'Fontsize', 12, 'Fontname', 'Times',...
'xlim', [0 cmin(end)/100], 'ylim', [0 .8], 'xlabel', 'Depth ( m )',...
'ylabel', '%wt Concentration', 'title', country, 'box', 'off' );
print('clorideProf.epsc')


figure
%set(gcf, 'Units','centimeters', 'Position',[0 0 13 11])
plot( [0:dt:te*12], D, 'linewidth', 1.5) 
set( gca, 'TickLength', [0 0], 'Fontsize', 12, 'Fontname', 'Times',...\
'linewidth', 1.5, 'xlabel', 'Month', 'ylabel', 'D   ( m^2/sec )',...
'title', country, 'box', 'off' );
print('diffusionCoef.epsc')

figure
plot( [0:dt:te*12], U(ceil(cmin./dx/100)+1,:)' )
legend(['c = ', num2str(cmin(1)), ' cm'], ['c = ', num2str(cmin(2)), ' cm'],...
['c = ', num2str(cmin(3)), ' cm'])
hold
plot([0 te*12], [ct ct], 'linestyle', '--', 'linewidth', 2, 'color', [0 0 0]) 
set( gca, 'TickLength', [0 0], 'Fontsize', 12, 'Fontname', 'Times',...
'linewidth', 1.5,'xlabel', 'Months', 'ylabel', '%wt Concentration' ,...
'title', country, 'box', 'off'  );

figure
plot(Temp)
set(gca, 'TickLength', [0 0], 'Fontsize', 12, 'Fontname', 'Times', 'xlabel',...
'Month', 'ylabel', 'Tempreture in C^0', 'linewidth', 1.5,...
'title', country, 'box', 'off');

figure
plot([0, maxSurfConcFirstYear, te]',[0, maxSurfConc, maxSurfConc]', 'linewidth', 1.5)
set(gca, 'TickLength', [0 0], 'Fontsize', 12, 'Fontname', 'Times', 'xlabel',...
'Years', 'ylabel', '%wt of concrete.', 'linewidth', 1.5, 'box', 'off',...
'ylim', [0 maxSurfConc*1.5]);
