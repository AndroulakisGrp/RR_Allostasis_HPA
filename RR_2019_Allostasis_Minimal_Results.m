%% Evaluate Allostasis Model
% To generate results with a specific point from the nominal (or stressed)
% parametric surface used to generate the Arnold Tongue in the Allostasis paper
load('RR_Allostasis_Minimal_workspace.mat')
k2 = nominal_grid(1,:); % If points from stressed surface are wanted use inter_grid or high_grid
lf = 1; % If light entrainment; lf = 0 when constructing Arnold tongue
entrain = 0; % Change to one if constructing Arnold Tongue
stress = 1; % stress = 1, for nominal; change to 1.5 for intermediate stress, 2 for high stress
[tc,yc]=ode45(@RR_2019_Allostasis_HPA,[0:0.02:6000], y0,[],h1, entrain, k, k2, ksq,1,1, 0,1);
figure, plot(tc, yc(:, 3));
ylabel('Cortisol')
xlabel('Time')
title('Cortisol Rhythm')
%% Points on (Female) Surfaces
% Plot of points from each of the female surfaces used to construct the
% Arnold Tongue. 9 points on each surface have been chosen. 
% Nominal Surface
figure, scatter3(nominal_space(:,1),nominal_space(:,2),nominal_space(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on, 
xlabel('K_p_1')
ylabel('K_p_2')
zlabel('k_p_3')
scatter3(nominal_grid(:,1),nominal_grid(:,2),nominal_grid(:,3),100, 'Filled')

% Intermediate Stress
figure, scatter3(inter_space(:,1),inter_space(:,2),inter_space(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on, 
xlabel('K_p_1')
ylabel('K_p_2')
zlabel('k_p_3')
scatter3(inter_grid(:,1),inter_grid(:,2),inter_grid(:,3),100, 'Filled')

% High Stress
figure, scatter3(high_space(:,1),high_space(:,2),high_space(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on, 
xlabel('K_p_1')
ylabel('K_p_2')
zlabel('k_p_3')
scatter3(high_grid(:,1),high_grid(:,2),high_grid(:,3),100, 'Filled')