%% Evaluate Allostasis Model
% To generate results with a specific point from the nominal (or stressed)
% female surface used to generate the Arnold Tongue in the Allostasis paper
load('RR_031418_Allostasis_Minimal_workspace.mat')
k2 = sob1f1v2_grid(1,:);
lf = 1; % If light entrainment; lf = 0 when constructing Arnold tongue
stress = 1; % Change to 1.5 for intermediate stress, 2 for high stress
[tc,yc]=ode45(@RR_021017l_SexDifferences_CST_Nocturnal_Arnold_v2,[0:0.02:6000], y0,[],h1, cp, k, k2, ksq,1,1, 0,1);
figure, plot(tc, yc(:, 3));
ylabel('Cortisol')
xlabel('Time')
title('Female Cortisol Rhythm')
%% Points on (Female) Surfaces
% Plot of points from each of the female surfaces used to construct the
% Arnold Tongue. 9 points on each surface have been chosen. 
% Nominal Surface
figure, scatter3(sob1f1v2(:,1),sob1f1v2(:,2),sob1f1v2(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on, 
xlabel('K_p_1')
ylabel('K_p_2')
zlabel('k_p_3')
scatter3(sob1f1v2_grid(:,1),sob1f1v2_grid(:,2),sob1f1v2_grid(:,3),100, 'Filled')

% Intermediate Stress
figure, scatter3(sob2f15v2(:,1),sob2f15v2(:,2),sob2f15v2(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on, 
xlabel('K_p_1')
ylabel('K_p_2')
zlabel('k_p_3')
scatter3(sob2f15v2_grid(:,1),sob2f15v2_grid(:,2),sob2f15v2_grid(:,3),100, 'Filled')

% High Stress
figure, scatter3(sob2f2v2(:,1),sob2f2v2(:,2),sob2f2v2(:,3),40,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on, 
xlabel('K_p_1')
ylabel('K_p_2')
zlabel('k_p_3')
scatter3(sob2f2v2_grid(:,1),sob2f2v2_grid(:,2),sob2f2v2_grid(:,3),100, 'Filled')

