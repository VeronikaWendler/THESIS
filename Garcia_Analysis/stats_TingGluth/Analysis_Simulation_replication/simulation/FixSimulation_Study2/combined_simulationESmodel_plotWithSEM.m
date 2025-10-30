function combined_simulationESmodel_plotWithSEM(...
    m_qRT, sem_qRT, m_qRT_S, sem_qRT_S, ...
    m_pS, sem_pS, m_RT, sem_RT, m_corrProb, sem_corrProb, ...
    s_qRT, s_qRT_sem, s_qRT_S, s_qRT_S_sem, ...
    s_pS, s_pS_sem, s_RT, s_RT_sem, ...
    s_corrProb, s_corrProb_sem)

figure; set(gcf,'Position',[100 100 800 800])

% Colours for plotting
color_E         = [0 0.749 1];   % deepskyblue
color_S         = [0.6 0 0.6];   % darkorchid
color_sim_E     = color_E;       % same for sim
color_sim_S     = color_S;
color_black     = [0 0 0];       % for lower plots

% Quantile Choice Prob plot
subplot(2,2,1); hold on

% Bar plot
b = bar([m_qRT; m_qRT_S]');
b(1).FaceColor = color_E;  % E empirical
b(2).FaceColor = color_S;  % S empirical
darken = @(c, f) max(0, c * f);
emp_E_err_col   = darken(color_E, 0.6);
emp_S_err_col   = darken(color_S, 0.6);
sim_E_err_col   = darken(color_sim_E, 0.6);
sim_S_err_col   = darken(color_sim_S, 0.6);
errs = [sem_qRT; sem_qRT_S]';
x1 = b(1).XEndPoints;
x2 = b(2).XEndPoints;
errorbar(x1, b(1).YEndPoints, errs(:,1), ...
    'Color', emp_E_err_col, 'CapSize', 8, 'LineStyle', 'none', 'LineWidth', 0.9);
errorbar(x2, b(2).YEndPoints, errs(:,2), ...
    'Color', emp_S_err_col, 'CapSize', 8, 'LineStyle', 'none', 'LineWidth', 0.9);

errorbar(1:5, s_qRT, s_qRT_sem, ':o', ...
    'Color', color_sim_E, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_sim_E, ...
    'MarkerFaceColor', color_sim_E, ...
    'CapSize', 6, 'LineStyle', 'none'); 
plot(1:5, s_qRT, ':o', 'Color', color_sim_E, 'LineWidth', 1.5); 

errorbar(1:5, s_qRT_S, s_qRT_S_sem, ':o', ...
    'Color', color_sim_S, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', color_sim_S, ...
    'MarkerFaceColor', color_sim_S, ...
    'CapSize', 6, 'LineStyle', 'none');
plot(1:5, s_qRT_S, ':o', 'Color', color_sim_S, 'LineWidth', 1.5);

hE_bar     = plot(nan, nan, '-', 'LineWidth', 8, 'Color', color_E);
hS_bar     = plot(nan, nan, '-', 'LineWidth', 8, 'Color', color_S);
hE_simline = plot(nan, nan, ':o', 'Color', color_sim_E, 'LineWidth', 1.5);
hS_simline = plot(nan, nan, ':o', 'Color', color_sim_S, 'LineWidth', 1.5);
legend([hE_bar hS_bar hE_simline hS_simline], ...
       {'E','S','E (sim)','S (sim)'}, ...
       'Box','off','Location','northeast');

set(gca, 'FontSize', 11)  

xlabel('RT Quintile', 'FontSize', 16)
ylabel('Choice prob', 'FontSize', 16)

ylim([0 0.9])
xlim([0.5,5.5])



% 2. Plot: P(choose S) by Dwell-time Advantage (Top Right)
subplot(2,2,2); hold on
errorbar(1:5, m_pS, sem_pS, '.-', 'MarkerSize',20, 'LineWidth',2, 'Color', color_S)
errorbar(1:5, s_pS, s_pS_sem, ':o', 'LineWidth', 1.5, 'Color', color_S)
set(gca,'XTick',1:5,'XTickLabel',{'E>>S','E>S','E~S','S>E','S>>E'})
xlabel('Dwell-time Advantage for S Quintile', 'FontSize', 16)
ylabel('P(choose S)','FontSize', 16)
title('')
xlim([0.5,5.5])
ylim([0.35 0.55])

% 3. Plot: RT by Dwell-time Advantage (Bottom Left)
subplot(2,2,3); hold on
errorbar(1:5, m_RT, sem_RT, '.-', 'MarkerSize',20, 'LineWidth',2, 'Color', color_black)
errorbar(1:5, s_RT, s_RT_sem, ':o', 'LineWidth', 1.5, 'Color', color_black)
set(gca,'XTick',1:5,'XTickLabel',{'E>>S','E>S','E~S','S>E','S>>E'})
xlabel('Dwell-time Advantage for S Quintile', 'FontSize', 14)
ylabel('RT', 'FontSize', 16)
title('')
xlim([0.5,5.5])

% 4. Plot: Accuracy by DwellAdvCorrect (Bottom Right)
subplot(2,2,4); hold on
errorbar(1:5, m_corrProb, sem_corrProb, '.-', 'MarkerSize',20, 'LineWidth',2, 'Color', color_black)
errorbar(1:5, s_corrProb, s_corrProb_sem, ':o', 'LineWidth', 1.5, 'Color', color_black)
xlabel('Dwell-time Advantage for Correct Option', 'FontSize', 14)
ylabel('P(correct)', 'FontSize', 16)
title('')
set(gca,'XTick',1:5,'XTickLabel',{'I>>C','I>C','C~I','C>I','C>>I'})
xlim([0.5,5.5])
ylim([0.75, 0.90]) 

end
