% Reaction-diffusion equation on a uniformly growing domain with perturbed kinetic parameters

clear all
close all
pwd

% global variables
global prev_nc final_nc final_ncS final_ncP mu_D mu_p mu_d mu_lambda D_mat p_mat d_mat CV k v Lf


% Properties of plot
LineWidth = 3;
FontSize = 18;
lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));


simulate = true;  % if false, read data from file
plotresults = false;
write = false; % save data to table
fit_lambda = false; % fit exponential decay to gradients, otherwise find lambda by calculating the expected x-coordidate
simulate_ss = false; % simulate steady state solution at Lf

% parameters
tol = 1e-10; % numerical tolerance
nruns = 1e2; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
L0 = 50; % initial domain size [µm]
Lf = 300; % final domain size [µm]
mu_D = 1e-3; % mean diffusion coefficient
mu_d = 4e-6; % mean degradation rate
mu_p = mu_d; % mean morphogen production rate
s_range = 0.2; % source length of domain is s_range*L
h_to_s = 3600; % hours to seconds conversion
T = 90 * h_to_s; % total time of domain growth [s]
v = (Lf - L0) / T; % domain growth rate [µm/s]
mu_lambda = sqrt(mu_D/(mu_d+v/Lf)); % mean gradient length at steady state with domain growth speed v and length Lf [µm]
CV = 0.3; % coefficient of variation in kinetic parameters D,d,p
m = 0; % symmetry parameter for pdepe
res = 1000; % resolution for numerical simulation in pdepe
t = linspace(0, T, res); % time steps
x_L = linspace(0, 1, res); % x mesh in Lagrangian framework


t_evaluate = 1:T/h_to_s;  % time points to find x_theta, sigma_x, lambda & C0
plot_gradients_at_t = linspace(10,T/h_to_s,5); % time points to plot pre-steady state gradients 
plot_results_at_t = 10:5:T/h_to_s; % time points to plot x_theta, sigma_x, lambda & C0 


LS = Lf * s_range; % final source length
LP = Lf - LS; % final patterning lenth

% readout positions for finding threshold concentrations
x_theta = [2 4 6] * mu_lambda;

% steady-state analytical deterministic solution at Lf
C_analytic = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));

% threshold concentrations at x_theta's
C_theta = C_analytic(x_theta);

% for each pre-steady state gradient, find the readout positions at theta * C0(t)
theta = C_theta / C_analytic(0);

% array of parameters for PDEPE
P = [L0, v, mu_p, mu_d, mu_D, diameter, s_range, CV];

CVfun = @(x) nanstd(x) ./ nanmean(x);
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

% loop through kinetic parameters
names = {'p', 'd', 'D', 'all'};

% for k = 1:numel(names)   
for k = 4
   

    if simulate

        % filename
        if k == 1
            param_name = 'Prod';
        elseif k == 2
            param_name = 'Decay';
        elseif k == 3
            param_name = 'Diff';
        elseif k == 4
            param_name = 'All';
        end
        title = 'GrowingDomain';
        prefix = ['data/' title '/Perturb' param_name '/'];

        % create folder to save files
        dir = fullfile(pwd, prefix);
        if isfolder(dir)
            disp('The directory already exists.');
        else
            mkdir(dir);
        end


        % absolute positional error
        sigma_x_abs = NaN(length(t_evaluate),length(theta));
        sigma_x_abs_SE = NaN(length(t_evaluate),length(theta));
        
        % positional error normalized by domain length (sigma_x(t) / L(t))
        sigma_x_rel = NaN(length(t_evaluate),length(theta));
        sigma_x_rel_SE = NaN(length(t_evaluate),length(theta));
        
        % average gradient length and amplitude
        avg_lambda = NaN(length(t_evaluate),1);
        avg_lambda_SE = NaN(length(t_evaluate),1);
        avg_C0 = NaN(length(t_evaluate),1);
        avg_C0_SE = NaN(length(t_evaluate),1);

        % variability in gradient length and amplitude
        CV_lambda = NaN(length(t_evaluate),1);
        CV_lambda_SE = NaN(length(t_evaluate),1);
        CV_C0 = NaN(length(t_evaluate),1);
        CV_C0_SE = NaN(length(t_evaluate),1);
        
        % absolute readout position
        x_theta_abs = NaN(nruns, length(t_evaluate), length(theta)); 
            
        % readout position normalized by domain length (x_theta(t) / L(t))
        x_theta_rel = NaN(nruns, length(t_evaluate), length(theta));       

        
        lambda = NaN(nruns,length(t_evaluate));
        C0 = NaN(nruns,length(t_evaluate));
        
        tic
        for j = 1:nruns

            j

            % initialize number of cells in L0
            prev_nc = round(L0 / diameter);
            p_mat = mu_p * ones(prev_nc, 1);
            d_mat = mu_d * ones(prev_nc, 1);
            D_mat = mu_D * ones(prev_nc, 1);
            
            %  draw random kinetic parameters for each cell
            if k == 1 || k == 4
                p_mat = random(logndist(mu_p, mu_p * CV), prev_nc, 1);
            end
            if k == 2 || k == 4
                d_mat = random(logndist(mu_d, mu_d * CV), prev_nc, 1);
            end
            if k == 3 || k == 4
                D_mat = random(logndist(mu_D, mu_D * CV), prev_nc, 1);
            end
            
            % solve the pde
            C = pdepe(m, @(x,t,u,dudx) UniGrowReacDiffPDEfun(x, t, u, dudx, P), ...
                @UniGrowReacDiffICfun, @(xl,ul,xr,ur,t) UniGrowReacDiffBCfun(xl, ul, xr, ur, t), x_L, t);
            
            % simulate steady state solution at Lf 
            if simulate_ss
                % solver initialization
                final_ncS = ceil(LS/diameter);
                final_ncP = ceil(LP/diameter);
                final_nc = final_ncS + final_ncP;
                x0 = (-final_ncS:final_ncP) * diameter;
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
                sol0 = bvpinit(x0, @y0);
                options = bvpset('Vectorized', 'on', 'NMax', 100*final_nc, 'RelTol', tol, 'AbsTol', tol);

                % solve the equation
                C_ss = bvp4c(@odefun, @bcfun, sol0, options);
            end
            

            % find the readout positions at all threshold concentrations
            % for pre-steady state gradients
            for p = 1:length(t_evaluate)

                
                [~, t_idx] = min(abs(t - t_evaluate(p) * h_to_s));
                C_t = C(t_idx,:); % concentration profile at current time 
                L = L0 + v * t(t_idx); % current domain length
                x = x_L * L; % x-coordinate in Eulerian framework

                LS = L * s_range; % current source length
                LP = L - LS; % current patterning length

                if fit_lambda
                    % fit an exponential to the pre-steady state gradients
                    x_high_res = linspace(LS,x(end),1000);
                    y_high_res = interp1(x,C_t,x_high_res);
                    init_C0 = interp1(x,C_t,LS);
                    init_lambda = 0.5 * mu_lambda;
                    expfit = fitnlm(x_high_res-LS, y_high_res, @(p,x)p(1)*exp(-(x/p(2))), [init_C0 init_lambda]);
                    C0(j, p) = expfit.Coefficients.Estimate(1);    
                    lambda(j, p) = expfit.Coefficients.Estimate(2);
                else
                    % instead of fitting, find lambda by calculating the expected x-coordidate
                    x_high_res = linspace(LS, x(end), 1000);
                    y_high_res = interp1(x, C_t, x_high_res);
                    lambda(j,p) = trapz(x_high_res - LS, (x_high_res - LS) .* y_high_res) / trapz(x_high_res - LS, y_high_res); % gradient length                
                    C0(j,p) = interp1(x, C_t, LS); % gradient amplitude
                end

                

                % find readout position at each threshold concentration
                for q = 1:length(theta)

                    C_theta = C0(j,p) * theta(q); % threshold concentration for readout
                    
                    % find all readout positions of C_theta
                    indices = find(diff(sign(C_t - C_theta)));
                    temp_x_theta = [];
                    for idx = indices
                        temp_x_theta = [temp_x_theta, interp1(C_t([idx idx+1]), x([idx idx+1]), C_theta) - LS];
                    end
        
                    % if the return array is empty, the index was out of scope,
                    % meaning the concentration was attained outside of the
                    % domain
                    if ~isempty(temp_x_theta)
                        x_theta_abs(j, p, q) = nanmean(temp_x_theta);
                        x_theta_rel(j, p, q) = nanmean(temp_x_theta) / LP;
                    end


                end

            end

        end
        toc


        % average lambda and C0 at all time points
        avg_lambda = nanmean(lambda)';
        avg_C0 = nanmean(C0)';
        avg_lambda_SE = SEfun(lambda)';   
        avg_C0_SE = SEfun(lambda)';
        for p = 1:length(t_evaluate)
            CV_lambda(p) = CVfun(lambda(:,p));
            CV_lambda_SE(p) = nanstd(bootstrp(nboot, CVfun, lambda(:,p)));
            CV_C0(p) = CVfun(C0(:,p));
            CV_C0_SE(p) = nanstd(bootstrp(nboot, CVfun, C0(:,p)));
        end
        
        % determine the positional error, over the independent runs and 
        % their standard errors from bootstrapping
        for p = 1:length(t_evaluate)
            for q = 1:length(theta)
                sigma_x_rel(p, q) = nanstd(x_theta_rel(:,p,q));
                sigma_x_rel_SE(p, q) = nanstd(bootstrp(nboot, @nanstd, x_theta_rel(:,p,q)));    

                sigma_x_abs(p, q) = nanstd(x_theta_abs(:,p,q));
                sigma_x_abs_SE(p, q) = nanstd(bootstrp(nboot, @nanstd, x_theta_abs(:,p,q)));
            end
        end   


        
        % write data to tables
        if write
            % avg_lambda, avg_C0, CV_lambda, CV_lambda_SE, CV_0, CV_0_SE
            T = table(avg_lambda', avg_lambda_SE', avg_C0', avg_C0_SE',CV_lambda', CV_lambda_SE', CV_C0', CV_C0_SE', ...
                'VariableNames', {'avg_lambda', 'avg_lambda_SE', 'avg_C0', 'avg_C0_SE', 'CV_lambda', 'CV_lambda_SE', 'CV_C0', 'CV_C0_SE'});
            writetable(T, [prefix title '_GradientVariability.csv']);
            
            % C
            [~, ncol] = size(C);
            colNames = cell(1, ncol);
            for i = 1:ncol
                colNames{i} = ['x=' num2str(i/res)];
            end
            T = array2table(C, 'VariableNames', colNames);
            writetable(T, [prefix title '_Concentration.csv']);


            % P
            T = table(P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8), ...
                'VariableNames',{'L0', 'v', 'mu_p', 'mu_d', 'mu_D', 'diameter', 's_range', 'CV'});
            writetable(table(P), [prefix title '_Param.csv']);


            % sigma_x_rel, sigma_x_rel_SE, sigma_x_abs, sigma_x_abs_SE
            for q = 1:length(theta)
                T = table(sigma_x_rel(:,q)',sigma_x_rel_SE(:,q)',sigma_x_abs(:,q)',sigma_x_abs_SE(:,q)', ...
                    'VariableNames', {'sigma_x_rel', 'sigma_x_rel_SE', 'sigma_x_abs', 'sigma_x_abs_SE'});
                writetable(T, [prefix title '_PosError_' num2str(theta(q)) 'theta.csv'])
            end

            % x_theta_rel, x_theta_abs
            for q = 1:length(theta)
                T = table(x_theta_rel(:,:,q), 'VariableNames', {'x_theta_relative'});
                writetable(T, [prefix title '_RelReadPos_' num2str(theta(q)) 'theta.csv']);
            
                T = table(x_theta_abs(:,:,q), 'VariableNames', {'x_theta_abs'});
                writetable(T, [prefix title '_AbsReadPos_' num2str(theta(q)) 'theta.csv']);
            end

        end
    
    else

        % read from table

        T = readtable([prefix title '_GradientVariability.csv']);
        avg_lambda = T.avg_lambda';
        avg_lambda_SE = T.avg_lambda_SE';
        avg_C0 = T.avg_C0';
        avg_C0_SE = T.avg_C0_SE';
        CV_lambda = T.CV_lambda';
        CV_lambda_SE = T.CV_lambda_SE';
        CV_C0 = T.CV_C0';
        CV_C0_SE = T.CV_C0_SE';
        
        T = readtable([prefix title '_Concentration.csv']);
        C = table2array(T);
        
        for q = 1:length(theta)
            T = readtable([prefix title '_PosError_' num2str(theta(q)) 'theta.csv']);
            sigma_x_rel(:,q) = T.sigma_x_rel';
            sigma_x_rel_SE(:,q) = T.sigma_x_rel_SE';
            sigma_x_abs(:,q) = T.sigma_x_abs';
            sigma_x_abs_SE(:,q) = T.sigma_x_abs_SE';            
        end

        for q = 1:length(theta)
            T = readtable([prefix title '_RelReadPos_' num2str(theta(q)) 'theta.csv']);
            x_theta_rel(:,:,q) = table2array(T);
    
            T = readtable([prefix title '_AbsReadPos_' num2str(theta(q)) 'theta.csv']);
            x_theta_abs(:,:,q) = table2array(T);
        end

    end
        
end


%% Plot x_theta and sigma_x

if plotresults

close all

% x_theta_rel
figure(1)
tile1 = tiledlayout(1,2);
ax1 = nexttile(tile1);
axis square
hold on
grid on
for q = 1:length(theta)
    plot(ax1, plot_results_at_t, mean(x_theta_rel(:,plot_results_at_t,q)), '-o', 'MarkerSize', 8, 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(q,length(theta)), 'DisplayName', ['C_\theta = ' num2str(theta(q)) 'C_0'])
end
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Normalized readout position x_\theta/LP(t)', 'FontSize', FontSize)
legend('FontSize', FontSize)

% sigma_x_rel
ax2 = nexttile(tile1);
axis square
hold on
grid on
for q = 1:length(theta)
    errorbar(ax2, plot_results_at_t, sigma_x_rel(plot_results_at_t, q), sigma_x_rel_SE(plot_results_at_t, q), 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(q,length(theta)), 'DisplayName', ['C_\theta = ' num2str(theta(q)) 'C_0'])
end
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Positional error of normalized readout position \sigma_{x/LP(t)}', 'FontSize', FontSize)
legend('FontSize', FontSize)

% x_theta_abs
figure(2)
tile2 = tiledlayout(1,2);
ax1 = nexttile(tile2);
axis square
hold on
grid on
for q = 1:length(theta)
    plot(ax1, plot_results_at_t, mean(x_theta_abs(:,plot_results_at_t,q)), '-o', 'MarkerSize', 8, 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(q,length(theta)), 'DisplayName', ['C_\theta = ' num2str(theta(q)) 'C_0'])
end
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Absolute readout position x_\theta [\mum]', 'FontSize', FontSize)
legend('FontSize', FontSize);


% sigma_x_abs / diameter
ax2 = nexttile(tile2);
axis square
hold on
grid on
for q = 1:length(theta)
    errorbar(ax2, plot_results_at_t, sigma_x_abs(plot_results_at_t, q)./diameter, sigma_x_abs_SE(plot_results_at_t, q)./diameter, 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(q,length(theta)),'DisplayName', ['C_\theta = ' num2str(theta(q)) 'C_0'])    
end
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Relative positional error \sigma_x/\delta [cells]', 'FontSize', FontSize)
legend('FontSize', FontSize)

end


%% Plot pre-steady state gradients

if plotresults


figure(3)
tile2 = tiledlayout(1,3);
ax1 = nexttile(tile2);
axis square;
hold on
grid on
for p = 1:length(plot_gradients_at_t)
    [~, t_idx] = min(abs(t - plot_gradients_at_t(p) * h_to_s));
    L = L0 + v * t(t_idx);
    x = x_L * L;
    [~, x_start_idx] = min(abs(x - s_range * L));
    [~, x_end_idx] = min(abs(x - L));

    plot(ax1, x(1:x_end_idx),C(t_idx, 1:x_end_idx), 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(p,length(plot_gradients_at_t)), 'DisplayName', ['t = ' num2str(plot_gradients_at_t(p)) 'h'])
    xlabel("Domain length L(t) [µm]", "FontSize", FontSize)
    ylabel("Morphogen concentration C(x) [a.u.]", "FontSize", FontSize)
    legend("FontSize", FontSize)
    xlim([0 Lf])
    % set(gca,'YScale','log')
end
LS = Lf * s_range;
LP = Lf-LS;
x_analytic = linspace(-LS, LP, 1000);
C_analytic = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));
plot(C_ss.x+LS, C_ss.y(1,:), 'r--', 'LineWidth', LineWidth, 'DisplayName', 'Noisy steady state solution')


ax2 = nexttile(tile2);
axis square;
hold on
grid on
plot(ax2, NaN, NaN, '--', 'LineWidth', LineWidth, 'Color', 'k', 'DisplayName', 'Source Boundary')
for p = 1:length(plot_gradients_at_t)
    [~, t_idx] = min(abs(t - plot_gradients_at_t(p) * h_to_s));
    L = L0 + v * t(t_idx);
    x = x_L * L;
    [~, x_start_idx] = min(abs(x - s_range * L));
    [~, x_end_idx] = min(abs(x - L));

    plot(ax2, x(1:x_end_idx),C(t_idx, 1:x_end_idx), 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(p,length(plot_gradients_at_t)), 'DisplayName', ['t = ' num2str(plot_gradients_at_t(p)) 'h'])
    plot(ax2, [x(x_start_idx) x(x_start_idx)], [0 C(t_idx, x_start_idx)], '--', 'LineWidth', 0.5*LineWidth, 'Color', 'k', 'HandleVisibility', 'off')
    xlabel("Domain length L(t) [µm]", "FontSize", FontSize)
    ylabel("Morphogen concentration C(x) [a.u.]", "FontSize", FontSize)
    legend("FontSize", FontSize)
    xlim([0 Lf])
    % set(gca,'YScale','log')
end

% normalized pre-steady state gradients
ax3 = nexttile(tile2);
axis square;
hold on
grid on
for q = 1:length(theta)
    plot(ax3, [0,1],[theta(q) theta(q)],'--','LineWidth',0.5*LineWidth,'Color','k','HandleVisibility','off')
end
plot(ax3, NaN,NaN,'--','LineWidth',0.5*LineWidth,'Color','k','DisplayName','C_{\theta}/C_0')
for p = 1:length(plot_gradients_at_t)  
    [~, t_idx] = min(abs(t - plot_gradients_at_t(p) * h_to_s));
    L(p) = L0 + v * t(t_idx);
    x = x_L * L(p);

    [~, x_start_idx] = min(abs(x - s_range * L(p)));
    [~, x_end_idx] = min(abs(x - L(p)));
    plot(ax3, (x(x_start_idx:x_end_idx)-x(x_start_idx)) ./ (x(x_end_idx)-x(x_start_idx)), C(t_idx,x_start_idx:x_end_idx) ./ C(t_idx,x_start_idx), 'LineWidth', LineWidth, 'Color', blueGRADIENTflexible(p,length(plot_gradients_at_t)), 'DisplayName', ['t = ' num2str(plot_gradients_at_t(p)) 'h'])
    xlabel("Normalized patterning domain length LP(t)/LP_{max}", "FontSize", FontSize)
    ylabel("Normalized morphogen concentration C(x)/C_0", "FontSize", FontSize)
    legend("FontSize", FontSize)
    set(gca,'YScale','log')
end


end


%% Plot lambda & C0

if plotresults

% average lambda 
figure(4)
tile3 = tiledlayout(1,3);
ax1 = nexttile(tile3);
axis square;
hold on
grid on
errorbar(ax1, plot_results_at_t, avg_lambda(plot_results_at_t), avg_lambda_SE(plot_results_at_t), 'LineWidth', LineWidth, 'HandleVisibility', 'off')
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Average gradient length \lambda [\mum]', 'FontSize', FontSize)
legend('FontSize', FontSize)

% lambda/L 
ax2= nexttile(tile3);
axis square;
hold on
grid on
for p = 1:length(t_evaluate)
    [~, t_idx] = min(abs(t - t_evaluate(p) * h_to_s));
    L(p) = L0 + v * t(t_idx);
end
plot(ax2, L(plot_results_at_t), avg_lambda(plot_results_at_t)./L(plot_results_at_t), '-o', 'MarkerSize', 8, 'LineWidth', LineWidth)
xlabel('Domain length [\mum]', 'FontSize', FontSize)
ylabel('Scaling factor \lambda/L', 'FontSize', FontSize)
ylim([0 1])


% average C0 VS L(t)
ax3 = nexttile(tile3);
axis square;
hold on
grid on
errorbar(ax3, plot_results_at_t, avg_C0(plot_results_at_t), avg_C0_SE(plot_results_at_t), 'LineWidth', LineWidth, 'HandleVisibility', 'off')
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Average gradient amplitude C_0 [a.u.]', 'FontSize', FontSize)
legend("FontSize", FontSize)

% CV_lambda 
figure(5)
tile4 = tiledlayout(1,2);
ax1 = nexttile(tile4);
axis square;
hold on
grid on
errorbar(ax1, plot_results_at_t, CV_lambda(plot_results_at_t), CV_lambda_SE(plot_results_at_t), 'LineWidth', LineWidth)
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Gradient length variability CV_{\lambda}', 'FontSize', FontSize)

% CV_0
ax2 = nexttile(tile4);
axis square;
hold on
grid on
errorbar(ax2, plot_results_at_t, CV_C0(plot_results_at_t), CV_C0_SE(plot_results_at_t), 'LineWidth', LineWidth)
xlabel('Growth time [h]', 'FontSize', FontSize)
ylabel('Gradient amplitude variability CV_0', 'FontSize', FontSize)

end


%% PDE of advection-dilution-reaction-diffusion equation on uniformly growing domain

function [c, f, s] = UniGrowReacDiffPDEfun(x_L, t, u, dudx, P)

% x_L = linspace(0, 1, res): x mesh in Lagrangian framework
% t = linspace(0, T, res): time steps
% u: concentration
% P: parameters for solving the pde

global prev_nc mu_D mu_p mu_d D_mat p_mat d_mat CV k Lf


L0 = P(1);
v = P(2);
diameter = P(6);
s_range = P(7);


L = L0 + v * t; % current domain length

nc = round(L / diameter); % total number of cells in L
ncS = round(L * s_range / diameter); % number of cells in the source
    
% check if number of cells increased, if yes, randomly pick one cell to
% divide and generate random kinetic parameters for the new cell
if nc > prev_nc 
    dividing_cell = randi([1, length(D_mat)], 1);
    p_cell = mu_p;
    d_cell = mu_d;
    D_cell = mu_D;
    if k == 1 || k == 4
        p_cell = random(logndist(mu_p, mu_p * CV), 1, 1);
    end
    if k == 2 || k == 4
        d_cell = random(logndist(mu_d, mu_d * CV), 1, 1);
    end
    if k == 3 || k == 4
        D_cell = random(logndist(mu_D, mu_D * CV), 1, 1);
    end
    p_mat = vertcat(p_mat(1:dividing_cell), p_cell, p_mat(dividing_cell+1:end));
    d_mat = vertcat(d_mat(1:dividing_cell), d_cell, d_mat(dividing_cell+1:end));
    D_mat = vertcat(D_mat(1:dividing_cell), D_cell, D_mat(dividing_cell+1:end));
    prev_nc = prev_nc + 1; % increment the current number of cells
end

% production rate is zero in the source
p_mat_temp = p_mat;
p_mat_temp(ncS+1:end) = 0;
   
% the current cell index: x_L is in [0,1], thus need x_L * L
cell_idx = ceil(x_L * L / diameter);

% ensure cell_idx does not exceed the size of D_mat
cell_idx = min(cell_idx, length(D_mat));

% kinetic parameters of current cell
D = D_mat(cell_idx);
p = p_mat_temp(cell_idx);
d = d_mat(cell_idx);


% pde: du/dt = D/L(t)^2 * d^u/dx_L^2 - v/L(t) * u + p*H(LS-x) - d*u
c = 1;
f = D / L^2 .* dudx;
s = -v / L * u  + p - d * u;


end

% Initial condition 
function u0 = UniGrowReacDiffICfun(x, c)

u0 = 0;

end

% Zero-flux boundary conditions at x=0 and x=L(t)
function [pl,ql,pr,qr] = UniGrowReacDiffBCfun(xl, ul, xr, ur, t)

pl = 0;
ql = 1;
pr = 0;
qr = 1;   

end


%% functions for the steady-state ODE

% reaction-diffusion equation
function dydx = odefun(x, y, c)

global D_mat p_mat d_mat final_ncS v Lf

dC = -y(2,:) / D_mat(c); % mass flux: j = -D*grad(C)
dj = p_mat(c) * (c <= final_ncS) - (v/Lf + d_mat(c)) * y(1,:); % conservation of mass: div(j) = p*H(-x) - (v/Lf+d)*C
dydx = [dC; dj];

end

% initial guess
function y = y0(x, c)

y = [0; 0];

end

% boundary & cell interface conditions
function res = bcfun(ya, yb)

global final_nc

res = ya(:);
res(1) = ya(2, 1); % zero flux at the left end of the source domain
res(2) = yb(2,final_nc); % zero flux at right end of the patterning domain
for c = 1:final_nc-1
    res(2*c+1) = ya(1,c+1) - yb(1,c); % concentration continuity
    res(2*c+2) = ya(2,c+1) - yb(2,c); % flux continuity
end

end



%% log-normal distribution with adjusted mean & stddev
function pd = logndist(mu, sigma)
    pd = makedist('Lognormal', 'mu', log(mu/sqrt(1+(sigma/mu)^2)), 'sigma', sqrt(log(1+(sigma/mu)^2)));
end
