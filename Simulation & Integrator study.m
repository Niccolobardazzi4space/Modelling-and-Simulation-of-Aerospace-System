% Modelling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: Niccolò Bardazzi 10800456

%% EX 1
clearvars; close all; clc;

% Comment: f(x) = cos(x)-x is a fucntion where -1 < cos(x) < 1, therefore a 
% good initial range would be from -1 since cos(-1) < cos (0) = 1 and 1 
% since cos(1) > cos(pi) = -1

% This can be seen also from the plot:
% x = linspace(-10,10);
% fx = cos(x)-x;
% figure()
% plot(x,fx)
% yline(0)

a = 0; b = 1; range = [a b];
tol = 1.e-8; 
f = @(x) (cos(x)-x);

options = optimset('TolFun',1e-16,'Display', 'off');
analytical_sol = fsolve(f, 1, options);

% Bisection method
BM = bisection_method(f,range,tol);
BM_obj = @() bisection_method(f,range,tol);
BM.time = timeit(BM_obj);

% Secant method
SM = secant_method(f,range,tol);
SM_obj = @() secant_method(f,range,tol);
SM.time = timeit(SM_obj);

% Regula Falsi method
RFM = regula_falsi_method(f,range,tol);
FRM_obj = @() regula_falsi_method(f,range,tol);
RFM.time = timeit(FRM_obj);

[bm_routine, sm_routine, frm_routine] = error_reorder(BM.routine, SM.routine, RFM.routine);

figure()
semilogy(max(eps,abs(bm_routine-analytical_sol)),'LineStyle','-', 'linewidth', 1,'DisplayName', 'Bisection'), hold on
semilogy(max(eps,abs(sm_routine-analytical_sol)),'LineStyle','-.', 'linewidth', 1,'DisplayName', 'Secant'), hold on
semilogy(max(eps,abs(frm_routine-analytical_sol)),'LineStyle','--', 'linewidth', 1,'DisplayName', 'Regula Falsi')
title(['Root finding algorithms with tol = ', num2str(tol)], 'FontSize', 12)
legend show, set(gca, 'XLimSpec', 'Tight'), legend('FontSize', 12)
xlabel('n° of steps', 'FontSize', 13);
ylh = ylabel('|err|','rotation',0,'VerticalAlignment','middle', 'FontSize', 13);
ylh.Position(1) = ylh.Position(1) - 0.8;
grid on

%% EX 2
clearvars; close all; clc; 

figure()
x = linspace(-1,2,100);
f3 = x.^4-2*x.^3+17/16*x.^2-1;
plot(x,f3), hold on, yline(0,'LineStyle','--')
title('Initial guesses', 'FontSize', 12)
ylh = ylabel('$f(x_1)$','Interpreter','latex','rotation',0,'VerticalAlignment','middle', 'FontSize', 13);
ylh.Position(1) = ylh.Position(1) - 0.1;
xlabel('$x_1$','Interpreter','latex', 'FontSize', 13)
legend('$x_1^4-2x_1^3+\frac{17}{16}x_1^2-1$', 'interpreter', 'latex', 'FontSize',14)
hold off

x = sym('x',[2 1]);
f = [x(1)^2-x(1)-x(2)
     x(1)^2/16+x(2)^2-1];

n = 4; % number of iterations

% 1st zero
x0_1 = [  1.5 
         0.75 ];
root_1 = NM_an_FD(f, x, x0_1, n);

% 2nd zero
x0_2 = [ -0.5
         0.75 ];
root_2 = NM_an_FD(f, x, x0_2, n);

%% EX 3
clearvars; close all; clc; 

h = [0.5 0.2 0.05 0.01];
tspan = [0 2];
y0 = 1/2;
y = @(t) t.^2+2*t+1-1/2.*exp(t);
f = @(x,t) x-t^2+1;

for i = 1:5
    a = [h h(end-1:end)./2];
end

% crowding more the interavls with divisors of tspan
rep = 5;
divisor = 2;   % must be divisor of tspan
h_crowded = [h zeros(1,rep*divisor)];
for i = (length(h)+1):length(h_crowded)
    h_crowded(i) = h_crowded(i-2)/divisor; 
end
w = linspace(tspan(1), tspan(2), 1000);
analytical_sol = y(w);

for i = 1:length(h_crowded)
    H.step{i} = h_crowded(i);
    RK_4.step{i} = h_crowded(i);
    [H.y{i}, H.t{i}] = Heun (f, y0, h_crowded(i), tspan);
    H_obj = @() Heun(f, y0, h_crowded(i), tspan);
    H.CPUtime{i} = timeit(H_obj);
    H.int_err_end{i} = H.y{i}(end)-analytical_sol(end);
    [RK_4.y{i}, RK_4.t{i}] = RK4 (f, y0, h_crowded(i), tspan);
    RK4_obj = @() RK4 (f, y0, h_crowded(i), tspan);
    RK_4.CPUtime{i} = timeit(RK4_obj);
    RK_4.int_err_end{i} = RK_4.y{i}(end)-analytical_sol(end);
end

wanted_title = 'Heun integration method';
wanted_legend = 'h';
plotting('', H.step(1:length(h)), H.t(1:length(h)), H.y(1:length(h)), wanted_title, wanted_legend), hold on
plot(w, analytical_sol,'LineWidth',0.7, 'DisplayName','exact')
xlabel('$t$','Interpreter','latex','FontSize', 12)
ylh = ylabel('$x(t)$','Interpreter','latex','rotation',0,'VerticalAlignment','middle','FontSize', 12);
ylh.Position(1) = ylh.Position(1) - 0.1;
Lgnd = legend('show', 'fontsize', 11);
legend('AutoUpdate','off')
Lgnd.Position(1) = 0.18;
Lgnd.Position(2) = 0.58;
xlim([tspan(1) tspan(2)]), ylim([y0 analytical_sol(end)])
zoomPlot(H.t{4}, H.y{4},[1.98 2],[5.27 5.307],[0.6 0.2 0.25 0.3],'-.',[0.4940, 0.1840, 0.5560],[1 3]);
hold on ,plot(H.t{3}, H.y{3},'LineStyle',':','Color', [0.9290, 0.6940, 0.1250]), hold on, plot(w, analytical_sol,'LineWidth',0.7,'Color',[0.4660, 0.6740, 0.1880]), grid on
hold off

wanted_title = 'RK4 integration method';
plotting('', RK_4.step(1:length(h)), RK_4.t(1:length(h)), RK_4.y(1:length(h)), wanted_title, wanted_legend), hold on
plot(w, analytical_sol,'LineWidth',0.7, 'DisplayName','exact')
xlabel('$t$','Interpreter','latex','FontSize', 12)
ylh = ylabel('$x(t)$','Interpreter','latex','rotation',0,'VerticalAlignment','middle','FontSize', 12);
ylh.Position(1) = ylh.Position(1) - 0.1;
Lgnd = legend('show', 'fontsize', 11);
legend('AutoUpdate','off')
Lgnd.Position(1) = 0.18;
Lgnd.Position(2) = 0.58;
xlim([tspan(1) tspan(2)]), ylim([y0 analytical_sol(end)])
zoomPlot(RK_4.t{4}, RK_4.y{4},[1.9999996 2],[5.3054708 5.3054718],[0.6 0.2 0.25 0.3],'-.',[0.4940, 0.1840, 0.5560],[1 3]);
hold on ,plot(RK_4.t{3}, RK_4.y{3},'LineStyle',':','Color', [0.9290, 0.6940, 0.1250])
hold on, plot(w, analytical_sol,'LineWidth',0.7,'Color',[0.4660, 0.6740, 0.1880])
hold on, plot(RK_4.t{2}, RK_4.y{2},'LineStyle','--','Color', [0.8500, 0.3250, 0.0980])
hold on, plot(RK_4.t{1}, RK_4.y{1})
grid on
hold off

figure('Renderer', 'painters', 'Position', [270 300 980 460])
subplot(1,2,1)
scatter(cell2mat(H.CPUtime), cell2mat(H.int_err_end), 'filled')
set(gca,'yscale','log')
title('Heun CPU time vs integration error', 'FontSize', 12)
ylh = ylabel('$err$','Interpreter','latex','rotation',0,'VerticalAlignment','middle', 'FontSize', 12);
ylh.Position(1) = ylh.Position(1) - 3.e-5;
xlabel('$t_{CPU}$', 'Interpreter','latex', 'FontSize', 12)

subplot(1,2,2)
scatter(cell2mat(RK_4.CPUtime), cell2mat(RK_4.int_err_end))
title('RK4 CPU time vs integration error', 'FontSize', 12)
set(gca,'yscale','log')
ylh = ylabel('$err$','Interpreter','latex','rotation',0,'VerticalAlignment','middle', 'FontSize', 12);
ylh.Position(1) = ylh.Position(1) - 3.e-5;
xlabel('$t_{CPU}$', 'Interpreter','latex', 'FontSize', 12)

%% EX 4
clearvars; close all; clc; 

N = 100;
alpha = linspace(0, pi, N);
eig_max1 = zeros(length(alpha),1);
eig_max2 = zeros(length(alpha),1);
hmax_RK = zeros(length(alpha),1);
tol = 1.e-7;
methods = {'RK2' 'RK4'};
range = [-1.e-1 3];
for k = 1:length(methods)
    for i = 1:length(alpha)
        A = A_mat(alpha(i));
        eig_max1(i) = cos(alpha(i)) + 1j*sin(alpha(i));
        eig_max2(i) = cos(alpha(i)+pi) + 1j*sin(alpha(i)+pi);
        hmax_RK(i) = RK_stability(A, tol, range, methods{k});
    end
    integrator.type{k} = methods{k};
    integrator.h_eig_max{k} = [eig_max1; eig_max2].*[hmax_RK; flipud(hmax_RK)];
    integrator.x{k} = 'none';
    integrator.tol{k} = tol;
end
wanted_title = 'Stability domain of RK4 and RK2';
wanted_legend = 'integrator';
plotting('', integrator.type, integrator.x, integrator.h_eig_max, wanted_title, wanted_legend), hold on, 
xlabel('$\Re\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
ylabel('$\Im\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
ylim([-4 4])
xlim([-4 4])
h = [0.5 0.2 0.05 0.01];
eig_A = 1;     % eigenvalues of 1 matrix
scatter_types = {'o', '*', '.', '^'};
for i = 1:length(h)
    scatter(real(h(i).*eig_A), imag(h(i).*eig_A), 'Marker', scatter_types{i},'DisplayName', ['h_', num2str(i), '=', num2str(h(i))]), hold on
end
hold on
Lgnd = legend('show');
legend('AutoUpdate','off')
yline(0), hold on
Lgnd.Position(1) = 0.15;
Lgnd.Position(2) = 0.58;
zoomPlot(real(integrator.h_eig_max{1}), imag(integrator.h_eig_max{1}), [-0.2 0.7],[-0.3 0.3],[0.6 0.2 0.3 0.3],'--',[0.8500, 0.3250, 0.0980],[2 4]);
hold on, grid on
plot(real(integrator.h_eig_max{2}), imag(integrator.h_eig_max{2}), 'Color',[0, 0.4470, 0.7410])
for i = 1:length(h)
    scatter(real(h(i).*eig_A), imag(h(i).*eig_A), 'Marker', scatter_types{i},'DisplayName', ['h_', num2str(i), '=', num2str(h(i))]), hold on
end
yline(0), hold off

%% EX 5 
clearvars; close all; clc; 

x0 = [ 1
       1 ];
tspan = [0 1]; % which will also be max and min stepsize for integration
alpha = linspace(0, pi, 100);
tol = [1.e-3 1.e-4 1.e-5 1.e-6];
maxerr = tol*1.e-2; % max 1% of global error
methods = {'RK1' 'RK2' 'RK4'};
autonomous = 1;
eig_max1 = zeros(length(alpha),1);
eig_max2 = zeros(length(alpha),1);
hmax_RK = zeros(length(alpha),1);
integrators = cell(length(methods),1);
for k = 1:length(methods)
    for m = 1:length(tol)
        for i = 1:length(alpha)
            A = A_mat(alpha(i));
            eig_max1(i) = cos(alpha(i)) + 1j*sin(alpha(i));
            eig_max2(i) = cos(alpha(i)+pi) + 1j*sin(alpha(i)+pi);
            hmax_RK(i) = accuracy_domain(A, x0, tspan, maxerr(m), tol(m), methods{k});
            if alpha(i) == pi
                if strcmp(methods{k}, 'RK1')
                    N = floor((tspan(2)-tspan(1))/hmax_RK(i)) + (rem(tspan(2)-tspan(1),hmax_RK(i))>0);                      
                    integrators{k}.fcn_evals{m} = 1*N;
                elseif strcmp(methods{k}, 'RK2')
                    N = floor((tspan(2)-tspan(1))/hmax_RK(i)) + (rem(tspan(2)-tspan(1),hmax_RK(i))>0);                       
                    integrators{k}.fcn_evals{m} = 2*N;
                elseif strcmp(methods{k}, 'RK4')
                    N = floor((tspan(2)-tspan(1))/hmax_RK(i)) + (rem(tspan(2)-tspan(1),hmax_RK(i))>0);                       
                    integrators{k}.fcn_evals{m} = 4*N;
                end
            end
        end  
        integrators{k}.type = methods{k};
        integrators{k}.tol{m} = tol(m);
        integrators{k}.h_eig_max{m} = [eig_max1; eig_max2].*[hmax_RK; flipud(hmax_RK) ];
        integrators{k}.x{m} = 'none';
    end
    wanted_title = 'Accuracy domain of';
    wanted_legend = 'tol';
    plotting(integrators{k}.type, integrators{k}.tol, integrators{k}.x, integrators{k}.h_eig_max, wanted_title, wanted_legend)
    xlabel('$\Re\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
    ylabel('$\Im\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
    legend show 
    legend('FontSize', 11)
    axis equal
    % Zoom for RK1
    if strcmp(methods{k}, 'RK1')
        lin = {':','-.'};
        colors = {[0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
        legend('AutoUpdate','off', 'FontSize', 11)
        zoomPlot(real(integrators{k}.h_eig_max{3}), imag(integrators{k}.h_eig_max{3}), [-3.e-5 3.e-5], [-3.e-5 3.e-5], [0.58 0.24 0.3 0.3], lin{1}, colors{1}, [2 4]);
        hold on
        plot(real(integrators{k}.h_eig_max{4}), imag(integrators{k}.h_eig_max{4}),'color', colors{2},'Linestyle','-.')
        grid on
    end
end

figure()
loglog(cell2mat(integrators{1}.tol), cell2mat(integrators{1}.fcn_evals), 'Color', [0.4660, 0.6740, 0.1880]), hold on
loglog(cell2mat(integrators{2}.tol), cell2mat(integrators{2}.fcn_evals), 'LineStyle', '--', 'Color', [0.3010, 0.7450, 0.9330]), hold on
loglog(cell2mat(integrators{3}.tol), cell2mat(integrators{3}.fcn_evals), 'LineStyle', '-.', 'Color', [0.6350, 0.0780, 0.1840]), hold off
legend('RK1','RK2','RK4', 'fontsize', 11)
xlabel('tolerance', 'FontSize', 12), ylabel('n° of function evaluations','FontSize', 12)
grid on
title('Price per accuracy for different RKs', 'FontSize', 12)

%% EX 6
clearvars; close all; clc;

N = 100;
alpha = linspace(0, pi, N);
eig_max1 = zeros(length(alpha),1);
eig_max2 = zeros(length(alpha),1);
hmax_BI2 = zeros(length(alpha),1);
tol = 1.e-7;
methods = {'BI2'};
range = [-1.e-1 12];
theta = [0.1 0.3 0.7 0.9 0.4];
method = str2double(erase(methods,'BI'));
for k = 1:length(theta)
    for i = 1:length(alpha)
        A = A_mat(alpha(i));
        eig_max1(i) = cos(alpha(i)) + 1j*sin(alpha(i));
        eig_max2(i) = cos(alpha(i)+pi) + 1j*sin(alpha(i)+pi);
        hmax_BI2(i) = BI_stability(A, theta(k), tol, range, method);
    end
    integratorBI2.type = methods{1};
    integratorBI2.theta{k} = theta(k);
    integratorBI2.tol = tol;
    integratorBI2.h_eig_max{k} = [eig_max1; eig_max2].*[hmax_BI2; flipud(hmax_BI2)];
    integratorBI2.x{k} = 'none';
end
wanted_title = 'Stability domain of';
wanted_legend = '\theta';
plotting(integratorBI2.type, integratorBI2.theta(1:4), integratorBI2.x(1:4), integratorBI2.h_eig_max(1:4), wanted_title, wanted_legend)
hold on
plot(integratorBI2.h_eig_max{5}, 'Marker', '.', 'color','k', 'DisplayName', [wanted_legend,'_5', ' = ', num2str(integratorBI2.theta{5})])
xlabel('$\Re\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
ylabel('$\Im\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
legend show, legend('FontSize', 11)
xlim([-6 12])
axis equal

%% EX 7
clearvars; close all; clc;

B = [-180.5   219.5
      179.5  -220.5];

x0 = [1
      1];
tspan = [0 5];
eigenvalues = eig(B);
h = 0.1;
method_RK = 'RK4';
method_BI = 2;
theta = 0.1;
[yRK4, tRK4] = RK_LTI (B, x0, h, tspan, method_RK);
[yBI2, tBI2] = BI_LTI(B, x0, h, theta, tspan, method_BI);

% Analytical solution 
sol = @(t) expm(B*t);
analytical_sol = zeros(2,length(tBI2));
t = linspace(tspan(1), tspan(2), length(tBI2));
for i = 1:length(tBI2)
    analytical_sol(:,i) = sol(tBI2(i))*x0;
end

figure()
subplot(2,2,1)
plot(tBI2, abs(yBI2(1,:)-analytical_sol(1,:)), 'LineStyle', '-', 'Marker', '.')
legend('$x_1$','Interpreter','latex', 'fontsize', 11)
% xlabel('$t$','Interpreter','latex','FontSize',12')
ylh = ylabel('$err$','Interpreter','latex','FontSize',12','rotation',0,'VerticalAlignment','middle');
ylh.Position(1) = ylh.Position(1) - 0.8;
title(['BI_', num2str(method_BI), 'method']), axis tight

subplot(2,2,2)
plot(tRK4, abs(yRK4(1,:)-analytical_sol(1,:)), 'LineStyle', '-')
% xlabel('$t$','Interpreter','latex','FontSize',12')
legend('$x_1$','Interpreter','latex', 'fontsize', 11)
title([method_RK(1:2),'_', method_RK(3), 'method']), axis tight

subplot(2,2,3)
plot(tBI2, abs(yBI2(2,:)-analytical_sol(2,:)), 'LineStyle', '-', 'Marker', '.','Color',[0.8500, 0.3250, 0.0980])
legend('$x_2$','Interpreter','latex', 'fontsize', 11)
xlabel('$t$','Interpreter','latex','FontSize',12')
ylh = ylabel('$err$','Interpreter','latex','FontSize',12','rotation',0,'VerticalAlignment','middle');
ylh.Position(1) = ylh.Position(1) - 0.6;
axis tight

subplot(2,2,4)
plot(tRK4, abs(yRK4(2,:)-analytical_sol(2,:)), 'LineStyle', '-','Color',[0.8500, 0.3250, 0.0980])
xlabel('$t$','Interpreter','latex','FontSize',12')
legend('$x_2$','Interpreter','latex', 'fontsize', 11), axis tight


N = 100;
alpha = linspace(0, pi, N);
eig_max1 = zeros(length(alpha),1);
eig_max2 = zeros(length(alpha),1);
hmax_BI2 = zeros(length(alpha),1);
hmax_RK4 = zeros(length(alpha),1);
tol = 1.e-7;
for i = 1:length(alpha)
    A = A_mat(alpha(i));
    eig_max1(i) = cos(alpha(i)) + 1j*sin(alpha(i));
    eig_max2(i) = cos(alpha(i)+pi) + 1j*sin(alpha(i)+pi);
    hmax_BI2(i) = BI_stability(A, theta, tol, [-1.e-1 12], method_BI);
    hmax_RK4(i) = RK_stability(A, tol, [-1.e-1 3], method_RK);
end
eig_max = [ eig_max1
            eig_max2 ];
hmax_BI2 = [     hmax_BI2
             flipud(hmax_BI2)];
hmax_RK4 = [     hmax_RK4
             flipud(hmax_RK4)];

figure()
plot(eig_max.*hmax_BI2, 'LineStyle', '--'), hold on, plot(eig_max.*hmax_RK4, 'LineStyle', '-'), hold on,
scatter(real(eigenvalues(1)), imag(eigenvalues(1)),'g*'), hold on
scatter(real(eigenvalues(2)), imag(eigenvalues(2)),'ro'), hold on, yline(0)
ylim([-20 120])
xlim([-420 10])

title('Stiff system integration')
xlabel('$\Re\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
ylabel('$\Im\{h\lambda\}$', 'Interpreter','latex','FontSize',12')
grid on
legend('RK4','BI2','1^{st} eigenvalue', '2^{nd} eigenvalue','AutoUpdate','off', 'Fontsize', 11)
Lgnd = legend('show');
Lgnd.Position(1) = 0.18;
Lgnd.Position(2) = 0.68;

% zoomPlot to highlight a portion of the major plot
zoomPlot(real(eig_max.*hmax_BI2),imag(eig_max.*hmax_BI2),[-4 4],[-4 4],[0.5 0.3 0.3 0.45], '--', [0, 0.4470, 0.7410], [2 4]);
hold on
plot(real(eig_max.*hmax_RK4),imag(eig_max.*hmax_RK4), 'LineStyle', '-', 'color', [0.8500, 0.3250, 0.0980])
hold on
scatter(real(eigenvalues(1)), imag(eigenvalues(1)),'g*')
hold on, grid on
yline(0)
legend hide

hold off


%% Functions
% EX 1
function BM = bisection_method(f,range,tol)
%
%     bisection_method - zero finding function with Bisection algorithm
%
%     DESCRIPTION:
%       Root finding algorithm based on halving the range around the zero
%
%     PROTOTYPE:
%        BM = bisection_method(f,range,tol)
%          
%     INPUT:
%       f = @(x) ()   Function that contains a zero in the domain taken
%       range [1,2]   Domain with one zero
%       tol           Tolerance to be satisfied 
%     
%     OUTPUT:
%       BM            struct with fields:
%                     - BM.fcn_evals   number of function evaluations   
%                     - BM.steps       number of steps needed for
%                                      convergence
%                     - BM.zero        zero found   
%                     - BM.routine     succesion of zero produced
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      10/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    a = range(1); b = range(2);
    fb = f(b);                             BM.fcn_evals = 1;
    % control:
    if f(a)*f(b) > 0
        error('Error. \nf(a) must have different sign with respect to f(b)!\n')
    end
    BM.steps = 1;
    BM.zero = (a+b)/2;
    fbm_0 = f(BM.zero);                    BM.fcn_evals = BM.fcn_evals+1;    
    BM.routine(BM.steps) = BM.zero;      % useful to display the error flow
    prev_zero = 1000;
    while abs(BM.zero-prev_zero)>tol
        if (fbm_0*fb)<0 
            a = BM.zero;
        else 
            b = BM.zero;
        end
        prev_zero = BM.zero;
        BM.zero = (a+b)/2;
        fbm_0 = f(BM.zero);                BM.fcn_evals = BM.fcn_evals+1;
        BM.steps = BM.steps+1;
        BM.routine(BM.steps) = BM.zero;
    end
end

function SM = secant_method(f,range,tol)
%
%     secant_method - zero finding function with Secant algorithm
%     
%     DESCRIPTION:
%       Root finding algorithm based on a succession of linear to better
%       approximate the function in the zero neighbourhood
%
%     PROTOTYPE:
%        SM = secant_method(f,range,tol)
%     
%     INPUT:
%       f = @(x) ()   Function that contains a zero in the domain taken
%       range [1,2]   Domain with one zero
%       tol           Tolerance to be satisfied 
%     
%     OUTPUT:
%       SM            struct with fields:
%                     - SM.fcn_evals   number of function evaluations   
%                     - SM.steps       number of steps needed for
%                                      convergence
%                     - SM.zero        zero found   
%                     - SM.routine     succesion of zero produced
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      10/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    a = range(1);  b = range(2);
    fa = f(a);                              SM.fcn_evals = 1;
    fb = f(b);                              SM.fcn_evals = SM.fcn_evals+1;
    SM.steps = 1;
    % control:
    if fa*fb > 0
        error('Error. \nf(a) must have different sign with respect to f(b)!\n')
    end
    SM.zero = b-fb/((fb-fa)/(b-a));
    SM.routine(SM.steps) = SM.zero;
    prev_zero = 1000;
    while abs(SM.zero-prev_zero)>tol
        a = b;
        b = SM.zero;
        SM.steps = SM.steps+1;
        fa = fb; 
        fb = f(b);                          SM.fcn_evals = SM.fcn_evals+1;
        prev_zero = SM.zero;
        SM.zero = b-fb/((fb-fa)/(b-a));
        SM.routine(SM.steps) = SM.zero;
    end
end

function RFM = regula_falsi_method(f,range,tol)
%
%     regula_falsi_method - zero finding function with Regula Falsi 
%                           algorithm
%     
%     DESCRIPTION:
%       Root finding algorithm based on a succession of linear to better
%       approximate the function in the zero neighbourhood. The range is
%       selected according to the sign of f(x)
%
%     PROTOTYPE:
%        RFM = regula_falsi_method(f,range,tol)
%     
%     INPUT:
%       f = @(x) ()   Function that contains a zero in the domain taken
%       range [1,2]   Domain with one zero
%       tol           Tolerance to be satisfied 
%     
%     OUTPUT:
%       RFM            struct with fields:
%                     - RFM.fcn_evals   number of function evaluations   
%                     - RFM.steps       number of steps needed for
%                                       convergence
%                     - RFM.zero        zero found   
%                     - RFM.routine     succesion of zero produced
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      10/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    a = range(1);   b = range(2);
    % control:
    fa = f(a);                              RFM.fcn_evals = 1;
    fb = f(b);                              RFM.fcn_evals = RFM.fcn_evals+1;
    if fa*fb > 0
        error('Error. \nf(a) must have different sign with respect to f(b)!\n')
    end
    RFM.steps = 1;
    RFM.zero = b-fb/((fb-fa)/(b-a));
    ffr_0 = f(RFM.zero);                    RFM.fcn_evals = RFM.fcn_evals+1;
    RFM.routine(RFM.steps) = RFM.zero;
    prev_zero = 1000;
    while abs(RFM.zero-prev_zero)>tol
        if (ffr_0*fb)<0
            a = RFM.zero;
            fa = f(a);                      RFM.fcn_evals = RFM.fcn_evals+1; 
        else
            b = RFM.zero;
            fb = f(b);                      RFM.fcn_evals = RFM.fcn_evals+1;
        end
        prev_zero = RFM.zero;
        RFM.zero = b-fb/((fb-fa)/(b-a));
        RFM.steps = RFM.steps+1;
        RFM.routine(RFM.steps) = RFM.zero;
    end
end

function varargout = error_reorder(varargin)
%
%     error_reorder - order arrays to give them to plot function with the
%                     same length
%     
%     DESCRIPTION:
%       This function makes the arrays given as input of the same length as
%       the longest one by repeating the last point as many times as needed
%
%     PROTOTYPE:
%        varargout = error_reorder(varargin)
%     
%     INPUT:
%       varargin      Arrays to be properly sized 
%     
%     OUTPUT:
%       varargout     Arrays with the same length
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      10/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

idx_max = 1;
max_steps = length(varargin{1,1});
for i = 1:nargin
        if length(varargin{1,i}) > max_steps
            max_steps = length(varargin{1,i});
            idx_max = i;
        end
end
nargout = nargin;
varargout = cell([1 nargout]);
for i = 1:nargout
    if i ~= idx_max
        varargout{1,i} = [varargin{1,i} ones(1,max_steps-length(varargin{1,i}))*varargin{1,i}(end)];
    else
        varargout{1,i} = [varargin{1,i}];
    end
end
end

% EX 2
function root = NM_an_FD(f, x, x0, n)
%
%     NM_an_FD - Newton Method with f'(x) analytical or approximated by
%                Finite Differences  
%     
%     DESCRIPTION:
%       This function gives the root of a symbolic function f(x) with a
%       prescribed number of iteration n with analytical, forward
%       differences and central differences jacobian.
%
%     PROTOTYPE:
%        root = NM_an_FD(f, x, x0, n)
%     
%     INPUT:
%       f = (x)       symbolic function whose jacobian must be found
%       x             symbolic variable used in f
%       x0 [n 1]      initial guess for Newoton Method
%       n             number of iteration used
%       
%     OUTPUT:
%       root                 struct with fields:
%                            - an         analytical
%                            - FD1        Forward Finite Differences 
%                            - FD2        Central Finite Differences
% 
%       an, FD1, FD2         struct with fields:
%                            - zero       root found
%                            - routine    succesion of roots found
%                            - steps      number of steps taken
%     
%     CALLED FUNCTIONS:
%       diffjac, dirder
% 
%     LAST UPDATED:
%      11/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

root.an.zero = x0;
root.FD1.zero = x0; % forward finite differences
root.FD2.zero = x0; % cebtral finite differences

root.an.routine = zeros(length(x0),n);
root.FD1.routine = zeros(length(x0),n);
root.FD2.routine = zeros(length(x0),n);

root.an.routine(:,1) = root.an.zero;
root.FD1.routine(:,1) = root.FD1.zero;
root.FD2.routine(:,1) = root.FD2.zero;

root.an.steps = n;
root.FD1.steps = n;
root.FD2.steps = n;

Jf = jacobian(f,x); % analytical computation
f_obj = {matlabFunction(f(1),'vars',{x}); matlabFunction(f(2),'vars',{x})};

for i = 1:n
    x1 = root.an.zero - (double(subs(Jf,x,root.an.zero)))\double(subs(f,x,root.an.zero));
    x1_FD1 = root.FD1.zero - diffjac(f_obj, root.FD1.zero, 'f')\double(subs(f,x,root.FD1.zero));
    x1_FD2 = root.FD2.zero - diffjac(f_obj, root.FD2.zero, 'c')\double(subs(f,x,root.FD2.zero));
    root.an.zero = x1;
    root.FD1.zero = x1_FD1;
    root.FD2.zero = x1_FD2;
    root.an.routine(:,i) = root.an.zero;
    root.FD1.routine(:,i) = root.FD1.zero;
    root.FD2.routine(:,i) = root.FD2.zero;
end
end

function Jf_FD = diffjac(f, x0, type)
%
%     diffjac - Finite Differences Jacobian
%     
%     DESCRIPTION:
%       This function computes the Jacobian with finite differences which 
%       can be central differences or forward differences according to the 
%       precision required
% 
%     PROTOTYPE:
%        Jf_FD = diffjac(f, x0, type)
%     
%     INPUT:
%       f = @(x) ()   objective function for the jacobian computation
%       x0 [n,1]      vector with initial guesses
%       type          char to be selected between 'c' and 'f'
%     
%     OUTPUT:
%       Jf_FD         Jacobian Finite Differences matrix
%     
%     CALLED FUNCTIONS:
%       dirder
% 
%     LAST UPDATED:
%      11/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

Jf_FD = zeros(length(x0), length(x0));
theta_f = max(sqrt(eps), sqrt(eps)*abs(x0));
theta_c = theta_f*4.e4;
for i = 1:length(f)
    for j = 1:length(x0)
        switch type
            case 'c'
                xp1 = x0;    xp1(j) = x0(j) + theta_c(j);
                xm1 = x0;    xm1(j)= x0(j) - theta_c(j);
                f_dot = (f{i}(xp1)-f{i}(xm1))/(2*theta_c(j));
            case 'f'
                xp1 = x0;    xp1(j) = x0(j) + theta_f(j);
                f_dot = (f{i}(xp1)-f{i}(x0))/(theta_f(j));
        end
        Jf_FD(i, j) = f_dot;
    end
end

end

% EX 3
function [yH, tH] = Heun (f, y0, h, tspan)
%
%     Heun - Heun integration method
%
%     DESCRIPTION:
%       This function integrates a function handle f with Heun's method 
%       (RK2)
%
%     PROTOTYPE:
%        [yH, tH] = Heun (f, y0, h, tspan)
%          
%     INPUT:
%       f = @(x) ()   Function to be integrated
%       y0 [n,1]      Initial integration array
%       h             Stepsize
%       tspan [1,2]   Domain with one zero
%     
%     OUTPUT:
%       yH [n,m]      Integrated points, m is the number of stepsize needed
%       tH [1,m]      Time vector associated to the integrated function
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

N = (tspan(2)-tspan(1))/h; % # of intervals
m = length(y0);
tH = linspace(tspan(1), tspan(2), N+1);
yH = zeros(m, N+1);
yH(:,1) = y0;
        for n = 1:N
            K1 = f(yH(:,n), tH(n));
            K2 = f(yH(:,n) + h*K1, tH(n) + h);
            yH(:,n+1) = yH(:,n) + h/2*(K1+K2);
        end
end

function [yRK4, tRK4] = RK4 (f, y0, h, tspan)
%
%     RK4 - Runge Kutta 4 integration method
%
%     DESCRIPTION:
%       This function integrates a function handle f with RK4 (Runge-Kutta
%       of 4th order)
%
%     PROTOTYPE:
%        [yRK4, tRK4] = RK4 (f, y0, h, tspan)
%          
%     INPUT:
%       f = @(x) ()   Function to be integrated
%       y0 [n,1]      Initial integration array
%       h             Stepsize
%       tspan [1,2]   Domain with one zero
%     
%     OUTPUT:
%       yRK4 [n,m]    Integrated points, m is the number of stepsize needed
%       tRK4 [1,m]    Time vector associated to the integrated function
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

N = (tspan(2)-tspan(1))/h;
m = length(y0);
tRK4 = linspace(tspan(1), tspan(2), N+1);
yRK4 = zeros(m,N+1);
yRK4(:,1) = y0;
        for n = 1:N
            K1 = f(yRK4(n), tRK4(n));
            K2 = f(yRK4(n) + h*K1/2, tRK4(n) + h/2);
            K3 = f(yRK4(n) + h*K2/2, tRK4(n) + h/2);
            K4 = f(yRK4(n) + h*K3, tRK4(n) + h);
            yRK4(:,n+1) = yRK4(:,n) + h*(K1+2*K2+2*K3+K4)/6;
        end 
end

% EX 4
function A = A_mat (alpha)
%
%     A_mat - matrix A creation
%
%     DESCRIPTION:
%       Creates a matrix A whose eigenvalues form the unitary circle in the
%       complex plane
%
%     PROTOTYPE:
%        A = A_mat (alpha)
%          
%     INPUT:
%       alpha         Variational parameter of A matrix
%     
%     OUTPUT:
%       A [2,2]       Matrix 
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

A = [ 0      1
     -1 2*cos(alpha) ];
end

function F = F_RK124(A, h, method)
%
%     F_RK124 - Forward operator of RK1, RK2 or RK4
%
%     DESCRIPTION:
%       This function gives as output the forward operator F useful for RKs
%       from step x_k to step x_k+1
%
%     PROTOTYPE:
%        F = F_RK124(A, h, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       h                  Stepsize
%       method             char selected between:
%                          - 'RK1'
%                          - 'RK2'
%                          - 'RK4'
%     
%     OUTPUT:
%       F [n,n]            Operator of the selected RK 
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

switch method
    case 'RK1'
        F = eye(length(A)) + h*A;
    case 'RK2'
        F = eye(length(A)) + h*A + (h*A)^2/2;
    case 'RK4'
        F = eye(length(A)) + h*A + (h*A)^2/2 + (h*A)^3/6 + (h*A)^4/24;
    otherwise
        error('Error. \nInput, method must be selected between "RK1", "RK2" or "RK4"')
end
end

function hmax = RK_stability(A, tol, range, method)
%
%     RK_stability - Stability region of RK1,RK2 or RK4
%
%     DESCRIPTION:
%       This function gives as output the maximum stepsize that can be used
%       as input of RKs methods to pass from step x_k to step x_k+1
%       according to the stability properties
%
%     PROTOTYPE:
%        hmax = RK_stability(A, tol, range, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       tol                Tolerance required wrt the exact hmax
%       range [1,2]        Interval where the stability region is located
%       method             char selected between:
%                          - 'RK1'
%                          - 'RK2'
%                          - 'RK4'
%     
%     OUTPUT:
%       hmax               Maximum stepsize for stable region
%     
%     CALLED FUNCTIONS:
%      F_RK124
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

h = (range(1) + range(2))/2;
F = F_RK124(A, h, method);
err = max(abs(eig(F))) - 1;
while abs(err) > tol
    h = (range(1) + range(2))/2;
    F = F_RK124(A, h, method);
    F_eig = eig(F);
    absolute_eig = abs(F_eig);
    eigenvalue = max(absolute_eig);
    err = eigenvalue - 1;
    if err > 0 
        range(2) = h;
    else
        range(1) = h;
    end
end
hmax = h;
end

function plotting(title_obj, legend_obj, x, y,varargin)
%
%     plotting - 4 plots together
%
%     DESCRIPTION:
%       This function plots 4 graphs in the same figure with different
%       linestyles and colors
%
%     PROTOTYPE:
%        plotting(title_obj, legend_obj, x, y,varargin)
%          
%     INPUT:
%       title_obj          char which appears in the title
%       legend_obj {1,n}   Variables that appear in the legend
%       x {1,n}            Variable of x in the plot, it can be set to
%                          'none' if it doesn't matter of x in the plot 
%                          function 
%       y {1,n}            Variable to be displayed in the plot
%       varargin {1,1}     char with the title wanted
%       varargin {1,2}     char with the legend wanted
%     
%     OUTPUT:
%       figure
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

figure()
if nargin > 5
    wanted_title = varargin{1};
    wanted_legend = varargin{2};
elseif nargin == 5
    wanted_legend = varargin{1};
    wanted_title = '';
else
    wanted_legend = 'data';
    wanted_title = '';
end
lin = {'-','--',':','-.'};
for i = 1:length(y)
    legends = [wanted_legend,'_', num2str(i), ' = ', num2str(legend_obj{i})];
    if strcmp(x{i},'none')
        plot(y{i}, 'linestyle', lin{i}, 'Displayname', legends), hold on
    else
        plot(x{i}, y{i}, 'linestyle', lin{i}, 'Displayname', legends), hold on
    end
end
title(sprintf('%s %s', wanted_title, title_obj),'FontSize',12)
grid on
hold off
end

function zoomPlot(x,y,xbounds,ybounds,pos,varargin)
% zoomPlot     add inlaid plot to current figure
%       [p,z]   = zoomPlot(x,y,xbounds,pos,vertex) 
%       zoomPlot(x,y,xbounds,ybounds,pos,varargin) where: 
%       x,y     = vectors being plotted
%       xbounds = [x1 x2] specifies the zoom indices, where x1 is the 
%               first x value displayed in the zoom plot and x2 is the last.
%       pos     = [left, bottom, width, height] specifies the location and
%               size of the side of the zoom box, relative to the lower-left
%               corner of the Figure window, in normalized units where (0,0)
%               is the lower-left corner and (1.0,1.0) is the upper-right.
% {opt} vertex  = toggles connecting lines corresponding to vertices, where 1 
%               corresponds to the top left vertex and continuing clockwise, 
%               4 corresponds to the bottom right vertex. All 4 vertices can 
%               be included.
% {opt} p       = axes handle for larger plot
% {opt} z       = axes handle for zoom plot
% 
% Note: place title, labels, and legend BEFORE placing zoom plot,
%     otherwise zoomPlot returns the handle of the original axes (p).
%     Change title using p.Title.String = 'insert title here'
% 
% Kelsey Bower (kelsey.bower@case.edu) 10/20/15
%
% Please retain the following:
% 
% Original Author: 
% Kelsey Bower, kelsey.bower@case.edu
%
% Last modification:
% Niccolò Bardazzi - 14/10/2021
%          
%     INPUT:
%       ybounds            Bounds for y variable
%       varargin {1,1}     char with the linstyle selected
%       varargin {1,2}     char with the color selected
%       varargin {1,3}     Angles of the zoom picture connected
%     
%     OUTPUT:
%       -
% 

if nargin > 8
    printf('Too many arguments. zoomPlot(x,y,xbounds,ybounds,pos,style,vertex)\n')
elseif nargin < 8
    vertex = [1 4];
elseif nargin == 8
    vertex = varargin{3};
end
% Get current axis position and limits
p = gca;

% Calculate x,y points of zoomPlot
x1 = (pos(1)-p.Position(1))/p.Position(3)*diff(p.XLim)+p.XLim(1);
x2 = ((pos(1)+pos(3)-p.Position(1))/p.Position(3))*diff(p.XLim)+(p.XLim(1));
y1 = (pos(2)-p.Position(2))/p.Position(4)*diff(p.YLim)+p.YLim(1);
y2 = ((pos(2)+pos(4)-p.Position(2))/p.Position(4))*diff(p.YLim)+(p.YLim(1));

% Plot lines connecting zoomPlot to original plot points
rectangle('Position',[xbounds(1) ybounds(1) diff(xbounds) diff(ybounds)], 'LineWidth', 0.01);
hold on
if any(vertex==1)
    plot([xbounds(1) x1], [ybounds(2) y2], 'k'); % Line to vertex 1
end
if any(vertex==2)
    plot([xbounds(2) x2], [ybounds(2) y2], 'k'); % Line to vertex 2
end
if any(vertex==3)
    plot([xbounds(2) x2], [ybounds(1) y1], 'k'); % Line to vertex 4
end
if any(vertex==4)
    plot([xbounds(1) x1], [ybounds(1) y1], 'k'); % Line to vertex 3
end

% Plot zoomPlot and change axis
axes('position',pos)
box on 
plot(x, y,"LineStyle",varargin{1},'Color',varargin{2})
axis([xbounds(1) xbounds(2) ybounds(1) ybounds(2)]);
legend hide
end

% EX 5
function xRK_end = RK_end (A, x0, h, tspan, method)
%
%     RK_end - End point integration value for RK1, RK2 or RK4
%
%     DESCRIPTION:
%       This function gives as output the value of the last point of a
%       numerical integration for a LTI system
%
%     PROTOTYPE:
%        xRK_end = RK_end (A, x0, h, tspan, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       x0 [1,n]           Initial integration point
%       h                  Stepsize
%       tspan [1,2]        Interval where the stability is located
%       method             char selected between:
%                          - 'RK1'
%                          - 'RK2'
%                          - 'RK4'
%     
%     OUTPUT:
%       xRK_end            Value at the last point of numerical integration
%     
%     CALLED FUNCTIONS:
%      F_RK124
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

N = floor((tspan(2)-tspan(1))/h);
last_h = rem((tspan(2)-tspan(1)),h);
xRK_end = x0;
F = F_RK124(A, h, method);
xRK_end = F^N*xRK_end;
if last_h ~= 0
    method = 'RK4';
    f = F_RK124(A, last_h, method);
    xRK_end = f*xRK_end;
end
end

function hmax = accuracy_domain(A, x0, tspan, maxerr, tol , method)
%
%     accuracy_domain - Maximum stepsize for the accuracy domain
%
%     DESCRIPTION:
%       This function gives as output the maximum stepsize for RKs methods
%       of a LTI system for an accurate numerical integration according to
%       the prescribed tolerance
%
%     PROTOTYPE:
%        hmax = accuracy_domain(A, x0, tspan, maxerr, tol , method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       x0 [1,n]           Initial integration point
%       tspan [1,2]        Interval where the stability is located
%       maxerr             Maximum admissible error on hmax (usually 1% of
%                          tol)
%       tol                Tolerance prescribed by the method
%       method             char selected between:
%                          - 'RK1'
%                          - 'RK2'
%                          - 'RK4'
%     
%     OUTPUT:
%       hmax               Maximum stepsize for the accuracy domain
%     
%     CALLED FUNCTIONS:
%      RK_end
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

range(2) = tspan(2);
range(1) = tspan(1);
h = (tspan(1) + tspan(2))/2;
x_anend = expm(A)*x0;
xRK_end = RK_end (A, x0, h, range, method);
eps_global = norm(x_anend - xRK_end, Inf);
err = eps_global-tol;
if err > 0 
    tspan(2) = h;
else
    tspan(1) = h;
end
while abs(err) > maxerr
    h = (tspan(1) + tspan(2))/2;
    xRK_end = RK_end (A, x0, h, range, method);
    eps_global = norm(x_anend - xRK_end, Inf);
    err = eps_global-tol;
    if err > 0
        tspan(2) = h;
    else
        tspan(1) = h;
    end
end
hmax = h;
end

% EX 6
function hmax = BI_stability(A, theta, tol, range, method)
%
%     BI_stability - Stability region of RK1,RK2 or RK4
%
%     DESCRIPTION:
%       This function gives as output the maximum stepsize that can be used
%       as input of BIs methods to pass from step x_k to step x_k+1
%       according to the stability properties
%
%     PROTOTYPE:
%        hmax = BI_stability(A, theta, tol, range, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       theta              Percentage of explicitely in BI methods
%       tol                Tolerance required wrt the exact hmax
%       range [1,2]        Interval where the stability is located
%       method             char ['Bi', num2string(#ofBI)]
%     
%     OUTPUT:
%       hmax               Maximum stepsize for stable region
%     
%     CALLED FUNCTIONS:
%      BI
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

h = (range(1) + range(2))/2;
F = BI(A, h, theta, method);
err = max(abs(eig(F))) - 1;
while abs(err) > tol
    h = (range(1) + range(2))/2;
    F = BI(A, h, theta, method);
    F_eig = eig(F);
    absolute_eig = abs(F_eig);
    eigenvalue = max(absolute_eig);
    err = eigenvalue - 1;
    if theta < 0.5
        err = -err;
    end
    if err > 0 
        range(2) = h;
    else
        range(1) = h;
    end
end
hmax = h;
end

function F = BI(A, h, theta, method)
%
%     BI - Back Interpolation operator
%
%     DESCRIPTION:
%       This function gives as output the operator used to pass from step
%       x_k to step x_k+1 according to the stability properties
%
%     PROTOTYPE:
%        F = BI(A, h, theta, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       h                  Stepsize
%       theta              Percentage of explicitely in BI methods
%       method             char ['BI', num2string(#ofBI)]
%     
%     OUTPUT:
%       F                  Operator from step x_k to step x_k+1
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

F_impl = eye(length(A));
F_espl = eye(length(A));
for i = 1:method
    F_impl = F_impl + (-1)^i*(A*(1-theta)*h)^i/factorial(i);
    F_espl = F_espl + (A*theta*h)^i/factorial(i);
end
F = F_impl\F_espl;
end

% EX 7
function [yBI, tBI] = BI_LTI(A, y0, h, theta, tspan, method)
%
%     BI_LTI - Back Interpolation integrator
%
%     DESCRIPTION:
%       This function gives as output the numerical solution to the
%       integration problem using a Back Interpolation method and the
%       respective time scale of an LTI system
%
%     PROTOTYPE:
%        [yBI2, tBI2] = BI_LTI(A, y0, h, theta, tspan, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       y0 [n,1]           Initial value array
%       h                  Stepsize
%       theta              Percentage of explicitely in BI methods
%       tspan [1,2]        Interval of integration
%       method             char ['BI', num2string(#ofBI)]
%     
%     OUTPUT:
%       yBI [n,m]          Solution of the integration, m: # of intervals+1 
%       tBI [1,m]          Time scale related to the solution
%     
%     CALLED FUNCTIONS:
%      BI
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

N = (tspan(2)-tspan(1))/h;
m = length(y0);
tBI = linspace(tspan(1), tspan(2), N+1);
yBI = zeros(m,N+1);
yBI(:,1) = y0;
F = BI(A, h, theta, method);
for n = 1:N
    yBI(:,n+1) = F*yBI(:,n);
end
end

function [yRK, tRK] = RK_LTI (A, y0, h, tspan, method)
%
%     RK_LTI - Runge Kutta integrator
%
%     DESCRIPTION:
%       This function gives as output the numerical solution to the
%       integration problem using a Runge Kutta method and the respective 
%       time scale of an LTI system
%
%     PROTOTYPE:
%        [yRK, tRK] = RK_LTI(A, y0, h, tspan, method)
%          
%     INPUT:
%       A [n,n]            Matrix which maps the dynamics of the system
%       y0 [n,1]           Initial value array
%       h                  Stepsize
%       tspan [1,2]        Interval of integration
%       method             char ['RK', num2string(#ofRK)]
%     
%     OUTPUT:
%       yRK [n,m]          Solution of the integration, m: # of intervals+1 
%       tRK [1,m]          Time scale related to the solution
%     
%     CALLED FUNCTIONS:
%      F_RK124
% 
%     LAST UPDATED:
%      20/10/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

N = (tspan(2)-tspan(1))/h;
m = length(y0);
tRK = linspace(tspan(1), tspan(2), N+1);
yRK = zeros(m,N+1);
yRK(:,1) = y0;
F = F_RK124(A, h, method);
for n = 1:N
    yRK(:,n+1) = F*yRK(:,n);
end    
end