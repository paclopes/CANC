% Careful Active Noise Control (CANC) algorithm simulation
%
% For the paper:
% Careful Active Noise Control Analysis
% Paulo A. C. Lopes
% to be published

global min_path_gain;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Parameters

N = 32;
Nw = N;         % control filter length
Ns = N;         % secondary path filter length
Np = N;         % primary path filter length
Nc = N;         % control filter update
Mx = 8;
M = Mx*N;          % memory
deltau = 1e-2;     % small number
deltay = deltau;   % small number
deltaw = deltau;   % small number
P = M;             % qvM calculation frames size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Simulation Setup Parameters

fs = 2000;                     % sampling frequency (Hz)
L = 2*fs;                      % simulation length (samples)
Q = 1;                         % backgound noise coloring filter size
min_path_gain = 0.1;           % mininum gain of the primary and secondary path at any frequency
qv0 = 0.3;                     % backgroung noise power
qr0 = qv0;                     % off-line auxiliary noise power
qr1 = 0;                       % on-line auxilixar noise power
kp = 2;                        % primary noise frequency is f = k/N*fs
qz = 0;                        % reference signal noise power

test_simulation = -1;          % use -1 to do large simuations
change_params = false;         % changes paths, frequencies, background noise and attenuation at every simulation
N_simulations = 100;            % number of simulations in a large simulation

change_at = L/2;               % sudden secondary path change time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Simulation Initializations

if test_simulation<0
    set_of_simulations = 1:N_simulations;
else
    N_simulations = 1;
    set_of_simulations = test_simulation;
end

anc_on = 3/2*M;                 % start of CANC algorithm

% loggers
primary_noise_power=zeros(N_simulations,1);      % primary noise power
residual_noise_power=zeros(N_simulations,1);     % residual noise power
min_residual_noise_power=zeros(N_simulations,1); % minimum residual noise power achivable
excess_noise_power=zeros(N_simulations,1);       % residual noise power

log_qv = zeros(L,N_simulations);        % estimated backgorund noise power
log_y = zeros(L,N_simulations);         % anti-noise signal
log_e = zeros(L,N_simulations);         % residual noise signal
log_d = zeros(L,N_simulations);         % primary noise signal
log_ex = zeros(L,N_simulations);        % excess residual noise signal
log_w = zeros(L,N_simulations);         % control filter at PN frequency, w_F(n+1)
log_p = zeros(L,N_simulations);         % primary path at PN frequency
log_s = zeros(L,N_simulations);         % secondary path at PN frequency
log_n = zeros(L,N_simulations);         % control filter update time

tic;
%warning('off', 'MATLAB:nearlySingularMatrix');

for simulation=1:N_simulations
    fprintf(1, 'Simulation %d\n', set_of_simulations(simulation));

    rng(637373 + set_of_simulations(simulation));

    % Individual Simulation Initializations
    if change_params || simulation==1
        % frequencies
        f = kp*fs/Nw;  % period is Nw/2

        % Primary Path
        p = calc_random_path(Np);

        % Secondary Path
        s1 = calc_random_path(Ns);
        s2 = - s1;

        % Feedback Path
        fb = calc_random_path(Ns)/10;

        % noise
        n=0:L+Np-2;
        theta=2*pi*f/fs*n+2*pi*rand;

        p_frequency_response = freqz(p, 1, [0, 2*pi*f/fs]);
        Au = sqrt(2)/abs(p_frequency_response(2));
        u = Au*cos(theta) + sqrt(qz)*randn(1, L+Np-1);
    end
    s = s1;

    v = conv(randn(1,L+Q-1), randn(1,Q), 'valid');
    v = sqrt(qv0)*v/std(v);

    % Variables
    w = zeros(Nw,1); w0 = w;
    wF = fft(w); wF = wF(3);
    pF = nan;
    sF = nan;
    qvM = nan;
    qr = qr0;
    n1 = nan;

    % Buffers
    yv = zeros(max(Ns,2*M-1),1);
    uv = zeros(M+2*max(Np,Ns)-2,1);
    ev = zeros(M,1);

    for n=1:L
        if n == change_at
            s = s2;
        end

        % physical simulation
        us = u(n) + fb'*yv(1:Ns);

        % w = p/s;
        % s = s/(1 - fb w)
        % w = p/s (1 - fb w)
            
        % proposed algorithm
        uv = [us; uv(1:end-1)];
        y = w'*uv(1:Nw) + sqrt(qr)*randn;

        % physical simulation
        yv = [y; yv(1:end-1)];
        d = p'*uv(1:Np) + v(n);
        e = d + s'*yv(1:Ns);
        ex = p'*uv(1:Np) + s'*yv(1:Ns); % to log

        % proposed algorithm
        ev = [e; ev(1:end-1)];

        if n>anc_on && mod(n,Nc)==0
            qr = qr1;
            U = uv((0:M-1)'+(1:Np));
            Y = yv((0:M-1)'+(1:Ns));
            H = [U, Y];
            % ex = H x
            Rx = H'*H + diag([deltau*ones(Np,1); deltay*ones(Ns,1)]);
            x = Rx\H'*ev;
            ph = x(1:Np); sh=x(Np+1:end);
            qvM = max(mean(reshape(abs(ev-H*x).^2, P, [])));
            Sx = qvM*inv(Rx);
            Sss = Sx(Np+1:end, Np+1:end);
            Spp = Sx(1:Np,1:Np);
            Sps = Sx(1:Np, Np+1:end);
            Ssp = Sx(Np+1:end,1:Np); 
 
            Qss = sh*sh'+Sss; Qsp = sh*ph'+Ssp; % careful
            % Qss = sh*sh'; Qsp = sh*ph'; % not careful
            Rc = zeros(Nw); Pv = 0; Qssx = zeros(Ns+Ns-1); Qspx = zeros(Ns+Ns-1, Np+Ns-1);
            for i=1:Ns
                Qssx(i:Ns+i-1, i:Ns+i-1) = Qssx(i:Ns+i-1, i:Ns+i-1) + Qss;
                Qspx(i:Ns+i-1, i:Np+i-1) = Qspx(i:Ns+i-1, i:Np+i-1) + Qsp;
            end
            for i=1:Ns:M
                Ui = uv((i-1:Nw-1+i-1)'+(1:Ns+Ns-1));
                Rc = Rc + Ui*Qssx*Ui'; % U(1:N,:)*sh = u'
                Pv = Pv + Ui*(Qspx*uv(i:Np+Ns-1+i-1));
            end
            deltaw1 = deltaw*trace(Rc)/Nw;
            Rc = Rc + deltaw1*eye(Nw);
        
            Rc = Rc/M; Pv = Pv/M;
            w0 = w;
            w = - Rc\Pv;
            % w = max(min(w, 1), -1);
        end

        view_plots = false;
        if view_plots
            cool_fig(10);
            [h, ws]=freqz(w);
            plot(ws,abs(h));
            grid on;
    
            cool_fig(11);
            plot(w);
            grid on;
            drawnow;
        end

        if mod(n, Nc) == 0
            n1 = n;
            sF = fft(s, Nw); sF = sF(3);
            pF = fft(p, Nw); pF = pF(3);
            wF = fft(w0); wF = wF(3);
        end

        % logs
        log_d(n, simulation) = d;
        log_e(n, simulation) = e;
        log_ex(n, simulation) = ex;
        log_qv(n, simulation) = qvM;
        log_y(n, simulation) = y;
        log_w(n, simulation) = wF;        
        log_p(n, simulation) = pF;        
        log_s(n, simulation) = sF;        
        log_n(n, simulation) = n1;
    end

    noise_power_i = mean(log_d(L/2+1:end, simulation).^2);
    residual_noise_power_i = mean(log_e(end-L/10:end, simulation).^2);
    min_residual_noise_power_i = qv0;
    excess_noise_power_i = mean(log_ex(end-L/10:end, simulation).^2);
    fprintf(1,'noise_power_dB: %f\n', 10*log10(noise_power_i));
    fprintf(1,'residual_noise_power_dB: %f\n', 10*log10(residual_noise_power_i));
    fprintf(1,'min_residual_noise_power_dB: %f\n', 10*log10(min_residual_noise_power_i));
    fprintf(1,'excess_noise_power_dB: %f\n', 10*log10(excess_noise_power_i));

    if test_simulation<0
        primary_noise_power(simulation)=noise_power_i;
        residual_noise_power(simulation)=residual_noise_power_i;
        excess_noise_power(simulation)=excess_noise_power_i;
        min_residual_noise_power(simulation)=min_residual_noise_power_i;
    end

end
%warning('on', 'MATLAB:nearlySingularMatrix');
toc;

if test_simulation<0
    mx = mean(excess_noise_power);
    sx = std(excess_noise_power);
    fprintf(1, 'mean excess noise power %f dB\n', 10*log10(mean(excess_noise_power(excess_noise_power<mx+10*sx))));
    fprintf(1, 'std excess noise power %f dB\n', 10*log10(std(excess_noise_power(excess_noise_power<mx+10*sx))));

    fprintf(1, '\nStats\n');
    fprintf(1, 'max excess noise: %f dB\n', 10*log10(max(excess_noise_power)));

    cool_fig(10);
    histogram(10*log10(excess_noise_power), -60:-0);
    xlabel('Excess Noise Power (dB)');
    ylabel(['Frequency (out of ', num2str(N_simulations), ' Simulations)']);
    grid on;
    set(gcf, 'name', 'histogram');
end

font_size = 12;
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextFontSize', font_size)
set(groot, 'DefaultAxesFontSize', font_size)
set(groot, 'DefaultLegendFontSize', font_size)

n = 1:L;
cool_fig(1);
plot_xy_p4((1:L)/fs, 10*log10(smooth(log_ex.^2,Nw)), 0.8*fs);
grid on;
xlabel('time (s)');
ylabel('Excess Residual Noise Power (dB)');
set(gcf, 'name', 'ex');
set(gcf, 'NumberTitle', 'off');
set(gca, 'Ylim', [-60   10]);
% save2pdf('CANC_ex', 1, 600);

cool_fig(2);
plot_xy_p3((1:L)/fs, 10*log10(smooth(log_e.^2,Nw)));
grid on;
xlabel('time (s)');
ylabel('Residual Noise Power (dB)');
set(gcf, 'name', 'e');
set(gcf, 'NumberTitle', 'off');

cool_fig(3);
spectrogram(log_e(:,1),256,[],[],fs,'yaxis');
caxis([-60 -20]);
set(gcf, 'name', 'spectrogram');
set(gcf, 'NumberTitle', 'off');

cool_fig(4);
plot_xy_p3((1:L)/fs, 10*log10(log_qv));
%set(gca, 'YLim', [-10, 15]);
grid on;
xlabel('time (s)');
ylabel('Estimated Backgroud Noise Power (dB)');
grid on;
set(gcf, 'name', 'qvM');
set(gcf, 'NumberTitle', 'off');

cool_fig(7);
plot_xy_p3((1:L)/fs, real(log_w));
grid on;
xlabel('time (s)');
ylabel('$\Re(w)$');
grid on;
set(gcf, 'name', 'w');
set(gcf, 'NumberTitle', 'off');

cool_fig(8);
sw = var1(log_w(1:Nc:L,:), Mx);
sw(1:anc_on/Nc,:) = nan;
plot_xy_p4((1:Nc:L)/fs, 10*log10(sw), 0.8*fs/Nc);
grid on;
xlabel('time (s)');
ylabel('$\sigma^2_{\mathrm{wF}}$ (dB)');
grid on;
set(gcf, 'name', 'sw');
set(gcf, 'NumberTitle', 'off');
% save2pdf('CANC_sw', 8, 600);

cool_fig(9);
eta1 = 1.23869;
eta2 = 0.864926;
%eta3 = 0.5485120;
eta3 = 1;
k = 0.33;
qu = Au^2/4; % u power per frequency bin
vsB = qv0./(M*qu*mean1(sw, Mx));
sb = mean1(log_s(1:Nc:L,:), Mx);
pb = mean1(log_p(1:Nc:L,:), Mx);
wb = mean1(log_w(1:Nc:L,:), Mx);
sab = abs(sb);
varBeta = k*vsB./(sab.^4 + eta3*8*vsB.*sab.^2 + 6*eta2*vsB.^2);
pd = sb.*wb + pb;
qp = qv0/(M*qu);
sw1T = max(varBeta.*(qp + mean1(abs(pd).^2, Mx)), 0);
gammaT = mean(sw1T./mean1(sw,Mx),2);
sw1 = [sw(2:end,:); nan*ones(1, N_simulations)];
gamma = mean(sw1./mean1(sw,Mx),2);
plot((1:Nc:L)/fs, 10*log10(gamma),'x');
hold on;
plot((1:Nc:L)/fs, 10*log10(gammaT),'o');
hold off;
grid on;
xlabel('time (s)');
ylabel('$\gamma(n)$ (dB)');
legend('simulation', 'theory');
grid on;
set(gcf, 'name', 'gamma');
set(gcf, 'NumberTitle', 'off');
