% Careful Active Noise Control proposed change 3 (CANC3) algorithm simulation.
%
% For the paper:
% Careful Active Noise Control
% Paulo A. C. Lopes and José A. B. Gerald
% to be published

global min_path_gain;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Parameters

Nw = 32;         % filter length
Ns = Nw;         % filter length
Np = Ns;         % filter length
Nc = Nw;         % control filter update
M = 8*Nw;        % memory
deltau = 1e-9;   % small number
deltay = 1e-9;   % small number
deltaw = 1e-9;   % small number
P = 16;          % qvM calculation frame size
alpha0 = 0.1;    % the reduction in caution
tau  = 100;      % relative changes time constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Simulation Setup Parameters

fs = 2000;                     % sampling frequency (Hz)
N_sins = 3;                    % number of sinusoids
L = 4*fs;                      % simulation length (samples)
change = round(L/2);           % sudden change time (samples)
Q = 8;                         % backgound noise coloring filter size
min_path_gain = 0.1;           % mininum gain of the primary and secondary path at any frequency
background_noise_interval_dB = [0 0];  % interval for the background noise level in dB
relative_primary_noise_dB = [0 0];     % interval for the primary noise power in dB relative to the background noise (equal to the maximum noise reduction)
auxiliare_noise_power_dB = -inf;

test_simulation = -1;          % use -1 to do large simuations
change_params = false;         % changes paths, frequencies, background noise and attenuation at every simulation
N_simulations = 100;           % number of simulations in a large simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Simulation Initializations

if test_simulation<0
    set_of_simulations = 1:N_simulations;
else
    N_simulations = 1;
    set_of_simulations = test_simulation;
end

auxiliare_noise_rms = 10^(auxiliare_noise_power_dB/20);

anc_on = M;

% loggers
primary_noise_power=zeros(N_simulations,1);      % primary noise power
residual_noise_power=zeros(N_simulations,1);     % residual noise power
min_residual_noise_power=zeros(N_simulations,1); % minimum residual noise power achivable
excess_noise_power=zeros(N_simulations,1);       % residual noise power

log_ss = zeros(L,N_simulations);                 % trace of Sss
log_qvM = zeros(L,N_simulations);                    
log_y = zeros(L,N_simulations);
log_e = zeros(L,N_simulations);
log_d = zeros(L,N_simulations);
log_eb = zeros(L,N_simulations);
log_w = zeros(L,N_simulations);

tic;
for simulation_i=1:N_simulations
    fprintf(1, 'Simulation %d\n', set_of_simulations(simulation_i));

    rng(637375+set_of_simulations(simulation_i));

    % Individual Simulation Initializations

    if change_params || simulation_i==1
        % frequencies
        fp = (300-60)*rand+60;
        f = (1:N_sins)'*fp;

        % Primary Path
        p = calc_random_path(Np);

        % Secondary Path
        s1 = calc_random_path(Ns); s = s1;
        s2 = -s1;

        % noise
        background_noise_rms = 10^((diff(background_noise_interval_dB)*rand+min(background_noise_interval_dB))/20);
        primary_noise_rms = 10^((diff(relative_primary_noise_dB)*rand+min(relative_primary_noise_dB))/20)*background_noise_rms;
        sin_amplitudes = [1,1.2,0.5];
        sin_amplitudes = sin_amplitudes*sqrt(2/sum(sin_amplitudes.^2))*primary_noise_rms;

        p_frequency_response = freqz(p, 1, 2*pi*f/fs);
        n=0:L+Np-2;
        theta=2*pi*f/fs*n+2*pi*rand(N_sins,1);
        u = sin_amplitudes./abs(p_frequency_response')*cos(theta-angle(p_frequency_response));
    end
    s = s1;
    
    v = conv(randn(1,L+Q-1), randn(1,Q), 'valid');
    v = background_noise_rms*v/std(v);

    % Variables
    w = zeros(Nw,1); w(1) = 0.1;
    sh = zeros(Ns,1);
    ph = zeros(Np,1);
    Sss = nan*eye(Ns);
    qvM = nan;
    alpha = alpha0;

    % Buffers
    yv = zeros(max(Ns,2*M-1),1);
    uv = zeros(M+2*max(Np,Ns)-2,1);
    ev = zeros(M,1);

    for n=1:L
        % proposed algorithm
        uv = [u(n); uv(1:end-1)];
        y = w'*uv(1:Nw) + auxiliare_noise_rms*randn;

        % physical simulation
        yv = [y; yv(1:end-1)];
        d = p'*uv(1:Np) + v(n);
        e = d + s'*yv(1:Ns);

        if n==change
            s = s2;
        end

        % proposed algorithm
        ev = [e; ev(1:end-1)];

        if n>anc_on && mod(n,Nc)==0
            U = uv((0:M-1)'+(1:Np));
            Y = yv((0:M-1)'+(1:Ns));
            H = [U, Y];
            % eb = H x
            Rx = H'*H + diag([deltau*ones(Np,1); deltay*ones(Ns,1)]);
            x = Rx\H'*ev;
            ph = x(1:Np); sh=x(Np+1:end);
            qvM = max(mean(reshape(abs(ev-H*x).^2, P, [])));
            Sx = alpha*qvM*inv(Rx);
            Spp = Sx(1:Np,1:Np); Sps = Sx(1:Np, Np+1:end);
            Ssp = Sx(Np+1:end,1:Np); Sss = Sx(Np+1:end, Np+1:end);
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

            dw =  w - w0;
            if (norm(dw)^2/norm(w)^2 > 1/tau)
               alpha = 1;
            else
               alpha = alpha0;
            end
        end

        % logs
        log_d(n, simulation_i) = d;
        log_e(n, simulation_i) = e;
        log_eb(n, simulation_i) = p'*uv(1:Np) + s'*yv(1:Ns);
        log_qvM(n, simulation_i) = qvM;
        log_ss(n, simulation_i) = trace(Sss);
        log_y(n, simulation_i) = y;
        log_w(n, simulation_i) = sum(w.*sin(2*pi*fp/fs*(1:Nw)'));
    end

    noise_power_i=mean(log_d(L/2+1:end, simulation_i).^2);
    residual_noise_power_i=mean(log_e(end-L/10:end, simulation_i).^2);
    min_residual_noise_power_i=background_noise_rms^2;
    excess_noise_power_i = mean(log_eb(end-L/10:end, simulation_i).^2);
    fprintf(1,'noise_power_dB: %f\n', 10*log10(noise_power_i));
    fprintf(1,'residual_noise_power_dB: %f\n', 10*log10(residual_noise_power_i));
    fprintf(1,'min_residual_noise_power_dB: %f\n', 10*log10(min_residual_noise_power_i));
    fprintf(1,'excess_noise_power_dB: %f\n', 10*log10(excess_noise_power_i));

    if test_simulation<0
        primary_noise_power(simulation_i)=noise_power_i;
        residual_noise_power(simulation_i)=residual_noise_power_i;
        excess_noise_power(simulation_i)=excess_noise_power_i;
        min_residual_noise_power(simulation_i)=min_residual_noise_power_i;
    end

end
toc;

if test_simulation<0
    cool_fig(1);
    histogram(10*log10(excess_noise_power), -60:0);
    xlabel('Excess Noise Power (dB)');
    ylabel(['Frequency (out of ', num2str(N_simulations), ' Simulations)']);
    set(gcf, 'name', 'histogram');
    set(gcf, 'NumberTitle', 'off');
    grid on;
    m = mean(excess_noise_power);
    s = std(excess_noise_power);
    fprintf(1, 'mean excess noise power %f dB\n', 10*log10(mean(excess_noise_power(excess_noise_power<m+10*s))));
    fprintf(1, 'std excess noise power %f dB\n', 10*log10(std(excess_noise_power(excess_noise_power<m+10*s))));

    fprintf(1, '\nStats\n');
    fprintf(1, 'max excess noise: %f dB\n', max(excess_noise_power));
end

cool_fig(3);
plot_xy_p3((1:L)/fs, 10*log10(smooth(log_e.^2,Nw)));
grid on;
xlabel('time (s)');
ylabel('Residual Noise Power (dB)');
set(gcf, 'name', 'e');
set(gcf, 'NumberTitle', 'off');

cool_fig(4);
spectrogram(log_e(:,1),256,[],[],fs,'yaxis');
caxis([-60 -20]);
set(gcf, 'name', 'spec');
set(gcf, 'NumberTitle', 'off');

cool_fig(5);
plot_xy_p3((1:L)/fs, 10*log10(log_qvM));
%set(gca, 'YLim', [-10, 15]);
grid on;
xlabel('time (s)');
ylabel('Estimated Backgroud Noise Power (dB)');
grid on;
set(gcf, 'name', 'qvM');
set(gcf, 'NumberTitle', 'off');

cool_fig(6);
plot_xy_p3((1:L)/fs, 10*log10(smooth(log_y.^2, Nw)));
grid on;
xlabel('time (s)');
ylabel('Antinoise Power (dB)');
grid on;
set(gcf, 'name', 'y');
set(gcf, 'NumberTitle', 'off');

cool_fig(7);
plot_xy_p3((1:L)/fs, 10*log10(log_ss));
grid on;
xlabel('time (s)');
ylabel('Trace $\Sigma_\mathrm{ss}$ (dB)');
grid on;
set(gcf, 'name', 'Sss');
set(gcf, 'NumberTitle', 'off');

cool_fig(8);
plot((1:L)/fs, log_w(:,1));
grid on;
xlabel('time (s)');
ylabel('w');
grid on;
set(gcf, 'name', 'w');
set(gcf, 'NumberTitle', 'off');

cool_fig(9);
plot_xy_p3((1:L)/fs, 10*log10(smooth(log_eb.^2,Nw)));
grid on;
xlabel('time (s)');
ylabel('Excess Residual Noise Power (dB)');
set(gcf, 'name', 'ex');
set(gcf, 'NumberTitle', 'off');
% save2pdf('CANCx_ex.pdf', 9, 600);
% save2pdf('CANC2_ex.pdf', 9, 600);