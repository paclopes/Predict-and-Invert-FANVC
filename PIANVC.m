% PIANVC 
% 
% Paulo Alexandre Crisóstomo Lopes, 21/4/2022
% Instituto de Engenharia de Sistemas e Computadores - Investigação e Desenvolvimento
% Instituto Superior Técnico
% Universidade de Lisboa

% simulation parameters
Nsim = 10;                  % number of simulations
fs = 2000;                  % sampling frequency
simulation_time = 20000;    % simulation samples
change_at = 2.5*fs;         % 2.5 s
change_time = fs;           % 1s fast
qn_steady = 1/400/fs;       % state noise
qn_change = 1/change_time;  % state noise at change

stability_margin = 0.1;     % distance of model poles and zeros from the unit circle
f = 200*(1:3)';             % primary noise sinusoids frequencies
amplitudes = [0.5, 1.2, 0.3]';      % primary noise sinusoids amplitudes
phases = [56, 170, -23]'*pi/180;    % primary noise sinusoids phases
frequency_noise = 1;        % rms Hz
Nx = 6;                     % model order (size-1)
on_id = 100;                % system identification start
on = 1000;                  % ANC start
qv0 = 0.01;                 % background noise power

% algorithm parameters
N = 16;                     % model order (size-1)
L = 64;                     % length of the filtering horizon
lambda = 0.999;             % forgetting factor of the RLS algorithm
alpha = 1.01;
Lu = 10;                    % saturation of the antinoise signal

% LOGS
log_e = zeros(simulation_time,Nsim);
log_u = zeros(simulation_time,Nsim);
log_z = zeros(simulation_time,Nsim);

warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

for n_sim = 1:Nsim
    tic
    rng(78245+n_sim);
   
    % simulation intializations
    frequency = f'+frequency_noise*randn(simulation_time+L, length(f));
    phase = 2*pi*cumsum(frequency)/fs + phases';
    d0 = sin(phase)*amplitudes;
    d = d0 + sqrt(qv0)*randn(simulation_time+L, 1);  % primary noise signal

    [a,b] = generate_plant(Nx, stability_margin);
    
    uv = zeros(L,1);        % anti-noise buffer
    e1v = zeros(Nx+1,1);    % residual noise minus background noise buffer

    % algorithm initialization
    u = 0;                  % anti-noise signal
    ev = zeros(N+1,1);      % residual noise buffer
    ah = zeros(N+1,1); ah(1) = 1; % \hat{a}: a estimate
    bh = zeros(N+1,1); bh(1) = 1; % \hat{b}: b estimate
    theta = [ah(2:end); bh];
    Pm = eye(2*N+1);         % RLS algorithm covariance matrix
    z = nan;                 % RLS algorithm error
    e0v = zeros(L,1);        % predicted error signal when u becomes 0

    for k = 1:simulation_time
        
        % simulation
        qn = qn_steady + qn_change*(abs(k-change_at)<=change_time/2);
        a = a + sqrt(qn)*randn(Nx+1,1);
        b = b + sqrt(qn)*randn(Nx+1,1);
        if qn > 0
            [a,b] = adjust_plant(a,b,stability_margin);
        end
        
        log_u(k,n_sim) = u;    % logs u(n) and not u(n+1)
              
        uv = [u; uv(1:end-1)]; % simulation and algorithm
        e1v = [0; e1v(1:end-1)];
        e1 = b'*uv(1:Nx+1) - a'*e1v;
        e1v(1) = e1;
        e = e1 + d(k);
        
        % algorithm
        ev = [e; ev(1:end-1)];
        if k >= on_id
            % id

            % RLS to calculate the parameters
            phi = [-ev(2:N+1); uv(1:N+1)]';
            Lm = Pm/lambda;
            Pm = Lm - Lm*phi'*(1+phi*Lm*phi')^-1*phi*Lm;
            z = e-phi*theta;
            theta = theta + Pm*phi'*z;
            ah(2:end) = theta(1:N);
            bh = theta(N+1:end);

            if k>=on
                % prediction
                evx = ev(1:N);
                uvx = uv(1:N+1);
                for n=1:L
                    uvx = [0; uvx(1:end-1)];
                    ex = bh'*uvx - ah(2:end)'*evx;
                    evx = [ex; evx(1:end-1)];
                    e0v(n) = ex;                
                end
                e_signal = [zeros(L,1); e0v];

                % inverse filter
                u_signal = non_causal_filter_v3(...
                    ah, ...
                    bh, ...
                    e_signal,...
                    alpha...
                );

                u = - u_signal(L+1);
                u = min(Lu, max(-Lu,u));
            else
                u = randn;
            end
            
        else
            u = 0;
        end

        log_e(k,n_sim) = e;
        log_z(k,n_sim) = z;
    end

    qe = mean(log_e(end-simulation_time/10:end,n_sim).^2);
    fprintf(1, 'sim: %d, residual noise power: %f\n', n_sim, qe);
    toc
end

warning('on', 'MATLAB:singularMatrix');
warning('on', 'MATLAB:nearlySingularMatrix');

figure(1)
plot_xy_p((0:size(log_e,1)-1)/fs, 10*log10(smooth(log_e.^2,min(round(simulation_time/200),200))));
set(gca, 'YLim', [-20,30]);
xlabel('time (s)');
ylabel('Noise (dB)');
grid on;
title('Residual noise versus time percentile plot');

i = 1:10;
figure(2)
plot((0:size(log_e,1)-1)/fs, 10*log10(smooth(log_e(:,i).^2,min(round(simulation_time/200),200))));
set(gca, 'YLim', [-20,30]);
xlabel('time (s)');
ylabel('Noise (dB)');
grid on;
title('Residual noise versus time plot of 10 simulations');

figure(3);
histogram(10*log10(mean(log_e(end-simulation_time/10:end,:).^2)),-25:10);
xlabel('Noise Power (dB)');
ylabel('Frequency');
grid on;
title('Final residual noise power histogram');

figure(4);
plot(roots(a),'x'); hold on;
plot(roots(b),'o');
plot(exp((0:0.1:2.1*pi)*1i));
%set(gca,'XLim',[-2,2]);
line([0,cos(2*pi*f(1)/fs)], [0,sin(2*pi*f(1)/fs)]);
hold off;
title('Actual model pole zero plot');

figure(5);
plot(roots(ah),'x'); hold on;
plot(roots(bh),'o');
plot(exp((0:0.1:2.1*pi)*1i));
% set(gca,'XLim',[-2,2]);
line([0,cos(2*pi*f(1)/fs)], [0,sin(2*pi*f(1)/fs)]);
hold off;
title('Estimated model pole zero plot');

figure(6);
[h, ~] = freqz(b, a, 1024, fs);
[hh, fx] = freqz(bh, ah, 1024, fs);
subplot(2,1,1);
plot(fx, 20*log10(abs(h))); hold on;
plot(fx, 20*log10(abs(hh))); hold off;
xlabel('frequency (Hz)');
ylabel('amplitude (dB)');
grid on;
title('Actual and Estimated model frequency response');
subplot(2,1,2);
plot(fx, unwrap(angle(h))/pi*180); hold on;
plot(fx, unwrap(angle(hh))/pi*180); hold off;
xlabel('frequency (Hz)');
ylabel('phase (deg)');
grid on;

figure(7);
plot_xy_p((0:size(log_e,1)-1)/fs, 10*log10(smooth(log_z.^2,100)));
grid on;
title('Model identification error signal');

figure(8);
plot(e0v); hold on;
plot(d(k+1:end)); hold off;
legend('predicted', 'actual');
title('Predicted and actual primary noise');
