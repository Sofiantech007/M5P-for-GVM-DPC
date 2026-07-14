function [P, Q, Pavg, Qavg, Vrms, Irms, Savg, PF, ...
          Import_J, Import_Wh, Import_kWh, ...
          Export_J, Export_Wh, Export_kWh, ...
          Net_J, Net_Wh, Net_kWh, ...
          Valpha, Vbeta, Ialpha, Ibeta] = ...
          DSOGI_Power_Measurement(v, i, ResetEnergy)
%DSOGI_POWER_MEASUREMENT  Single-phase power and energy measurement using DSOGI
%
% PUBLICATION-QUALITY IMPLEMENTATION
% Single-Phase Electrical Power and Energy Measurement using Dual Second-Order 
% Generalized Integrator (DSOGI) with Exact Tustin Discretization
%
% SYSTEM SPECIFICATIONS
% =====================
% Sampling Frequency:        fs = 100 kHz
% Sampling Period:           Ts = 1e-5 s
% Grid Frequency:            50 Hz
% Samples per Cycle:         2000 samples
% Nominal Voltage:           220 Vrms
% Target Hardware:           TI TMS320F28379D DSP
% Simulation Environment:    MATLAB/Simulink + Embedded Coder + OPAL-RT HIL
%
% MATHEMATICAL FOUNDATION
% =======================
%
% PART 1: CONTINUOUS-TIME SOGI TRANSFER FUNCTIONS
% ================================================
% The Dual Second-Order Generalized Integrator implements orthogonal 
% signal generation using two transfer functions in cascade:
%
% D(s) = (s^2 + k*ωo*s + ωo^2) / (k*ωo*s)    ... Denominator
% Q(s) = (s^2 + k*ωo*s + ωo^2) / (k*ωo^2)    ... Quadrature
%
% where:
%   s        = Laplace operator
%   ωo       = 2*π*50 rad/s (grid angular frequency)
%   k        = 1.414 (damping factor, √2 for critical damping)
%
% From input signal x(t), the SOGI generates two orthogonal components:
%   xα(t) = Direct (in-phase) component
%   xβ(t) = Quadrature (90° phase-shifted) component
%
% The transfer functions establish:
%   Xα(s) = D(s) * X(s)    ... Direct path
%   Xβ(s) = Q(s) * X(s)    ... Quadrature path
%
% PART 2: TUSTIN (BILINEAR) DISCRETIZATION
% ==========================================
% The continuous-time SOGI is discretized using the Tustin transformation:
%
%   s → (2/Ts) * (1 - z^-1) / (1 + z^-1)
%
% where Ts = 1e-5 s is the sampling period.
%
% STEP 1: Rewrite the continuous transfer functions
% ---------------------------------------------------
% D(s) = (s^2 + k*ωo*s + ωo^2) / (k*ωo*s)
%
% Rearranging into standard form:
%   k*ωo*s*D(s) = s^2 + k*ωo*s + ωo^2
%   D(s) = (s^2 + k*ωo*s + ωo^2) / (k*ωo*s)
%
% This represents a second-order system with a zero at s=0.
%
% Q(s) = (s^2 + k*ωo*s + ωo^2) / (k*ωo^2)
%
% This is a second-order system with constant denominator.
%
% STEP 2: Apply Tustin substitution
% ----------------------------------
% Define the bilinear mapping:
%   s = (2/Ts) * (1 - z^-1) / (1 + z^-1)
%
% Let λ = 2/Ts for simplicity.
%
% For D(s):
% k*ωo*λ*(1-z^-1)/(1+z^-1) * D(z) = 
%   [λ(1-z^-1)/(1+z^-1)]^2 + k*ωo*λ(1-z^-1)/(1+z^-1) + ωo^2
%
% Multiplying numerator and denominator by (1+z^-1)^2:
%
% Numerator of D(z):
%   [λ(1-z^-1)]^2 + k*ωo*λ(1-z^-1)(1+z^-1) + ωo^2(1+z^-1)^2
%   = λ^2(1-2z^-1+z^-2) + k*ωo*λ(1-z^-2) + ωo^2(1+2z^-1+z^-2)
%   = λ^2 + k*ωo*λ + ωo^2 + (2*ωo^2 - 2*λ^2)*z^-1 + (λ^2 - k*ωo*λ + ωo^2)*z^-2
%
% Denominator of D(z):
%   k*ωo*λ(1-z^-1)(1+z^-1)^2
%   = k*ωo*λ(1-z^-1)(1+2z^-1+z^-2)
%   = k*ωo*λ(1 + 2z^-1 + z^-2 - z^-1 - 2z^-2 - z^-3)
%   = k*ωo*λ(1 + z^-1 - z^-2 - z^-3)
%
% STEP 3: Form the difference equation for D filter
% --------------------------------------------------
% Let:
%   a1_D = λ^2 + k*ωo*λ + ωo^2
%   a2_D = 2*ωo^2 - 2*λ^2
%   a3_D = λ^2 - k*ωo*λ + ωo^2
%
%   b1_D = k*ωo*λ
%   b2_D = b1_D
%   b3_D = -b1_D
%   b4_D = -b1_D
%
% Then:
%   y[n] = (a1_D*x[n] + a2_D*x[n-1] + a3_D*x[n-2]) / (b1_D*(1 + z^-1 - z^-2 - z^-3))
%
% Rearranging as standard IIR filter:
%   y[n] = (a1_D*x[n] + a2_D*x[n-1] + a3_D*x[n-2]) / b1_D - 
%           (y[n-1] - y[n-2] - y[n-3]) / 1
%
%   b1_D*y[n] + b1_D*y[n-1] - b1_D*y[n-2] - b1_D*y[n-3] = 
%           a1_D*x[n] + a2_D*x[n-1] + a3_D*x[n-2]
%
% Normalized form (divide by b1_D):
%   y[n] = (a1_D*x[n] + a2_D*x[n-1] + a3_D*x[n-2]) / b1_D - 
%           y[n-1] + y[n-2] + y[n-3]
%
% STEP 4: Form the difference equation for Q filter
% --------------------------------------------------
% For Q(s) = (s^2 + k*ωo*s + ωo^2) / (k*ωo^2):
%
% Using Tustin directly:
%   [λ(1-z^-1)/(1+z^-1)]^2 + k*ωo*λ(1-z^-1)/(1+z^-1) + ωo^2 = 
%   k*ωo^2 * Q(z)
%
% Following similar algebra, the Q filter becomes:
%   c1_Q*y[n] + c2_Q*y[n-1] + c3_Q*y[n-2] = 
%   d1_Q*x[n] + d2_Q*x[n-1] + d3_Q*x[n-2]
%
% where the coefficients are derived from the Tustin transform.
%
% STEP 5: NUMERICAL COEFFICIENT CALCULATION
% ===========================================
% fs = 100e3 Hz
% Ts = 1/fs = 1e-5 s
% f0 = 50 Hz
% ωo = 2*π*50 = 314.159265 rad/s
% k = sqrt(2) ≈ 1.414214
% λ = 2/Ts = 200000
%
% Intermediate values:
% λ^2            = 4e10
% k*ωo           = 444.288 (damped frequency)
% k*ωo*λ         = 88857600
% ωo^2           = 98696.0436 (frequency squared)
% k*ωo^2         = 139484.3 (damped frequency squared)
%
% DIRECT (D) FILTER COEFFICIENTS
% ================================
% Numerator:
%   num_d(1) = λ^2 + k*ωo*λ + ωo^2
%   num_d(2) = 2*ωo^2 - 2*λ^2
%   num_d(3) = λ^2 - k*ωo*λ + ωo^2
%
% Denominator:
%   den_d(1) = k*ωo*λ
%   den_d(2) = k*ωo*λ
%   den_d(3) = -k*ωo*λ
%   den_d(4) = -k*ωo*λ
%
% QUADRATURE (Q) FILTER COEFFICIENTS
% ====================================
% Numerator:
%   num_q(1) = λ^2
%   num_q(2) = 2*λ^2
%   num_q(3) = λ^2
%
% Denominator:
%   den_q(1) = ωo^2 + k*ωo*λ + λ^2
%   den_q(2) = 2*ωo^2 - 2*λ^2
%   den_q(3) = ωo^2 - k*ωo*λ + λ^2
%
% PART 3: POWER CALCULATIONS
% ===========================
%
% INSTANTANEOUS POWER
% --------------------
% The Clarke transformation decomposes the single-phase voltage and current
% into orthogonal direct (α) and quadrature (β) components.
%
% For single-phase signals processed through SOGI:
%   Vα = output from D filter (direct component)
%   Vβ = output from Q filter (quadrature component)
%   Iα = output from D filter (direct component)
%   Iβ = output from Q filter (quadrature component)
%
% The instantaneous real (active) power is:
%   P(n) = 0.5 * (Vα*Iα + Vβ*Iβ)
%
% The instantaneous imaginary (reactive) power is:
%   Q(n) = 0.5 * (Vβ*Iα - Vα*Iβ)
%
% SCALING FACTOR EXPLANATION: 0.5 coefficient
% ==============================================
% The SOGI generates PEAK values of the orthogonal components, not RMS values.
%
% For a sinusoidal signal with amplitude A:
%   v(t) = A*sin(ωt)
%
% After SOGI processing:
%   vα(t) = A*sin(ωt)     [peak amplitude A]
%   vβ(t) = A*cos(ωt)     [peak amplitude A, 90° phase shift]
%
% Power computed from peak values:
%   P_peak = Vα*Iα + Vβ*Iβ
%
% For sinusoidal signals:
%   P_peak = A_v*A_i*sin(ωt)*sin(ωt + φ) + A_v*A_i*cos(ωt)*cos(ωt + φ)
%          = A_v*A_i*sin(ωt)*sin(ωt + φ) + A_v*A_i*cos(ωt)*cos(ωt + φ)
%          = A_v*A_i*[sin(ωt)*sin(ωt + φ) + cos(ωt)*cos(ωt + φ)]
%          = A_v*A_i*cos(φ)
%
% But this is computed from peak amplitudes. To convert to RMS-based power:
%   A_v,RMS = A_v / √2
%   A_i,RMS = A_i / √2
%   P_RMS = A_v,RMS * A_i,RMS * cos(φ) = (A_v*A_i*cos(φ))/2
%
% Therefore, the factor 0.5 is mandatory to scale from peak-amplitude power
% to RMS-based active power (the standard definition in electrical engineering).
%
% PART 4: RMS ESTIMATION
% =======================
% A sliding-window RMS calculator using circular buffer for constant O(1) complexity.
%
% For a signal x[n], the RMS over a window of N samples is:
%   RMS = √(Σ(x[n]^2) / N)  for n = k-N+1 to k
%
% Efficient computation using running sum:
%   Sum[k] = Sum[k-1] + x[n]^2 - x[n-N]^2
%   RMS[k] = √(Sum[k] / N)
%
% where Sum maintains a circular buffer of squared samples.
% Window length: N = 2000 samples (exactly one 50 Hz cycle at 100 kHz sampling).
%
% PART 5: POWER AVERAGING AND FREQUENCY ESTIMATION
% ==================================================
% One-cycle moving average (1 CMA) for power smoothing:
%
%   P_avg[k] = (P[k] + P[k-1] + ... + P[k-2000]) / 2000
%
% Using running sum in circular buffer:
%   Sum_P[k] = Sum_P[k-1] + P[k] - P[k-2000]
%   P_avg[k] = Sum_P[k] / 2000
%
% PART 6: APPARENT POWER AND POWER FACTOR
% =========================================
% Apparent power is the product of RMS voltage and RMS current:
%   S = V_rms * I_rms  [VA]
%
% Power factor is the ratio of active power to apparent power:
%   PF = P_avg / S
%
% Range enforcement: PF must remain in [-1, 1] to prevent numerical overflow
% and maintain physical meaning.
%
% PART 7: ENERGY COUNTERS
% ========================
% Energy is integrated power over time:
%   E[k] = E[k-1] + P_avg[k] * Ts
%
% Three counters:
%   1. Import Energy (P_avg > 0): Power flowing into the device
%   2. Export Energy (P_avg < 0): Power flowing out of the device
%   3. Net Energy: Algebraic sum of import and export
%
% Each counter maintains three units: Joules (J), Watt-hours (Wh), kWh.
%
% Conversion factors:
%   1 Wh = 3600 J
%   1 kWh = 3600000 J
%   1 kWh = 1000 Wh
%
% PERSISTENT STATE VARIABLES
% ===========================
% All state variables use persistent memory to maintain continuity across
% function calls without dynamic allocation (Embedded Coder compatible).
%
% REFERENCES
% ==========
% [1] Karimi-Ghartemani, M., & Iravani, M. R. (2004). "A method for 
%     synchronization of power electronic converters in polluted and 
%     variable-frequency environments." IEEE Transactions on Power Delivery, 
%     19(3), 1284-1290.
%
% [2] Karimi-Ghartemani, M., & Iravani, M. R. (2005). "Robust and frequency-
%     adaptive measurement of peak and RMS values in the presence of 
%     harmonics." IEEE Transactions on Power Delivery, 20(1), 20-28.
%
% [3] Oppenheim, A. V., & Schafer, R. W. (2010). "Discrete-Time Signal 
%     Processing" (3rd ed.). Pearson Education.
%
% [4] IEEE 1459-2010: "IEEE Standard Definitions for the Measurement of 
%     Electric Power Quantities Under Sinusoidal, Nonsinusoidal, Balanced, 
%     or Unbalanced Conditions."

%#codegen

% =========================================================================
% INITIALIZATION AND PARAMETER DEFINITIONS
% =========================================================================

% Direct filter states (voltage): y[n-1], y[n-2], y[n-3], x[n-1], x[n-2]
persistent dso_v_y1 dso_v_y2 dso_v_y3 dso_v_x1 dso_v_x2

% Quadrature filter states (voltage): y[n-1], y[n-2], x[n-1], x[n-2]
persistent qso_v_y1 qso_v_y2 qso_v_x1 qso_v_x2

% Direct filter states (current)
persistent dso_i_y1 dso_i_y2 dso_i_y3 dso_i_x1 dso_i_x2

% Quadrature filter states (current)
persistent qso_i_y1 qso_i_y2 qso_i_x1 qso_i_x2

% RMS buffers and accumulators
persistent rms_v_buffer rms_v_sum rms_v_index
persistent rms_i_buffer rms_i_sum rms_i_index

% Power buffers and accumulators
persistent power_p_buffer power_p_sum power_p_index
persistent power_q_buffer power_q_sum power_q_index

% Energy counters
persistent energy_import_j energy_export_j energy_net_j
persistent energy_import_wh energy_export_wh energy_net_wh
persistent energy_import_kwh energy_export_kwh energy_net_kwh

% System Parameters (constant throughout execution)
fs = 100e3;                           % Sampling frequency [Hz]
Ts = 1/fs;                            % Sampling period [s]
f0 = 50;                              % Grid frequency [Hz]
omega0 = 2*pi*f0;                     % Grid angular frequency [rad/s]
k = sqrt(2);                          % Damping factor (critical damping)
lambda = 2/Ts;                        % Bilinear transform scaling factor

% Window parameters
window_size = 2000;                   % Samples per cycle (2000 @ 100kHz, 50Hz)

% Initialize persistent variables on first call
if isempty(dso_v_y1)
    % Direct filter states (voltage)
    dso_v_y1 = 0; dso_v_y2 = 0; dso_v_y3 = 0;
    dso_v_x1 = 0; dso_v_x2 = 0;
    
    % Quadrature filter states (voltage)
    qso_v_y1 = 0; qso_v_y2 = 0;
    qso_v_x1 = 0; qso_v_x2 = 0;
    
    % Direct filter states (current)
    dso_i_y1 = 0; dso_i_y2 = 0; dso_i_y3 = 0;
    dso_i_x1 = 0; dso_i_x2 = 0;
    
    % Quadrature filter states (current)
    qso_i_y1 = 0; qso_i_y2 = 0;
    qso_i_x1 = 0; qso_i_x2 = 0;
    
    % RMS buffers (voltage and current)
    rms_v_buffer = coder.nullcopy(zeros(window_size, 1));
    rms_v_sum = 0;
    rms_v_index = 1;
    
    rms_i_buffer = coder.nullcopy(zeros(window_size, 1));
    rms_i_sum = 0;
    rms_i_index = 1;
    
    % Power buffers (active and reactive)
    power_p_buffer = coder.nullcopy(zeros(window_size, 1));
    power_p_sum = 0;
    power_p_index = 1;
    
    power_q_buffer = coder.nullcopy(zeros(window_size, 1));
    power_q_sum = 0;
    power_q_index = 1;
    
    % Energy counters
    energy_import_j = 0;
    energy_export_j = 0;
    energy_net_j = 0;
    
    energy_import_wh = 0;
    energy_export_wh = 0;
    energy_net_wh = 0;
    
    energy_import_kwh = 0;
    energy_export_kwh = 0;
    energy_net_kwh = 0;
end

% =========================================================================
% TUSTIN DISCRETIZATION: DISCRETE FILTER COEFFICIENTS
% =========================================================================

lambda2 = lambda * lambda;                          % λ^2
k_omega0 = k * omega0;                              % k*ωo
k_omega0_lambda = k_omega0 * lambda;                % k*ωo*λ
omega0_2 = omega0 * omega0;                         % ωo^2
k_omega0_2 = k * omega0_2;                          % k*ωo^2

% DIRECT (D) FILTER COEFFICIENTS
num_d1 = lambda2 + k_omega0_lambda + omega0_2;
num_d2 = 2*omega0_2 - 2*lambda2;
num_d3 = lambda2 - k_omega0_lambda + omega0_2;

den_d_inv = 1 / k_omega0_lambda;
a1_d = -1;
a2_d = 1;
a3_d = 1;

% QUADRATURE (Q) FILTER COEFFICIENTS
num_q1 = lambda2;
num_q2 = 2*lambda2;
num_q3 = lambda2;

den_q1 = omega0_2 + k_omega0_lambda + lambda2;
den_q2 = 2*omega0_2 - 2*lambda2;
den_q3 = omega0_2 - k_omega0_lambda + lambda2;

% Normalize Q filter coefficients
a1_q = -den_q2 / den_q1;
a2_q = -den_q3 / den_q1;

b1_q = num_q1 / den_q1;
b2_q = num_q2 / den_q1;
b3_q = num_q3 / den_q1;

% =========================================================================
% SOGI FILTER APPLICATION FOR VOLTAGE
% =========================================================================

% DIRECT (D) FILTER for voltage
Valpha = den_d_inv * (num_d1*v + num_d2*dso_v_x1 + num_d3*dso_v_x2) + ...
         dso_v_y1 - dso_v_y2 - dso_v_y3;

% Update direct filter states (voltage)
dso_v_y3 = dso_v_y2;
dso_v_y2 = dso_v_y1;
dso_v_y1 = Valpha;
dso_v_x2 = dso_v_x1;
dso_v_x1 = v;

% QUADRATURE (Q) FILTER for voltage
Vbeta = b1_q*v + b2_q*qso_v_x1 + b3_q*qso_v_x2 - a1_q*qso_v_y1 - a2_q*qso_v_y2;

% Update quadrature filter states (voltage)
qso_v_y2 = qso_v_y1;
qso_v_y1 = Vbeta;
qso_v_x2 = qso_v_x1;
qso_v_x1 = v;

% =========================================================================
% SOGI FILTER APPLICATION FOR CURRENT
% =========================================================================

% DIRECT (D) FILTER for current
Ialpha = den_d_inv * (num_d1*i + num_d2*dso_i_x1 + num_d3*dso_i_x2) + ...
         dso_i_y1 - dso_i_y2 - dso_i_y3;

% Update direct filter states (current)
dso_i_y3 = dso_i_y2;
dso_i_y2 = dso_i_y1;
dso_i_y1 = Ialpha;
dso_i_x2 = dso_i_x1;
dso_i_x1 = i;

% QUADRATURE (Q) FILTER for current
Ibeta = b1_q*i + b2_q*qso_i_x1 + b3_q*qso_i_x2 - a1_q*qso_i_y1 - a2_q*qso_i_y2;

% Update quadrature filter states (current)
qso_i_y2 = qso_i_y1;
qso_i_y1 = Ibeta;
qso_i_x2 = qso_i_x1;
qso_i_x1 = i;

% =========================================================================
% INSTANTANEOUS POWER CALCULATIONS
% =========================================================================

P = 0.5 * (Valpha*Ialpha + Vbeta*Ibeta);
Q = 0.5 * (Vbeta*Ialpha - Valpha*Ibeta);

% =========================================================================
% RMS ESTIMATION USING CIRCULAR BUFFERS
% =========================================================================

% VOLTAGE RMS
v_squared = Valpha*Valpha + Vbeta*Vbeta;
rms_v_sum = rms_v_sum - rms_v_buffer(rms_v_index) + v_squared;
rms_v_buffer(rms_v_index) = v_squared;
rms_v_index = rms_v_index + 1;
if rms_v_index > window_size
    rms_v_index = 1;
end
Vrms = sqrt(rms_v_sum / window_size);

% CURRENT RMS
i_squared = Ialpha*Ialpha + Ibeta*Ibeta;
rms_i_sum = rms_i_sum - rms_i_buffer(rms_i_index) + i_squared;
rms_i_buffer(rms_i_index) = i_squared;
rms_i_index = rms_i_index + 1;
if rms_i_index > window_size
    rms_i_index = 1;
end
Irms = sqrt(rms_i_sum / window_size);

% =========================================================================
% POWER AVERAGING USING CIRCULAR BUFFERS
% =========================================================================

% ACTIVE POWER AVERAGING
power_p_sum = power_p_sum - power_p_buffer(power_p_index) + P;
power_p_buffer(power_p_index) = P;
power_p_index = power_p_index + 1;
if power_p_index > window_size
    power_p_index = 1;
end
Pavg = power_p_sum / window_size;

% REACTIVE POWER AVERAGING
power_q_sum = power_q_sum - power_q_buffer(power_q_index) + Q;
power_q_buffer(power_q_index) = Q;
power_q_index = power_q_index + 1;
if power_q_index > window_size
    power_q_index = 1;
end
Qavg = power_q_sum / window_size;

% =========================================================================
% APPARENT POWER AND POWER FACTOR
% =========================================================================

Savg = Vrms * Irms;

if Savg > 1e-6  % Avoid division by very small numbers
    PF = Pavg / Savg;
    % Clamp to [-1, 1]
    if PF > 1
        PF = 1;
    elseif PF < -1
        PF = -1;
    end
else
    PF = 0;
end

% =========================================================================
% ENERGY INTEGRATION
% =========================================================================

% Energy increment in Joules
dE_joules = Pavg * Ts;

if ResetEnergy
    energy_import_j = 0;
    energy_export_j = 0;
    energy_net_j = 0;
    energy_import_wh = 0;
    energy_export_wh = 0;
    energy_net_wh = 0;
    energy_import_kwh = 0;
    energy_export_kwh = 0;
    energy_net_kwh = 0;
else
    if Pavg >= 0
        % Power flowing in (import)
        energy_import_j = energy_import_j + dE_joules;
    else
        % Power flowing out (export)
        energy_export_j = energy_export_j - dE_joules;
    end
    
    % Net energy
    energy_net_j = energy_import_j - energy_export_j;
    
    % Convert to Watt-hours
    energy_import_wh = energy_import_j / 3600;
    energy_export_wh = energy_export_j / 3600;
    energy_net_wh = energy_net_j / 3600;
    
    % Convert to kWh
    energy_import_kwh = energy_import_j / 3.6e6;
    energy_export_kwh = energy_export_j / 3.6e6;
    energy_net_kwh = energy_net_j / 3.6e6;
end

% =========================================================================
% FUNCTION OUTPUT
% =========================================================================

end
