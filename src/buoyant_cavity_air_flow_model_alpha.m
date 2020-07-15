function [Ra_c, V_dot] = ...
    buoyant_cavity_air_flow_model_alpha(...
        T_ext,...
        T_0,...
        alpha_0,...
        RH,...
        theta,...
        L,...
        b,...
        h) 

    % Specific heat capacity:
    c_pa = 1006;                        % [J/kg*K]

    % Air density:
    rho = density_a(T_ext, RH);                  % [kg/m^3]
    
    %%% Data processing, Geometry
    % Total roof height: 
    H = sin(theta*pi/180)*L;                        % [m]      

    % Hydraulic diameter:
    D_h = 2*h*b / (h+b);                            % [m]

    syms L_0 real
    eqn = rayleigh(T_ext,H,T_0,b,h,L,L_0,alpha_0,D_h,rho,c_pa)...
        *(1 - L_0/L + L_0/L * exp(-L/L_0))...
        - L_0 == 0;
    %eqn = ((g*beta*rho*H*(T_0-T_ext)/(nu*(S_f + rho * L_0 * alpha_0 / (rho * c_pa * h) * (1 + 10.59*(rho * L_0 * alpha_0 / (rho * c_pa * h) * D_h / mu)^(-0.374)) / (D_h * b) + 2*1500 * (2 * L_0 * alpha_0 / (rho * c_pa * h) * b * h)^0.8))) * Pr_c) * (1 - L_0/L + L_0/L * exp(-L/L_0))-L_0 == 0;
    
    L_0sol = vpasolve(eqn, L_0, [-Inf, Inf]);
        
    [Ra_c, S, Pr_c, Gr_c] = rayleigh_b(T_ext,H,T_0,b,h,L,L_0sol,alpha_0,D_h,rho,c_pa);
    
    u = L_0sol * alpha_0 / (rho * c_pa * h);
    V_dot = u * h * b;    
end

function Ra_c = rayleigh(T_ext,H,T_0,b,h,L,L_0,alpha_0,D_h,rho,c_pa)
    % Gravity constant 
    g = 9.82;                           % [m/s^2] 
    
    % Air temperature
    T_ext_K = T_ext + 273.15;                    % [K]

    % Coefficient of thermal expansion: 
    beta = 1/T_ext_K;                               % [1/K]

    % Dynamic viscosity [Pa*s]:
    C1 = 1.458e-6;
    C2 = 110.4;
    mu = C1 * T_ext_K^(3/2) / (T_ext_K + C2); 
    
    % Frictional flow resistance
    % Formula modified from Kronvall (1980)
    phi = 2/3 + 11/24 * h/b * (2 - h/b);         % Dimensional factor
    S_f = 32 * mu * L / (phi * D_h^2 * b * h);   % [Pa/(m^3/s)]
    
    % Kinematic viscosity:
    nu = mu/rho;                                    % [m^2/s]
    
    u = L_0 * alpha_0 / (rho * c_pa * h);
    Re = rho * u * D_h / mu;
    
    %if Re < 1000
        K_ca = 0.98*Re^(-0.03);
    %else
        K_cb= 10.59*Re^(-0.374);
    %end
    K_c = (K_ca+K_cb)/2;
    
    S_xi = rho * u * (1+K_c) / (D_h * b);
    %S_vent = 1500 * (2*u*b*h)^0.8;
    S = S_f + S_xi;% + 2*S_vent;

    Gr_c = g * beta * rho * H * (T_0 - T_ext) / (nu * S * b);
    Pr_c = nu * rho * c_pa * b / alpha_0;
    Ra_c = Gr_c * Pr_c;
    
end

function [Ra_c, S, Pr_c, Gr_c] = rayleigh_b(T_ext,H,T_0,b,h,L,L_0,alpha_0,D_h,rho,c_pa)
    % Gravity constant 
    g = 9.82;                           % [m/s^2] 
    
    % Air temperature
    T_ext_K = T_ext + 273.15;                    % [K]

    % Coefficient of thermal expansion: 
    beta = 1/T_ext_K;                               % [1/K]

    % Dynamic viscosity [Pa*s]:
    C1 = 1.458e-6;
    C2 = 110.4;
    mu = C1 * T_ext_K^(3/2) / (T_ext_K + C2); 
    
    % Frictional flow resistance
    % Formula modified from Kronvall (1980)
    phi = 2/3 + 11/24 * h/b * (2 - h/b);         % Dimensional factor
    S_f = 32 * mu * L / (phi * D_h^2 * b * h);   % [Pa/(m^3/s)]
    
    % Kinematic viscosity:
    nu = mu/rho;                                    % [m^2/s]
    
    u = L_0 * alpha_0 / (rho * c_pa * h);
    Re = rho * u * D_h / mu;
    
    %if Re < 1000
        K_ca = 0.98*Re^(-0.03);
    %else
        K_cb= 10.59*Re^(-0.374);
    %end
    K_c = (K_ca+K_cb)/2;
    
    S_xi = rho * u * (1+K_c) / (D_h * b);
    S_vent = 1500 * (2*u*b*h)^0.8;
    S = S_f + S_xi;% + 2*S_vent;

    Gr_c = g * beta * rho * H * (T_0 - T_ext) / (nu * S * b);
    Pr_c = nu * rho * c_pa * b / alpha_0;
    Ra_c = Gr_c * Pr_c;
    
end

% Calculates density of moist air
% Input: Air temperature [degC], RH 0.0-1.0 [-]
function rho_a = density_a(T_a, RH)
    T_a_K = T_a+273.15;                % Temperature [K]
    
    p_d = 101325;                      % Partial pressure of dry air
    R_d = 287.058;                     % Specific gas constant of dry air
    
    p_sat = 6.1078*10^(7.5*T_a/T_a_K); % Saturation vapor pressure 
    p_v = RH*p_sat;                    % Vapor pressure
    R_v = 461.495;                     % Specific gas constant for water vapor
    
    rho_a = p_d./(R_d*T_a_K)+p_v./(R_v*T_a_K); % Air density of inlet air
end
