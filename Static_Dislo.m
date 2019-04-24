% Simple model to get the formation of LAGB during static annealing
%adopted from my work-hardening model
strainRate  	    = inputs(1);  % Strain rate [1/s]
%maxStrainStep	    = inputs(3);  % Maximum allowed strain step [1]
Temperature 	    = inputs(2);  % Temperature [K]
%E    			    = inputs(3);  % Youngs modulus [Pa]
G   			    = inputs(3);  % Shear modulus [Pa], 
                                  %can include the temperature dependence as in thesis of Peter van Liempt
taylorFactor	    = inputs(4);  % Taylor factor [1]
D      			    = inputs(5);  % Grain size [m]
slipSystems 	    = inputs(6);  % Number of active slip systems [1]
sfe                 = inputs(7); % Stacking fault energy [J/m^2]
alphaPass		    = inputs(8); % Passing coefficient for shear [1] (used in forest hardening here)
burgers			    = inputs(9); % Burgers vector [m]
EdgeDipMinDistance  = inputs(10); % Min displacement in slip dir. [m] ##in dislotwin this is EdgeDipMinDistance
                                  % different values in 3IVM (1 nm) and
                                  % Dislotwin (1 burgers) taken 4.96E-10
                                  % originally
dAnnihilClimb	    = inputs(11); % Min displacement in climb dir. [m]
Q_slip			    = inputs(12); % Actvtn E for cutting process slip [eV], originally 4E-19
Q_climb			    = inputs(13); % Actvtn E for climb during annihil [eV], taken 4.4e-19 (Self diffusion coeff. in iron)
                                  % refer to Graham, Tomlin (1963) or
                                  % Buffington,Hirano,Cohen (1961)
rho_m(1,1)          = inputs(14); % Mobile dislocation density [m^-2]
rho_c(1,1)          = inputs(15); % Immobile dislocation density in the cells [m^-2]
rho_w               = inputs(16); % Immobile dislocation density in the walls [m^-2]
r_forest(1,1)       = inputs(17); % Ratio to decide the separation between cell and wall forest dislocations
K_f                 = inputs(18); % rho_f contribution factor to mean free path, currently 0.25
K_w                 = inputs(19); % rho_w contribution factor to mean free path, currently 0.15
K_d                 = inputs(20); % subgrain size contribution to mean free path, currently 0.25
v_0                 = inputs(21); % velocity pre-factor in the orowan equation, higher values as Karo? or low like dislotwin?originally 1000
c_1                 = inputs(22); % constant for the size effect term in flow stress equation
c_2                 = inputs(23); % constant for the locked dislocation term in flow stress equation
c_4                 = inputs(24); % constant for lock formation, need to check the value properly, originally 0.3 or 0.25
K                   = inputs(25); % similitude constant for initial cell size, ideally 10
theta_HAGB          = inputs(26); % misorientation of HAGB (radians)
nu                  = inputs(27); % poissons ratio
E_crss              = inputs(28); % Activation energy of cross slip, currently also includes the stress contribution, originally 4E-19
atom_freq           = inputs(29); % vibration frequency of atoms
tau_crit            = inputs(30); % critical shear stress required for dislocation movement, originally 1.5E8 (seems okay because austenite has more C than ferrite. Karo takes it 50E6 for ferrite)
p_s                 = inputs(31); % p-exponent in glide velocity
q_s                 = inputs(32); % q-exponent in glide velocity
K_B                 = inputs(33); % Boltzmann constant
M_HAGB              = inputs(34); % pre-exponential factor from Roberts (1978)
slipinteraction     = inputs(35); %interaction coefficient between dislocations on diff slip planes (lattice_interactionslipslip)
climb_freq          = inputs(36); %the climb frequency factor defined as bv by Markus based on Argon Moffat,
                                  %this will be in order of 1000

%________SHEAR MODULUS TEMPERATURE DEPENDENCE___________________________
G = G*(1 - 7.9921E-7*(Temperature^2.0) + 3.3171E-10*(Temperature^3.0));

%________M_HAGB TEMPERATURE DEPENDENCE___________________________
%## mobility temperature effect needs to be formulated like mentioned by
%Humphreys in the unified model. Currently using the method from Wang Shan
%(2011)

M_HAGB = M_HAGB*exp(-158992/(8.314*Temperature))/(8.314*Temperature);

%_________________________________________________________________________

%increment related data
final_strain = 0.1;       %for SRX there is no strain, but this can be time
increments = 10000;       %defined by the user
time_end   = final_strain/strainRate;
delta_t    = time_end/increments;
t = 0;
count_1 = 1;

%------INPUTS FROM HOT DEFORMATION MODEL------------------------------------
rho_w_l(1,1) = 1E12;
rho_w_dip(1,1) = rho_w;          
d_s(1,1) = D;
diff_dipole_distance = 1E-09;
%______TOTAL DISLOCATION DENSITY FOR FOREST________________________________
rho_total(1,1) = rho_m(1,1) + rho_w_dip(1,1);
%______CALCULATING FOREST DISLOCATIONS_____________________________________
rho_f(1,1) = rho_total(1,1)*(slipSystems - 1)/(slipSystems);
%----------Initialization--------------------------------------------------
inv_lambda(1,1) = K_f*sqrt(rho_f(1,1)) + K_w*sqrt(rho_w_l(1,1)) + 1.0/D + K_d/d_s(1,1) ;
lambda(1,1) = 1.0/inv_lambda(1,1);
rel_rho_decrease = 1e9;

%----------GENERATION TERM------------------------------------------------
% No generation  terms for this model. As it is a static RV model
%applying the simplistic model consisting only of mobile, wall and dipoles
%assuming that initially there are just homogeneously distributed
%dislocations

%%-------------LOOPING OVER TIME STEPS------------------------------------- 
while t <= time_end
    delta_t = 0.00003;
    if count_1 > 1
        delta_t = abs(rel_rho_decrease/(drho(count_1 - 1,1)));
        if strain(count_1 - 1,1) > 0.05 
            if log10(strainRate) > -2
                delta_t = 1.0*(10^(-log10(strainRate)))/(10^2.0);
            else
                delta_t = 1;
            end
        end 
    end
    shear_rate = rho_m(count_1,1)*burgers*v_0*exp((-Q_slip/(K_B*Temperature)));
                                         % no tau terms for static cases
    
%=======================================================================
%---------Misorientation behaviour--------------------------------------
M_LAGB(count_1,1) = M_HAGB;
    %% _____________DISLOCATION DENSITY UPDATE AND STORAGE________________________________________
    count_1 = count_1 + 1
    tspan = [t t+delta_t];
    %Mobile dislocations
%     if isreal(DotRhoLockFormation) == 0
%         break
%     end
%     rho_m(count_1,1) = rho_m(count_1 - 1,1) + (DotRhoMultiplication(count_1 - 1,1)...
%                      ...- DotRhoWallDipStress(count_1 - 1,1) ...
%                      - DotRhoLockFormation(count_1 - 1,1) ...
%                      ...- DotRhoAnnihilation(count_1 - 1,1)...
%                      - DotRhoDipFormation(count_1 - 1,1))*(delta_t);
    %---------------Mobile dislocation evolution----------------------------             
%     [t,rho] = ode45(@(t,rho) rate_m(rho_m(count_1 - 1,1),inputs,lambda(count_1 - 1), ...
%                              rho_f(count_1 - 1),diff_dipole_distance),...
%                              tspan,rho_m(count_1 - 1,1));
%     % ode45: ode45(@variables, function to solve, time span, initial value)
%     rho_m(count_1,1) = rho(end:end);
    drho_m_dt        = rate_m(rho_m(count_1 - 1,1),inputs,lambda(count_1 - 1), ...
                             shear_rate,rho_f(count_1 - 1),diff_dipole_distance);
    rho_m(count_1,1) = rho_m(count_1 - 1,1) + drho_m_dt*delta_t;
    drho_m(count_1 - 1,1) = drho_m_dt;
    
%     if rho_m(count_1,1) < rho_m(count_1 - 1,1)
%         rho_m(count_1,1) = rho_m(count_1 - 1,1);
%     end     %avoiding reduction in the dislocation density
    
    %---------------Immobile cell dislocation evolution---------------------------- 
    %Immobile Dislocations in the cells
%     rho_c(count_1,1) = rho_c(count_1 - 1,1) + (DotRhoCellMulti(count_1 - 1,1)...
%                      - DotRhoCellLAGB(count_1 - 1,1)...
%                      - DotRhoCellCross(count_1 - 1,1) ...
%                      - DotRhoCellClimb(count_1 - 1,1))*(delta_t);
%     drho_c_l_dt          = rate_c_l(inputs,rho_f(count_1 - 1),rho_c(count_1 - 1));
%     rho_c(count_1,1)     = rho_c(count_1 - 1,1) + drho_c_l_dt*delta_t;
%     drho_c(count_1 - 1,1)= drho_c_l_dt;
    %Immobile Dislocations in the walls
%     rho_w_l(count_1,1)   = rho_w_l(count_1 - 1,1) + (4.0*DotRhoWallLock(count_1 - 1,1) ...
%                         ... + DotRhoCellLAGB(count_1 - 1,1) ...
%                          - DotRhoWallCross(count_1 -1,1) ...
%                          - DotRhoWallClimb(count_1 - 1,1))*(delta_t); %because even locked dislocations get concentrated in a small area
   
%     d_rho_w_l_dt         = (4.0*DotRhoWallLock(count_1 - 1,1) ...
%                          ... + DotRhoCellLAGB(count_1 - 1,1) ...
%                          - DotRhoWallCross(count_1 -1,1) ...
%                          - DotRhoWallClimb(count_1 - 1,1));

    %---------------Locked dislocations in wall evolution----------------------------                 
%     [t,rhowl] = ode45(@(t,rhowl) rate_w_l(inputs, rho_f(count_1 - 1),rho_w_l(count_1 - 1)),...
%                              tspan,rho_w_l(count_1 - 1,1));
%     rho_w_l(count_1,1) = rhowl(end:end);
    drho_w_l_dt          = rate_w_l(inputs, rho_f(count_1 - 1),rho_w_l(count_1 - 1));
    rho_w_l(count_1,1)   = rho_w_l(count_1 - 1,1) + drho_w_l_dt*delta_t;
    drho(count_1 - 1,1)                 = drho_w_l_dt;                                 
                     
%     rho_w_dip(count_1,1) = rho_w_dip(count_1 - 1,1) + (4.0*DotRhoDipFormation(count_1 - 1,1) ...
%                            ...+ DotRhoWallDipStress(count_1 - 1,1)...
%                            - DotRhoWallDipClimb(count_1 - 1,1) - ...
%                            DotRhoWallDipCross(count_1 - 1,1))*(delta_t) ;%assuming 25% volume of cell walls leading to factor 4

    %---------------Dipole evolution--------------------------------------- 
%     [t,rhodip] = ode45(@(t,rhodip) rate_dip(inputs,rho_w_l(count_1 - 1,1),...
%                                             tau_applied(count_1 - 1,1),...
%                                    diff_dipole_distance,rho_w_dip(count_1 - 1,1),...
%                                    rho_m(count_1 - 1,1)),...
%                              tspan,rho_w_dip(count_1 - 1,1));
%     rho_w_dip(count_1,1) = rhodip(end:end);
    drho_dip_dt          = rate_dip(inputs,rho_w_l(count_1 - 1,1),...
                                            tau_applied(count_1 - 1,1),...
                                   diff_dipole_distance,rho_w_dip(count_1 - 1,1),...
                                   rho_m(count_1 - 1,1));
    rho_w_dip(count_1,1) = rho_w_dip(count_1 - 1,1) + drho_dip_dt*delta_t;
    drho_dip(count_1 - 1,1) = drho_dip_dt;
    
    if rho_w_dip(count_1,1) < rho_w_dip(count_1 - 1,1)
        rho_w_dip(count_1,1) = rho_w_dip(count_1 - 1,1);
    end
    
    %Total dislocations contributing to forest dislocations excluding
    rho_total(count_1,1) = rho_m(count_1,1) + (rho_w_dip(count_1,1)/4.0); 
    %division to account for dipoles being in specific part of volume
    
    %Forest dislocations
    rho_f(count_1,1) = rho_total(count_1,1)*(slipSystems - 1)/(slipSystems);
    
    %FOrest dislocation ratio update
    %r_forest(count_1,1) = rho_c(count_1,1)/(rho_w_l(count_1,1) + rho_w_dip(count_1,1))
    r_forest(count_1,1) = inputs(17);
    
    %=======================================================================
    
    %-------------CELL/SUBGRAIN SIZE UPDATE-----------------------------------------
    %%___________SUBGRAIN SIZE CHANGE______________________________________
%     if rho_w_l(count_1 - 1, 1) ~= 0.0
%         %d_s_reduction(count_1,1) = -0.5*K*(rho_w_l(count_1,1)^(-1.5))*...
%            %                       (rho_w_l(count_1,1) - rho_w_l(count_1 - 1,1));%/delta_t;
%         %d_s_reduction(count_1,1) = K*(rho_w_l(count_1,1))^(-0.5) - K*(rho_w_l(count_1 - 1,1))^(-0.5);
%         d_s_reduction(count_1,1) = -0.5*K*((rho_w_l(count_1,1))^(-1.5))*(rho_w_l(count_1,1) - rho_w_l(count_1 - 1,1));
%         %d_s_reduction(count_1,1) = -0.5*K*((rho_w_l(count_1,1))^(-1.5))*d_rho_w_l_dt;
%     else
%         d_s_reduction(count_1,1) = 0.0
%     end         
%     d_s_increase(count_1,1)  = M_LAGB(count_1 - 1,1)*...
%                     (4.0*gamma(count_1 - 1,1)/d_s(count_1 - 1,1)); %exponential term of temperature not needed as M_HAGB
%     %takes that into account
%     
%     
%     %d_s(count_1,1)  = d_s(count_1 - 1,1) + (d_s_reduction(count_1,1)) + (d_s_increase(count_1,1))*(delta_t);
%     d_s(count_1,1)  = d_s(count_1 - 1,1) + d_s_reduction(count_1,1) ...
%                                          + (d_s_increase(count_1,1))*(delta_t);
%     v_LAGB(count_1,1) = d_s_increase(count_1,1);
%   
    if rho_w_l(count_1 - 1, 1) ~= 0.0
%         [t,dsub] = ode45(@(t,dsub) rate_grain(inputs,rho_w_l(count_1,1),...
%                                               M_LAGB(count_1 - 1,1),...
%                                               gamma(count_1 - 1,1),...
%                                               d_s(count_1 - 1,1),...
%                                               drho_w_l_dt),...
%                                               tspan,d_s(count_1 - 1,1));
%         d_s(count_1,1) = dsub(end:end);
        
        dd_dt          = rate_grain(inputs,rho_w_l(count_1,1),...
                                    M_LAGB(count_1 - 1,1),...
                                    gamma(count_1 - 1,1),...
                                    d_s(count_1 - 1,1),...
                                    drho_w_l_dt);
        d_s(count_1,1) = d_s(count_1 - 1,1) + dd_dt*delta_t;
        d_d_s(count_1,1) = dd_dt;
    else
        d_s(count_1,1) = d_s(count_1 - 1);
    end

    %MEAN FREE PATH UPDATE
    inv_lambda(count_1,1) = K_f*sqrt(rho_f(count_1,1)) + K_w*sqrt(rho_w_l(count_1,1)/4.0) + 1.0/D + K_d/d_s(count_1,1) ;
    %inv_lambda(count_1,1) = K_f*sqrt(rho_f(count_1,1)) + K_d/d_s(count_1,1) ;
    %dividing rho_w_l by fraction because here again we assume homogeneous 
    %distribution
    lambda(count_1,1) = 1.0/inv_lambda(count_1,1);
    t = tspan(end:end);    
 
end
 %%Plot Section
 
 
 
 %%========================================================================
 % ------------------FUNCTIONS---------------------------------------------
 %========================================================================
 %-------------------MOBILE DISLOCATIONS-----------------------------------
 function drho_m_dt = rate_m(rho_m,inputs,lambda, shear_rate,rho_f,diff_dipole_distance)
 %declaring variables
 burgers         = inputs(9);
 slipSystems     = inputs(6);
 c_4             = inputs(24);
 EdgeDipMinDistance = inputs(10);
 slipSystems     = inputs(6);
 
 %calculations
 multiplication  = shear_rate/(burgers*lambda);
 annihilation    = 2.0*(EdgeDipMinDistance)*shear_rate*rho_m/(burgers*slipSystems);
 lockformation   = (c_4*shear_rate*sqrt(rho_f)*(slipSystems - 1)...
     .../(r_forest(count_1,1) + 1.0))...
     /(burgers*slipSystems));
 if diff_dipole_distance > 0.0
     dipole_formation = 2.0*(diff_dipole_distance)* ...
         shear_rate*rho_m/...
         (burgers*slipSystems);
 else
     dipole_formation = 0.0;
 end
 
 drho_m_dt        = multiplication - lockformation - dipole_formation ...
                    - annihilation;
 end
 
 
 %------------------WALL-LOCKED DISLOCATIONS-------------------------------
function drho_w_l_dt = rate_w_l(inputs, rho_f,rho_w_l)
%declaring variables
c_4             = inputs(24);
shear_rate      = inputs(1)*inputs(4);
burgers         = inputs(9);
G               = inputs(3);
nu              = inputs(27);
climb_freq      = inputs(36);
sfe             = inputs(7);
Q_climb			= inputs(13);
K_B             = 1.38e-23;
Temperature 	= inputs(2); 
EdgeDipMinDistance  = inputs(10);
atom_freq           = inputs(29);
E_crss              = inputs(28);
slipSystems     = inputs(6);
%r_forest        = inputs(17);

G = G*(1 - 7.9921E-7*(Temperature^2.0) + 3.3171E-10*(Temperature^3.0));
%rho_fw                = (1 - r_forest)*rho_f;
%calculations
b_d                    = 24.0*pi*((1.0 - nu)/(2.0 + nu))*(sfe/(G*burgers));

generation             = c_4*shear_rate*sqrt(rho_f)*(slipSystems - 1)...
                                      .../(r_forest(count_1,1) + 1.0))...
                                      /(burgers*slipSystems) ;

sigma_wall_climb       = G*burgers*sqrt(rho_w_l)/(2.0*pi*(1 - nu)) ;
v_wallclimb            = 2.0*climb_freq*11.0*(b_d^2.0)*...
                                     exp( -Q_climb/(K_B*Temperature))*...
                                    (exp(abs(sigma_wall_climb)*(burgers^3.0)/(K_B*Temperature)) - 1.0);
climb                  = 2.0*v_wallclimb*(5.0*EdgeDipMinDistance)...
                            *((rho_w_l)^2.0);
cross                  = 2.0*(burgers^2)*atom_freq*...
                         exp(-(E_crss - abs(sigma_wall_climb)*(burgers^3.0))/(K_B*Temperature))*...
                          ((rho_w_l)^2)/4.0;
drho_w_l_dt            = 4.0*generation - climb - cross;
end
 
 %------------------WALL-DIPOLES-------------------------------------------
 function drho_dip_dt = rate_dip(inputs,rho_w_l,tau_applied,diff_dipole_distance,rho_w_dip,rho_m)
%declaring variables
climb_freq          = inputs(36);
Q_climb			    = inputs(13);
sfe                 = inputs(7);
G   			    = inputs(3);
burgers			    = inputs(9);
K_B                 = 1.38e-23;
Temperature 	    = inputs(2); 
EdgeDipMinDistance  = inputs(10);
atom_freq           = inputs(29); 
E_crss              = inputs(28);
slipSystems 	    = inputs(6);
shear_rate      = inputs(1)*inputs(4);
nu              = inputs(27);

%calculations
if diff_dipole_distance > 0.0
    dipole_formation = 2.0*(diff_dipole_distance)* ...
                        shear_rate*rho_m/...
                        (burgers*slipSystems);
else 
    dipole_formation = 0.0;
end

b_d               = 24.0*pi*((1.0 - nu)/(2.0 + nu))*(sfe/(G*burgers));
sigma_climb       = tau_applied*1.0 ;
vclimb            = 2.0*climb_freq*11.0*(b_d^2.0)*... 
                      exp( -(Q_climb)/(K_B*Temperature))*... 
                      (exp(abs(sigma_climb)*(burgers^3.0)/(K_B*Temperature)) - 1.0)  ;

                  % v_0*11 is the freq factor for climb
    if (diff_dipole_distance) > 0.0 
        climb  = rho_w_dip*4.0*vclimb/...
                                         (diff_dipole_distance);
    else
        climb  = 2.0*vclimb*(5.0*EdgeDipMinDistance)...
                                        *((rho_w_l)^2.0); 
    end

    cross  = 2.0*(burgers^2)*atom_freq*...
             exp(-(E_crss - (abs(sigma_climb)*(burgers^3.0)))/(K_B*Temperature))*...
                          ((rho_w_dip)^2)/4.0;

drho_dip_dt = 4.0*dipole_formation - climb -cross;
 end

 %--------------------------------------------------------------------------
%------------------Subgrain diameters--------------------------------------
function dd_dt = rate_grain(inputs,rho_w_l,M_LAGB,gamma,d_s,d_rho_w_l_dt)
%declaring variables
K                   = inputs(25);

%calculations
increase = M_LAGB*(4.0*gamma/d_s);

decrease = -0.5*(10)*((rho_w_l)^(-1.5))*d_rho_w_l_dt;
%decrease = d_lambda_dt;

dd_dt    = increase + decrease;

end