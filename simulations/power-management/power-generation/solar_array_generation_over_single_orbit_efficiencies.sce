clear
CL_init()

// date/time of orbital elements (TREF)
cjd0 = CL_dat_cal2cjd(2019,6,3,0,0,0); 
// keplerian mean orbital elements, frame = ECI
sma = 35786.e3;                                 // semi major axis
ecc = 0;                                        // eccentricity
inc = 0 * %pi/180;                              // inclination
pom = %pi/2;                                    // Argument of perigee
mlh = 0;                                        // MLTAN (hours)
gom = CL_op_locTime(cjd0, "mlh", mlh, "ra");    // RAAN
anm = 0;                                        // Mean anomaly

// form the keplerian orbital elements vector
kep0 = [%CL_eqRad + sma; ecc; inc; pom; gom; anm]; 

// Simulation dates/times (duration = 1 day)
one_second = 1/86400;
time_step = 30*one_second;          // every thirty seconds 
sim_length = 1                      // days
cjd = cjd0 + (0 : time_step : sim_length); 

// Propagate with "lydlp" model (output = osculating elements)
//kep = CL_ex_propagate("lydlp", "kep", cjd0, kep0, cjd, "o"); 
kep = CL_ex_propagate("lydlp", "kep", cjd0, kep0, cjd, "o"); 

// Position and velocity in ECI
[pos_eci, vel_eci] = CL_oe_kep2car(kep);
sun_eci = CL_eph_sun(cjd);

// Calculate the vectors to the earths center and sun from spacecraft
earth_vector = -pos_eci
sun_vector = sun_eci - pos_eci

// create the six satellite vectors
for n = 1:size(pos_eci)(2)
    
    // find orthogonal vectors to the vector from satellite to earth
    perp_angles = kernel(pos_eci(:,n).')

    z_neg(:,n) = pos_eci(:,n)
    z_pos(:,n) = -pos_eci(:,n)
    y_neg(:,n) = perp_angles(:,1)
    y_pos(:,n) = -perp_angles(:,1)   
    
    // compensate for negative earth angle
    if pos_eci(1,n) > 0
        x_neg(:,n) = perp_angles(:,2)
        x_pos(:,n) = -perp_angles(:,2)
    else
        x_neg(:,n) = -perp_angles(:,2)
        x_pos(:,n) = perp_angles(:,2)
    end
end



// calculate the angle between the vectors
z_neg_angle = CL_vectAngle(z_neg, sun_eci);
z_pos_angle = CL_vectAngle(z_pos, sun_eci);
x_neg_angle = CL_vectAngle(x_neg, sun_eci);
x_pos_angle = CL_vectAngle(x_pos, sun_eci);
y_neg_angle = CL_vectAngle(y_neg, sun_eci);
y_pos_angle = CL_vectAngle(y_pos, sun_eci);

// find when we're in eclipse
interv = CL_ev_eclipse(cjd, pos_eci, sun_eci, typ = "umb");

// if no eclipse cleanse the data
if interv == [] then
    eclipse_mask = ones(size(pos_eci)(2),1);
else
    
    // calculate time steps when the spacecraft is in eclipse
    pre_eclipse = int((interv(1)-cjd0)/time_step);
    eclipse_duration = int(2*(interv(2) - interv(1))/time_step + 1);
    post_eclipse = int(1/time_step - pre_eclipse - eclipse_duration + 1)
    eclipse_mask = [ones(pre_eclipse,1); zeros(eclipse_duration,1); ones(post_eclipse,1)];
end

// calculate the solar array factors
//  negative angles means the sun is behind the panel and no power should be generated
z_neg_factor = max(0,cos(z_neg_angle(:)) .* eclipse_mask);
z_pos_factor = max(0,cos(z_pos_angle(:)) .* eclipse_mask);
x_neg_factor = max(0,cos(x_neg_angle(:)) .* eclipse_mask);
x_pos_factor = max(0,cos(x_pos_angle(:)) .* eclipse_mask);
y_neg_factor = max(0,cos(y_neg_angle(:)) .* eclipse_mask);
y_pos_factor = max(0,cos(y_pos_angle(:)) .* eclipse_mask);


// plot
scf();


// create a time vector over a day
time = 0:time_step:sim_length;

// plot the results
plot(time', [z_neg_factor, z_pos_factor, x_neg_factor, x_pos_factor], "thickness", 2);
xtitle("Solar Panel Efficiency", "Days", "Efficiency (%)");
legend(["Z Negative", "Z Positive", "X Negative", "X Positive", "Z Negative"]);

// print average efficiencies
printf("Average z negative efficiency = %0.2f\n", mean(z_neg_factor));
printf("Average z positive efficiency = %0.2f\n", mean(z_pos_factor));
printf("Average x negative efficiency = %0.2f\n", mean(x_neg_factor));
printf("Average x positive efficiency = %0.2f\n", mean(x_pos_factor));
printf("Average y negative efficiency = %0.2f\n", mean(y_neg_factor));
printf("Average y positive efficiency = %0.2f\n", mean(y_pos_factor));
