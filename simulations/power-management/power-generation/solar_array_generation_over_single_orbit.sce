// date/time of orbital elements (TREF)
cjd0 = CL_dat_cal2cjd(2020,3,2,0,0,0); 

// keplerian mean orbital elements, frame = ECI
sma = 35786.e3;                                 // semi major axis
ecc = 0;                                        // eccentricity
inc = 0 * %pi/180;                              // inclination
pom = %pi/2;                                    // Argument of perigee
mlh = 0;                                        // MLTAN (hours)
gom = CL_op_locTime(cjd0, "mlh", mlh, "ra");    // RAAN
anm = 0;                                        // Mean anomaly

// form the keplerian orbital elements vector
kep0 = [sma; ecc; inc; pom; gom; anm]; 

// Simulation dates/times (duration = 1 day)
one_second = 1/86400;
time_step = 30*one_second; // every thirty seconds
cjd = cjd0 + (0 : time_step : 1); 

// Propagate with "lydlp" model (output = osculating elements)
kep = CL_ex_propagate("lydlp", "kep", cjd0, kep0, cjd, "o"); 

// Position and velocity in ECI
[pos_eci, vel_eci] = CL_oe_kep2car(kep);
sun_eci = CL_eph_sun(cjd);

// Calculate the vectors to the earths center and sun from spacecraft
earth_vector = -pos_eci
sun_vector = sun_eci - pos_eci

// calculate the angle between the vectors
sc_sun_angle = CL_vectAngle(pos_eci, sun_eci)

// find when we're in eclipse
interv = CL_ev_eclipse(cjd, pos_eci, sun_eci, typ = "umb");

// calculate time steps when the spacecraft is in eclipse
pre_eclipse = (interv(1,:)-cjd0)/time_step;
eclipse_duration = 2*(interv(2,:) - interv(1,:))/time_step + 1;
post_eclipse = 1/time_step - pre_eclipse - eclipse_duration + 2;
eclipse_mask = [ones(pre_eclipse,1); zeros(eclipse_duration,1); ones(post_eclipse,1)];

// calculate the solar array factors
//  negative angles means the sun is behind the panel and no power should be generated
nadir_factor = max(0,cos(sc_sun_angle(:)) .* eclipse_mask)

// plot
scf();


// create a time vector over a day
time = 0:time_step:1

// plot the results
plot(time, [100*nadir_factor'], "thickness", 2); 
xtitle("Solar Panel Efficiency", "Days", "Efficiency (%)"); 
CL_g_stdaxes();
