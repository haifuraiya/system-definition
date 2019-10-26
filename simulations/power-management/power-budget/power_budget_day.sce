clear

//// define the simulation parameters

// keplerian mean orbital elements, frame = ECI
sma = 35786.e3;                                 // semi major axis
ecc = 0;                                        // eccentricity
inc = 0 * %pi/180;                              // inclination
pom = %pi/2;                                    // Argument of perigee
mlh = 0;                                        // MLTAN (hours)
anm = 0;                                        // Mean anomaly

// solar cell information
//  based on SpectroLab XTJ Prime cells
//  using a linear yearly radiation loss of efficiency based on a 
//  15 year GEO (ECSS) drop to 87% of BOL efficiency
cell_area = 0.0027;              // m^2
cell_bol_efficiency = 0.307;
cell_daily_efficiency_drop = nthroot(0.87,15*365);

// solar panel configuration
//solar_cells_x_pos = 20 * 3;        // double deployed panel plus body panel (20 x 30 cm sides)
//solar_cells_x_neg = 20 * 3;        // double deployed panel plus body panel (20 x 30 cm sides)
//solar_cells_y_pos = 0 * 0;         // no panel
//solar_cells_y_neg = 0 * 0;         // no panel
//solar_cells_z_pos = 7 * 1;         // body panel (10 x 30 cm sides)
//solar_cells_z_neg = 0 * 0;         // no panel

solar_cells_x_pos = 20 * 2 + 7 * 1;     // double deployed panel plus body panel (20 x 30 cm sides)
solar_cells_x_neg = 20 * 2 + 7 * 1;     // double deployed panel plus body panel (20 x 30 cm sides)
solar_cells_y_pos = 0 * 0;              // no panel
solar_cells_y_neg = 0 * 0;              // no panel
solar_cells_z_pos = 20 * 2;             // body panel (10 x 30 cm sides)
solar_cells_z_neg = 20 * 2;             // no panel

// solar panel regulator efficiency
solar_array_regulator_efficiency = 0.9;

// battery configuration
battery_start_capacity = 120;        // Watt hour
battery_cycle_life = 500*5;          // based on 100% DoD (Depth of Discharge)
battery_total_charge = 0;            // running total of all battery charge to use in capacity degradation

// load configuration
load_power = 40;                    // orbit average power usage


// import the CelestLab constants
CL_init;
total_solar_irradiance = CL_dataGet("totalSolarIrradiance");

// create empty matrices
sim_length = 1;
nbpts = (2*60*24) * sim_length + 1;

// start with full battery
battery_charge(1) = battery_start_capacity;
battery_current_capacity(1) = battery_start_capacity;


// date/time of orbital elements (TREF)
day = 1;
cjd0 = CL_dat_cal2cjd(2019,3,1,0,0,0) + day;
gom = CL_op_locTime(cjd0, "mlh", mlh, "ra");    // RAAN

// form the keplerian orbital elements vector
kep0 = [%CL_eqRad + sma; ecc; inc; pom; gom; anm]; 

// Simulation dates/times (duration = 1 day)
one_second = 1/86400;
time_step = 30*one_second; // every thirty seconds

// form the new data
cjd = cjd0 + (0 : time_step : sim_length); 

// Propagate with "lydlp" model (output = osculating elements)
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
    
    // compensate for neg earth angle
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
    post_eclipse = int(sim_length/time_step - pre_eclipse - eclipse_duration + 1)
    eclipse_mask = [ones(pre_eclipse,1); zeros(eclipse_duration,1); ones(post_eclipse,1)];
end

// calculate the solar array powers
//  neg angles means the sun is behind the panel and no power should be generated
z_neg_power = max(0,cos(z_neg_angle(:)) .* eclipse_mask) * solar_cells_z_neg * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
z_pos_power = max(0,cos(z_pos_angle(:)) .* eclipse_mask) * solar_cells_z_pos * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
x_neg_power = max(0,cos(x_neg_angle(:)) .* eclipse_mask) * solar_cells_x_neg * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
x_pos_power = max(0,cos(x_pos_angle(:)) .* eclipse_mask) * solar_cells_x_pos * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
y_neg_power = max(0,cos(y_neg_angle(:)) .* eclipse_mask) * solar_cells_y_neg * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
y_pos_power = max(0,cos(y_pos_angle(:)) .* eclipse_mask) * solar_cells_y_pos * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;

// add all the power from each panel
total_power = z_neg_power + z_pos_power + x_neg_power + x_pos_power + y_neg_power + y_pos_power;

// include power convertor efficiencies
total_power = total_power * solar_array_regulator_efficiency;

// loop through finding the power balance
battery_charge(1) = battery_charge($);
for t=2:length(total_power)+1
    
    // calculate the power differential
    power_step = total_power(t-1) - load_power;
    power_integral = power_step * time_step*24;      // Watt hours
    
    // update battery charge, limit to full or zero capacity
    battery_charge(t) = battery_charge(t-1) + power_integral;
    if battery_charge(t) > battery_current_capacity(t-1) then
        battery_charge(t) = battery_current_capacity(t-1);
    elseif battery_charge(t) < 0 then
        battery_charge(t) = 0;
    end
    
    // update the total battery charge sum
    battery_total_charge = battery_total_charge + max(0, battery_charge(t)-battery_charge(t-1));
    
    // update the new battery capacity
    battery_current_capacity(t) = battery_start_capacity * max(0,(battery_cycle_life - battery_total_charge/battery_start_capacity)/battery_cycle_life);
    
    // record the current depth of discharge
    //   this is the ratio of charge present in the battery to the current total capacity available
    battery_dod(t) = battery_charge(t)/battery_current_capacity(t);
end
//previous_battery_charge = battery_current_capacity($);

// save the battery capacity at the end of the day
battery_capacity_day(day) = battery_current_capacity($);

// save the min, mean and max battery charge state for the day
battery_dod_min(day) = min(battery_dod(2:$));
battery_dod_mean(day) = mean(battery_dod(2:$));
battery_dod_max(day) = max(battery_dod(2:$));

previous_battery_charge = battery_current_capacity($);

// plot the generated power
subplot(211)
plot(linspace(0, 24, length(total_power)), total_power, "thickness", 2);
xtitle("Battery State of Charge", "Hours", "Generated Power (W)");

// plot the battery state of charge
subplot(212)
plot(linspace(0, 24, length(battery_charge)), 100*(1 - battery_charge/battery_start_capacity), "thickness", 2);
xtitle("Battery State of Charge", "Hours", "Depth of Discharge (%)");
