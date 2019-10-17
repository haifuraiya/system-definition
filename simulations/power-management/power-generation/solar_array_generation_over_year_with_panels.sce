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


// import the CelestLab constants
CL_init;
total_solar_irradiance = CL_dataGet("totalSolarIrradiance");

// create empty matrices
sim_length = 1;
nbpts = (2*60*24) * sim_length + 1;
z_neg_factor_matrix = zeros(nbpts, 100);
z_pos_factor_matrix = zeros(nbpts, 100);
y_neg_factor_matrix = zeros(nbpts, 100);
y_pos_factor_matrix = zeros(nbpts, 100);
x_neg_factor_matrix = zeros(nbpts, 100);
x_pos_factor_matrix = zeros(nbpts, 100);

// loop through each day in the day
for day = 1:365

    // date/time of orbital elements (TREF)
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

    // calculate the solar array factors
    //  neg angles means the sun is behind the panel and no power should be generated
    z_neg_power = max(0,cos(z_neg_angle(:)) .* eclipse_mask) * solar_cells_z_neg * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
    z_pos_power = max(0,cos(z_pos_angle(:)) .* eclipse_mask) * solar_cells_z_pos * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
    x_neg_power = max(0,cos(x_neg_angle(:)) .* eclipse_mask) * solar_cells_x_neg * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
    x_pos_power = max(0,cos(x_pos_angle(:)) .* eclipse_mask) * solar_cells_x_pos * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
    y_neg_power = max(0,cos(y_neg_angle(:)) .* eclipse_mask) * solar_cells_y_neg * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;
    y_pos_power = max(0,cos(y_pos_angle(:)) .* eclipse_mask) * solar_cells_y_pos * cell_area * total_solar_irradiance * cell_daily_efficiency_drop^day;

    // calculate daily average efficiencies
    z_neg_power_average(day) = mean(z_neg_power);
    z_pos_power_average(day) = mean(z_pos_power);
    y_neg_power_average(day) = mean(y_neg_power);
    y_pos_power_average(day) = mean(y_pos_power);
    x_neg_power_average(day) = mean(x_neg_power);
    x_pos_power_average(day) = mean(x_pos_power);
end

// plot the daily average power generated for each panel
plot(linspace(0, 365, length(z_neg_power_average)), [z_neg_power_average, z_pos_power_average, y_neg_power_average, y_pos_power_average, x_neg_power_average, x_pos_power_average], "thickness", 2);
xtitle("Individual Solar Panel Generated Power", "Days", "Daily Average Power (W)");
legend(["Z neg", "Z pos", "Y neg", "Y pos", "X neg", "X pos"], legend_location="in_lower_right");

// plot the total daily average power generated
total_power_average = (z_neg_power_average + z_pos_power_average + y_neg_power_average + y_pos_power_average + x_neg_power_average + x_pos_power_average)
plot(linspace(0, 365, length(z_neg_power_average)), total_power_average, "thickness", 2);
xtitle("Total Solar Panel Generated Power", "Days", "Daily Average Power (W)");
