clear

// sum the vectors onto a matrix
function matrix = populate_matrix(values)

  // create matrix of zeros to populate data
  matrix = zeros(length(values), 100)

  // loop through each value marking where an efficiency value is
  for i = 1:size(values)(1)
      matrix(i, int(99*values(i))+1) = 1
  end
  
endfunction


// plot the matrix
function plot_matrix(matrix_var, name)

    // find the size of the input and create vectors for each axis
    nbpts = size(matrix_var)(1)
    hlocv = linspace(0,24,nbpts);
    decl = linspace(0,100,100);    
    [Hlocv, Decl] = ndgrid(hlocv, decl);

    // create figure window
    f=scf();
    a=gca();
    f.visible="on";
    f.immediate_drawing="off";

    // form a colour map
    Nmap = 64;
    f.color_map = 0.4+0.6*jetcolormap(Nmap);
    Noir = addcolor([1,1,1]*0);
    Gris = addcolor([1,1,1]*0.3);
    GrisF = addcolor([1,1,1]*0.2);

    // arrange the bounds
    matrix_var(:,1) = matrix_var(:,1)/365;
    zmin = min(abs(matrix_var));
    zmax = max(abs(matrix_var));
    Sgrayplot(hlocv, decl, abs(matrix_var), colminmax=[1,Nmap], zminmax=[zmin,zmax]);

    // create the x axis settings
    x_ticks = 0:2:24;
    a.x_ticks = tlist(['ticks','locations','labels'], x_ticks, string(x_ticks));
    
    // general setting
    CL_g_stdaxes(a)
    a.data_bounds = [0,min(decl);24,max(decl)];
    a.tight_limits="on";
    a.title.text = ("Solar Power Generation Over a Year - " + name);
    a.x_label.text = "Time (h)";
    a.y_label.text = "Array Efficiency (%)";
    
    // adjustments
    h = CL_g_select(a, "Text");
    h.font_style=8;
    h.font_size=2;
    h.font_foreground=GrisF;
    f.immediate_drawing="on";
    f.visible = "on";
    
    // plot
    scf();
    
endfunction


// import the CelestLab constants
CL_init;

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
        post_eclipse = int(sim_length/time_step - pre_eclipse - eclipse_duration + 1)
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
    
    // sum the factors onto the matrix
    z_neg_factor_matrix = z_neg_factor_matrix + populate_matrix(z_neg_factor);
    z_pos_factor_matrix = z_pos_factor_matrix + populate_matrix(z_pos_factor);
    y_neg_factor_matrix = y_neg_factor_matrix + populate_matrix(y_neg_factor);
    y_pos_factor_matrix = y_pos_factor_matrix + populate_matrix(y_pos_factor);
    x_neg_factor_matrix = x_neg_factor_matrix + populate_matrix(x_neg_factor);
    x_pos_factor_matrix = x_pos_factor_matrix + populate_matrix(x_pos_factor);
    
    // calculate daily average efficiencies
    z_neg_factor_average(day) = mean(z_neg_factor);
    z_pos_factor_average(day) = mean(z_pos_factor);
    y_neg_factor_average(day) = mean(y_neg_factor);
    y_pos_factor_average(day) = mean(y_pos_factor);
    x_neg_factor_average(day) = mean(x_neg_factor);
    x_pos_factor_average(day) = mean(x_pos_factor);

end

// plot the heat map of each direction
plot_matrix(z_neg_factor_matrix, "Z Negative")
plot_matrix(z_pos_factor_matrix, "Z Positive")
plot_matrix(y_neg_factor_matrix, "Y Negative")
plot_matrix(y_pos_factor_matrix, "Y Positive")
plot_matrix(x_neg_factor_matrix, "X Negative")
plot_matrix(x_pos_factor_matrix, "X Positive")


// plot the daily average efficiency
plot(linspace(0, 365, length(z_neg_factor_average)), [100*z_neg_factor_average, 100*z_pos_factor_average, 100*y_neg_factor_average, 100*y_pos_factor_average, 100*x_neg_factor_average, 100*x_pos_factor_average], "thickness", 2);
xtitle("Solar Panel Efficiency", "Days", "Efficiency (%)");
legend(["Z Negative", "Z Positive", "Y Negative", "Y Positive", "X Negative", "X Positive"], legend_location="in_lower_right");



// plot average efficiencies
