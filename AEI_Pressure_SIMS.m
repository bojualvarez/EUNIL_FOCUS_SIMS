function RecArray_SIM_forRuss_v2(TransParams)
%% 
%This script will execute a simulation of a phased transducer array for use in
%acoustoelectric cardiac imaging (ACI)
%Arguments:
%   TransParams - A 1x1 structural array made up of the following structures:
%       Shape      (string): String of values indicating the shape of the array
%       ElemX      (double): Number of elements in the x direction
%       ElemY      (double): Number of elements in the y direction
%       KerfX      (double): end-to-end spacing of the elements in the x direction (in m)
%       KerfY      (double): end-to-end spacing of the elements in the y direction (in m)
%       Width      (double): width of each element in the array (in m) 
%       Height     (double): height of each element in the array (in m)
%       FocX       (double): lateral (x) focus of the array (in m)
%       FocY       (double): elevational (y) focus of the array (in m)
%       FocZ       (double): axial (z) focus of the array (in m)
%       Rad        (double): radius of curvature of the array (in m)
%       Pressure   (double): peak pressure of the ultrasound array (in MPa)
%       CenFreq    (double): center frequency of the ultrasound array (in MHz)
%       Power      (double): input power into the array for simulation (in V)
%       Efficiency (double): efficiency of power conversion (between 0 and 1)
%Modified from Focus_Array_Pulse2.m written by Yexian Qin/Andy Tseng
%Written by Alexander Alvarez 12/18/17;

%% Define transducer array parameters from the input Structural Array

Shape   = TransParams.Shape; % Define the variable Shape as the corresponding element in the TransParams array
elem_x  = TransParams.ElemX; % Define the variable elem_x as the corresponding element in the TransParams array
elem_y  = TransParams.ElemY; % Define the variable elem_y as the corresponding element in the TransParams array
kerf_x  = TransParams.KerfX; % Define the variable kerf_x as the corresponding element in the TransParams array
kerf_y  = TransParams.KerfY; % Define the variable kerf_y as the corresponding element in the TransParams array
width   = TransParams.Width; % Define the variable width as the corresponding element in the TransParams array
height  = TransParams.Height; % Define the variable height as the corresponding element in the TransParams array
focus_x = TransParams.FocX; % Define the variable focus_x as the corresponding element in the TransParams array
focus_y = TransParams.FocY; % Define the variable focus_y as the corresponding element in the TransParams array
focus_z = TransParams.FocZ; % Define the variable focuz_z as the corresponding element in the TransParams array
r_curv  = TransParams.Rad; % Define the variable r_curv as the corresponding element in the TransParams array
P0 = TransParams.Pressure; %Define a variable P0 as the corresponding Pressure value in the TransParams structural array
f0 = TransParams.CenFreq; %Define a variable f0 as the corresponding CenFreq value in the TransParams structural array
elec_pwr_in = TransParams.Power; %Define a variable elec_pwr_in as the corresponding Power value in the TransParams structural array
conv_eff = TransParams.Efficiency; %Define a variable conv_eff as the corresponding Efficiency value in the TransParams structural array

%% Use create_rect functions in FOCUS to create a transducer of the shape defined in the TransParams structural array
% Inputs for the create_rect functions: 
%   elem_x (double): the number of elements in the x direction
%   elem_y (double): the number of elements in the y direction
%   width  (double): with of each element in the array (in m)
%   height (double): height of each element in the array (in m)
%   kerf_x (double): the edge-to-edge spacing between elements in the x direction (in m)
%   kerf_y (double): the edge-to-edge spacing between elements in the y directon (in m)
%   r_curv (double): radius of curvature of the array (in m)
% Outputs for the create_rect functions:
%   xdcr_array: an 8x8 array of FOCUS Transducer Structures

%Conditional definitions depending on array shape (and creation of
%arrays for simulations)
if strcmp(Shape,'Cylindrical Section') == 1
    xdcr_array = create_rect_csa(elem_x, elem_y, width, height, kerf_x, kerf_y, r_curv); %Create the array using the create_rect_csa function included in FOCUS
elseif strcmp(Shape,'Curved Strip') == 1
    xdcr_array = create_rect_curved_strip_array(elem_x, elem_y, width, height, kerf_x, r_curv); %Create the array using the create_rect_curved_strip_array function included in FOCUS
elseif strcmp(Shape,'Enclosed Cylindrical Section') == 1
    xdcr_array = create_rect_enclosed_csa(elem_x, elem_y, width, height, kerf_y, r_curv); %Create the array using the create_rect_enclosed_csa function included in FOCUS
elseif strcmp(Shape,'Planar') == 1
    xdcr_array = create_rect_planar_array(elem_x, elem_y, width, height, kerf_x, kerf_y); %Create the array using the create_rect_planar_array function included in FOCUS
end

%% Display the transducer array design

%Open figure and clear figure
figure();clf;
%Change the background color of the figure
set(gcf,'color','w');
%Use the FOCUS function draw_array to create a map of the transducer
%elements in 3D space
draw_array(xdcr_array);
%Set the color map
set(gca,'fontsize',12);
grid off
%Set the labels for the array map
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

%% Define Parameters for the simulation
fs = 50e6; %sampling frequency for the simulation (in Hz)
dt = 1/fs; %Find the steps in t for the sampling regime
ncycles = 1; %number of cycles for the ultrasound parameters to run through
prf = 1e3;  % pulse repetition frequency
height = height/elem_y; % Change the height definition from the one given by the curved strip
Plane_cell = {'xy','yz', 'xz'}; %TransParams.Plane; % Define the plane to simulate in based on the corresponding function

% Define the length of the array in the x and y direction as the number of
% elements times the width of each element + the number of elements minus 1
% times the kerf spacing 
x_len = elem_x * width + (elem_x - 1) * kerf_x; % Define x_len (in m)
y_len = elem_y * height + (elem_y - 1) * kerf_y; % Define y_len (in m)
fprintf('X_len = %f\n',x_len);
fprintf('Y_len = %f\n',y_len);


% Define the media for the simulation to take place in using the 
% define_media function included with the FOCUS software which will allow
% population with pre-defined media expressions. The following specific
% parameters can also be defined in the define_media expression:
%   cb: specific heat of blood,  units are J/kg/C
%   wb: blood perfusion (optional), units are kg/m^3/second
%   rho: density of the tissue, units are kg/m^3
%   c: speed of sound, units are m/s
%   b: power law coefficient for attenuation
%   atten: calculated attenuation
%   ct: specific heat of tissue, units are J/kg/C
%   kappa: thermal conductivity, units are W/m/C
%   beta: nonlinear parameter 
%   atten_coeff: attenuation coefficient, db/Cm/MHz
define_media();

%Define the medium to simulate in  
medium = lossless; % Define the medium to simulate in

% Also available are the following types of pre-defined acoustic media
%   medium = lossless;
%   medium = water;
%   medium = attenuated; 
%   medium = fat;
%   medium = skin; 
%   medium = muscle;
%   medium = liver

lambda = medium.soundspeed/f0; % Define the wavelength of the ultrasound in the mediumn

%% Coordinate Grid Setup
% Use a case-insensitive string comparison (strcmpi) to determine what 
% plane the user wants to simulate in; these variables will be used to
% set up our coordinate grid to cover the width of the transducer
% array on the x axis and from 0 to 15 wavelengths on the z axis. 
for i = 1:size(Plane_cell,2)
    plane = Plane_cell{i};
    if(strcmpi(plane,'xy')) % If the xy plane is chosen

        xmin = -(width + kerf_x) * (elem_x/2 + 1); % Define variable xmin as the lowest x plane (in m)
        xmax = (width + kerf_x) * (elem_x/2 + 1); % Define variable xmax as the highest x plane (in m)
        ymin = -(height + kerf_y) * (elem_y/2 + 1); % Define variable ymin as the lowest y plane (in m)
        ymax = (height + kerf_y) * (elem_y/2 + 1); % Define variable ymax as the highest y plane (in m)
        zmin = focus_z; % Define variable zmin as the lowest z plane (in m); in this case, it will be the focus of the transducer (since we are operating in the xy plane)
        zmax = focus_z; % Define variable zmax as the highest z plane (in m); in this case, it will be the focus of the transducer (since we are operating in the xy plane)

        % or why xpoints and ypoints are defined as this many

        xpoints = 400; % Define variable xpoints as the number of samples in the x direction
        ypoints = 400; % Define variable ypoints as the number of samples in the y direction

        dx = (xmax-xmin)/xpoints; % Define the step size in the x direction as the range divided by the number of samples
        dy = (ymax-ymin)/ypoints; % Define the step size in the y direction as the range divided by the number of samples
        dz = 1; % Define the step size in the z direction as 1 because we are operating in the xy plane

        x_lim = [-1,1]; % Define the limit in the x direction for later plotting as the x limit
        y_lim = [-1,1]; % Define the limit in the y direction for later plotting as the y limit

        x = xmin:dx:xmax; % Define the array x as dx-spaced points between the minimum and maximum x values defined in the coordinate grid above
        y = ymin:dy:ymax; % Define the array y as dy-spaced points between the minimum and maximum y values defined in the coordinate grid above
        z = zmin:dz:zmax; % Define the array z as dz-spaced points between the minimum and maximum z values defined in the coordinate grid above
        disp(strcat('dx = ', num2str(dx)));
        disp(strcat('dy = ', num2str(dy)));
    elseif (strcmpi(plane,'xz'))

        % I am not sure why the min and max values are defined as they are

        xmin = -(width + kerf_x) * (elem_x/2 + 1); % Define variable xmin as the lowest x plane (in m)
        xmax = (width + kerf_x) * (elem_x/2 + 1); % Define variable xmax as the highest x plane (in m)
        ymin = 0; % Define variable ymin as the lowest y plane (in m); in this case, it will be 0 since we are in the xz plane
        ymax = 0; % Define variable ymax as the highest y plane (in m); in this case, it will be 0 since we are in the xz plane   
        zmin = 40 * lambda; % Define variable zmin as the lowest z plane (in m)
        zmax = 60 * lambda; % Define variable zmax as the highest z plane (in m)

        % or why xpoints and zpoints are defined as this many

        xpoints = 200; % Define variable xpoints as the number of samples in the x direction
        zpoints = 200; % Define variable zpoints as the number of samples in the z direction

        dx = (xmax-xmin)/xpoints; % Define the step size in the x direction as the range divided by the number of samples
        dy = 1; % Define the step size in the y direction as 1 because we are operating in the xz plane
        dz = (zmax-zmin)/zpoints; % Define the step size in the z direction as the range divided by the number of samples

        x_lim = [-1,1]; % Define the limit in the x direction for later plotting as the x limit
        y_lim = [4,6]; % Define the limit in the z direction for later plotting as the y limit

        x = xmin:dx:xmax; % Define the array x as dx-spaced points between the minimum and maximum x values defined in the coordinate grid above
        y = ymin:dy:ymax; % Define the array y as dy-spaced points between the minimum and maximum y values defined in the coordinate grid above
        z = zmin:dz:zmax; % Define the array z as dz-spaced points between the minimum and maximum z values defined in the coordinate grid above
        disp(strcat('dz = ', num2str(dz)));
    elseif (strcmpi(plane,'yz'))
        
        xmin = 0; % Define variable xmin as the lowest x plane (in m); in this case, it will be 0 since we are in the yz plane
        xmax = 0; % Define variable xmax as the highest x plane (in m); in this case, it will be 0 since we are in the yz plane
        ymin = -(height + kerf_y) * (elem_y/2 + 1); % Define variable ymin as the lowest y plane (in m)
        ymax = (height + kerf_y) * (elem_y/2 + 1); % Define variable ymax as the highest y plane (in m)
        zmin = 40 * lambda; % Define variable zmin as the lowest z plane (in m)
        zmax = 60 * lambda; % Define variable zmax as the highest z plane (in m)

        % or why ypoints and zpoints are defined as this many

        ypoints = 200; % Define variable xpoints as the number of samples in the y direction
        zpoints = 200; % Define variable ypoints as the number of samples in the z direction

        dx = 1; % Define the step size in the x direction as 1 because we are operating in the yz plane
        dy = (ymax-ymin)/ypoints; % Define the step size in the y direction as the range divided by the number of samples    
        dz = (zmax-zmin)/zpoints; % Define the step size in the z direction as the range divided by the number of samples

        x_lim = [-1,1]; % Define the limit in the x direction for later plotting as the x limit
        y_lim = [4,6]; % Define the limit in the z direction for later plotting as the y limit
        
        x = xmin:dx:xmax; % Define the array x as dx-spaced points between the minimum and maximum x values defined in the coordinate grid above
        y = ymin:dy:ymax; % Define the array y as dy-spaced points between the minimum and maximum y values defined in the coordinate grid above
        z = zmin:dz:zmax; % Define the array z as dz-spaced points between the minimum and maximum z values defined in the coordinate grid above

    end

    coord_grid=set_coordinate_grid([dx,dy,dz], xmin, xmax, ymin, ymax, zmin, zmax);

    % Set up a time sampling structure that tells fnm_tsd to sample from
    % t=0 to t=2 periods at intervals of the period of the sampling frequency
    % (5MHz). 

    tmin = 1/f0;
    tmax = 100*tmin;
    t = tmin:dt:tmax; % Define variable t as the time structure between tmin and tmax with steps of 1/sampling rate
    T = max(t)-min(t); % Define variable T as the total time in the simulation
    time_struct = set_time_samples(dt, tmin, tmax);
    % This sets our excitation function to be a Hanning weighted tone burst
    % with an amplitude of 1 and center frequency of f0 as defined in TransParams. See the
    % documentation for set_excitation_function for details on how FOCUS
    % handles excitation functions. The next step is to focus the array.
    disp(strcat('dt = ' , num2str(dt)));
    input_func = set_excitation_function(2, f0, ncycles/f0, 0);

    disp(['Focusing array at (', num2str(focus_x), ', ', num2str(focus_y), ...
        ', ', num2str(focus_z), ')']);

    %% Time Delays

    ndiv = find_ndiv(xdcr_array,coord_grid,medium,f0,0.2);
% 
%     % Use the FOCUS function set_time delays to determine the time delays 
%     % required to focus a transducer array at a given point. The output
%     % transducer array - xdcr_array - will have these adjusted time delays 
%     xdcr_array = set_time_delays(xdcr_array, focus_x, focus_y, focus_z, medium, fs, 1);
% 
%     % Use the FOCUS function find_single_focus_phase will produce an array 
%     % of complex_weight shifts to focus an array defined by xdcr_array at a 
%     % given point (x,y,z). The output transducer array - xdcr_array - will
%     % have these adjusted complex shifts
    xdcr_array = find_single_focus_phase(xdcr_array, focus_x, focus_y, focus_z, medium, f0, ndiv, 1);
    focal_pt = [focus_x, focus_y, focus_z]; % Define a vector focal_pt as the focus in the x, y, and z directions
    tdelay = zeros(size(xdcr_array)); % Set a matrix tdelay as a matrix of zeros the size of the xdcr_array structural array

    for ni=1:numel(tdelay) % iterate over the number of elements in tdelay using ni as the iterator 
        pos = xdcr_array(ni).center; % Define the temperary variable pos as the center value of xdcr_array at the element number defined by the iteration
        tdel = sqrt(sum((focal_pt - pos).^2))/medium.soundspeed; % Define the temporary variable tdel as the distance of the pos value from the focus [using the distance formula sqrt(x^2 + y^2 + z^2)] divided by the speed of sound in the previously defined medium
        tdelay(ni) = tdel; % Set the element in the matrix tdelay defined by the iteration as tdel
    end

    tdelay = max(tdelay(:)) - tdelay; % Scale the tdelay by subtracting each value in the tdelay matrix from the max value of the tdelay matrix 

    % Iterate over the length of the xdcr_array using the iterator ni
    for ni=1:length(xdcr_array(:))
        xdcr_array(ni).time_delay = tdelay(ni); % Set the time_delay value for each ni element in xdcr_array to the corresponding value in the time_delay array
    end
    % Define a matrix of time delays - tdelay - using the FOCUS function
    % get_time_delays from the redefined xdcr_array
    tdelay = get_time_delays(xdcr_array);

    % Open figure and clear it
    figure();clf;

    % Set the current figure color to white
    set(gcf,'color','w');

    % Display the time delays as a colormap on the xy plane (colormap gives
    % time delays in microseconds)
    imagesc(tdelay'*1e6);

    % Set the current axis value fontsize to 12
    set(gca,'fontsize',12);

    % Set the titles and x/y labels for the figure 
    title('Time Delay')
    xlabel('x');
    ylabel('y');
    colorbar

%     tdelay = max(tdelay(:)) - tdelay; % Scale the tdelay by subtracting each value in the tdelay matrix from the max value of the tdelay matrix 
% 
%     % Iterate over the length of the xdcr_array using the iterator ni
%     for ni=1:length(xdcr_array(:))
%         xdcr_array(ni).time_delay = tdelay(ni); % Set the time_delay value for each ni element in xdcr_array to the corresponding value in the time_delay array
%     end
    
    fprintf('f-number: %.2f x %.2f\n',focus_z/x_len,focus_z/y_len); % print in the command window the f-number in both directions by dividing the z focus (imaging depth) by the length of the aperture in each direction (x_len and height) 
    %% Calculate the Pressure Field

    tic(); % Start a time counter 
    disp('Calculating pressure field...'); % Display the string in the command window

    % Use the Fast Nearfield Method with time sequence decomposition (fsn_tmd) 
    % to find the pressure of the transducer array using the coordinate grid,
    % time structure, number of divisions, and excitation function defined
    % above and define it as the variable p_tsd
    p_tsd = transient_pressure(xdcr_array, coord_grid, medium, time_struct, ndiv, input_func);

    % Display the simulation time using toc() and the strings below
    disp(['Simulation complete in ', num2str(toc()), ' seconds.'])

    p_size = size(p_tsd); % Determine the variable p_size as the size of the transient pressure array defined above
    max_p0 = max(p_tsd(:)); % Define the variable max_p0 as the maximum of the transient pressure array defined above

    % Scale the pressure value in different methods depending on the plane of
    % simulation

    % Define a variable Tp as the total number of cycles divided by the
    % center frequency of the array (i.e., the time in transmitting the
    % beam the number of cycles speccified
    Tp = ncycles/f0;

    imp_z = medium.density * medium.soundspeed; % Define a variable imp_z as the acoustic impedance (in kg/m^3 * m/s = kg/m^2/s)

    % Define a variable I_avg as the average intensity of the ultrasound
    % beam by taking the sum of the squared complex magnitude of the 4D
    % pressure matrix (p_tsd)
    % along the time dimension divided by the impedance (imp_z) (multiplied by 2 because
    % the wave has to travel in both directions) multiplied by the time step
    % (dt) and divided by the Total range of times (T) (to give an
    % average)
    I_avg = sum(abs(p_tsd).^2,4)/(2*imp_z)*dt/T; %(kg/m/s^2)^2/(kg/m^2/s)*s = kg/s^2 

    % Define the variable pwr_avg as the average power by taking the sum of
    % all intensity values and multiplying by the step size in x and y and
    % by the ratio of the time range in the simulation divided by the time
    % for the transducer to complete its cycling (Tp)
    pwr_avg = sum(I_avg(:))*dx*dy*T/Tp; % kg/s^2 * m * m = kg*m^2/s^2 = J/s = W 

    % Define the variable elec_pwr as the electric power by dividing the 
    % pwr_avg by the conversion efficiency defined by the user
    elec_pwr = pwr_avg/conv_eff;  

    % Scale the pressure by finding the sqrt of the ratio of the input
    % electric power input defined by the user and the electric power
    % calculated from the simulation pressures
    p_scale = sqrt(elec_pwr_in/(elec_pwr)); % times 2 for 50 Ohm impedance match
    p_tsd = p_tsd * p_scale; % Re-define p_tsd with scaling
    pnp = max(-p_tsd(:));
    disp(strcat('Peak Negative Pressure is: ', num2str(pnp*1e-6), 'MPa'));

    %% Find the envelope of the pressure wave

    % Define array p_env as the envelope of the pressure file by taking the
    % Hilbert transform on the reshaped array of pressure values (scaled as
    % above)
    p_env = abs(hilbert(reshape(p_tsd,[prod(p_size(1:3)),p_size(4)])'));
    p_env = reshape(p_env',p_size); % Re-define p_env as the p_env array reshaped to the original array size
    [max_p, maxidx] = max(abs(p_tsd(:))); % Define variable max_p as the maximum pressure in the p_tsd array and maxidx as its index 
    
    [indx, indy, indz, indt] = ind2sub(size(p_tsd),maxidx); % Define the variables indx, indy, indz, indt as the subscript equivalents of the maxidx index in the p_tsd array

    nt = size(p_tsd, 4); % Define the variable nt as the number of time points in the p_tsd array's fourth dimension
    %% Display the main and side lobes of the output pressure wave

    % If the plane is in the xy direction, then execute the following
    if(strcmpi(plane,'xy'))
        figure();clf; % Open figure and clear it
        subplot(121); % Open subplot 121
        x_p = p_env(:, indy, indz, indt); % Define the vector x_p as all p_env (envelope pressure) values at all x values and at the index for the maximum value in the y, z, and t directions
        plot(x*100,mag2db(x_p/max(x_p)),'linewidth',2); % Plot the scaled x_p pressure values converted to dB and scaled by the maximum pressure against the x values multiplied by 100 with an increased linewidth
        ylim([-60,5]); % Set the y limit as -60, 5
        % Set the title, xlabel, and ylabel of the plot
        title('Beam Profile')
        xlabel('Azimuth (cm)');
        ylabel('dB');
        grid on; % turn gridlines on

        % Repeat above for the y (Elevation) direction
        subplot(122); % Open subplot 122
        y_p = p_env(indx,:,indz, indt);
        plot(y*100,mag2db(y_p/max(y_p)),'linewidth',2);
        grid on;
        ylim([-60,5]);
        title('Beam Profile')
        xlabel('Elevation (cm)');
        set(gcf,'color','white')
    end
    %% Map the Transmit Pulse

    % Define an array TxPulse as p_tsd array with all singleton dimension removed
    % at all the time values and at the spatial indices corresponding to the 
    % maximum pressure divided by 1e6 (to put it into Pa)
    TxPulse = squeeze(p_tsd(indx,indy,indz,:))/1e6;
    % Define the array TxPulseEnv as the envelope of the TxPulse array by
    % taking the hilbert transform and finding the complex magnitude of that
    % TxPulsse array
    TxPulseEnv = abs(hilbert(TxPulse));

    % Define an array n1 as the indices in the TxPulsEnv where TxPulseEnv is 
    % greater than or equal to the maximum of the array divided by 2 
    n1 = find(TxPulseEnv >= max(TxPulseEnv)/2);
    % Define the variable PulseFWHM as the FWHM of the envelope (vs. time) by determining the length of the vector of indices
    % where TxPulseEnv is greater than half the maximum
    PulseFWHM = length(n1)*dt;  %units in us;
    fprintf('Pulse width = %.3f (us)\n',PulseFWHM*1e6); % Print the pulse width as the PulseFWHM

    if strcmpi(plane,'xy')
        image_pl = squeeze(p_env(:,:,indz,indt))';

        lt_px = squeeze(p_tsd(:,indy,indz,indt))/1e6;
        lt_px_env = abs(hilbert(lt_px));

        lt_py = squeeze(p_tsd(indx,:,indz,indt))/1e6;
        lt_py_env = abs(hilbert(lt_py));

        ax_x = x*100;
        ax_y = y*100;
        lbl_x = 'Azimuth (cm)';
        lbl_y = 'Elevation (cm)';
        lbl_z = 'Pressure (MPa)';
        title_s = sprintf('Pressure at z=%.1f (cm)',focus_z*100);
        
        figure();clf;
        set(gcf,'color','w');
        [c,h] = contour(ax_x,ax_y,squeeze(p_env(:,:,:,indt))'/max_p,[1,1]/2);
        c_x = c(1,2:end);
        c_y = c(2,2:end);
        plot(c_x,c_y);
        axis equal;
        title('Focal size');
        xlabel(lbl_x);
        ylabel(lbl_y);

        xFWHM = max(c_x) - min(c_x);
        yFWHM = max(c_y) - min(c_y);
        FocalSize = polyarea(c_x,c_y);
        FocalGain = 20*log10(width*100*elem_x*height*100/FocalSize);

        fprintf('FWHM = %.2f (X) x %.2f (Y) mm\n',xFWHM*10,yFWHM*10);
        fprintf('Focal Size = %.3f (mm^2) \n',FocalSize*100);
        fprintf('Focal Gain = %.3f (dB) \n',FocalGain);
    elseif strcmpi(plane,'xz')
        image_pl = squeeze(p_env(:,indy,:,indt))';
        ax_x = x*100;
        ax_y = z*100;

        lt_px = squeeze(p_tsd(:,indy,indz,indt))/1e6;
        lt_px_env = abs(hilbert(lt_px));

        lt_py = squeeze(p_tsd(indx,indy,:,indt))/1e6;
        lt_py_env = abs(hilbert(lt_py));

        lbl_x = 'Azimuth (cm)';%'x (cm)';
        lbl_y = 'Depth (cm)';%'z (cm)';
        lbl_z = 'Pressure (MPa)';
        title_s = sprintf('Pressure at y=%.1f (mm)',y(indy)*100);

        figure();clf;
        set(gcf,'color','w');
        [c,h] = contour(ax_x,ax_y,squeeze(p_env(:,:,:,indt))'/max_p,[1,1]/2);
        c_x = c(1,2:end);
        c_y = c(2,2:end);
        plot(c_x,c_y);
        axis equal;
        title('Focal size');
        xlabel(lbl_x);
        ylabel(lbl_y);

        xFWHM = max(c_x) - min(c_x);
        zFWHM = max(c_y) - min(c_y);
        FocalSize = polyarea(c_x,c_y);
        FocalGain = 20*log10(width*100*elem_x*height*100/FocalSize);

        fprintf('FWHM = %.2f (X) x %.2f (Z) mm\n',xFWHM*10,zFWHM*10);
        fprintf('Focal Size = %.3f (mm^2) \n',FocalSize*100);
        fprintf('Focal Gain = %.3f (dB) \n',FocalGain);
    elseif strcmpi(plane,'yz')
        image_pl = squeeze(p_env(indx,:,:,indt))';
        ax_x = y*100;
        ax_y = z*100;

        lt_px = squeeze(p_tsd(indx,:,indz,indt))/1e6;
        lt_px_env = abs(hilbert(lt_px));

        lt_py = squeeze(p_tsd(indx,indy,:,indt))/1e6;
        lt_py_env = abs(hilbert(lt_py));

        lbl_x = 'Elevation (cm)';%'y (cm)';
        lbl_y = 'Depth (cm)';%'z (cm)';
        lbl_z = 'Pressure (MPa)';
        title_s = sprintf('Pressure at x=%.1f (mm)',x(indx)*1e3);
        
        figure();clf;
        set(gcf,'color','w');
        [c,h] = contour(ax_x,ax_y,squeeze(p_env(:,:,:,indt))'/max_p,[1,1]/2);
        c_x = c(1,2:end);
        c_y = c(2,2:end);
        plot(c_x,c_y);
        axis equal;
        title('Focal size');
        xlabel(lbl_x);
        ylabel(lbl_y);

        yFWHM = max(c_x) - min(c_x);
        zFWHM = max(c_y) - min(c_y);
        FocalSize = polyarea(c_x,c_y);
        FocalGain = 20*log10(width*100*elem_x*height*100/FocalSize);

        fprintf('FWHM = %.2f (Y) x %.2f (Z) mm\n',yFWHM*10,zFWHM*10);
        fprintf('Focal Size = %.3f (mm^2) \n',FocalSize*100);
        fprintf('Focal Gain = %.3f (dB) \n',FocalGain);        
    end

    figure();clf;
    set(gcf,'color','w');
    subplot(131);
    plot(t*1e6,TxPulse)
    xlabel('Time (\mus)');
    ylabel('Pressure (MPa)');
    % title(['Tx pressure: x=',num2str(x(indx)),'cm, y=',num2str(y(indy)),'cm']);
    % xlim([20,40]);

    subplot(132);
    plot(ax_x,lt_px_env);
    title('Beam Profile')
    xlabel(lbl_x);
    ylabel('Pressure (MPa)');
    xlim(round([min(ax_x),max(ax_x)]));

    subplot(133);
    plot(ax_y,lt_py_env);
    title('Beam Profile')
    xlabel(lbl_y);
    ylabel('Pressure (MPa)');
    xlim(round([min(ax_y),max(ax_y)]));

    figure();clf;
    set(gcf,'color','w');
    % colormap default;
    colormap(jet(256));
    % subplot(211);
    imagesc(ax_x,ax_y,image_pl/1e6);
    set(gca,'fontsize',12);
    xlabel(lbl_x);
    ylabel(lbl_y);colorbar;
    axis equal tight;
    % title(['Beam pattern in ',plane,' plane']);
    title(title_s);
    xlim(x_lim);
    ylim(y_lim);

    fprintf('\nMax pressure:\t%.4f (MPa)\n',max_p/1e6);
      
    Tprf = 1/prf;   %
    Tp = ncycles/f0;

    imp_z = medium.density * medium.soundspeed;
    I_avg = sum(abs(p_tsd).^2,4)/(2*imp_z)*dt/T;

    fprintf('Average Intensity @ focus = %.3f (W/cm^2) \n',I_avg(indx,indy)/100^2);

    pwr_avg = sum(I_avg(:))*dx*dy*T/Tp;        % average power
    elec_pwr = pwr_avg/conv_eff;

    fprintf('Average Acoustic Power = %.3f (W) \n',pwr_avg);
    fprintf('Average Electric Power = %.3f (W) \n',elec_pwr);

    t1 = 0:1/fs:Tp;

    s = sin(2*pi*input_func.f0*t1);
    s = s(:).*hann(length(t1));
    P_in = sum(s.^2)/(2*50)*dt/Tp;

    V = sqrt(elec_pwr/P_in);
    fprintf('Input voltage: \t%.2f V\n',V);
    
    %%
    %
    % figure();
    % for it = 1:nt,
    %     mesh(z*100, x*100, squeeze(p_tsd(:, :, :, it)))
    %     title(['FNM TSD Example, t = ', sprintf('%0.3f',t(it) * 1e6), '\mu','s'])
    %     zlabel('pressure (Pa)')
    %     xlabel('z (cm)')
    %     ylabel('x (cm)')
    %     temp=axis();
    %     temp(5)=-max_p;
    %     temp(6)=max_p;
    %     axis(temp);
    %     drawnow
    % end
end

PressureField = struct;
PressureField.freq = f0;
PressureField.Shape = Shape;
PressureField.units = {'MPa','mm','dB'};
PressureField.kerf_x = kerf_x*1e3;
PressureField.r_curve = r_curv*1e3;
PressureField.focus = focus_z*1e3;
PressureField.pk_pressure = max_p/1e6;
if strcmpi(plane,'xy')
   PressureField.az_fwhm = xFWHM*10;
   PressureField.el_fwhm = yFWHM*10;
   PressureField.focal_size = FocalSize*100;
   PressureField.focal_gain = FocalGain;
end

disp(PressureField);
%%

d_path=fullfile(pwd,'Mat_Files');
if (~isdir(d_path))
    mkdir(d_path);
end
infmt = 'yyyy-MM-dd_HH-mm-ss';
date_time = datetime(datetime,'Format',infmt);
fname = sprintf('%s_PressureField_%d_%s_%dx%d_R%d_F%d.mat',date_time,f0*1e-6,Shape,elem_x,elem_y,r_curv*1e3,focus_z*1e3);
save(fullfile(d_path,fname),'PressureField','-mat');

end