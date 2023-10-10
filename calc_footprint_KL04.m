function [PHI, PHIalong, PHIcross, Xcen, Ycen, Xrot, Yrot,xr,yr]=calc_footprint_KL04(angle, Csize, sigV, sigW, ustar, zmeas, z0, h, Psi)

%License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
%Flux footprint after Kljun et al. (2004), Metzger et al. (2012).
%Ported from the eddy4R family of R-packages (https://github.com/NEONScience/eddy4R) eddy4R.turb::def.foot.k04.R function to Matlab by Qindan Zhu and Eva Pfannerstill, UC Berkeley.
%Reference: Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P.: A simple parameterisation for flux footprint predictions, Boundary Layer Meteorol., 112, 503-523, doi:10.1023/B:BOUN.0000030653.71031.96, 2004.
%Reference: Metzger, S., Junkermann, W., Mauder, M., Beyrich, F., Butterbach-Bahl, K., Schmid, H. P., and Foken, T.: Eddy-covariance flux measurements with a weight-shift microlight aircraft, Atmos. Meas. Tech., 5, 1699-1717, doi:10.5194/amt-5-1699-2012, 2012.

%INPUT VARIABLES
  %angle	%wind direction [degrees] to rotate the inertial footprint matrix; used for both models
  %Csize	%cell size [m] of result grid; used for both models
  %hor		%horizontal wind speed [m/s]; used in KM01, and to determine z0 for K04
  %sigV		%crosswind fluctuations [m/s]; used for both models to calculate the crosswind dispersion
  %sigW		%vertical wind fluctuations [m/s]; used for K04 to characterize the vertical transport
  %ustar	%friction velocity [m/s]; used a) for both models to characterize the shear stress, and b) to determine z0 for K04
  %Lvirt	%Obukhov length from local friction velocity and buoyancy flux [m]; used a) for KM01 to characterize the vertical transport, and b) for the universal functions to determine z0 for K04
  %zmeas	%height of measurement - displacement [m]
  %z0       %roughness height (?)
  %h        %PBL Height
  %Psi		%set to 0

  %%%%%%PREPARATION OF INPUT PARAMETERS
  %constant parameters from pages 507 and 516
    alpha1 = 0.8;
    Ac = 4.28;	%+-0.13
    Ad = 1.68;	%+-0.11
    Af = 0.18;	%+-0.01
    Ax = 2.59;	%+-0.17
    Bup = 3.42;	%+-0.35

  %calculation of fitted parameters Eqs. (13) - (16)
    a = Af / (Bup - log(z0));	%maximum value of the distribution
    b = 3.70;	%+-0.30
    c = Ac * (Bup - log(z0));
    d = Ad * (Bup - log(z0));	%maximum upwind extend (non-dimensional)
  
   %%%%%%SIZE ESTIMATION OF FOOTPRINT MATRIX

  %scaling from non-dimensional to dimensional framework, Eqs. (17) - (18)
    scal = zmeas * (sigW / ustar)^(-alpha1);

  %alongwind (crosswind integrald) density distribution

    %position of maximum along x
      %non-dimensional, Eq. (11)
	Xmax = Ax * (Bup - log(z0));
      %dimensional, Eq. (17)
	xmax = Xmax * scal;

    %maximum value, Eqs. (12) - (13)
      fmax = a;
      
    %alongwind footprint extend

      %in lee (negative) direction, Eq. (10)
	whrxn = ceil((-(d * scal + Csize/2)) / Csize);

      %in luv (postitiv) direction until contribution falls below 1 % of fmax
	whri = xmax; whro = fmax;	%start from distribution peak
	while(whro > fmax / 100) 
	  whri = whri + Csize;	%use step width of landuse matrix
	  whro = FFPalong(whri);	%calculate
    end
	whrxp = ceil(whri / Csize);	%cell length necessay in X direction

    ymax = 0;
    fmax = FFPcross(whrxp * Csize, ymax, sigV, sigW, ustar, zmeas, z0, h, Psi);
    
    %crosswind footprint extend until contribution falls below 1 % fmax
     if fmax >= 0
          while whro > (fmax / 100)
              whri = whri + Csize;	%use step width of landuse matrix
              whro = FFPcross(whrxp * Csize, whri, sigV, sigW, ustar, zmeas, z0, h, Psi); 
          end
      elseif fmax <0 %EP added this part
          while whro < (fmax / 100)
              whri = whri + Csize;	%use step width of landuse matrix
              whro = FFPcross(whrxp * Csize, whri, sigV, sigW, ustar, zmeas, z0, h, Psi); 
          end
      end
     whry = ceil(whri / Csize);	%cell length necessay in Y direction 
    if whry <=2
        whry = 3; %otherwise too small a disctribution, code will stop
    end
     
     
     
     %%%%%%%CELL ALLOCATION AND INTEGRATION

  %place aircraft in cell center around zero

    %alongwind integration boundaries with aicraft centered in 0
    if(whrxn < 0) 
        XRng = [((whrxn:(-1)) + 0.5), (1:whrxp) - 0.5] * Csize;
    else 
        XRng = [0, (1:whrxp) - 0.5]* Csize;
    end
    
    %crosswind integration boundaries with aicraft centered in 0
    YRng = [0, 1:whry - 0.5] * Csize;
    
    %alongwind cell center coordinates

    Xcen = arrayfun(@(x) mean(XRng(x:(x+1))), 1:(length(XRng)-1));
    
    %crosswind cell center coordinates
    Ycen = [0, arrayfun(@(y) mean(YRng(y:(y+1))), 2:(length(YRng)-1))];

    %integration of alongwind footprint

    %function to integral over, Eq. (A10) - (A11)
      gam = @(t) t.^b .* exp(-t);

    %auxilary dimensionless distance
      Lhat = (XRng / scal + d) / c;

    %integral
      Gam = arrayfun(@(xwhr) integral(gam, b*Lhat(xwhr), b*Lhat(xwhr+1)), 1:(length(Lhat)-1));

    %cellwise alongwind footprint, Eq. (A10)
      PHIalong = a * c * exp(b) * b^(-b) / b *Gam;

    %integral over the entire footprint
      INTall = a * c * exp(b) * b^(-b) * gamma(b);

    %percentage of alongwind footprint covered
      cover = sum(PHIalong) / INTall * 100;

    %normalisation to unity
      PHIalong = PHIalong / sum(PHIalong);
      PHIalong = real(PHIalong); %only take real part
    
    %%%%%% #integration of crosswind footprint
    % #integration, output: top -> bottom == upwind -> downwind, left -> right == alongwind axis -> outside
    
    PHIcross = [];
    for i=1:length(Xcen)
        thiscross = FFPcrossXY(Xcen(i),YRng, sigV, sigW, ustar, zmeas, z0, h, Psi);
        PHIcross = [PHIcross; thiscross];
    end
 
    do_plot = false;
    if do_plot
        figure;
        subplot(2,1,1);
        yyaxis left;
        line(Xcen/1e3, PHIalong,'linewidth',2);
        ylabel('Covered Footprint');
        yyaxis right
        line(Xcen/1e3, cumsum(PHIalong),'linewidth',2);
        ylabel('Cum. Covered Footprint');
        title('Along Wind');
        set(gca,'linewidth',2);
        set(gca, 'fontsize',12);
        
        
    end
    
    %%%%%% COMBINE ALONG- AND CROSSWIND DENSITY DISTRIBUTIONS AND ROTATE INTO MEAN WIND
    %combine crosswind contributions on alongwind axis; will always yield uneven column number; rows sum to unity
    ncol = size(PHIcross, 2);
    nrow = size(PHIcross, 1);
    
    %PHIcross = [fliplr(PHIcross(:,2:ncol)), 2*PHIcross(:,1), PHIcross(:,2:ncol)];
    PHIcross = [fliplr(PHIcross(:,2:ncol)), PHIcross(:,1), PHIcross(:,2:ncol)];
    for i=1:size(PHIcross,1)
        PHIcross(i,:) = PHIcross(i,:)/sum(PHIcross(i,:));
    end
    PHIcross = real(PHIcross);
    

    %combination of along- and cross component; along wind from up to down, crosswind from left to right; sums to unity
    PHI = [];
    for x=1:nrow
        this_phi =  PHIalong(x) * PHIcross(x,:);
        PHI = real([PHI; this_phi]);
    end
    
    if do_plot
        figure;
        mesh(meshx/1e3, meshy/1e3, PHI);
        ylabel('Cross wind (km)');
        xlabel('Along wind (km)');
        h = colorbar;
        title('footprint (%)','fontsize',14);
        set(h,{'linew'},{2});
        set(gca,'linewidth',2);
        set(gca, 'fontsize',12);
    end
%     contour(meshx/1e3, meshy/1e3, PHI);
% if do_plot
%         figure;
%         newYcen = [fliplr(-Ycen), Ycen(2:end)];
%         [meshy, meshx] = meshgrid(newYcen, Xcen);
%         contour(meshx/1e3, meshy/1e3, PHI);
%         ylabel('Cross wind (km)');
%         xlabel('Along wind (km)');
%         h = colorbar;
%         title('footprint (%)','fontsize',14);
%         set(h,{'linew'},{2});
%         set(gca,'linewidth',2);
%         set(gca, 'fontsize',12);
% end

    Ycen = [fliplr(-Ycen), Ycen(2:end)];

%--------------------------------------------------------------------
% Derive footprint ellipsoid incorporating R% of the flux
% starting at peak value, if requested
%--------------------------------------------------------------------
r = 0.9;
   if ~isnan(r)
       rs = r;
   else
       if crop==1
           rs = 0.8;
       else
           rs = NaN;
       end
   end
   %dy = Csize;

   if ~isnan(rs)
       % Calculate integral of f_2d starting at peak value until R% are reached
       f_array = reshape(PHI,1,size(PHI,1)*size(PHI,2));
       f_sort  = sort(f_array,'descend');
       f_sort  = f_sort(~isnan(f_sort));
       f_cum   = cumsum(f_sort); %.*Csize.*dy;
       for i = 1:length(rs)
           f_diff     = abs(f_cum-rs(i));
           [~, ind_r] = min(f_diff);
           fr         = real(f_sort(ind_r));
           contour_r  = contourc(real(Ycen),real(Xcen),real(PHI),[fr fr]);
    
           % Contourc adds info on level and number of vertices - replace with NaN
           ind_nan = find(contour_r(1,:)==fr);
           contour_r(:,ind_nan) = NaN;
           
           % Decrease number of digits and sort/unique
           try
               % matlab release 2014b or newer only
               contour_r = round(contour_r,1);
           catch
               contour_r = round(10.*contour_r)./10;
           end
       
           if ~isnan(r)               
               % Fill output structure
               xr  = contour_r(1,:);
               yr  = contour_r(2,:);
           end
       end %for i
   end

    % start the rotation of the footprint based on the wind direction (angle)
    theta = 90-angle;
    R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
    
    [meshy, meshx] = meshgrid(Ycen, Xcen);
    
    Xrot = nan(size(meshx));
    Yrot = nan(size(meshy));
    for corner=1:numel(meshx)
        out = R * [meshx(corner); meshy(corner)];
        Xrot(corner) = out(1) ;
        Yrot(corner) = out(2) ;
    end
    
      if ~isnan(rs)
        [meshyr, meshxr] = meshgrid(yr, xr);
        
        Xrotr = nan(size(meshxr));
        Yrotr = nan(size(meshyr));
        for corner=1:numel(meshxr)
            out = R * [meshxr(corner); meshyr(corner)];
            Xrotr(corner) = out(1) ;
            Yrotr(corner) = out(2) ;
        end
    end

%     if do_plot
%         figure;
%         contour(Xrot/1e3,Yrot/1e3, PHI);
%         ylabel('Y (km)');
%         xlabel('X (km)');
%         h = colorbar;
%         title('rotated footprint (%)','fontsize',14);
%         set(h,{'linew'},{2});
%         set(gca,'linewidth',2);
%         set(gca, 'fontsize',12);
%     end
%--------------------------------------------------------------------
% Rotate footprint for plot in 2d
%--------------------------------------------------------------------
   if ~isnan(angle)
       wind_dir_rad = (angle).*pi./180;
       dist         = sqrt(meshx.^2 + meshy.^2);
       angl        = atan2(meshy, meshx);
       x_2d_rot     = dist.*sin(wind_dir_rad - angl);
       y_2d_rot     = dist.*cos(wind_dir_rad - angl);
       
       if ~isnan(r)
           for i = 1:length(r)
               dist      = sqrt(xr.^2 + yr.^2);
               angl     = atan2(xr, yr);
               % Fill output structure
               xr = dist.*sin(wind_dir_rad - angl);
               yr = dist.*cos(wind_dir_rad - angl);
           end 
       end
   end


    if do_plot
    figure;
    colormap jet;
    surf( (Xrot-1*Csize)/1e3, (Yrot-1*Csize)/1e3,PHI);shading flat;view(2);
    hold all;
    if ~isnan(rs)
            z = fr.*ones(size(yr)).*100;
            plot3(xr/1e3,yr/1e3,z,'r');
    end 
    xlabel('x [km]'); ylabel('y [km]');
    %axis([min(Xrot./1e3,'all') max(Xrot./1e3,'all') min(Yrot./1e3,'all') max(Yrot./1e3,'all')]);
    axis tight
    axis image
    set(gca,'TickDir','out');
    colorbar;
    end
% 
%     if do_plot
%         figure;
%         mesh(Yrot/1e3, Xrot/1e3, PHI);
%         ylabel('Northing (km)');
%         xlabel('Easting (km)');
%         h = colorbar;
%         title('footprint (%)','fontsize',14);
%         set(h,{'linew'},{2});
%         set(gca,'linewidth',2);
%         set(gca, 'fontsize',12);
%     end
    
    % build in functions
    %crosswind integrald flux footprint, Eq. (7)
    function Fstar = FFPalong(x)
        Xstar = ((sigW / ustar).^alpha1) * x / zmeas;
        %seperated term 
          Lprime = (Xstar + d) / c;
        %crosswind integrald flux footprint, Eq. (7)
          Fstar = a * (Lprime.^b) * exp(b * (1 - Lprime));
    end

    %crosswind distribution of footprint (Heidbach, 2010)
    function Dy = FFPcross(x, y, sigV, sigW, ustar, zmeas, z0, h, Psi)
        %describing the friction within the air layer / column [s]
          Tly = 0.08 * h^2 / (h - zmeas) / ustar;
        %column average transport velocity [s]
          Ulog = ustar / 0.4 * (log(zmeas / z0) - (zmeas - z0) / zmeas - Psi);
        %average travel time of a particle
          tau = sqrt((x / Ulog)^2 + ((zmeas - z0) / sigW)^2);
        %scaled crosswind fluctuations
          sigma = tau / (1 + sqrt(tau / (2 * Tly))) * tau / Tly * sigV;
        %crosswind distribution
          Dy = (1 / (sqrt(2 * pi) * sigma)) * exp((-y^2) / (2 * (sigma^2)));
    end

%     %%function for crosswind dispersion
%     function Dy = FFPcrossY(y, sigma)
%         Dy = (1 / (sqrt(2 * pi) * sigma)) * exp((-y^2) / (2 * (sigma^2)));
%     end

    %alongwind distance dependence of crosswind dispersion
    function PHIcross = FFPcrossXY(x, y, sigV, sigW, ustar, zmeas, z0, h, Psi )
        %describing the friction within the air layer / column [s]
          Tly = 0.08 * h^2 / (h - zmeas) / ustar;
        %column average transport velocity [s]
          Ulog = ustar / 0.4 * (log(zmeas / z0) - (zmeas - z0) / zmeas - Psi);
        %average travel time of a particle
          tau = sqrt((x / Ulog)^2 + ((zmeas - z0) / sigW)^2);
        %scaled crosswind fluctuations
          sigma = tau / (1 + sqrt(tau / (2 * Tly))) * tau / Tly * sigV;
        %%function for crosswind dispersion
          FFPcrossY = @(y) (1 / (sqrt(2 * pi) * sigma)) .* exp((-y.^2) ./ (2 * (sigma^2)));
        %call function for crosswind dispersion (integration slightly increases density towards the outside)
          %PHIcross = FFPcrossY(Ycen, sigma)
          PHIcross = arrayfun(@(ywhr) integral(FFPcrossY, y(ywhr), y(ywhr+1)), 1:(length(y)-1));
        %normalisation to 0.5
          PHIcross = PHIcross / (2 * sum(PHIcross));
    end
    
end


