classdef misc_flux_calculation
    properties(Constant)
        campaign_filepath = '/Users/monicazhu/Box/LA_Flux/Flux campaign files/';
        highway_filepath = '/Users/monicazhu/Box/LA_Flux/inputs/stanford-xc453kn9742-shapefile';
        input_filepath = '/Users/monicazhu/Box/LA_Flux/inputs/';
        
    end
    
    methods(Static)
        function assemble_flux_measurements(thisday_filepath)
            %example: thisday_filepath June15_SJVFresno_Tuesday_11am
            close all;
%             clear all;
            filepath = fullfile(misc_flux_calculation.campaign_filepath, thisday_filepath);
%             filepath = '/Users/monicazhu/Box/LA_Flux/Flux campaign files/June15_SJVFresno_Tuesday_11am';


            data = load(fullfile(filepath, 'NOx_ppb.mat'));

            %Convert time to seconds, I think?
            LIF_time = abs(rem(data.ppbday_ave,1)).*86400;
            dtime = data.ppbday_ave+datenum('2021-01-01')-1-7/24;

            %Load in all airplane/meteorology measurements
            Table = readtable(fullfile(filepath, 'met_aircraft_data.txt'));

            %Find array of plane_time measurements - all done at 0.1s
            if ismember('UTCTime_secs_',Table.Properties.VariableNames) == 1
                plane_time = Table.UTCTime_secs_;
            else
                plane_time = Table.UTCTime_HH_hhhhh_;
            end

            %Assign temporary plane measurements
            A = Table.Lat;
            B = Table.Long;
            C = Table.Vert_Wind_m_s_;
            D = Table.Tamb_C_;
            E = Table.Pitch_deg_;
            F = Table.Ps_mb_;
            G = Table.WindSpeed_m_s_;
            H = Table.WindDir_Deg_;
            I = Table.NovAtelAlt_m_;
            J = Table.TAS_m_s_;
            K = Table.Roll_deg_;
            L = Table.Pitch_deg_;

            %Filter out NaNs
            apple = isnan(plane_time);
            plane_time(apple)=[]; 
            A(apple) = [];
            B(apple) = [];
            C(apple) = [];
            D(apple) = [];
            E(apple) = [];
            F(apple) = [];
            G(apple) = [];
            H(apple) = [];
            I(apple) = [];
            J(apple) = [];
            K(apple) = [];
            L(apple) = [];

            %This code is to deal with flights that went over midnight UTC that reset
            %the time of day
            for k = 1:length(plane_time)
                if plane_time(k)<40000
                    plane_time(k) = plane_time(k) + 86400;
                end
            end

            for l = 1:length(LIF_time)
                if LIF_time(l)<40000
                    LIF_time(l) = LIF_time(l) + 86400;
                end
            end

            savedata.lat = misc_flux_calculation.interp_to_lif_time(plane_time,A,LIF_time);
            savedata.lon = misc_flux_calculation.interp_to_lif_time(plane_time,-B,LIF_time);
            savedata.v_wind = misc_flux_calculation.interp_to_lif_time(plane_time,C,LIF_time);
            savedata.temp = misc_flux_calculation.interp_to_lif_time(plane_time,D,LIF_time);
            savedata.pot_temp = misc_flux_calculation.interp_to_lif_time(plane_time,E,LIF_time);
            savedata.pressure = misc_flux_calculation.interp_to_lif_time(plane_time,F,LIF_time);
            savedata.windspeed = misc_flux_calculation.interp_to_lif_time(plane_time,G,LIF_time);
            savedata.wind_dir = misc_flux_calculation.interp_to_lif_time(plane_time,H,LIF_time);
            savedata.alt = misc_flux_calculation.interp_to_lif_time(plane_time,I,LIF_time);
            savedata.airspeed = misc_flux_calculation.interp_to_lif_time(plane_time,J,LIF_time);
            savedata.roll = misc_flux_calculation.interp_to_lif_time(plane_time,K,LIF_time);
            savedata.pitch = misc_flux_calculation.interp_to_lif_time(plane_time,L,LIF_time);
            savedata.NOmgm3=(transpose(data.aveNO).*32.*12.187)./(273.15+savedata.temp)./1000; %calculate mg/m3 from ppb assuming 1 atm
            savedata.NO2mgm3=(transpose(data.aveNO2).*46.*12.187)./(273.15+savedata.temp)./1000; %calculate mg/m3 from ppb assuming 1 atm
            savedata.LIF_time = LIF_time;
            no2 = data.aveNO2;
            nox = data.aveNO + data.aveNO2;
            
            
%             misc_flux_calculation.plot_flight_track(savedata.lon, savedata.lat, savedata.alt, no2, nox, savedata.LIF_time)
%             [savedata.no2_lag_time, lag_indx] = misc_flux_calculation.calc_lag_time(LIF_time, nox, savedata.alt, savedata.temp);
%             savedata.nox_lag_time = misc_flux_calculation.calc_lag_time(LIF_time, no2, savedata.alt, savedata.temp);
            
            save(fullfile(filepath, 'merge_obs'), 'savedata');
            clear all;
        end
        
        function make_flux_calculation(dayofmonth)
            close all
            f_filepath = misc_flux_calculation.campaign_filepath;
            h_filepath = dir(f_filepath);
            for i =1:numel(h_filepath)
                if contains(h_filepath(i).name, sprintf('June%d',dayofmonth))
                    filepath = fullfile(f_filepath, h_filepath(i).name);
                    break;
                end
            end
            data = load(fullfile(filepath, 'merge_obs.mat'));
            lon = data.savedata.lon;
            lat = data.savedata.lat;
            v_wind = data.savedata.v_wind;
            temp = data.savedata.temp;
            pot_temp = data.savedata.pot_temp;
            pressure = data.savedata.pressure;
            windspeed = data.savedata.windspeed;
            wind_dir = data.savedata.wind_dir;
            alt = data.savedata.alt;
            airspeed = data.savedata.airspeed;
            roll = data.savedata.roll;
            pitch = data.savedata.pitch;
            NOmgm3 = data.savedata.NOmgm3;
            NO2mgm3 = data.savedata.NO2mgm3;
            LIF_time = data.savedata.LIF_time;
            
            % add time filtering
            time_window = misc_flux_calculation.make_time_window(dayofmonth);
            
            % Remove spikes in roll and filter data for them

            %This is weird, but what I decide to do is say that whenever the NOx instrument is
            %running a calibration or zero (ie. NO and NO2 are NaN) say that the
            %airplane is doing a high roll/turn, which will get filtered out next step
            roll(isnan(NOmgm3)) = 50;
            roll(isnan(NO2mgm3)) = 50;

            %Set the highest degree of roll allowable and use nanwindow to filter out data
            %using 8 from previous Goldstein flux flights
            rollmax = 8; 
            nanwindow = 10;
            roll_spikes = find(abs(roll) >= rollmax);
            roll_spikes = roll_spikes+10;
            for i=-nanwindow:nanwindow
                if (roll_spikes)<=length(NO2mgm3)
                    NO2mgm3(roll_spikes+i) = NaN;
                    NOmgm3(roll_spikes+i) = NaN;
                    pitch(roll_spikes+i) = NaN;
                    v_wind(roll_spikes+i) = NaN;
                    temp(roll_spikes+i) = NaN;
                    airspeed(roll_spikes+i) = NaN;
                    pot_temp(roll_spikes+i) = NaN;
                    windspeed(roll_spikes+i) = NaN;
                    wind_dir(roll_spikes+i) = NaN;
                    alt(roll_spikes+i) = NaN;
                    pressure(roll_spikes+i) = NaN;
                    lat(roll_spikes+i) = NaN;
                    lon(roll_spikes+i) = NaN;
                end
            end 
            roll_spikes = roll_spikes-10;

            %Select flight segments between all roll spikes and zeros/calibrations
            n_spec = 1;
            segments = struct();
            
            for j = 1:length(roll_spikes)-1
            time_beg = LIF_time(roll_spikes(j));
            time_end = LIF_time(roll_spikes(j+1));
            
            %Calculate length and time of segment
            segmentlengthkm = ((time_end-time_beg)*mean(airspeed(roll_spikes(j):roll_spikes(j+1)),'omitnan'))/1000;
            
            % determine if it is with in the time window
            in_time_window = misc_flux_calculation.is_in_time_window(time_beg - LIF_time(1), time_end - LIF_time(1), time_window);

            %Need segment to be greater than 10km (or maybe even 15km) for good wavelet analysis otherwise ignore 
            if segmentlengthkm > 10 && in_time_window
                
                %Select segments of data    
                NO_seg=NOmgm3(roll_spikes(j)+1:roll_spikes(j+1)-1);  
                n_obs = numel(NO_seg);
                NO2_seg=NO2mgm3(roll_spikes(j)+1:roll_spikes(j+1)-1);
                time_seg=LIF_time(roll_spikes(j)+1:roll_spikes(j+1)-1);
                temp_seg=temp(roll_spikes(j)+1:roll_spikes(j+1)-1);
                vw = v_wind(roll_spikes(j)+1:roll_spikes(j+1)-1);
                lat_seg = lat(roll_spikes(j)+1:roll_spikes(j+1)-1);
                lon_seg = lon(roll_spikes(j)+1:roll_spikes(j+1)-1);
                pitch_seg = pitch(roll_spikes(j)+1:roll_spikes(j+1)-1);
                roll_seg = roll(roll_spikes(j)+1:roll_spikes(j+1)-1);
                airspeed_seg = airspeed(roll_spikes(j)+1:roll_spikes(j+1)-1);
                pot_temp_seg = pot_temp(roll_spikes(j)+1:roll_spikes(j+1)-1);
                windspeed_seg = windspeed(roll_spikes(j)+1:roll_spikes(j+1)-1);
                wind_dir_seg = wind_dir(roll_spikes(j)+1:roll_spikes(j+1)-1);
                alt_seg = alt(roll_spikes(j)+1:roll_spikes(j+1)-1);
                pressure_seg = pressure(roll_spikes(j)+1:roll_spikes(j+1)-1);
                enhance = (time_seg(2:n_obs) - time_seg(1:n_obs-1)).* airspeed_seg(1:n_obs-1);
                length_seg = cat(1,0,cumsum(enhance));
%                 
                %%median despiking of NO, NO2 and vertical wind following technique
                %%proposed in London flux paper
                NO2_seg = misc_flux_calculation.despike_lag(NO2_seg);
                NOx_seg = misc_flux_calculation.despike_lag(NO_seg + NO2_seg);
                vw = misc_flux_calculation.despike_lag(vw);
                
                temp_waves(:,1) = NOx_seg;
                temp_waves(:,2) = NO2_seg;
                n_obs = numel(NO2_seg);
                
                for i = 1:2
                    [flux, lagc] = xcov(temp_waves(:,i),vw,2000,'normalized');
                    center = ceil(length(flux)./2); %define center to find local maximum of covariance
                    peaks = find(max(flux) == flux) - center +1;
                    fprintf('this peaks diff: %d \n', peaks);
                    ajpeaks(i) = peaks;
                    if peaks >0
                        
                        waves_adjust_index = peaks:n_obs;
                        vw_adjust_index = 1:n_obs - peaks+1;
                        ajwaves{i} = temp_waves(waves_adjust_index,i);
                        ajtime{i} = time_seg(vw_adjust_index);

                    else
                        waves_adjust_index = 1:n_obs+peaks;
                        vw_adjust_index = -peaks+1:n_obs;
                        ajwaves{i} = temp_waves(waves_adjust_index,i);
                        ajtime{i} = time_seg(vw_adjust_index);
                    end
                end
                
                if all(abs(ajpeaks)<40)
                    fprintf('Meet the criteria\n');
                    [ctime, ctime_indx1, ctime_indx2] = intersect(ajtime{1}, ajtime{2});

                    waves(:,1) = ajwaves{1}(ctime_indx1);
                    waves(:,2) = ajwaves{2}(ctime_indx2);

                    [~, ~, ctime_indx2] = intersect(ctime, time_seg);
                    vw = vw(ctime_indx2);
                    segments(n_spec).waves = waves;
                    segments(n_spec).vw = vw;
                    segments(n_spec).alt = alt_seg(ctime_indx2);
                    segments(n_spec).lon = lon_seg(ctime_indx2);
                    segments(n_spec).lat = lat_seg(ctime_indx2);
                    segments(n_spec).airspeed = airspeed_seg(ctime_indx2);
                    segments(n_spec).windspeed = windspeed_seg(ctime_indx2);
                    segments(n_spec).wind_dir = wind_dir_seg(ctime_indx2);
                    segments(n_spec).temp = temp_seg(ctime_indx2);
                    segments(n_spec).time = ctime;
                    segments(n_spec).xkm = length_seg(ctime_indx2)/1e3;
                    
                for i=1:2
                    sst69 = waves(:,i); 
                    sst71 = vw;  
                    variance1 = std(sst69)^2;
                    sst69 = (sst69 - mean(sst69))/sqrt(variance1) ;
                    variance2 = std(sst71)^2;
                    sst71 = (sst71 - mean(sst71))/sqrt(variance2) ;
                    
                    %pbl=1; %1.5; %km 
%                     zi=mean(alt_seg(vw_adjust_index),'omitnan')/1000; %km
                    n = length(sst69);
                    dt=0.2; %Sampling interval NOx
%                     time = [0:length(sst69)-1]*dt + 0 ;  % construct time array
%                     xlim = [0,max(time)*xkm];  % plotting range
                    pad = 1;      % pad the time series with zeroes (recommended)
                    s0 = 0.4;    % 2 x sampling time
                    dj = 0.25; %octaves
                    %j1 = round(log2(round(size(sst69,1))/32))/dj;    % e.g. use log2 of half the sample size with dj sub-octaves each 64 = scaling factor to choose the right wavelet scales
                    j1 = round(log2(round(size(sst69,1))*dt/s0))/dj; 
                    lag1 = 0.72;  % lag-1 autocorrelation for red noise background
                    mother = 'Morlet';

                    % Wavelet transform:
                    [wave69,period69,scale69,coi69] = WAVELET(sst69,dt,pad,dj,s0,j1,mother);
                    [wave71,period71,scale71,coi71] = WAVELET(sst71,dt,pad,dj,s0,j1,mother);
                    %power = (abs(wave)).^2 ;        % compute wavelet power spectrum
                    spec(i).coivoc=coi69;
                    spec(i).scalevoc=scale69;

                    Cdelta = 0.776;   % this is for the MORLET wavelet

                    %crossspectrum
                    wc6971=real(wave69).*real(wave71)+imag(wave69).*imag(wave71);
                    spec(i).wavew=(variance1*variance2)*sum(wc6971,2)/n;
                    spec(i).fqwave=1./period69;
                    spec(i).crosspec = wc6971;

                    scale_avg = (scale69')*(ones(1,n));  % expand scale --> (J+1)x(N) array
                    scale_avg = wc6971./ scale_avg;   % [Eqn(24)]
                    scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg);   % [Eqn(24)]
                    spec(i).flux = scale_avg;

                end
                
                segments(n_spec).spec = spec;
                n_spec = n_spec + 1;
                end
            end
            temp_waves = [];
            waves = [];
            end
            fprintf('%d segments are found.\n', n_spec - 1);
            save(fullfile(filepath, 'wavelet_decomp'), 'segments');
        end
        
        function [var] = despike_lag(var)
%             figure;
%             hold on;
%             plot(var,'color','blue');
            med = movmedian(var,7);
            differ = var - med;
            [N, edges] = histcounts(differ,25);
            zs=find(N==0);
            [val,idx] = min(abs(zs-12));
            DT = abs(edges(zs(idx)));
            for i = 1:length(differ)
                if abs(differ(i)) > DT
                    var(i) = med(i);
                end
            end
%             plot(var,'color','red');
        end
        
        function time_window = make_time_window(dayofmonth)
            %SJV
            if dayofmonth == 3
                time_window = [4780, 11184; 13301, 13778; 14167, 14660;...
                    14700, 15158; 15221, 15583; 15730, 16147; 16229, 17562; 17635, 18475];
            end
            if dayofmonth == 8
                time_window = [6.8e3, 1.33e4; 1.5432e4, 1.5948e4; 1.6459e4, 1.6984e4; ...
                    1.7099e4, 1.7539e4;  1.7662e4, 1.8123e4; 1.82e4, 1.8625e4;...
                    1.8719e4, 1.9161e4; 1.9232e4, 2.0984e4; 2.1028e4, 2.1959;];
            end
            
            if dayofmonth == 16
                time_window = [6.288e3, 8.146e3; 8.447e3, 8.935e3; 8.999e3, 9.573e3;...
                    9598, 10155; 10262, 10761; 10809, 14300; 16513, 18633; 18710, 19639];
            end
            
            if dayofmonth == 9
                time_window = 1.0e+04 *[0.6747, 1.3344; 1.5775, 1.9752; 1.9878, 2.0636];
            end
            
            if dayofmonth == 13
                time_window = 1.0e+04 *[0.6578, 1.2491; 1.5085, 1.9718];
            end
            
            if dayofmonth == 15
                 time_window = 1.0e+04 *[0.6374, 1.1955; 1.4222, 1.8234];
            end
            
            function in_time_window = is_in_time_window(time_beg, time_end, time_window)
            in_time_window = false;
            for i = 1:size(time_window, 1)
                this_window = time_beg >= time_window(i, 1) & time_end <= time_window(i, 2);
                in_time_window = in_time_window | this_window;
            end
        end
            
        end
        
        function in_time_window = is_in_time_window(time_beg, time_end, time_window)
            in_time_window = false;
            for i = 1:size(time_window, 1)
                this_window = time_beg >= time_window(i, 1) & time_end <= time_window(i, 2);
                in_time_window = in_time_window | this_window;
            end
        end
        
        function q_int = interp_to_lif_time(planetime, q, liftime)
            x = planetime;
            [x, index] = unique(x); 
            q_int = interp1(x,q(index), liftime);
        end
        
    end
end