classdef misc_flux_calculation
    properties(Constant)
        highway_filepath = './inputs/highway_national';
        input_filepath = './inputs/';
        output_filepath = './output';
        figure_filepath = './output/figures';
        emission_filepath = './output/emissions';
        la_filepath = './Flux campaign files';
        hrrr_filepath = '../Database/HRRR';
        
        % emissions
        carb_filepath = '../Database/CARB';
        iowa_filepath = '../Database/Soil_IOWA/daily';
        five_filepath = '../Database/Katelyn/Anthropogenic_emissions';
        bm_filepath = '../Database/BM';
        dayofmonths = [3,8,9,13,15,16,22];

%         mLat = 34.5:0.005:37.5;
%         mLon = -122:0.005:-116;

        % 500 m spatial resolution
        mLat = 34.5:0.0045:37.5;
        mLon = -122:0.0055:-116;
        dXX = 0.0055/2;
        dYY = 0.0045/2;
    end
    
    methods(Static)
        %%%%MAKE FUCNTIONS%%%%
        
        % Combine Met and NOx
        function merge_nox_met(dayofmonth, do_plot)
            %Merge NOx measurement with meterological data and allign them
            %with the LIF time series.
            %OUTPUT : Merge_obs.mat
            
            close all
            f_filepath = misc_flux_calculation.la_filepath;
            h_filepath = dir(f_filepath);
            for i =1:numel(h_filepath)
                if contains(h_filepath(i).name, sprintf('June%d',dayofmonth))
                    filepath = fullfile(f_filepath, h_filepath(i).name);
                    break;
                end
            end

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
                plane_time = Table.UTCTime_swcs_;
            end

            %Assign temporary plane measurements
            A = Table.Lat;
            B = Table.Long;
            C = Table.Vert_Wind_m_s_;
            D = Table.Tamb_C_;
            E = Table.Theta_K_;
            F = Table.Ps_mb_;
            G = Table.WindSpeed_m_s_;
            H = Table.WindDir_Deg_;
            I = Table.NovAtelAlt_m_;
            J = Table.TAS_m_s_;
            K = Table.Roll_deg_;
            L = Table.Pitch_deg_;
            M = Table.RadAlt__m_;

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
            M(apple) = [];

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

            savedata.lat = misc_footprint_analysis.interp_to_lif_time(plane_time,A,LIF_time);
            savedata.lon = misc_footprint_analysis.interp_to_lif_time(plane_time,-B,LIF_time);
            savedata.v_wind = misc_footprint_analysis.interp_to_lif_time(plane_time,C,LIF_time);
            savedata.temp = misc_footprint_analysis.interp_to_lif_time(plane_time,D,LIF_time);
            savedata.pot_temp = misc_footprint_analysis.interp_to_lif_time(plane_time,E,LIF_time);
            savedata.pressure = misc_footprint_analysis.interp_to_lif_time(plane_time,F,LIF_time);
            savedata.windspeed = misc_footprint_analysis.interp_to_lif_time(plane_time,G,LIF_time);
            savedata.wind_dir = misc_footprint_analysis.interp_to_lif_time(plane_time,H,LIF_time);
            savedata.alt = misc_footprint_analysis.interp_to_lif_time(plane_time,I,LIF_time);
            savedata.airspeed = misc_footprint_analysis.interp_to_lif_time(plane_time,J,LIF_time);
            savedata.roll = misc_footprint_analysis.interp_to_lif_time(plane_time,K,LIF_time);
            savedata.pitch = misc_footprint_analysis.interp_to_lif_time(plane_time,L,LIF_time);
            savedata.radalt = misc_footprint_analysis.interp_to_lif_time(plane_time,M,LIF_time);
            savedata.galt = savedata.alt - savedata.radalt;
            

            savedata.NOmgm3=(transpose(data.aveNO).*30.*12.187)./(273.15+savedata.temp)./nanmedian(savedata.pressure); %calculate mg/m3 from ppb
            savedata.NO2mgm3=(transpose(data.aveNO2).*46.*12.187)./(273.15+savedata.temp)./nanmedian(savedata.pressure); %calculate mg/m3 from ppb
%             savedata.NOzmgm3=(transpose(data.aveNOz).*46.*12.187)./(273.15+savedata.temp)./1000; %calculate mg/m3 from ppb
            savedata.LIF_time = LIF_time;
           % disp(nanmedian(savedata.pressure));
            no2 = data.aveNO2;
            nox = data.aveNO + data.aveNO2;
            
            % pblh filtering
            pblh_filtering();

            line(savedata.LIF_time - savedata.LIF_time(1), savedata.NOmgm3);
            % plot the altitude vs time, manually select the time window
            if do_plot 
                indx = savedata.LIF_time<86400;
                LIF_time = savedata.LIF_time(indx);
                alt = savedata.alt(indx);
                time = LIF_time - LIF_time(1);
                
                figure;
                line(time, alt, 'linestyle', 'none', 'marker', '.');
                ylabel('Altitude (m)');
                xlabel('Time (s)');
            end
            save(fullfile(filepath, 'merge_obs'), 'savedata');
           
            clear all;

            function pblh_filtering()
                hrrr = load('./inputs/hrrr_inputs_regrid_final');
                hrrr = hrrr.hrrr;
                dayofmonths = misc_flux_calculation.dayofmonths;
                indx = dayofmonth == dayofmonths;
                mpblh = hrrr(indx).mpblh;
                hours = hrrr(indx).hours;
                
                for i=1:numel(savedata.lon)
                    this_hour = round(savedata.LIF_time(i)/86400*24);
                    hour_indx = this_hour == hours;
                    if ~isnan(savedata.lon(i)) && ~isnan(savedata.lat(i))
                        lon_indx = round((savedata.lon(i) - misc_flux_calculation.mLon(1))/0.0055);
                        lat_indx = round((savedata.lat(i) - misc_flux_calculation.mLat(1))/0.0045);
                        if lon_indx >= 1 && lon_indx <= numel(misc_flux_calculation.mLon) && ...
                                lat_indx >=1 && lat_indx <= numel(misc_flux_calculation.mLat)
                            this_pblh = mpblh(lon_indx, lat_indx, hour_indx);
                            if savedata.radalt(i) >= this_pblh
                                savedata.NOmgm3(i) = nan;
                                savedata.NO2mgm3(i) = nan;
                            end
                        end
                    end
                end

            end
        end
        
        % Make flux calculation
        function make_flux_calculation_final(dayofmonth, lag_corr_opt)
            % READ DATA INPUT
            close all
            f_filepath = misc_flux_calculation.la_filepath;
            h_filepath = dir(f_filepath);
            for i =1:numel(h_filepath)
                if contains(h_filepath(i).name, sprintf('June%d',dayofmonth))
                    filepath = fullfile(f_filepath, h_filepath(i).name);
                    break;
                end
            end
            data = load(fullfile(filepath, 'merge_obs.mat'));
            
            % Addtional steps to make data ready to use
            data_ready = preprocessing_obs(data);
            
            % Remove spikes in roll
            %Set the highest degree of roll allowable and use nanwindow to filter out data
            %using 8 from previous Goldstein flux flights
            rollmax = 8; 
            nanwindow = 10;
            roll_spikes = find(abs(data_ready.roll) >= rollmax);
            for i=-nanwindow:nanwindow
                if (roll_spikes+i)<=length(data_ready.NO2mgm3) & (roll_spikes+i)>0
                    data_ready.NO2mgm3(roll_spikes+i) = NaN;
                    data_ready.NOmgm3(roll_spikes+i) = NaN;
                end
            end 
            
            % add time filtering
            time_window = misc_flux_calculation.make_time_window(dayofmonth);
            time_window_forest = misc_flux_calculation.make_deposition_time_window(dayofmonth);
            time_window = [time_window; time_window_forest];

            % find all break points, which in between are possible segments
            roll_spikes = find(isnan(data_ready.NO2mgm3) | isnan(data_ready.NOmgm3));
            
            % Select flight segment
            %Select flight segments between all roll spikes and zeros/calibrations
            n_spec = 1;
            segments = struct();
            
            for j = 1:length(roll_spikes)-1
                
                time_beg = data_ready.LIF_time(roll_spikes(j));
                time_end = data_ready.LIF_time(roll_spikes(j+1));

                %Calculate length and time of segment
                segmentlengthkm = ((time_end-time_beg)*mean(data_ready.airspeed(roll_spikes(j):roll_spikes(j+1)),'omitnan'))/1000;

                % determine if it is with in the time window
                in_time_window = misc_flux_calculation.is_in_time_window(time_beg - data_ready.LIF_time(1), time_end - data_ready.LIF_time(1), time_window);

                %Need segment to be greater than 10km (or maybe even 15km) for good wavelet analysis otherwise ignore 
                if segmentlengthkm > 10 && in_time_window
                    fprintf('segmentlengh:%f \n ',segmentlengthkm);
                    %Select segments of data
                    seg_indx = roll_spikes(j)+1:roll_spikes(j+1)-1;
                    seg = select_seg(seg_indx, data_ready);
                    
                    %lag_corr_opt = 'const_lag';
                    switch lag_corr_opt
                        case 'seg_lag'
                            seg_corr = corr_lag_seg(seg);
                        case 'const_lag'
                            seg_corr = corr_lag_seg_const(seg);
                    end
                    
                    for i=1:size(seg_corr.waves, 2)
                        sst69 = seg_corr.waves(:,i); 
                        sst71 = seg_corr.vw;  
                        variance1 = std(sst69)^2;
                        sst69 = (sst69 - mean(sst69))/sqrt(variance1) ;
                        variance2 = std(sst71)^2;
                        sst71 = (sst71 - mean(sst71))/sqrt(variance2) ;
                        n = length(sst69);
                        dt=0.2; %Sampling interval NOx
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
                        spec(i).aflux = misc_flux_calculation.make_moveavg(scale_avg, seg_corr.airspeed);
                        
                        %lod (limit of detection) calculation
                        spec(i).lod = calc_lod_seg(seg_corr, i);
                        spec(i).std = zeros(size(spec(i).flux)) + spec(i).lod;
                    end
                    
                    if size(seg_corr.waves, 2)
                        seg_corr.spec = spec;
                        % copy seg_sorr to segments(n_spec)
                        fds = fieldnames(seg_corr);
                        for i_fd = 1:numel(fds)
                            segments(n_spec).(fds{i_fd}) = seg_corr.(fds{i_fd});
                        end
                        n_spec = n_spec + 1;
                    end
                end

            end

            fprintf('%d segments are found.\n', n_spec - 1);
            switch lag_corr_opt
                case 'seg_lag'
                    save(fullfile(filepath, 'wavelet_decomp_seg_lag_final'), 'segments');
                case 'const_lag'
                    save(fullfile(filepath, 'wavelet_decomp_const_lag_final'), 'segments');
            end
            
               
            function merge = preprocessing_obs(data)
                % 1. trim data if LIF_time > 86400, which corresponds the measurements from prevous day
                LIF_time = data.savedata.LIF_time;
                indx = LIF_time < 86400;
                merge.LIF_time = LIF_time(indx);
                merge.lon = data.savedata.lon(indx);
                merge.lat = data.savedata.lat(indx);
                merge.v_wind = data.savedata.v_wind(indx);
                merge.temp = data.savedata.temp(indx);
                merge.pot_temp = data.savedata.pot_temp(indx);
                merge.pressure = data.savedata.pressure(indx);
                merge.windspeed = data.savedata.windspeed(indx);
                merge.wind_dir = data.savedata.wind_dir(indx);
                % Note: Altitude above sea level - surface altitude =
                % Altitude above ground level.
                merge.alt = data.savedata.alt(indx)- data.savedata.galt(indx);
                merge.airspeed = data.savedata.airspeed(indx);
                merge.roll = data.savedata.roll(indx);
                merge.pitch = data.savedata.pitch(indx);
                merge.NOmgm3 = data.savedata.NOmgm3(indx);
                merge.NO2mgm3 = data.savedata.NO2mgm3(indx);
            end
            
            function seg = select_seg(seg_indx, data_ready)
                seg.NO = data_ready.NOmgm3(seg_indx);  
                seg.NO2 = data_ready.NO2mgm3(seg_indx);
                seg.time = data_ready.LIF_time(seg_indx);
                seg.temp = data_ready.temp(seg_indx);
                seg.vw = data_ready.v_wind(seg_indx);
                seg.lat = data_ready.lat(seg_indx);
                seg.lon = data_ready.lon(seg_indx);
                seg.pitch = data_ready.pitch(seg_indx);
                seg.roll = data_ready.roll(seg_indx);
                seg.airspeed = data_ready.airspeed(seg_indx);
                seg.pot_temp = data_ready.pot_temp(seg_indx);
                seg.windspeed = data_ready.windspeed(seg_indx);
                seg.wind_dir = data_ready.wind_dir(seg_indx);
                seg.alt = data_ready.alt(seg_indx);
                seg.pressure = data_ready.pressure(seg_indx);
                %%median despiking of NO, NO2 and vertical wind following technique
                    %%proposed in London flux paper
                seg.NO2 = misc_flux_calculation.despike_lag(seg.NO2);
                seg.NOx = misc_flux_calculation.despike_lag(seg.NO + seg.NO2);
                seg.vw = misc_flux_calculation.despike_lag(seg.vw);
            end
            
            function seg_corr = corr_lag_seg(seg)
                this_waves = nan(size(seg.NOx,1),3);
                this_waves(:,1) = seg.NOx;
                this_waves(:,2) = seg.NO;
                this_waves(:,3) = seg.pot_temp;

                n_obs = numel(seg.NOx);

                for i = 1:3
                    [flux, lagc] = xcov(this_waves(:,i),seg.vw,2000,'normalized');
                    center = ceil(length(flux)./2); %define center to find local maximum of covariance
                    peaks = find(max((flux)) == (flux)) - center +1;
                    fprintf('this peaks diff: %d \n', peaks);
                    ajpeaks(i) = peaks;
                    if peaks >0
                        waves_adjust_index = peaks:n_obs;
                        vw_adjust_index = 1:n_obs - peaks+1;
                        ajwaves{i} = this_waves(waves_adjust_index,i);
                        ajtime{i} = seg.time(vw_adjust_index);

                    else
                        waves_adjust_index = 1:n_obs+peaks;
                        vw_adjust_index = -peaks+1:n_obs;
                        ajwaves{i} = this_waves(waves_adjust_index,i);
                        ajtime{i} = seg.time(vw_adjust_index);
                    end
                end

                if all(abs(ajpeaks(1:2))<40)
                    fprintf('Meet the criteria\n');
                    [ctime, ctime_indx1, ctime_indx2] = intersect(ajtime{1}, ajtime{2});

                    waves(:,1) = ajwaves{1}(ctime_indx1);
                    waves(:,2) = ajwaves{2}(ctime_indx2);

                    [~, ~, ctime_indx] = intersect(ctime, seg.time);
                    waves(:,3) = seg.pot_temp(ctime_indx);
                    
                    seg_corr.waves = waves;
                    seg_corr.vw =  seg.vw(ctime_indx);
                    seg_corr.alt = seg.alt(ctime_indx);
                    seg_corr.lon = seg.lon(ctime_indx);
                    seg_corr.lat = seg.lat(ctime_indx);
                    seg_corr.airspeed = seg.airspeed(ctime_indx);
                    seg_corr.windspeed = seg.windspeed(ctime_indx);
                    seg_corr.wind_dir = seg.wind_dir(ctime_indx);
                    seg_corr.temp = seg.temp(ctime_indx);
                    seg_corr.time = ctime;
                    seg_corr.lag_peaks = ajpeaks;
%                     segments(n_spec).xkm = length_seg(ctime_indx2)/1e3;
                else
                    seg_corr.waves = [];
                    fprintf('Does not meet the criteria\n');
                end
            end
            
            function seg_corr = corr_lag_seg_const(seg)
                data = load(fullfile(filepath, 'wavelet_decomp_seg_lag_final'));
                orig_segments = data.segments;
                nox_peaks = zeros(numel(orig_segments), 3);
                n_seg = 1;
                for i = 1:numel(orig_segments)
                    if ~isempty(orig_segments(i).lag_peaks)
                        nox_peaks(n_seg,:) = orig_segments(i).lag_peaks;
                        n_seg = n_seg + 1;
                    end
                end
                const_peaks= median(nox_peaks(1:n_seg-1,:),1);
                
                this_waves = nan(size(seg.NOx,1),3);
                this_waves(:,1) = seg.NOx;
                this_waves(:,2) = seg.NO;
                this_waves(:,3) = seg.pot_temp;

                n_obs = numel(seg.NOx);
                for i = 1:3
                    peaks = round(const_peaks(i));
                    if peaks >0
                        waves_adjust_index = peaks:n_obs;
                        vw_adjust_index = 1:n_obs - peaks+1;
                        ajwaves{i} = this_waves(waves_adjust_index,i);
                        ajtime{i} = seg.time(vw_adjust_index);

                    else
                        waves_adjust_index = 1:n_obs+peaks;
                        vw_adjust_index = -peaks+1:n_obs;
                        ajwaves{i} = this_waves(waves_adjust_index,i);
                        ajtime{i} = seg.time(vw_adjust_index);
                    end
                end
    
                [ctime, ctime_indx1, ctime_indx2] = intersect(ajtime{1}, ajtime{2});

                waves(:,1) = ajwaves{1}(ctime_indx1);
                waves(:,2) = ajwaves{2}(ctime_indx2);

                [~, ~, ctime_indx] = intersect(ctime, seg.time);
                waves(:,3) = seg.pot_temp(ctime_indx);

                seg_corr.waves = waves;
                seg_corr.vw =  seg.vw(ctime_indx);
                seg_corr.alt = seg.alt(ctime_indx);
                seg_corr.lon = seg.lon(ctime_indx);
                seg_corr.lat = seg.lat(ctime_indx);
                seg_corr.airspeed = seg.airspeed(ctime_indx);
                seg_corr.windspeed = seg.windspeed(ctime_indx);
                seg_corr.wind_dir = seg.wind_dir(ctime_indx);
                seg_corr.temp = seg.temp(ctime_indx);
                seg_corr.time = ctime;
                seg_corr.lag_peaks = const_peaks;
    %                     segments(n_spec).xkm = length_seg(ctime_indx2)/1e3;
                
            end
            
            function lod = calc_lod_seg(seg, i_wave)
                    delaylod = 1:5;
                    this_aflux = [];
                    for j=1:numel(delaylod)
                        this_waves = seg.waves;
                        waves_adjust_index = randperm(size(this_waves,1));
                        vw_adjust_index = 1:size(this_waves,1);
                        this_waves = this_waves(waves_adjust_index,i_wave);
                        this_vw = seg.vw(vw_adjust_index);
                        sst69 = this_waves; 
                        sst71 = this_vw;  
                        variance1 = std(sst69)^2;
                        sst69 = (sst69 - mean(sst69))/sqrt(variance1) ;
                        variance2 = std(sst71)^2;
                        sst71 = (sst71 - mean(sst71))/sqrt(variance2) ;
                        n = length(sst69);
                        dt=0.2; %Sampling interval NOx
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
                        Cdelta = 0.776;   % this is for the MORLET wavelet
                        wc6971=real(wave69).*real(wave71)+imag(wave69).*imag(wave71);
                        scale_avg = (scale69')*(ones(1,n));  % expand scale --> (J+1)x(N) array
                        scale_avg = wc6971./ scale_avg;   % [Eqn(24)]
                        scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg);   % [Eqn(24)]
                        
                        this_aflux = [this_aflux, std(misc_flux_calculation.make_moveavg(scale_avg, seg_corr.airspeed))];
                    end
                    lod = mean(this_aflux);
            end
            
        end   
        
        function make_flux_calculation_assemble()
            dayofmonths= misc_flux_calculation.dayofmonths;
            for i = 1:numel(dayofmonths)
                misc_flux_calculation.merge_nox_met(dayofmonths(i), false);
                misc_flux_calculation.make_flux_calculation_final(dayofmonths(i), 'seg_lag');
                misc_flux_calculation.make_flux_calculation_final(dayofmonths(i), 'const_lag');
            end
        end
        
        % HRRR 
        function make_hrrr_inputs()
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';
            
            hrrr = struct();
            
            dayofmonths = misc_flux_calculation.dayofmonths;
            for i_day=1:numel(dayofmonths)
                dayofmonth = dayofmonths(i_day);
                hrrr_filename = sprintf('hrrr_assemble_%02d.nc', dayofmonth);
                hrrr_filename = fullfile(misc_flux_calculation.hrrr_filepath, hrrr_filename);

                data = ncinfo(hrrr_filename);
                pblh = ncread(data.Filename, 'hpbl');
                temp = ncread(data.Filename, 'temp');
                fricv = ncread(data.Filename, 'fricv');
                roughness = ncread(data.Filename, 'roughness');

                hrrr_lon = ncread(data.Filename, 'lon');
                hrrr_lat = ncread(data.Filename, 'lat');
                days = ncread(data.Filename, 'day');
                hours = ncread(data.Filename, 'hour');
                
                mpblh = nan(size(mlon,1),size(mlon,2),size(pblh,2));
                mtemp = nan(size(mlon,1),size(mlon,2),size(pblh,2));
                mfricv = nan(size(mlon,1),size(mlon,2),size(pblh,2));
                mrs = nan(size(mlon,1),size(mlon,2),size(pblh,2));
                
                for i=1:size(pblh,2)
                    mpblh(:,:,i) = regrid_vars(pblh);
                    mtemp(:,:,i) = regrid_vars(temp);
                    mfricv(:,:,i) = regrid_vars(fricv);
                    mrs(:,:,i) = regrid_vars(roughness);
                end
                
                hrrr(i_day).days = days;
                hrrr(i_day).hours = hours;
                hrrr(i_day).mrs = mrs;
                hrrr(i_day).mfricv = mfricv;
                hrrr(i_day).mtemp = mtemp;
                hrrr(i_day).mpblh = mpblh;
                
            end
            
            save('./inputs/hrrr_inputs_regrid_final','hrrr');
            
            function this_var = regrid_vars(var)
                this_var = griddata(double(hrrr_lon), double(hrrr_lat), double(var(:,i)), mlon(:), mlat(:));
                
                this_var = reshape(this_var, size(mlon));
            end
            
        end
        
        % Land cover
        function make_crops_landcover_ing()
            %%% =======================================================================

            %%% CropScape (from here: "https://nassgeodata.gmu.edu/CropScape/")

            %%% =======================================================================



            %%% Get the crops
            input_filepath = misc_flux_calculation.input_filepath;
            overwrite = true;
            if (exist(fullfile(input_filepath,'CropData-500m.mat'), 'file') == 2) & ~overwrite

                load(fullfile(input_filepath,'CropData-500m.mat'))

            else

                [crops,~,R] = geotiffread('/Users/monicazhu/Dropbox/LA_Flux/inputs/CropScape_CA_2018.tif');

                crops = flipud(crops)'; % Transform to the same type as other plots

                CROPlat = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1))';

                CROPlon = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2))';

                % Make the crop IDs (from here: "https://www.nass.usda.gov/Research_and_Science/Cropland/metadata/metadata_ca18.htm")

                [~,txt,raw] = xlsread('/Users/monicazhu/Dropbox/LA_Flux/inputs/CropScapeAttributes.xls');

                CropNames = txt(3:end,1);

                CropIDs   = cell2mat(raw(3:end,2));

                CropType  = txt(3:end,3);



                %%% Build the native resolution crop fraction

                allCrops = nan(size(crops));

                for i = 1:length(CropIDs)

                    inds = crops == CropIDs(i);

                    if strcmp(CropType(i),'CROPS')

                        allCrops(inds) = 1;

                    else

                        allCrops(inds) = 0;

                    end

                end
                

                % Other land types
                allDelpOpen                = zeros(size(crops));

                allDelpOpen(crops == 121)  = 1;
                
                allDelpLow                 = zeros(size(crops));

                allDelpLow(crops == 122)   = 1;
                
                allDelpMed                 = zeros(size(crops));

                allDelpMed(crops == 123)   = 1;
                
                allDelpHig                 = zeros(size(crops));

                allDelpHig(crops == 124)   = 1;
                
                allGrass                   = zeros(size(crops));

                allGrass(crops == 176)     = 1;

                allMixed                   = zeros(size(crops));

                allMixed(crops == 143)     = 1;

                allEvergreen               = zeros(size(crops));

                allEvergreen(crops == 142) = 1;

                allDeciduous               = zeros(size(crops));

                allDeciduous(crops == 141) = 1;

                allShrubland               = zeros(size(crops));

                allShrubland(crops == 152)    = 1;
                
                mLat = misc_flux_calculation.mLat;
                mLon = misc_flux_calculation.mLon;
                %%% Regrid the crops to the MODIS grid

                dXX        = abs(mLon(2) - mLon(1))/2;

                dYY        = abs(mLat(2) - mLat(1))/2;

                cropFrac   = nan(length(mLon),length(mLat));

                grassFrac  = nan(length(mLon),length(mLat));

                mixedFrac  = nan(length(mLon),length(mLat));

                evergFrac  = nan(length(mLon),length(mLat));

                decidFrac  = nan(length(mLon),length(mLat));

                shrubFrac = nan(length(mLon),length(mLat));
                
                delpopenFrac = nan(length(mLon),length(mLat));
                
                delplowFrac = nan(length(mLon),length(mLat));
                
                delpmedFrac = nan(length(mLon),length(mLat));
                
                delphighFrac = nan(length(mLon),length(mLat));
                
                for i = 1:length(mLon)

                    iX = mLon(i)-dXX < CROPlon & CROPlon < mLon(i)+dXX;

                    for j = 1:length(mLat)

                        iY    = mLat(j)-dYY < CROPlat & CROPlat < mLat(j)+dYY;
                        
                        % Develop Openland

                        vals  = allDelpOpen(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            delpopenFrac(i,j) = nansum(vals(:))/nVals;

                        end
                        
                        % Develop Low 
                        
                        vals  = allDelpLow(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            delplowFrac(i,j) = nansum(vals(:))/nVals;

                        end
                        
                        % Develop Med
                        
                        vals  = allDelpMed(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            delpmedFrac(i,j) = nansum(vals(:))/nVals;

                        end
                        
                        % Develop High 
                        vals  = allDelpHig(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            delphighFrac(i,j) = nansum(vals(:))/nVals;

                        end
                        
                        % Crops

                        vals  = allCrops(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            cropFrac(i,j) = nansum(vals(:))/nVals;

                        end

                        % Grassland

                        vals  = allGrass(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            grassFrac(i,j) = nansum(vals(:))/nVals;

                        end

                        % Mixed

                        vals  = allMixed(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            mixedFrac(i,j) = nansum(vals(:))/nVals;

                        end

                        % Evergreen

                        vals  = allEvergreen(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            evergFrac(i,j) = nansum(vals(:))/nVals;

                        end

                        % Deciduous

                        vals  = allDeciduous(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            decidFrac(i,j) = nansum(vals(:))/nVals;

                        end

                        % Shrublands

                        vals  = allShrubland(iX,iY);

                        nVals = sum(~isnan(vals(:)));

                        if nVals > 0

                            shrubFrac(i,j) = nansum(vals(:))/nVals;

                        end

                    end

                end

                save('./inputs/CropData-500m.mat','-v7.3','crops','CROPlat','CROPlon','CropNames','CropIDs','CropType',...
                    'allCrops','cropFrac','grassFrac','mixedFrac','evergFrac','decidFrac','shrubFrac',...
                    'delpopenFrac','delplowFrac','delpmedFrac','delphighFrac');

            end
        end
        
        function make_crops_landcover_sim()
            clear all;
            input_filepath = misc_flux_calculation.input_filepath;
            cdata = load(fullfile(input_filepath, 'CropData-500m.mat'));
            
            [mhw, mlon, mlat] = misc_flux_calculation.make_highway_landcover;

            shrubFrac = cdata.shrubFrac;
            cropsFrac = cdata.cropFrac;
            decidFrac = cdata.decidFrac;
            evergFrac = cdata.evergFrac;
            grassFrac = cdata.grassFrac;
            mixedFrac = cdata.mixedFrac;
            openFrac = cdata.delpopenFrac;
            lowFrac = cdata.delplowFrac;
            medFrac = cdata.delpmedFrac;
            highFrac = cdata.delphighFrac;
            urbFrac = openFrac + lowFrac + medFrac + highFrac;

            clc = nan(size(mlon));

            
            for i=1:size(mlon, 1)
                for j=1:size(mlon, 2)
                    % reconcile the classification
                    % 1 developed 2 cropsland 3 grassland 4 forest 5 barren
                    % 6 water
                    %cropscape
                    if urbFrac(i, j) > 0.5
                        clc(i, j) = 1;
                    elseif cropsFrac(i, j) > 0.5
                        clc(i, j) = 2;
                    elseif grassFrac(i, j) > 0.5
                        clc(i, j) = 3;
                    elseif decidFrac(i, j)+evergFrac(i, j)+mixedFrac(i, j)>0.5
                        clc(i, j) = 4;
                    elseif shrubFrac(i, j) > 0.5
                        clc(i, j) = 5;
                    end  
                   
                    % determine which points contains highway
                    if mhw(i,j) == 1
                        clc(i, j) = 0;
                    end
                end
            end
            
            save('./inputs/landcover-500m.mat','-v7.3', 'mlon','mlat','clc');
          
        end
        
        % Make spatial averaged flux over the Central Valley
        function make_merge_avg()
            % conduct spatial average to 500m
            % match to hrrr inputs
            % match to landcover inputs
            
            mLon = misc_flux_calculation.mLon;
            mLat = misc_flux_calculation.mLat;
            dXX = misc_flux_calculation.dXX;
            dYY = misc_flux_calculation.dYY;
            
            dayofmonths = misc_flux_calculation.dayofmonths;
            data = load('./inputs/hrrr_inputs_regrid_final');
            hrrr = data.hrrr;
            hrrr_fds = {'mrs';'mfricv';'mtemp';'mpblh'};
            
            data = load('./inputs/landcover-500m.mat');
            clc = data.clc;
            
            lc_fds = {'lc'; 'dayofmonth'};
            fds = {'lon';'lat';'vw';'alt';'airspeed';'windspeed';'wind_dir';'temp';'time';'aflux';'err'};
                
            merge = misc_flux_calculation.merge_sjv_segments();
            for i_seg = 1:numel(merge)
                segment = merge(i_seg);
                
                segment_avg(i_seg) = make_empty_struct_from_cell([fds; hrrr_fds; lc_fds]);
                
                i_day = segment.dayofmonth == dayofmonths;
                this_hrrr = hrrr(i_day);
                
                for i = 1:length(mLon)

                    iX = mLon(i)-dXX < segment.lon & segment.lon < mLon(i)+dXX;
                    if sum(iX)
                        for j = 1:length(mLat)
                            iY    = mLat(j)-dYY < segment.lat & segment.lat < mLat(j)+dYY;
                            if sum(iX & iY)

                                segment_avg(i_seg).lon = [segment_avg(i_seg).lon; mLon(i)];
                                segment_avg(i_seg).lat = [segment_avg(i_seg).lat; mLat(j)];
                                % regrid flux to 500 m
                                for i_fd = 3:numel(fds)
                                    if any(contains(fieldnames(segment), fds(i_fd)))
                                        this_fd_val = mean(segment.(fds{i_fd})(iX & iY));
                                    elseif strcmpi(fds{i_fd}, 'aflux')
                                        this_fd_val = mean(segment.spec(1).(fds{i_fd})(iX & iY));
                                    elseif strcmpi(fds{i_fd}, 'err')    
                                        this_fd_val = segment.spec(1).lod/sqrt(sum(iX & iY));
                                    end
                                    segment_avg(i_seg).(fds{i_fd}) = [segment_avg(i_seg).(fds{i_fd}); this_fd_val];
                                end

                                % match it to hrrr
                                for i_fd = 1:numel(hrrr_fds)
                                    this_hour = round(segment_avg(i_seg).time(end)/86400*24);
                                    hour_indx = this_hour == this_hrrr.hours;
                                    this_fd_val = this_hrrr.(hrrr_fds{i_fd})(i, j, hour_indx);
                                    segment_avg(i_seg).(hrrr_fds{i_fd}) = [segment_avg(i_seg).(hrrr_fds{i_fd}); this_fd_val];
                                end

                                % match it to clc (land cover)
                                this_fd_val = clc(i,j);
                                segment_avg(i_seg).lc = [segment_avg(i_seg).lc; this_fd_val];

                                %record dayofmonth
                                segment_avg(i_seg).dayofmonth = segment.dayofmonth;
                            end
                        end
                    end
                end  
            end
            
            save(fullfile(misc_flux_calculation.output_filepath,'recap-sjv'), 'segment_avg');
        end
        
        % Make footprints
        function mphi = make_footprints_per_seg(segment)
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';
            
            sigW = std(segment.vw);
            sigU = std(segment.windspeed);
            
            mphi = nan(numel(segment.lon), size(mlon, 1), size(mlon,2));
            
            for i = 1:numel(segment.lon)
                % calculate footprint
                wd = segment.wind_dir(i);
                alt = segment.alt(i);
                rs = segment.mrs(i);
                pblh = segment.mpblh(i);
                ustar = segment.windspeed(i).*0.41./log(alt./rs);
                try
                    [phi, xrot, yrot, ~, ~, ~, ~] = calc_footprint_KL04(wd, 1000, sigU, sigW, ustar, alt, rs, pblh, 0);
                    xrot = segment.lon(i) + xrot/91e3;
                    yrot = segment.lat(i) + yrot/111e3;
                    mphi(i,:,:) = reshape(griddata(xrot, yrot, abs(phi), mlon(:), mlat(:)), size(mlon));
                catch err
                    this_mphi = nan(size(mlon));
                    indx = mlon == segment.lon(i) & mlat == segment.lat(i);
                    this_mphi(indx) = 1;
                    mphi(i,:,:) = this_mphi;
                end

            end
        end
        
        % I decided not the save the footprint since the output is very
        % large, ignore 'recap-sjv-fp', use 'recap-sjv'
        function make_footprints()

            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv'));
            segs = data.segment_avg;
            
            for i=1:numel(segs)
                try
                    mphi = misc_flux_calculation.make_footprints_per_seg(segs(i));
                catch err
                    mphi = nan;
                end
                segs(i).mphi = mphi;
                save(fullfile(misc_flux_calculation.output_filepath,'recap-sjv-fp'), '-v7.3', 'segs');
            end
            
            save(fullfile(misc_flux_calculation.output_filepath,'recap-sjv-fp'), '-v7.3','segs');
        end
        
        % Make emission evaluation
        function carb = read_carb_emissions_perday(dayofmonth)
            %filepath = '/Users/monicazhu/Documents/Database/CARB';
            filepath = misc_flux_calculation.carb_filepath;
            dayofyear = datenum('2021-06-01') + dayofmonth - 1 - datenum('2021-01-01') +1;
            filepattern = sprintf('st_4k_%d', dayofyear);
            filedir = dir(fullfile(filepath, filepattern));
            filename = fullfile(filepath, filedir.name);
            
            data = ncinfo(filename);
            no = ncread(data.Filename, 'NO');
            no2 = ncread(data.Filename, 'NO2');
            %vertical sumup
            no2 = squeeze(nansum(no2, 3));
            no = squeeze(nansum(no, 3));
            
            %unit conversion: moles/s to mg/h/m²
            no2 = no2 .* 3600 .* 46 .*1./(16000000).* 1000;
            no = no .* 3600 .* 30 .*1./(16000000).* 1000;
            nox = no2 + no;
            
            %read geospatial coordinates
            met = readtable(fullfile(filepath, '4km.csv'));
            ic = met.I_CELL;
            jc = met.J_CELL;
            lonc = met.Centroid_X;
            latc = met.Centroid_Y;
            
            lon = nan(size(no,1),size(no,2));
            lat = nan(size(no,1),size(no,2));
            
            for i=1:numel(ic)
                lon(ic(i), jc(i)) = lonc(i);
                lat(ic(i), jc(i)) = latc(i);
            end

            [tlon,tlat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = tlon';
            mlat = tlat';
            
            carb = nan(size(mlon,1),size(mlon,2),size(nox,3));
            for i=1:size(nox,3)
                this_nox = double(nox(:,:,i));
                carb(:,:,i) = reshape(griddata(lon(:), lat(:), this_nox(:), mlon(:), mlat(:)),size(mlon));
            end
            
        end
        
        function make_carb_emissions()
            dayofmonths = misc_flux_calculation.dayofmonths;
            carb = [];
            for i=1:numel(dayofmonths)
                carb_daily = misc_flux_calculation.read_carb_emissions_perday(dayofmonths(i));
                carb = cat(4, carb, carb_daily);
            end
            filepath = misc_flux_calculation.carb_filepath;
            save(fullfile(filepath, 'carb_ready'),'carb');
        end
        
        function no = read_five_emissions_pertype(filepath, filename, dayofmonth)
            [tlon,tlat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = tlon';
            mlat = tlat';
            daynum = weekday(datenum('2021-06-01') + dayofmonth -1);
            if daynum == 1
                filepath = fullfile(filepath, 'sundy');
            elseif daynum == 7
                filepath = fullfile(filepath, 'satdy');
            else
                filepath = fullfile(filepath, 'weekdy');
            end
            data = ncinfo(fullfile(filepath, filename));
            if contains(filename, 'chemi')
                eno = ncread(data.Filename, 'E_NO');
                eno2 = ncread(data.Filename, 'E_NO2');
                eno = double(squeeze(nansum(eno, 3)));
                eno2 = double(squeeze(nansum(eno2,3)));
                % 'mol km^-2 hr^-1' to mg/m2/hr
                eno = eno*30*1e3/1e6;
                eno2 = eno2*46*1e3/1e6;
                eno = eno + eno2;
            elseif contains(filename, 'road')
                % 'metric_Ton(NO2 equiv) hr^-1' to mg/m2/hr
                eno = ncread(data.Filename, 'NOX');
                eno = eno*1e9/1e6;
            end
            

            elon = ncread(data.Filename,'XLONG');
            elat = ncread(data.Filename, 'XLAT');
            elon = double(elon);
            elat = double(elat);
            
            figure;
            pcolor(elon, elat, squeeze(eno(:,:,7)));
            shading flat;
            colorbar;
            no = nan(size(mlon,1),size(mlon,2),size(eno,3));
            for j=1:size(eno,3)
                this_nox = double(eno(:,:,j));
                no(:,:,j) = reshape(griddata(elon(:), elat(:), this_nox(:), mlon(:), mlat(:)),size(mlon));
            end
        end

        function make_five_emissions_4km()
            dayofmonths = misc_flux_calculation.dayofmonths;
            five = [];
            for i=1:numel(dayofmonths)
                dayofmonth = dayofmonths(i);
                filepath = misc_flux_calculation.five_filepath;
                filename = 'wrfchemi_12z_d01';
                no = misc_flux_calculation.read_five_emissions_pertype(filepath, filename, dayofmonth);
                five = cat(4, five, no);
            end
            
            save(fullfile(misc_flux_calculation.five_filepath, 'five_ready'),'five');
        end
        
        function make_five_emissions_4km_sectors()
            % type: offroad, onroad-gas, onroad-dsl
            
            dayofmonths = misc_flux_calculation.dayofmonths;
            five = [];
            for i=1:numel(dayofmonths)
                dayofmonth = dayofmonths(i);
                types = {'onroad-gas', 'offroad', 'onroad-dsl'};
                filenames = {'onroad_12to24Z.nc', 'offroad_12to24Z.nc', 'onroad_12to24Z.nc'};
                t_emis = [];
                for i_type = 1:numel(types)
                    filepath = fullfile(misc_flux_calculation.bm_filepath, types{i_type});
                    filename = filenames{i_type};
                    no = misc_flux_calculation.read_five_emissions_pertype(filepath, filename, dayofmonth);
                    t_emis = t_emis + no;
                end
                five = cat(4, five, t_emis);
            end
            
            save(fullfile(misc_flux_calculation.bm_filepath, 'five_ready'),'five');
        end
        
        function make_five_emissions_onroad_1km()
            [tlon,tlat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = tlon';
            mlat = tlat';
            dayofmonths = misc_flux_calculation.dayofmonths;
            
            filepath = misc_flux_calculation.bm_filepath;
            data = load(fullfile(filepath, 'NOx_mobile_emissions_McDonald'));
            % unit: grams per day to mg/(h m2)
            dies = data.dailyavg_june_NOx_dies *1e3/(1e6);
            gas = data.dailyavg_june_NOx_gas*1e3/(1e6);
            flon = data.lon;
            flat = data.lat;

            utc_times = 17:24;
            local_times = utc_times - 7;
            
            five = nan(size(mlon, 1), size(mlon, 2), numel(local_times), numel(dayofmonths));
            for i_day = 1:numel(dayofmonths)
                for i_time = 1:numel(local_times)
                    dayofweek = weekday(datenum('2021-06-01') + dayofmonths(i_day) - 1);
                    [gasf, diesf] = em_emission_scaling(dayofweek, local_times(i_time));
                    this_emis = dies * diesf + gas * gasf;
                    this_emis = reshape(griddata(flon(:), flat(:), this_emis(:), mlon(:), mlat(:)),size(mlon));
                    five(:,:,i_time, i_day) = this_emis;
                end
            end

            save(fullfile(misc_flux_calculation.bm_filepath, 'five_1km_onroad_ready'),'five', 'utc_times');


            function [gasf, dieself] = em_emission_scaling(dayofweek, this_hour)
            
                ref_hour = 10:18;
                %sunday through saturday
                gas_factor = [0.052113,0.058116,0.062292,0.063423,0.063423,0.06264,0.061944,0.058812,0.052461;...
                    0.049995,0.05148,0.053856,0.05643,0.062568,0.069399,0.073656,0.073557,0.05841;...
                    0.051005,0.05252,0.054944,0.05757,0.063832,0.070801,0.075144,0.075043,0.05959;...
                    0.05151,0.05304,0.055488,0.05814,0.064464,0.071502,0.075888,0.075786,0.06018;...
                    0.053025,0.0546,0.05712,0.05985,0.06636,0.073605,0.07812,0.078015,0.06195;...
                    0.05489,0.05874,0.06237,0.06611,0.07194,0.07711,0.07953,0.07854,0.0671;...
                    0.059584,0.064386,0.06664,0.067032,0.06713,0.066934,0.066542,0.063504,0.05635;...
                    ];
    
                diesel_factor= [0.015341,0.015515,0.015312,0.014645,0.014065,0.013572,0.013601,0.013456,0.013514;...
                    0.081995,0.08211,0.079235,0.07498,0.068425,0.05842,0.0483,0.04002,0.03634;...
                    0.09269,0.09282,0.08957,0.08476,0.07735,0.06604,0.0546,0.04524,0.04108;...
                    0.091977,0.092106,0.088881,0.084108,0.076755,0.065532,0.05418,0.044892,0.040764;...
                    0.091977,0.092106,0.088881,0.084108,0.076755,0.065532,0.05418,0.044892,0.040764;...
                    0.0882,0.08712,0.08268,0.0768,0.06828,0.05748,0.04716,0.0396,0.03504;...
                    0.032781,0.031752,0.029155,0.026166,0.023275,0.02009,0.01764,0.015631,0.014161;...
                    ];
                
                indx = ref_hour == this_hour;
                gasf = gas_factor(dayofweek, indx);
                dieself = diesel_factor(dayofweek, indx);
                
            end
        

        end
        
        function [emis_agri, emis_natu] = read_dbsnp_emissions_perday(dayofmonth)
            filepath = misc_flux_calculation.iowa_filepath;
            filepattern = sprintf('wrfout_d02_2021-06-%02d.nc', dayofmonth);
            filedir = dir(fullfile(filepath, filepattern));
            filename = fullfile(filepath, filedir.name);
            
            data = ncinfo(filename);
            
            lon = ncread(data.Filename, 'XLONG');
            lat = ncread(data.Filename, 'XLAT');
            
            agri = ncread(data.Filename, 'EBIO_NO_AG');
            natu = ncread(data.Filename, 'EBIO_NO_NA');

            % only select utc 12:23
            lon = lon(:,:,13:end);
            lat = lat(:,:,13:end);
            agri = agri(:,:,13:end);
            natu = natu(:,:,13:end);

            smois = ncread(data.Filename, 'SMOIS');
            ts = ncread(data.Filename, 'TSLB');
            
            %unit conversion: mol/h/km2 to mg/h/m²
            agri = agri*30*1e3/1e6;
            natu = natu*30*1e3/1e6;

            [tlon,tlat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = tlon';
            mlat = tlat';

            emis_agri = nan(size(mlon,1),size(mlon,2),size(agri,3));
            emis_natu = nan(size(mlon,1),size(mlon,2),size(natu,3));
            for i=1:size(agri,3)
                this_agri = double(agri(:,:,i));
                this_natu = double(natu(:,:,i));
                this_lon = double(squeeze(lon(:,:,i)));
                this_lat = double(squeeze(lat(:,:,i)));
                emis_agri(:,:,i) = reshape(griddata(this_lon(:), this_lat(:), this_agri(:), mlon(:), mlat(:)),size(mlon));
                emis_natu(:,:,i) = reshape(griddata(this_lon(:), this_lat(:), this_natu(:), mlon(:), mlat(:)),size(mlon));
            end

        end

        function make_dbsnp_emissions_2km()
            dayofmonths = misc_flux_calculation.dayofmonths;
            emis_agri = [];
            emis_natu = [];
            for i=1:numel(dayofmonths)
                [emis_agri_daily, emis_natu_daily] = misc_flux_calculation.read_dbsnp_emissions_perday(dayofmonths(i));
                emis_agri = cat(4, emis_agri, emis_agri_daily);
                emis_natu = cat(4, emis_natu, emis_natu_daily);
            end
            filepath = misc_flux_calculation.iowa_filepath;
            save(fullfile(filepath, 'dbsnp_ready'), 'emis_agri', 'emis_natu');
        end

        function make_merge_full_emissions(opt)
            % list of emissions:
            % 1 CARB
            % 2 FIVE 4km + MEGAN = F4M
            % 3 FIVE 1km (on road) + MEGAN = F1M
            % 4 FIVE 4km + BDSNP = F4B
            % 5 FIVE 1km (on road) + BDSNP = F1B
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';

            switch opt
                case 1
                    data = load(fullfile(misc_flux_calculation.carb_filepath, 'carb_ready'));
                    emis = data.carb;
                    emis = emis(:,:,13:24,:);
                    utc_times = 12:23;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'carb'), 'emis','utc_times');
                case 2
                    data = load(fullfile(misc_flux_calculation.carb_filepath, 'carb_ready'));
                    emis = data.carb;
                    emis = emis(:,:,13:24,:);
                    utc_times = 12:23;
                    data = load(fullfile(misc_flux_calculation.five_filepath, 'five_ready'));
                    five = data.five;
                    indx = five > 0;
                    emis(indx) = five(indx);
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F4M'), 'emis','utc_times');
                case 3
                    
                    data = load(fullfile(misc_flux_calculation.bm_filepath, 'five_1km_onroad_ready'));
                    five = data.five;
                    utc_times = data.utc_times;
                    
                    data = load(fullfile(misc_flux_calculation.carb_filepath, 'carb_ready'));
                    emis = data.carb;
                    emis = emis(:,:,18:25,:);

                    indx = five > 0;
                    emis(indx) = five(indx);
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F1M'), 'emis','utc_times');

                case 4
                    data = load(fullfile(misc_flux_calculation.five_filepath, 'five_ready'));
                    five = data.five;
                    utc_times = 12:23;
                    data = load(fullfile(misc_flux_calculation.iowa_filepath, 'dbsnp_ready'));
                    emis = data.emis_agri + data.emis_natu+ five;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F4B'), 'emis','utc_times');
                case 5
                    data = load(fullfile(misc_flux_calculation.iowa_filepath, 'dbsnp_ready'));
                    emis = data.emis_agri + data.emis_natu;
                    utc_steps = 13:24;

                    data = load(fullfile(misc_flux_calculation.bm_filepath, 'five_1km_onroad_ready'));
                    five = data.five;
                    utc_times = data.utc_times;
                    indx = find(utc_steps == utc_times(1));
                    emis = emis(:,:,indx:end,:)+ five;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F1B'), 'emis','utc_times');
            end
            
        end

        function make_flux_match_emissions()
            % list of emissions:
            % 1 CARB
            % 2 FIVE 4km + MEGAN = F4M
            % 3 FIVE 1km (on road) + MEGAN = F1M
            % 4 FIVE 4km + BDSNP = F4B
            % 5 FIVE 1km (on road) + BDSNP = F1B
            emis_str = {'CARB','F4M','F1M','F4B','F1B'};

%             data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-fp-emis'));
%             segs = data.segs;
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv'));
            segs = data.segment_avg;
            
            for i_seg=42:numel(segs)
                fprintf('Seg: %d\n', i_seg);
                mphi = misc_flux_calculation.make_footprints_per_seg(segs(i_seg));
                for i_emis = 1:numel(emis_str)
                    emis_name = emis_str{i_emis};
                    data = load(fullfile(misc_flux_calculation.emission_filepath, emis_name));
                    emis = data.emis;
                    utc_times = data.utc_times;
                    dayofmonths = misc_flux_calculation.dayofmonths;
                    memis = [];
                    for i = 1:numel(segs(i_seg).lon)
                        % match the time
                        this_time = segs(i_seg).time(i);
                        this_day = segs(i_seg).dayofmonth;
                        this_hour = round(this_time/86400*24);
    
                        indx_hour = this_hour == utc_times;
                        indx_day = this_day == dayofmonths;
    
                        this_mphi = mphi(i, :, :);

                        this_emis = emis(:,:,indx_hour, indx_day);
                        this_emis = nansum(this_emis(:).*this_mphi(:));
                        memis = [memis, this_emis];
                    end
                    segs(i_seg).(emis_name) = memis;
                end
                segs(i_seg).fp = mphi;
                save(fullfile(misc_flux_calculation.output_filepath,'recap-sjv-fp-emis'), '-v7.3', 'segs');
            end
            save(fullfile(misc_flux_calculation.output_filepath,'recap-sjv-fp-emis'), '-v7.3', 'segs');
        end
        
        function make_emissions_csv()
            total_segs = misc_flux_calculation.combine_segs_emissions();
            %combine shrublands and forest
            total_segs.lc(total_segs.lc == 5) = 4;
            recap = repmat(total_segs.aflux*5*3600, [3,1]);
            emis = [total_segs.CARB; total_segs.F4M; total_segs.F4B];
            gp_names = { 'CARB', 'FIVE+MEGAN', 'FIVE+BDSNP'};
            groups = [];
            nobs = numel(total_segs.aflux);
            for i=1:numel(gp_names)
                this_group = repmat(gp_names(i),nobs,1);
                groups = [groups; this_group];
            end

            data.recap = recap;
            data.emis = emis;
            data.groups = groups;

            T = struct2table(data);
            writetable(T,'./output/emis_eval.csv','Delimiter',',')  
        end

        %%%% UTILITY %%%%
        function merge = merge_sjv_segments()
            dayofmonths = misc_flux_calculation.dayofmonths;

            for i=1:numel(dayofmonths)
                segments = misc_flux_calculation.read_segments(dayofmonths(i));
                if i==1
                    merge = segments;
                else
                    merge = [merge, segments];
                end
            end
           
        end
        
        function segments = read_segments(dayofmonth)
            close all
            f_filepath = misc_flux_calculation.la_filepath;
            h_filepath = dir(f_filepath);
            for i =1:numel(h_filepath)
                if contains(h_filepath(i).name, sprintf('June%d',dayofmonth))
                    filepath = fullfile(f_filepath, h_filepath(i).name);
                    break;
                end
            end

            data = load(fullfile(filepath, 'wavelet_decomp_const_lag_final.mat'));
            segments = data.segments;
            for i=1:numel(segments)
                segments(i).dayofmonth = dayofmonth;
            end
        end
        
        function [merged_segments] = read_merge_segments(dayofmonth)
            close all
            f_filepath = '/Users/monicazhu/Dropbox/LA_Flux/Flux campaign files';
            h_filepath = dir(f_filepath);
            for i =1:numel(h_filepath)
                if contains(h_filepath(i).name, sprintf('June%d',dayofmonth))
                    filepath = fullfile(f_filepath, h_filepath(i).name);
                    break;
                end
            end

            data = load(fullfile(filepath, 'wavelet_decomp_const_lag_final.mat'));
            merged_filepath = filepath;
            segments = data.segments;
            
            seg_fds = fieldnames(segments);
            spec_fds = {'flux';'aflux';'lod'};%fieldnames(segments(1).spec);
            fds = [seg_fds; spec_fds];
            merged_segments = make_empty_struct_from_cell(fds);
            for i=1:numel(segments)
                for j=1:numel(fds)
                    if contains(fds{j}, spec_fds)
                        this_spec = segments(i).spec;
                        this_flux = cat(1, this_spec(1).(fds{j}), this_spec(2).(fds{j}), this_spec(3).(fds{j}));
                        merged_segments.(fds{j}) = cat(1, merged_segments.(fds{j}), this_flux');
                    else
                        merged_segments.(fds{j}) = cat(1, merged_segments.(fds{j}), segments(i).(fds{j}));
                    end
                    
                end
            end
        end
        
        function flux_movavg = make_moveavg(flux, airspeed)
            avgspeed = mean(airspeed);
            %2km moving average
            n_avg = round((2000/avgspeed)/0.2);
            flux_movavg = movmean(flux, n_avg);

        end
        
        function [hw_lon, hw_lat, hw_selected] = read_highway()
            hw = shaperead(fullfile(misc_flux_calculation.highway_filepath, 'xc453kn9742.shp'));
            latlim = [34,38];
            lonlim = [-121,-118.5];
            routes = [99, 43, 5, 46, 41, 198, 180, 145, 137,65];
            in_domain = [];
            for i=1:numel(hw)
                if any(routes == hw(i).route)
                    bb = hw(i).BoundingBox;
                    bb_lons = bb(:,1);
                    bb_lats = bb(:,2);
                    in_domain_lim = bb_lons(1) >= lonlim(1) && bb_lons(2) <= lonlim(2) && ...
                        bb_lats(1) >= latlim(1) && bb_lats(2) <= latlim(2);
                    if in_domain_lim
                        in_domain = cat(1, in_domain, i);
                    end
                end
            end
            hw_selected = hw(in_domain);
            hw_lon = [];
            hw_lat = [];
            
            for i=1:numel(hw_selected)
                hw_lon = cat(2, hw_lon, hw_selected(:).X);
                hw_lat = cat(2, hw_lat, hw_selected(:).Y);
            end
%             figure;
%             geoshow(hw);
        end
        
        function [mhw, mlon, mlat] = make_highway_landcover()
            [hw_lon, hw_lat] = misc_flux_calculation.read_highway();
            hw_lon = [hw_lon, movmean(hw_lon,2), movmean(movmean(hw_lon,2),2)];
            hw_lat = [hw_lat, movmean(hw_lat,2), movmean(movmean(hw_lat,2),2)];
            
            mhw_lon = round((hw_lon + 122)/0.0055)+1;
            mhw_lat = round((hw_lat -34.5)/0.0045)+1;
            
            
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            
            mlon = lon';
            mlat = lat';
           
            mhw = zeros(size(mlon));
            
            indx = mhw_lat >0 & mhw_lon>0 & mhw_lat <= size(mhw,2) & mhw_lon <= size(mhw, 1);
            mhw_lon = mhw_lon(indx);
            mhw_lat = mhw_lat(indx);
            
            for i=1:numel(mhw_lon)
                mhw(mhw_lon(i), mhw_lat(i)) = 1;
            end
            
        end
        
        function [daily_segs, total_segs] = combine_segs()
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv'));
            segs = data.segment_avg;
            fds = fieldnames(segs(1));
            dayofmonths = misc_flux_calculation.dayofmonths;
            
            % daily flux distribution
            seg_days = [];
            for i=1:numel(segs)
                seg_days = [seg_days, segs(i).dayofmonth];
            end
            
            for i=1:numel(dayofmonths)
                indx = dayofmonths(i) == seg_days;
                this_segs = segs(indx);
                daily_segs(i) = make_empty_struct_from_cell(fds);
                for i_seg=1:numel(this_segs)
                    for i_fd = 1:numel(fds)
                        daily_segs(i).(fds{i_fd}) = [daily_segs(i).(fds{i_fd}); this_segs(i_seg).(fds{i_fd})];
                    end
                end
            end
            
            % total flux distribution
            total_segs = make_empty_struct_from_cell(fds);
            for i_seg=1:numel(daily_segs)
                for i_fd = 1:numel(fds)
                    total_segs.(fds{i_fd}) = [total_segs.(fds{i_fd}); daily_segs(i_seg).(fds{i_fd})];
                end
            end
        end
        
        function [total_segs, daily_segs] = combine_segs_emissions()
            % fp is too large, won't read in;
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-fp-emis'));
            segs = data.segs;
            fds = fieldnames(segs(1));
            fds = setdiff(fds, {'fp'});
            dayofmonths = misc_flux_calculation.dayofmonths;
            
            % total flux distribution
            total_segs = make_empty_struct_from_cell(fds);
            for i_seg=1:numel(segs)
                for i_fd = 1:numel(fds)
                    if size(segs(i_seg).(fds{i_fd}),1) == 1
                        total_segs.(fds{i_fd}) = [total_segs.(fds{i_fd}); segs(i_seg).(fds{i_fd})'];
                    else
                        total_segs.(fds{i_fd}) = [total_segs.(fds{i_fd}); segs(i_seg).(fds{i_fd})];
                    end
                    
                end
            end

            % daily flux distribution
            seg_days = [];
            for i=1:numel(segs)
                seg_days = [seg_days, segs(i).dayofmonth];
            end
            
            for i=1:numel(dayofmonths)
                indx = dayofmonths(i) == seg_days;
                this_segs = segs(indx);
                daily_segs(i) = make_empty_struct_from_cell(fds);
                for i_seg=1:numel(this_segs)
                    for i_fd = 1:numel(fds)
                        if size(segs(i_seg).(fds{i_fd}),1) == 1
                            daily_segs(i).(fds{i_fd}) = [daily_segs(i).(fds{i_fd}); segs(i_seg).(fds{i_fd})'];
                        else
                            daily_segs(i).(fds{i_fd}) = [daily_segs(i).(fds{i_fd}); segs(i_seg).(fds{i_fd})];
                        end
                    end
                end
            end
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
        
        function time_window = make_deposition_time_window(dayofmonth)
            %SJV
            if dayofmonth == 3
                time_window = [11184,13301];
            end
            if dayofmonth == 8
                time_window = [1.33e4,1.5432e4];
            end
            
            if dayofmonth == 16
                time_window = [8.146e3, 8.447e3];
            end
            
            if dayofmonth == 9
                time_window = 1.0e+04 *[1.3344, 1.5775];
            end
            
            if dayofmonth == 13
                time_window = 1.0e+04 *[1.2491, 1.5085];
            end
            
            if dayofmonth == 15
                 time_window = 1.0e+04 *[ 1.1955, 1.4222];
            end
            
            if dayofmonth == 22
                time_window = 1.0e+04 *[1.1636, 1.4457];
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