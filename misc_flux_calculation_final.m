classdef misc_flux_calculation_final
    properties(Constant)
        highway_filepath = './inputs/highway_national';
        input_filepath = './inputs/';
        output_filepath = './output';
        figure_filepath = '/Users/zhu/Dropbox/Publications/recap-nox/Figures/';
        emission_filepath = './output/emissions';
        la_filepath = './Flux campaign files';
        hrrr_filepath = '../Database/HRRR';
        vista_filepath = '/Users/zhu/Dropbox/Database/NACP-Vista/data_unzip';

        % emissions
        carb_filepath = '../Database/CARB';
        carb_2k_filepath = '../Database/CARB_2021'
        iowa_filepath = '../Database/Soil_IOWA/daily';
        five_filepath = '../Database/Katelyn/Anthropogenic_emissions';
        beis_filepath = '../Database/BEIS/';
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
            %remove negative values in pressure
            savedata.pressure(savedata.pressure<0) = nan;
            %savedata.NOmgm3=(transpose(data.aveNO).*30.*12.187)./(273.15+savedata.temp)./savedata.pressure; %calculate mg/m3 from ppb
            %savedata.NO2mgm3=(transpose(data.aveNO2).*46.*12.187)./(273.15+savedata.temp)./savedata.pressure; %calculate mg/m3 from ppb
            savedata.NOmgm3 = savedata.pressure*1e2.*transpose(data.aveNO)*1e-3*30./((273.15+savedata.temp)*28.97*287.05);
            savedata.NO2mgm3 = savedata.pressure*1e2.*transpose(data.aveNO2)*1e-3*46./((273.15+savedata.temp)*28.97*287.05);
            %savedata.NOmgm3=(transpose(data.aveNO).*30.*12.187)./(273.15+savedata.temp)./nanmedian(savedata.pressure); %calculate mg/m3 from ppb
            %savedata.NO2mgm3=(transpose(data.aveNO2).*46.*12.187)./(273.15+savedata.temp)./nanmedian(savedata.pressure); %calculate mg/m3 from ppb
%             savedata.NOzmgm3=(transpose(data.aveNOz).*46.*12.187)./(273.15+savedata.temp)./1000; %calculate mg/m3 from ppb
            savedata.LIF_time = LIF_time;
           % disp(nanmedian(savedata.pressure));
            no2 = data.aveNO2;
            nox = data.aveNO + data.aveNO2;
            
            % pblh filtering
            pblh_filtering();

            %line(savedata.LIF_time - savedata.LIF_time(1), savedata.NOmgm3);
            % plot the altitude vs time, manually select the time window
            if do_plot 
                indx = savedata.LIF_time<86400;
                LIF_time = savedata.LIF_time(indx);
                alt = savedata.alt(indx);
                time = LIF_time - LIF_time(1);
                
                figure;
                subplot(1,2,1);
                line(time, alt, 'linestyle', 'none', 'marker', '.');
                ylabel('Altitude (m)');
                xlabel('Time (s)');
                subplot(1,2,2);
                scatter(savedata.NOmgm3 + savedata.NO2mgm3, savedata.alt);
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
                n = 1;
                for i=1:numel(savedata.lon)
                    this_hour = round(savedata.LIF_time(i)/86400*24);
                    hour_indx = this_hour == hours;
                    if ~isnan(savedata.lon(i)) && ~isnan(savedata.lat(i))
                        lon_indx = round((savedata.lon(i) - misc_flux_calculation.mLon(1))/0.0055);
                        lat_indx = round((savedata.lat(i) - misc_flux_calculation.mLat(1))/0.0045);
                        if lon_indx >= 1 && lon_indx <= numel(misc_flux_calculation.mLon) && ...
                                lat_indx >=1 && lat_indx <= numel(misc_flux_calculation.mLat)
                            if numel(mpblh(lon_indx, lat_indx, hour_indx))==0
                                this_pblh = nan;
                            else
                                this_pblh = mpblh(lon_indx, lat_indx, hour_indx);
                            end

                            savedata.hrrr_pblh(i) = this_pblh;
                            if savedata.radalt(i) >= this_pblh
                                savedata.NOmgm3(i) = nan;
                                savedata.NO2mgm3(i) = nan;
                                n = n+1;
                            end
                        end
                    end
                end

            end
        end
        
        % Make flux calculation
        function make_flux_calculation_final(dayofmonth, lag_corr_opt)
            % READ DATA INPUT
            %close all
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
            %data_ready = data;
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
            % set nox to be nan if it is not in the time window
            data_ready = filter_time_window(data_ready);

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
                % calculate the max difference in altitudes

                alt_diff = nanmax(data_ready.alt(roll_spikes(j):roll_spikes(j+1))) - nanmin(data_ready.alt(roll_spikes(j):roll_spikes(j+1)));
                % determine if it is with in the time window
                %in_time_window = misc_flux_calculation.is_in_time_window(time_beg - data_ready.LIF_time(1), time_end - data_ready.LIF_time(1), time_window);
                
                %Need segment to be greater than 10km (or maybe even 15km) for good wavelet analysis otherwise ignore 
                if segmentlengthkm > 10 && alt_diff < 200
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
                        j1 = round(log2(size(sst69,1)*dt/s0))/dj; 
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
                        
                        avg = scale69 < 128;
                        scale_avg = (scale69')*(ones(1,n));  % expand scale --> (J+1)x(N) array
                        scale_avg = wc6971./ scale_avg;   % [Eqn(24)]
                        %scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
                        scale_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(scale_avg);
                        spec(i).scale_avg = scale_avg;
                        % COI influence
                        period_big = repmat(period69',1,n);
                        coi_big    = repmat(coi69,length(scale69),1);
                        icoi = (coi_big>period_big);
                        crosspeccoi = (wc6971.*~icoi);
                        ratiocoi = sum(crosspeccoi)./sum(wc6971);
                        coi_avg = (scale69')*(ones(1,n));
                        coi_avg = crosspeccoi./ coi_avg;   % [Eqn(24)]
                        coi_avg = sqrt(variance1*variance2)*dj*dt/Cdelta*sum(coi_avg(avg,:));   % [Eqn(24)]
                        qcoi = abs(coi_avg./scale_avg) ; %fraction of flux under COI
                        spec(i).qcoi = qcoi;
                        scale_avg(qcoi>0.8) = nan;
                        spec(i).flux = scale_avg;
                        aflux = misc_flux_calculation.make_moveavg(scale_avg, seg_corr.airspeed);
                        spec(i).aqcoi = misc_flux_calculation.make_moveavg(qcoi, seg_corr.airspeed);
                        indx = spec(i).aqcoi>0.2;
                        aflux(indx) = nan;
                        spec(i).aflux = aflux;
                        %lod (limit of detection) calculation
                        spec(i).lod = calc_lod_seg(seg_corr, i);
                        spec(i).lod_wn = calc_lod_seg_white_noise(seg_corr, i);
                        spec(i).std = zeros(size(spec(i).flux)) + spec(i).lod_wn;
                        spec(i).std_rdm = spec(i).lod_wn./nanmean(spec(i).aflux);
                        std_lowfreq = 2.2*(seg_corr.alt./seg_corr.hrrr_pblh).^0.5.*seg_corr.hrrr_pblh./segmentlengthkm/1e3;
                        spec(i).std_lowfreq = std_lowfreq;
% 
%                         spec(i).se=2.2*pbl*sqrt(zi/pbl)/(xkm*(vw(size(vw,1),1)-vw(1,1))*3600*24);
                        spec(i).re= 1.75*(seg_corr.alt./seg_corr.hrrr_pblh).^(0.25).*sqrt(seg.hrrr_pblh./segmentlengthkm/1000);
                    end
                    
                    if size(seg_corr.waves, 2)
                        seg_corr.segmentlengthkm = zeros(size(seg_corr.time)) + segmentlengthkm;
                        seg_corr.spec = spec;
                        % copy seg_sorr to segments(n_spec)
                        fds = fieldnames(seg_corr);
                        for i_fd = 1:numel(fds)
                            segments(n_spec).(fds{i_fd}) = seg_corr.(fds{i_fd});
                        end
                        n_spec = n_spec + 1;
                        hold on;
                        plot(seg_corr.time);
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
                merge.hrrr_pblh = nanmedian(data.savedata.hrrr_pblh);
                
                time_filter = merge.LIF_time;
                radalt = data.savedata.alt(indx)- data.savedata.galt(indx);
                galt = data.savedata.galt(indx);
                radalt(radalt>700) = nan;
                indx_a = ~isnan(radalt);
                galt = interp1(time_filter(indx_a), galt(indx_a), time_filter);

                % Note: Altitude above sea level - surface altitude =
                % Altitude above ground level.
                merge.alt = data.savedata.alt(indx)- galt;
                merge.airspeed = data.savedata.airspeed(indx);
                merge.roll = data.savedata.roll(indx);
                merge.pitch = data.savedata.pitch(indx);
                merge.NOmgm3 = data.savedata.NOmgm3(indx);
                merge.NO2mgm3 = data.savedata.NO2mgm3(indx);
            end
            
            function data_ready = filter_time_window(data_ready)

                LIF_time = data_ready.LIF_time;
                LIF_time = LIF_time - LIF_time(1);
                indx = zeros(size(LIF_time));

                for i_time=1:numel(LIF_time)
                    time_window = misc_flux_calculation.make_time_window(dayofmonth);
                    if misc_flux_calculation.is_in_time_window(LIF_time(i_time), LIF_time(i_time), time_window)
                        indx(i_time) = 1;
                    end
                    
                    time_window = misc_flux_calculation.make_deposition_time_window(dayofmonth);
                    if misc_flux_calculation.is_in_time_window(LIF_time(i_time), LIF_time(i_time), time_window)
                        indx(i_time) = 2;
                    end
                    
                    time_window = misc_flux_calculation.make_racetrack_time_window(dayofmonth);
                    if misc_flux_calculation.is_in_time_window(LIF_time(i_time), LIF_time(i_time), time_window)
                        indx(i_time) = 3;
                    end
                end
                data_ready.time_filter = indx;
                data_ready.NOmgm3(indx==0) = nan;
                data_ready.NO2mgm3(indx==0) = nan;
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
                seg.hrrr_pblh = data_ready.hrrr_pblh;
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
                    seg_corr.hrrr_pblh = seg.hrrr_pblh;
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
                seg_corr.hrrr_pblh = seg.hrrr_pblh;
                seg_corr.nox = seg.NOx;
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
            
            function lod = calc_lod_seg_white_noise(seg, i_wave)
                    delaylod = 1:5;
                    this_aflux = [];
                    for j=1:numel(delaylod)
                        this_waves = seg.waves;
                        this_waves = this_waves(:,i_wave);
                        this_vw = seg.vw;
                        sst69 = this_waves; 
                        sst71 = this_vw;  
                        variance1 = std(sst69)^2;
                        sst69 = (sst69 - mean(sst69))/sqrt(variance1) ;
                        variance2 = std(sst71)^2;
                        sst71 = (sst71 - mean(sst71))/sqrt(variance2) ;
                        n = length(sst69);
                        %generate the white noise
                        sst71 = 1.*randn(n,1);
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
%                 misc_flux_calculation.merge_nox_met(dayofmonths(i), true);
                misc_flux_calculation.make_flux_calculation_final(dayofmonths(i), 'seg_lag');
                misc_flux_calculation.make_flux_calculation_final(dayofmonths(i), 'const_lag');
            end
        end

        % correct the vertical divergence
        function make_vertical_correction()
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-fpcontour'));
            segs = data.segs;
            s = [];
            re = [];
            w = [];
            u = [];
            for i = 1:numel(segs)
                ratio = segs(i).alt./segs(i).mpblh;
                
                %scaling based on nox vertical gradient
                %scaling = 1.7./(-1.53*ratio+1.7);
                %segs(i).scaling = scaling;
                %segs(i).aflux_adj = segs(i).aflux + ratio/0.4635/3600;
                segs(i).aflux_adj = segs(i).aflux .* 1./(1-0.5613.*ratio);
%                 segs(i).aflux_adj = segs(i).aflux./(1-0.64/0.4635*ratio);
%                 indx = segs(i).aflux < 0;
%                 segs(i).aflux_adj(indx) = -segs(i).aflux(indx)./(1-0.64/0.4635*ratio(indx)) + 2*segs(i).aflux(indx);
                %scaling based on nox vertical gradient
                scaling = segs(i).aflux + ratio/0.4635/3600;
                segs(i).scaling = scaling;
                segs(i).aflux_hadj = segs(i).aflux.*scaling;
                
                fwstar = segs(i).aflux_temp;
                temp = segs(i).temp;
                pblh = segs(i).mpblh;
                wstar=(9.81.*pblh.*abs(fwstar./(temp+273))).^(1/3);
                segs(i).wstar = wstar;
                segs(i).tstar = pblh./wstar;
                % the time is messy?
                time_diff = max(segs(i).time) - min(segs(i).time);
                xkm = time_diff*mean(segs(i).airspeed);

                segs(i).re= 1.75*(segs(i).alt./segs(i).mpblh).^(0.25).*sqrt(segs(i).mpblh./xkm);

                segs(i).ustar = segs(i).windspeed.*0.41./log(segs(i).alt./segs(i).mrs);

                s = [s; segs(i).scaling];
                re = [re; segs(i).re];
                w = [w; segs(i).wstar];
                u = [u; segs(i).ustar];
            end

            figure;
            histogram(u, 40);
            xlabel('u^*');

            figure;
            histogram(w, 40);
            xlabel('w^*');

            figure;
            subplot(1,2,1);
            histogram(re, 40);
            xlabel('Random Error');

            subplot(1,2,2);
            scatter(re, s, 3);
            xlabel('Random Error');
            ylabel('Scaling Factor');
            hold on;
            line([0.35,0.35],[0,10]);

            output_filename = 'recap-sjv-vert-fpcontour';
            save(fullfile(misc_flux_calculation.output_filepath,output_filename), '-v7.3','segs');

        end

        % make footprint
        function [rcontours_x, rcontours_y, phialong, dx_half] = make_footprints_contour_per_seg(segment)
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';
            
            sigW = std(segment.vw);
            sigU = std(segment.windspeed);
            
            rcontours_x = cell(size(segment.lon));
            rcontours_y = cell(size(segment.lon));

            phialong = zeros(size(segment.lon));
            dx_half = zeros(size(segment.lon));
            for i = 1:numel(segment.lon)
                % calculate footprint
                wd = segment.wind_dir(i);
                alt = segment.alt(i);
                rs = segment.mrs(i);
                pblh = segment.mpblh(i);
%                 fwstar = segment.aflux_temp(i);
                ws = segment.windspeed(i);
                temp = segment.temp(i);
%                 wstar=(9.81.*pblh.*abs(fwstar/(temp+273))).^(1/3);

%                 dx_half(i) = 0.9.*ws.*alt.^(2/3).*pblh.^(1/3)./wstar/1e3;
                ustar = segment.windspeed(i).*0.41./log(alt./rs);
                try
                    [phi, xrot, yrot, xr_rots, yr_rots, PHIalong, ~,Xcen,~] = calc_footprint_KL04(wd, 250, sigU, sigW, ustar, alt, rs, pblh, 0,0.1:0.1:0.9);
                    %[phi, xrot, yrot, xr_rots, yr_rots, PHIalong, ~,Xcen,~] = calc_footprint_KL04(wd, 500, sigU, sigW, ustar, alt, rs, pblh, 0, [0.9,0.95,0.99]);
                    for j=1:numel(xr_rots)
                        xr_rots{j} = segment.lon(i) + xr_rots{j}/91e3;
                        yr_rots{j} = segment.lat(i) + yr_rots{j}/111e3;
                    end
                    rcontours_x{i} = xr_rots;
                    rcontours_y{i} = yr_rots;
                    cumphi = cumsum(abs(PHIalong));
                    phialong(i) = interp1(cumphi, Xcen/1e3, 0.9);
                catch err
                    rcontours_x{i} = {nan};
                    rcontours_y{i} = {nan};
                    phialong(i)  = nan;
                end
            end
        end
        
        function make_footprints_contour(lag_corr_opt)
            switch lag_corr_opt
                case 'seg_lag'
                    filename = 'recap-sjv-seg';
                case 'const_lag'
                    filename = 'recap-sjv';
            end
            data = load(fullfile(misc_flux_calculation.output_filepath, filename));
            
            segs = data.segment_avg;
            
            for i=1:numel(segs)
                fprintf('Segs index: %d\n', i);
                [xr_rots, yr_rots, phialong, dx_half] = misc_flux_calculation.make_footprints_contour_per_seg(segs(i));
                segs(i).xr_rots = xr_rots;
                segs(i).yr_rots = yr_rots;
                segs(i).phialong = phialong;
                segs(i).dx_half = dx_half;
            end
            output_filename = sprintf('%s-fpcontour', filename);
            save(fullfile(misc_flux_calculation.output_filepath,output_filename), '-v7.3','segs');
        end
        
        % make source disaggregation: linear regression, use area
        % weights
        function make_source_attribution_contour_linear_reg_prep()
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';

            % read fps
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-vert-fpcontour'));
            segs = data.segs;

            data = load('./inputs/landcover-500m.mat');
            clc = data.clc;

           
            for i=1:numel(segs)
                this_fp = segs(i);

                fclcs = cell(size(this_fp.lon));

                for j=1:numel(this_fp.lon)
                    try
                        xr_rots = this_fp.xr_rots{j}{9};
                        yr_rots = this_fp.yr_rots{j}{9};
                        indx = ~isnan(xr_rots) & ~isnan(yr_rots);
                        xr_rots = xr_rots(indx);
                        yr_rots = yr_rots(indx);

                        yrange = [min(yr_rots)-0.1, max(yr_rots)+0.1];
                        xrange = [min(xr_rots)-0.1, max(xr_rots)+0.1];
                        indx = mlon >= xrange(1) & mlon < xrange(2) & mlat >= yrange(1) & mlat < yrange(2);

                        this_mlon = mlon(indx);
                        this_mlat = mlat(indx);
                        this_clc = clc(indx);
                        
                        in = inpolygon(this_mlon,this_mlat,xr_rots,yr_rots);
                        clc_incontour = this_clc(in);


                        this_clc = zeros(4,1);
                        this_clc(1) = nansum(clc_incontour==0)/numel(clc_incontour);
                        this_clc(2) = nansum(clc_incontour==1)/numel(clc_incontour);
                        this_clc(3) = nansum(clc_incontour == 2)/numel(clc_incontour);
                       
                        this_clc(4) = 1- sum(this_clc(1:3));
                        fclcs{j} = this_clc;
                    catch err
                        this_clc = nan(4,1);
                        fclcs{j} = this_clc;
                    end
                end
                
                segs(i).clc_fp = fclcs;
            end
            save(fullfile(misc_flux_calculation.output_filepath,'recap-sjv-vert-source-reg-contour'), 'segs');

        end
        
        function make_source_attribution_linear_contour()
            clear all;
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-vert-source-reg-contour'));
            segs = data.segs;
            regs = [];
            flux = [];
            flux_adj = [];
            flux_hadj = [];
            re = [];
            scaling = [];
            zzi = [];

            for i=1:numel(segs)
                fclcs = segs(i).clc_fp;
                flux = [flux; segs(i).aflux*3600];
                flux_adj = [flux_adj; segs(i).aflux_adj*3600];
                flux_hadj = [flux_hadj; segs(i).aflux_hadj*3600];
                re = [re; segs(i).re];
                scaling = [scaling; segs(i).scaling];
                zzi = [zzi; segs(i).alt./segs(i).mpblh];
                for j=1:numel(segs(i).lon)
                    this_clc = fclcs{j};
%                     this_clc(4) = sum(this_clc(5));
%                     this_clc(5) = 1 - sum(this_clc(1:4));
                    this_clc = this_clc([1,2,3,4]);
                    regs = cat(1, regs, this_clc');
                end
            end  
            %indx = zzi > 0.6;
            %flux(indx) = nan;
            %flux_adj(indx) = nan;
            
            slopes = zeros(3,3);
            se = zeros(3,3);

            mdl = fitlm(regs, flux,'Intercept',false);
            slopes(1,:) = mdl.Coefficients.Estimate(1:3);
            se(1,:) = mdl.Coefficients.SE(1:3);

            mdl = fitlm(regs, flux_adj,'Intercept',false);
            slopes(2,:) = mdl.Coefficients.Estimate(1:3);
            se(2,:) = mdl.Coefficients.SE(1:3);

            mdl = fitlm(regs, flux_hadj,'Intercept',false);
            slopes(3,:) = mdl.Coefficients.Estimate(1:3);
            se(3,:) = mdl.Coefficients.SE(1:3);

            figure;
            hold on;
            bar(1:3, slopes',0.8);
%             errorbar(1:3, slopes(1:3), se(1:3),'linestyle','-','linewidth',2,'color','k');
            xticks([1 2 3]);
            xticklabels({'Highway','Urban','Soil'});
            ylabel('NO_x flux (mg/(m^2 hr))');
            xlim([0.5,3.5]);
            ax = gca;
            ax.FontSize = 14; 
            %export_fig(fullfile(misc_flux_calculation.figure_filepath,'plot_source'),'-png','-m2','-transparent','-painters','-cmyk');

            mdl = fitlm(regs, flux_adj,'Intercept',false);
            slopes= mdl.Coefficients.Estimate;
            se = mdl.Coefficients.SE;
            figure;
            hold on;
            bar(1:3, slopes(1:3),0.3);
            errorbar(1:3, slopes(1:3), se(1:3),'linestyle','-','linewidth',2,'color','k');
            xticks([1 2 3]);
            xticklabels({'Highway','Urban','Soil'});
            ylabel('NO_x flux (mg/(m^2 hr))');
            xlim([0.5,3.5]);
            ax = gca;
            ax.FontSize = 14; 

        end
        
        function make_source_attribution_linear_contour_uncertainty()
            clear all;
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-source-reg-contour'));
            segs = data.segs;
            regs = [];
            flux = [];
            
            for i=1:numel(segs)
                fclcs = segs(i).clc_fp;
                flux = [flux; segs(i).aflux*3600];
                for j=1:numel(segs(i).lon)
                    this_clc = fclcs{j};
                    this_clc(3) = this_clc(3)+this_clc(4);
                    this_clc = this_clc([1,2,3,5]);
                    regs = cat(1, regs, this_clc');
                end
            end  

            
            mdl = fitlm(regs(:,1:3), flux,'Intercept',false);
            slopes_orig = mdl.Coefficients.Estimate;
            
            indx = sum(regs(:,1:3),2) ==1;
            regs = regs(indx,:);
            flux = flux(indx,:);
            slopes = [];
            for i=1:100
                y = datasample(1:numel(flux),numel(flux),'Replace',true);
                mdl = fitlm(regs(y,1:3), flux(y),'Intercept',false);
                slope = mdl.Coefficients.Estimate;
                slopes = cat(2,slopes,slope);
            end

            slopes = [];
            for i=1:100
                err = randn(numel(flux),1).*flux*0.5;
                mdl = fitlm(regs(:,1:3), flux+err,'Intercept',false);
                slope = mdl.Coefficients.Estimate;
                slopes = cat(2,slopes,slope);
            end

        end
        
    end
end