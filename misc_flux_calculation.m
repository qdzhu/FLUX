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
            

            savedata.NOmgm3=(transpose(data.aveNO).*30.*12.187)./(273.15+savedata.temp)./prctile(savedata.pressure, 50); %calculate mg/m3 from ppb
            savedata.NO2mgm3=(transpose(data.aveNO2).*46.*12.187)./(273.15+savedata.temp)./prctile(savedata.pressure, 50); %calculate mg/m3 from ppb
%             savedata.NOzmgm3=(transpose(data.aveNOz).*46.*12.187)./(273.15+savedata.temp)./1000; %calculate mg/m3 from ppb
            savedata.LIF_time = LIF_time;
           % disp(nanmedian(savedata.pressure));
            no2 = data.aveNO2;
            nox = data.aveNO + data.aveNO2;
            save(fullfile(filepath, 'merge_obs'), 'savedata');
            clear all;
        end
        
        function make_flux_calculation_final(dayofmonth)
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
                    
                    lag_corr_opt = 'const_lag';
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

                    end
                    seg_corr.spec = spec;
                    % copy seg_sorr to segments(n_spec)
                    fds = fieldnames(seg_corr);
                    for i_fd = 1:numel(fds)
                        segments(n_spec).(fds{i_fd}) = seg_corr.(fds{i_fd});
                    end
                    n_spec = n_spec + 1;
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