classdef misc_flux_emis_comp
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
        
        function convert_carb_4k_emissions()
            %filepath = '/Users/monicazhu/Documents/Database/CARB';
            filepath = misc_flux_calculation.carb_2k_filepath;
            data = ncinfo(fullfile(misc_flux_calculation.carb_2k_filepath,'RECAP_emis_st4k_gid649_comb3d_June2021.nc'));

            no = ncread(data.Filename, 'NO');
            no = squeeze(nansum(no, 3));
            %unit conversion: moles/s to mg/h/m²
            
            no = no .* 3600 .* 30 .*1./(16000000).* 1000;
            
            no2 = ncread(data.Filename, 'NO2');
            %vertical sumup
            no2 = squeeze(nansum(no2, 3));
            no2 = no2 .* 3600 .* 46 .*1./(16000000).* 1000;
            nox = no2 + no;

            
            %read geospatial coordinates
            met = readtable(fullfile(misc_flux_calculation.carb_filepath, '4km.csv'));
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


            dayofyear = datenum('2021-06-01') + 3 - 1 - datenum('2021-01-01') +1;
            anthro_filepattern = sprintf('anthro/st_2k*%d*', dayofyear);
            anthro = read_nox_emissions(anthro_filepattern);

            met = shaperead(fullfile(filepath, 'state2k1/state2k1.shp'));
            xlon = zeros(size(anthro,1),size(anthro,2));
            xlat = zeros(size(anthro,1),size(anthro,2));
            for i=1:numel(met)
                ic = met(i).I_CELL;
                jc = met(i).J_CELL;
                xlon(ic, jc) = nanmean(met(i).X(1:4));
                xlat(ic, jc) = nanmean(met(i).Y(1:4));
            end

            
            carb = nan(size(xlon,1),size(xlon,2),size(nox,3));
            for i=1:size(nox,3)
                fprintf('Working on %d \n',i);
                this_nox = double(nox(:,:,i));
                carb(:,:,i) = reshape(griddata(lon(:), lat(:), this_nox(:), xlon(:), xlat(:)),size(xlon));
            end
            
            save(fullfile(filepath,'total_carb_nox_emis_2k'),'carb','-v7.3');
            
            function nox = read_nox_emissions(filepattern)
                filedir = dir(fullfile(filepath, filepattern));
                filename = fullfile(filedir.folder, filedir.name);
                
                data = ncinfo(filename);
                no = ncread(data.Filename, 'NO');
                no = squeeze(nansum(no, 3));
                %unit conversion: moles/s to mg/h/m²
                
                no = no .* 3600 .* 30 .*1./(4000000).* 1000;

                if contains('NO2', data.Filename)
                    no2 = ncread(data.Filename, 'NO2');
                    %vertical sumup
                    no2 = squeeze(nansum(no2, 3));
                    no2 = no2 .* 3600 .* 46 .*1./(4000000).* 1000;
                    nox = no2 + no;
                else
                    nox = no;
                end
                    
            end
        end

        function [carb_anthro, carb_anthro_total, carb_bio] = read_carb_2k_emissions_perday(dayofmonth)
            %filepath = '/Users/monicazhu/Documents/Database/CARB';
            filepath = misc_flux_calculation.carb_2k_filepath;
            dayofyear = datenum('2021-06-01') + dayofmonth - 1 - datenum('2021-01-01') +1;

            bio_filepattern = sprintf('biogenic/MEGAN*%d*', dayofyear);
            bio = read_nox_emissions(bio_filepattern);
            
            anthro_filepattern = sprintf('anthro/st_2k*%d*', dayofyear);
            anthro = read_nox_emissions(anthro_filepattern);
            utc_times = 18:23;
            local_time = utc_times-7;
            anthro = anthro(:,:,local_time);
            bio = bio(:,:,local_time);

            data = load(fullfile(filepath,'total_carb_nox_emis_2k'));
            carb = data.carb;
            carb = carb(:,:,local_time+24*(dayofmonth-1));
            anthro_total = carb - bio;
            
            %read geospatial coordinates
            lon = nan(size(anthro,1),size(anthro,2));
            lat = nan(size(anthro,1),size(anthro,2));

            met = shaperead(fullfile(filepath, 'state2k1/state2k1.shp'));
            xloncorn = zeros(4, size(anthro,1),size(anthro,2));
            xlatcorn = zeros(4, size(anthro,1),size(anthro,2));
            for i=1:numel(met)
                ic = met(i).I_CELL;
                jc = met(i).J_CELL;
                xloncorn(:, ic, jc) = met(i).X(1:4);
                xlatcorn(:, ic, jc) = met(i).Y(1:4);
            end

            data = load('./inputs/landcover-500m.mat');
            mlon = data.mlon;
            mlat = data.mlat;
            clc = data.clc;

            
            carb_anthro = nan(size(mlon,1),size(mlon,2),size(anthro,3));
            carb_bio = nan(size(mlon,1),size(mlon,2),size(bio,3));
            carb_anthro_total = nan(size(mlon,1),size(mlon,2),size(anthro,3));
            
            for i=1:size(carb_anthro,3)
                fprintf('Working on hour %d\n', i);
                this_anthro = double(anthro(:,:,i));
                this_anthro_total = double(anthro_total(:,:,i));
                this_bio = double(bio(:,:,i));

                this_emis_anthro = nan(size(mlon,1),size(mlon,2));
                this_emis_anthro_total = nan(size(mlon,1),size(mlon,2));
                this_emis_bio = nan(size(mlon,1),size(mlon,2));
                for j=1:size(anthro,1)
                    for k=1:size(anthro,2)
                        in = inpolygon(mlon,mlat,xloncorn(:,j,k),xlatcorn(:,j,k));
                        if sum(in(:)) > 0
%                             this_emis_bio(in) = this_bio(j,k);
%                             this_emis_anthro(in) = (this_anthro(j,k));
%                             this_emis_anthro_total(in) = (this_anthro_total(j,k));
                            bioratio = sum(clc(in) == 2)/numel(clc(in)); 
                            this_emis_bio(in) = 0;
                            if bioratio > 0
                                indx = in & clc == 2;
                                this_emis_bio(indx) = this_bio(j,k)/bioratio;
                            end
                            
                            anthroratio = sum(clc(in) < 2)/numel(clc(in)); 
                            this_emis_anthro(in) = 0;
                            this_emis_anthro_total(in) = 0;
                            if anthroratio > 0
                                indx = in & clc <2;
                                this_emis_anthro(indx) = (this_anthro(j,k))/anthroratio;
                                this_emis_anthro_total(indx) = (this_anthro_total(j,k))/anthroratio;

                            end

                        end
                    end
                end
                carb_anthro(:,:,i) = this_emis_anthro;
                carb_anthro_total(:,:,i) = this_emis_anthro_total;
                carb_bio(:,:,i) = this_emis_bio;
            end

            outputfilename = sprintf('carb_ready_%d', dayofmonth);
            save(fullfile(filepath, outputfilename),'carb_anthro','carb_anthro_total','carb_bio');

            function nox = read_nox_emissions(filepattern)
                filedir = dir(fullfile(filepath, filepattern));
                filename = fullfile(filedir.folder, filedir.name);
                
                data = ncinfo(filename);
                no = ncread(data.Filename, 'NO');
                no = squeeze(nansum(no, 3));
                %unit conversion: moles/s to mg/h/m²
                
                no = no .* 3600 .* 30 .*1./(4000000).* 1000;

                if contains('NO2', data.Filename)
                    no2 = ncread(data.Filename, 'NO2');
                    %vertical sumup
                    no2 = squeeze(nansum(no2, 3));
                    no2 = no2 .* 3600 .* 46 .*1./(4000000).* 1000;
                    nox = no2 + no;
                else
                    nox = no;
                end
                    
            end
            
        end
        
        function make_carb_2k_emissions()
            dayofmonths = misc_flux_calculation.dayofmonths;

            carb_anthro = [];
            carb_anthro_total = [];
            carb_bio = [];

            for i=1:numel(dayofmonths)
                fprintf('Working on day %d \n',i);
%                 dayofmonth = dayofmonths(i);
%                 data = load(fullfile(misc_flux_calculation.carb_2k_filepath, sprintf('carb_ready_%d', dayofmonth)));
%                 carb_anthro_daily = data.carb_anthro;
%                 carb_anthro_total_daily = data.carb_anthro_total;
%                 carb_bio_daily = data.carb_bio;
                [carb_anthro_daily, carb_anthro_total_daily, carb_bio_daily] = misc_flux_emis_comp.read_carb_2k_emissions_perday(dayofmonths(i));
                carb_anthro = cat(4, carb_anthro, carb_anthro_daily);
                carb_anthro_total = cat(4, carb_anthro_total, carb_anthro_total_daily);
                carb_bio = cat(4, carb_bio, carb_bio_daily);

            end
            filepath = misc_flux_calculation.carb_2k_filepath;
            save(fullfile(filepath, 'carb_ready'),'carb_anthro','carb_anthro_total','carb_bio');
        end
        
        % Make emission evaluation
        
        function no = read_five_emissions_pertype(filepath, filename, dayofmonth)
            data = load('./inputs/landcover-500m.mat');
            mlon = data.mlon;
            mlat = data.mlat;
            clc = data.clc;
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

                %eno = eno(:,:,13:end);
            elseif contains(filename, 'road')
                % 'metric_Ton(NO2 equiv) hr^-1' to mg/m2/hr
                eno = ncread(data.Filename, 'NOX');
                eno = eno*1e9/1e6;
            end
            

            elon = ncread(data.Filename,'XLONG');
            elat = ncread(data.Filename, 'XLAT');
            elon = double(elon);
            elat = double(elat);
            [ xloncorn, xlatcorn ] = wrf_grid_corners( elon(:,:,1), elat(:,:,1) );

%             figure;
%             subplot(1,2,1);
%             pcolor(elon, elat, squeeze(eno(:,:,1)));
%             shading flat;
%             colorbar;

            no = nan(size(mlon,1),size(mlon,2),size(eno,3));
            for i=1:size(eno,3)
                this_nox = double(eno(:,:,i));
                this_emis = nan(size(mlon,1),size(mlon,2));
                for j=1:size(elon,1)
                    for k=1:size(elon,2)
                        in = inpolygon(mlon,mlat,xloncorn(:,j,k),xlatcorn(:,j,k));
%                         if sum(in(:)) > 0
%                             this_emis(in) = (this_nox(j,k));
%                         end
                        if sum(in(:)) > 0
                            anthroratio = sum(clc(in) < 2)/numel(clc(in)); 
                            this_emis(in) = 0;
                            if anthroratio > 0
                                indx = in & clc <2;
                                this_emis(indx) = (this_nox(j,k))/anthroratio;
                            end
                        end
                    end
                end
                no(:,:,i) = this_emis; 
%                 subplot(1,2,2);
%                 pcolor(mlon, mlat, no(:,:,1));
%                 shading flat;
%                 colorbar;
            end
        end

        function make_five_emissions_4km()
            dayofmonths = misc_flux_calculation.dayofmonths;
            five = [];
            for i=1:numel(dayofmonths)
                dayofmonth = dayofmonths(i);
                filepath = misc_flux_calculation.five_filepath;
                filename = 'wrfchemi_12z_d01';
                no = misc_flux_emis_comp.read_five_emissions_pertype(filepath, filename, dayofmonth);
                five = cat(4, five, no);
            end
            
            save(fullfile(misc_flux_calculation.five_filepath, 'five_ready'),'five');
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
                case -1
                    data = load(fullfile(misc_flux_calculation.carb_2k_filepath, 'carb_ready'));
                    emis_bio = data.carb_bio;
                    utc_times = 18:23;
                    local_time = utc_times-7;
                    emis = emis_bio(:,:,:,:);
                    save(fullfile(misc_flux_calculation.emission_filepath, 'bio_carb'), 'emis','utc_times');

                case 0
                    data = load(fullfile(misc_flux_calculation.carb_2k_filepath, 'carb_ready'));
                    emis_anthro = data.carb_anthro;
                    emis_anthro_total = data.carb_anthro_total;
                    utc_times = 18:23;
                    local_time = utc_times-7;
                    emis = emis_anthro(:,:,:,:);
                    emis_total = emis_anthro_total(:,:,:,:);
                    save(fullfile(misc_flux_calculation.emission_filepath, 'anthro_carb'), 'emis','emis_total','utc_times');

                case 1
                    data = load(fullfile(misc_flux_calculation.carb_2k_filepath, 'carb_ready'));
                    emis_anthro = data.carb_anthro;
                    emis_anthro_total = data.carb_anthro_total;
                    emis_bio = data.carb_bio;
                    utc_times = 18:23;
                    local_time = utc_times-7;
                    %emis = emis_anthro + emis_bio;
                    emis = emis_anthro_total + emis_bio;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'carb'), 'emis','utc_times');
                case 2
                    data = load(fullfile(misc_flux_calculation.carb_2k_filepath, 'carb_ready'));
                    emis_bio = data.carb_bio;
                    utc_times = 18:23;
                    local_time = utc_times-7;

                    data = load(fullfile(misc_flux_calculation.five_filepath, 'five_ready'));
                    five = data.five;
                    %utc_times = 12:23;
                    emis = five(:,:,7:end,:) + emis_bio;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F4M'), 'emis','utc_times');
                case 3
                    
                    data = load(fullfile(misc_flux_calculation.bm_filepath, 'five_1km_onroad_ready'));
                    five = data.five;
                    utc_times = data.utc_times;
                    
                    data = load(fullfile(misc_flux_calculation.carb_filepath, 'carb_ready'));
                    emis = data.carb;
                    local_time = utc_times-7;
                    emis = emis(:,:,local_time,:);

                    indx = five > 0;
                    emis(indx) = five(indx);
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F1M'), 'emis','utc_times');

                case 4
                    data = load(fullfile(misc_flux_calculation.five_filepath, 'five_ready'));
                    five = data.five;
                    utc_times = 12:23;
                    data = load(fullfile(misc_flux_calculation.iowa_filepath, 'dbsnp_ready'));
                    emis_bdsnp = data.emis_bio;
                    emis_five = five;
                    emis = data.emis_bio + five;
                    
                    soil_temp = data.bd_ts;
                    
                    data = load('./inputs/hrrr_inputs_regrid_final');
                    hrrr = data.hrrr;
                    mtemp = [];
                    for i=1:7
                        mtemp = cat(4, mtemp, hrrr(i).mtemp);
                    end
                    hrrr_times = 17:23;

                    save(fullfile(misc_flux_calculation.emission_filepath, 'F4B'), 'emis','emis_bdsnp','emis_five','soil_temp','utc_times','mtemp','hrrr_times');


                    emis = five;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'five'), 'emis','utc_times');


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
                case 6
                    data = load(fullfile(misc_flux_calculation.five_filepath, 'five_ready'));
                    five = data.five;
                    utc_times = 12:23;
                    data = load(fullfile(misc_flux_calculation.beis_filepath, 'beis_ready'));
                    emis = data.emis_agri + five;
                    emis_beis = data.emis_agri;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'F4BE'), 'emis','emis_beis','utc_times');

                case 6.5
         
                    utc_times = 12:23;
                    data = load(fullfile(misc_flux_calculation.beis_filepath, 'beis_ready'));
                    emis = data.emis_agri;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'beis_only'), 'emis','utc_times');

                case 7
                    data = load(fullfile(misc_flux_calculation.iowa_filepath, 'dbsnp_ready'));
                    emis_bio = data.emis_bio;
                    utc_times = 13:24;

                    data = load(fullfile(misc_flux_calculation.carb_2k_filepath, 'carb_ready'));
                    emis_anthro = data.carb_anthro;
                    utc_times = 18:23;
                    local_time = utc_times-7;

                    emis = emis_bio(:,:,6:end-1,:)+ emis_anthro;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'C2B'), 'emis','utc_times');

                    emis = emis_bio(:,:,6:end-1,:);
                    save(fullfile(misc_flux_calculation.emission_filepath, 'bdsnp_only'), 'emis','utc_times');

                case 8
                    data = load(fullfile(misc_flux_calculation.beis_filepath, 'beis_ready'));
                    emis_bio = data.emis_agri;

                    data = load(fullfile(misc_flux_calculation.carb_2k_filepath, 'carb_ready'));
                    emis_anthro = data.carb_anthro;
                    utc_times = 18:23;
                    local_time = utc_times-7;

                    emis = emis_bio(:,:,7:end,:)+ emis_anthro;
                    save(fullfile(misc_flux_calculation.emission_filepath, 'C2BE'), 'emis','utc_times');
            end
            
        end

        function make_flux_emis_map_comp_median_perday(dayofmonth)
            data = load(fullfile(misc_flux_calculation.output_filepath, 'recap-sjv-vert-fpcontour'));
            segs = data.segs;
            
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';
            sphi = [];
            sphi2flux = [];%nan(size(mlon));
%             sphi = nan(numel(segs), size(mlon,1),size(mlon,2));
%             sphi2flux = zeros(size(mlon));

            for i_seg=1:numel(segs)
                this_fp = segs(i_seg);
                %fprintf('Seg: %d\n', i_seg);
                if segs(i_seg).dayofmonth == dayofmonth
                for j = 1:numel(segs(i_seg).aflux)
%                     try
                    if numel(this_fp.xr_rots{j}) > 1 && ~isnan(segs(i_seg).aflux(j))
                        xr_rots = this_fp.xr_rots{j}{9};
                        yr_rots = this_fp.yr_rots{j}{9};
%                         xr_rots = this_fp.xr_rots{j}{9};
%                         yr_rots = this_fp.yr_rots{j}{9};
                        indx = ~isnan(xr_rots) & ~isnan(yr_rots);
                        xr_rots = xr_rots(indx);
                        yr_rots = yr_rots(indx);

                        in = inpolygon(mlon,mlat,xr_rots,yr_rots);
                        this_mphi = nan(size(mlon));
                        this_mphi(in) = 1;
                        %sphi = sphi + this_mphi;

                        this_mphiflux = this_mphi.*segs(i_seg).aflux_adj(j)*3600;
                        sphi2flux = cat(3, sphi2flux, this_mphiflux);
                        sphi = cat(3, sphi, this_mphi);
                    end
%                     catch err
%                     end
                end
                end
            end
            medflux = nanmedian(sphi2flux,3);
            avgflux = nanmean(sphi2flux, 3);
            sphi = nansum(sphi, 3);

%             indx  = sphi < 5;
%             flux_emis_map(indx) = nan;
            emis_comp.flux_emis_avg = avgflux;
            emis_comp.flux_emis_med = medflux;
            emis_comp.flux_sphi = sphi;

            emis_str = {'bio_carb','anthro_carb','beis_only','bdsnp_only','five','CARB','C2B','C2BE','F4M','F4B','F4BE'};
           
            [lon,lat] = meshgrid(misc_flux_calculation.mLon, misc_flux_calculation.mLat);
            mlon = lon';
            mlat = lat';

            for i_emis = 1:numel(emis_str)
                emis_name = emis_str{i_emis};
                fprintf('Working on: %s\n', emis_name);
                data = load(fullfile(misc_flux_calculation.emission_filepath, emis_name));
                emis = data.emis;
                if strcmpi(emis_name, 'F4B')
                    soil_temp = data.soil_temp;
                    hrrr_times = data.hrrr_times;
                    hrrr_mtemp = data.mtemp;
                    msoil_temp = zeros(size(mlon));
                    mtemp = zeros(size(mlon));
                    emis_five = data.emis_five;
                    five = zeros(size(mlon));
                end
                utc_times = data.utc_times;
                dayofmonths = misc_flux_calculation.dayofmonths;

                memis = [];
                n = 0;
                for i_seg=1:numel(segs)
                    this_fp = segs(i_seg);
                    if segs(i_seg).dayofmonth == dayofmonth
                    for i = 1:numel(segs(i_seg).lon)
                        % match the time
                        if numel(this_fp.xr_rots{i}) > 1 && ~isnan(segs(i_seg).aflux(i))
                            xr_rots = this_fp.xr_rots{i}{9};
                            yr_rots = this_fp.yr_rots{i}{9};
    %                         xr_rots = this_fp.xr_rots{j}{9};
    %                         yr_rots = this_fp.yr_rots{j}{9};
                            indx = ~isnan(xr_rots) & ~isnan(yr_rots);
                            xr_rots = xr_rots(indx);
                            yr_rots = yr_rots(indx);
    
                            in = inpolygon(mlon,mlat,xr_rots,yr_rots);
                            this_mphi = nan(size(mlon));
                            this_mphi(in) = 1;

                            this_time = segs(i_seg).time(i);
                            this_day = segs(i_seg).dayofmonth;
                            this_hour = round(this_time/86400*24);
        
                            indx_hour = this_hour == utc_times;
                            indx_day = this_day == dayofmonths;
    
                            this_emis = squeeze(emis(:,:,indx_hour, indx_day)).*this_mphi;
                            
                            %
                            if strcmpi(emis_name, 'F4B')
                                this_soil_temp = soil_temp(:,:,indx_hour, indx_day);
                                msoil_temp = msoil_temp + squeeze(this_soil_temp);
                                this_emis_five = emis_five(:,:,indx_hour, indx_day);
                                five = five + squeeze(this_emis_five);
    
                                indx_hrrr_hour = this_hour == hrrr_times;
                                this_mtemp = hrrr_mtemp(:,:,indx_hrrr_hour,indx_day);
                                mtemp = mtemp + squeeze(this_mtemp);
                                n = n +1;
                            end
                            
                            memis = cat(3, memis, this_emis);
                        end
                    end
                    end
                end
                emis_comp.(emis_name) = nanmean(memis, 3);
                if strcmpi(emis_name, 'F4B')
                    emis_comp.soil_temp = msoil_temp/n;
                    emis_comp.mtemp = mtemp/n;
                    emis_comp.five = five/n;
                end
            end

            emis_comp.mlon = mlon;
            emis_comp.mlat = mlat;
            data = load('./inputs/landcover-500m.mat');
            clc = data.clc;
            emis_comp.clc = clc;
            
         
            output_filename = sprintf('./output/emissions/emis_eval_%02d.mat', dayofmonth);
            save(output_filename,'emis_comp'); 

        end
        
        function regrid_emis_comp_2k(dayofmonth)
            output_filename = sprintf('./output/emissions/emis_eval_%02d.mat', dayofmonth);
            data = load(output_filename);
            emis_comp = data.emis_comp;
            mlon = emis_comp.mlon;
            mlat = emis_comp.mlat;

            filepath = misc_flux_calculation.carb_2k_filepath;
            dayofyear = datenum('2021-06-01') + dayofmonth - 1 - datenum('2021-01-01') +1;
            filepattern = sprintf('anthro/st_2k*%d*', dayofyear);
            filedir = dir(fullfile(filepath, filepattern));
            filename = fullfile(filedir.folder, filedir.name);
            
            data = ncinfo(filename);
            no = ncread(data.Filename, 'NO');
            no = squeeze(no(:,:,1,1));
            lon = zeros(size(no));

            met = shaperead(fullfile(filepath, 'state2k1/state2k1.shp'));
            clon = zeros(4, size(lon,1),size(lon,2));
            clat = zeros(4, size(lon,1),size(lon,2));
            for i=1:numel(met)
                ic = met(i).I_CELL;
                jc = met(i).J_CELL;
                clon(:, ic, jc) = met(i).X(1:4);
                clat(:, ic, jc) = met(i).Y(1:4);
            end
            data_4k.flux_emis = nan(numel(lon),1);
            data_4k.flux_emis_contour = nan(numel(lon),1);
            data_4k.carb_emis = nan(numel(lon),1);
            data_4k.carb_anthro_emis = nan(numel(lon),1);
            data_4k.bio_carb = nan(numel(lon),1);
            data_4k.beis_only = nan(numel(lon),1);
            data_4k.bdsnp_only = nan(numel(lon),1);
            data_4k.five = nan(numel(lon),1);

            data_4k.c2b_emis = nan(numel(lon),1);
            data_4k.c2be_emis = nan(numel(lon),1);
            data_4k.f4m_emis = nan(numel(lon),1);
            data_4k.f4b_emis = nan(numel(lon),1);
            data_4k.f4be_emis = nan(numel(lon),1);
            data_4k.flux_sphi = nan(numel(lon),1);
            data_4k.flux_sphi_contour = nan(numel(lon),1);
            data_4k.highway_clc = nan(numel(lon),1);
            data_4k.urban_clc = nan(numel(lon),1);
            data_4k.bio_clc = nan(numel(lon),1);
            data_4k.diaries_clc = nan(numel(lon),1);
            data_4k.res_clc = nan(numel(lon),1);
            data_4k.clc = nan(numel(lon),1);
            data_4k.soil_temp = nan(numel(lon),1);

            data_4k.mtemp = nan(numel(lon),1);
            for i=1:numel(lon)
                in = inpolygon(mlon,mlat,clon(:,i),clat(:,i));
                if sum(in(:)) > 0 && any(~isnan(emis_comp.flux_emis_avg(in)))
                    data_4k.flux_emis(i) = nanmean(emis_comp.flux_emis_avg(in));
                    data_4k.carb_emis(i) = nanmean(emis_comp.CARB(in));
                    data_4k.carb_anthro_emis(i) = nanmean(emis_comp.anthro_carb(in));

                    data_4k.bio_carb(i) = nanmedian(emis_comp.bio_carb(in));
                    data_4k.beis_only(i) = nanmedian(emis_comp.beis_only(in));
                    data_4k.bdsnp_only(i) = nanmedian(emis_comp.bdsnp_only(in));
                    data_4k.five(i) = nanmedian(emis_comp.five(in));

                    data_4k.c2b_emis(i) = nanmean(emis_comp.C2B(in));
                    data_4k.c2be_emis(i) = nanmean(emis_comp.C2BE(in));
                    data_4k.f4m_emis(i) = nanmean(emis_comp.F4M(in));
                    data_4k.f4b_emis(i) = nanmean(emis_comp.F4B(in));
                    data_4k.f4be_emis(i) = nanmean(emis_comp.F4BE(in));
                    data_4k.flux_sphi(i) = nansum(emis_comp.flux_sphi(in));
                    data_4k.soil_temp(i) = nanmean(emis_comp.soil_temp(in));
                    data_4k.mtemp(i) = nanmean(emis_comp.mtemp(in));
                    clcs = emis_comp.clc(in);
                    data_4k.highway_clc(i) = nansum(clcs == 0)/numel(clcs);
                    data_4k.urban_clc(i) = nansum(clcs == 1)/numel(clcs);
                    data_4k.bio_clc(i) = (nansum(clcs == 2))/numel(clcs);
                    data_4k.diaries_clc(i) = nansum(clcs == 12) /numel(clcs);
                    data_4k.res_clc(i) = 1 - data_4k.highway_clc(i) - data_4k.urban_clc(i) - data_4k.bio_clc(i) - data_4k.diaries_clc(i);

                    if sum(clcs == 0) >0
                        data_4k.clc(i) = 0;
                    elseif sum(clcs ==1) >0
                        data_4k.clc(i) = 1;
                    elseif sum(~isnan(clcs)) == 0
                        data_4k.clc(i) = 3;
                    else
                        data_4k.clc(i) = 2;
                    end
                end
            end

   
            data_4k.clon1 = clon(1,:)';
            data_4k.clon2 = clon(2,:)';
            data_4k.clon3 = clon(3,:)';
            data_4k.clon4 = clon(4,:)';
  
            data_4k.clat1 = clat(1,:)';
            data_4k.clat2 = clat(2,:)';
            data_4k.clat3 = clat(3,:)';
            data_4k.clat4 = clat(4,:)';
            
            
            T = struct2table(data_4k);
            output_filename = sprintf('./output/emissions/emis_eval_2km_%02d.csv', dayofmonth);
            writetable(T,output_filename,'Delimiter',',') 

        end

        function regrid_emis_comp_4k(dayofmonth)
            output_filename = sprintf('./output/emissions/emis_eval_%02d.mat', dayofmonth);
            data = load(output_filename);
            emis_comp = data.emis_comp;
            mlon = emis_comp.mlon;
            mlat = emis_comp.mlat;

            met = readtable(fullfile(misc_flux_calculation.carb_filepath, '4km.csv'));
            ic = met.I_CELL;
            jc = met.J_CELL;
            lon = met.Centroid_X;
            lat = met.Centroid_Y;
            %corner
            data = load(fullfile(misc_flux_calculation.carb_filepath,'carb_corner.mat'));
            corner_lat = data.corner_lat;
            corner_lon = data.corner_lon;
            ix = data.ix;
            iy = data.iy;
            
            clon = nan(4, numel(lon));
            clat = nan(4, numel(lon));
            
            for i=1:numel(ix)
                indx = ix(i) == ic & iy(i) == jc;
                clon(:, indx) = corner_lon(i,1:4);
                clat(:, indx) = corner_lat(i,1:4);
            end

            data_4k.flux_emis = nan(numel(lon),1);
            data_4k.flux_emis_contour = nan(numel(lon),1);
            data_4k.carb_emis = nan(numel(lon),1);
            data_4k.carb_anthro_emis = nan(numel(lon),1);
            data_4k.bio_carb = nan(numel(lon),1);
            data_4k.beis_only = nan(numel(lon),1);
            data_4k.bdsnp_only = nan(numel(lon),1);
            data_4k.five = nan(numel(lon),1);

            data_4k.c2b_emis = nan(numel(lon),1);
            data_4k.c2be_emis = nan(numel(lon),1);
            data_4k.f4m_emis = nan(numel(lon),1);
            data_4k.f4b_emis = nan(numel(lon),1);
            data_4k.f4be_emis = nan(numel(lon),1);
            data_4k.flux_sphi = nan(numel(lon),1);
            data_4k.flux_sphi_contour = nan(numel(lon),1);
            data_4k.highway_clc = nan(numel(lon),1);
            data_4k.urban_clc = nan(numel(lon),1);
            data_4k.bio_clc = nan(numel(lon),1);
            data_4k.diaries_clc = nan(numel(lon),1);
            data_4k.res_clc = nan(numel(lon),1);
            data_4k.clc = nan(numel(lon),1);
            data_4k.soil_temp = nan(numel(lon),1);

            data_4k.mtemp = nan(numel(lon),1);
            for i=1:numel(lon)
                in = inpolygon(mlon,mlat,clon(:,i),clat(:,i));
                if sum(in(:)) > 0 && any(~isnan(emis_comp.flux_emis_avg(in)))
                    data_4k.flux_emis(i) = nanmean(emis_comp.flux_emis_avg(in));
                    data_4k.carb_emis(i) = nanmean(emis_comp.CARB(in));
                    data_4k.carb_anthro_emis(i) = nanmean(emis_comp.anthro_carb(in));

                    data_4k.bio_carb(i) = nanmean(emis_comp.bio_carb(in));
                    data_4k.beis_only(i) = nanmean(emis_comp.beis_only(in));
                    data_4k.bdsnp_only(i) = nanmean(emis_comp.bdsnp_only(in));
                    data_4k.five(i) = nanmean(emis_comp.five(in));

                    data_4k.c2b_emis(i) = nanmean(emis_comp.C2B(in));
                    data_4k.c2be_emis(i) = nanmean(emis_comp.C2BE(in));
                    data_4k.f4m_emis(i) = nanmean(emis_comp.F4M(in));
                    data_4k.f4b_emis(i) = nanmean(emis_comp.F4B(in));
                    data_4k.f4be_emis(i) = nanmean(emis_comp.F4BE(in));
                    data_4k.flux_sphi(i) = nansum(emis_comp.flux_sphi(in));
                    data_4k.soil_temp(i) = nanmean(emis_comp.soil_temp(in));
                    data_4k.mtemp(i) = nanmean(emis_comp.mtemp(in));
                    clcs = emis_comp.clc(in);
                    data_4k.highway_clc(i) = nansum(clcs == 0)/numel(clcs);
                    data_4k.urban_clc(i) = nansum(clcs == 1)/numel(clcs);
                    data_4k.bio_clc(i) = (nansum(clcs == 2))/numel(clcs);
                    data_4k.diaries_clc(i) = nansum(clcs == 12) /numel(clcs);
                    data_4k.res_clc(i) = 1 - data_4k.highway_clc(i) - data_4k.urban_clc(i) - data_4k.bio_clc(i) - data_4k.diaries_clc(i);

                    if sum(clcs == 0) >0
                        data_4k.clc(i) = 0;
                    elseif sum(clcs ==1) >0
                        data_4k.clc(i) = 1;
                    elseif sum(~isnan(clcs)) == 0
                        data_4k.clc(i) = 3;
                    else
                        data_4k.clc(i) = 2;
                    end
                end
            end

   
            data_4k.clon1 = clon(1,:)';
            data_4k.clon2 = clon(2,:)';
            data_4k.clon3 = clon(3,:)';
            data_4k.clon4 = clon(4,:)';
  
            data_4k.clat1 = clat(1,:)';
            data_4k.clat2 = clat(2,:)';
            data_4k.clat3 = clat(3,:)';
            data_4k.clat4 = clat(4,:)';
            
            
            T = struct2table(data_4k);
            output_filename = sprintf('./output/emissions/emis_eval_4km_%02d_mean.csv', dayofmonth);
            writetable(T,output_filename,'Delimiter',',') 

        end

        function make_flux_emis_map_comp_assemble()
            dayofmonths = [3,8,9,13,15,16,22];
            for i =1:numel(dayofmonths)
                fprintf("working on %d\n", dayofmonths(i))
                misc_flux_emis_comp.make_flux_emis_map_comp_median_perday(dayofmonths(i));
                %misc_flux_emis_comp.regrid_emis_comp_2k(dayofmonths(i));
                misc_flux_emis_comp.regrid_emis_comp_4k(dayofmonths(i));
            end
        end

    end
end