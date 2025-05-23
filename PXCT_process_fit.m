% Visualize microfossil PXCT data
% Ted Present, Feb 2024

data_folder = "E:\PXCT\PXCT_data\S4_2" ;
tif_file_slice_list = dir(fullfile(data_folder,"*.tif"));

%%
firstTiff = Tiff(fullfile(data_folder,tif_file_slice_list(1).name));
firstTiffImg = read(firstTiff);

%%
metadata = readMetadata("E:\PXCT\PXCT_data\metadata\S4_2\TIFF_delta_FBP_ram-lak_freqscl_1.00_cutoffs.txt");
metadata.photon_energy = 6.2; % keV       TODO: append these experimental info to the metadata structure
metadata.photon_flux = 8e5; % photons/s, for 5.7mm working distance
metadata.field_of_view = [80 83]; % um
metadata.step_size = 2.5; % um
metadata.exposure_time = 0.05; % s
metadata.projections = 710; %TODO: pull from dose estimate file
metadata.mu = NaN; % TODO: pull from Henke 1993 for quartz
metadata.rho = 2.65; % g/cc, bulk density of specimen assuming quartz; quartz GLASS is 2.2g/cc, and this is what the rendering picked up, too; also what Diaz uses for their capillary
metadata.N0 = 3.461e6; % photons/um^2 % TODO: calculate from flux, FOV, step size, exposure time, projections
metadata.dose = metadata.mu*metadata.projections*metadata.N0*metadata.photon_energy/metadata.rho; % Gy, total dose imparted to specimen, Diaz et al., J. Struct. Bio., 2015
        % Scan = 
        % Measured photons per frame = 2.249e+07
        % N_0 used for imaging for one projection (I*1e-12) = 3.480e+06 photons/micron^2
        % num_proj = 646
        % hv = 9.93478e-16
        % 
        % D = mu*I*hv*num_proj/rho 
        % D_with_gas = D/setup_transmission 
        % D_with_overhead = D_with_gas*(1+overhead) 
        % 
        % If mu = 1.618e+04 m^-1, rho = 1.000e+03 kg/m^3, setup_transmission = 0.980, then:
        % Accounting for experiment transmission, D_with_gas = 3.687e+07 Gy = 36.866 MGy
        % 
        % Not accounting for experiment transmission, Absorbed photons = 6.839e+12 
        % 
        % Not accounting for experiment transmission, Absorbed dose = 1.267e+07 
        % 
        % Accounting for experiment transmission, Absorbed dose = 1.293e+07 
        % 
        % Accounting for experiment transmission and overhead, Absorbed dose = 1.293e+07 
        % 
        % asize = 350, pixel size = 39.7073 nm

%%

edensityImg = edensity(firstTiffImg,metadata);
densityImg = density(edensityImg, 60.08, 30); % calculate material density assuming quartz

xc = 1025;
yc = 975;
r = 950;

figure;
    subplot(1,2,1)
        histogram(edensityImg)
    subplot(1,2,2)
        imshow(edensityImg)
        viscircles([xc yc],r,'EdgeColor','r', 'LineWidth', 1);

%%
[xDim,yDim] = size(edensityImg);
[xx,yy] = meshgrid(1:yDim,1:xDim); % convert to single?
mask = false(xDim,yDim);
mask = mask | hypot(xx - xc, yy - yc) < r;
clear xx yy

croppedImg = edensityImg; %.*mask;
croppedImg = imcrop(croppedImg,[xc-r yc-r 2.*r 2.*r]);

figure;
    subplot(1,2,1)
        histogram(croppedImg)
    subplot(1,2,2)
        imshow(croppedImg)

%%
tag_struct.ImageLength = firstTiff.getTag("ImageLength"); %todo: loop through Tiff.getTagNames
tag_struct.ImageWidth = firstTiff.getTag("ImageWidth");
tag_struct.SampleFormat = firstTiff.getTag("SampleFormat");
tag_struct.Photometric = firstTiff.getTag("Photometric");
tag_struct.BitsPerSample = firstTiff.getTag("BitsPerSample");
tag_struct.SamplesPerPixel = firstTiff.getTag("SamplesPerPixel");
tag_struct.Compression = firstTiff.getTag("Compression");
tag_struct.PlanarConfiguration = firstTiff.getTag("PlanarConfiguration");
tag_struct.ResolutionUnit = firstTiff.getTag("ResolutionUnit");
tag_struct.XResolution = firstTiff.getTag("XResolution");
tag_struct.YResolution = firstTiff.getTag("YResolution");
tag_struct.Orientation = firstTiff.getTag("Orientation");

% overwrite tag_struct with file-specific metadata
% tag_struct.SampleFormat = Tiff.SampleFormat.IEEEFP; tag_struct.BitsPerSample = 32;
tag_struct.SampleFormat = Tiff.SampleFormat.UInt; tag_struct.BitsPerSample = 16;
tag_struct.SamplesPerPixel = 1; 
tag_struct.Compression = Tiff.Compression.None;

tag_struct.ImageLength = size(croppedImg,1);
tag_struct.ImageWidth = size(croppedImg,2);
tag_struct.ResolutionUnit = 3; % cm
tag_struct.XResolution = metadata.pixel_size;
tag_struct.YResolution = metadata.pixel_size;

tiff_vol = Tiff('S4_2_edensity.tif', 'w8'); % prepare a BigTIFF for writing
for i = 1:length(tif_file_slice_list)
    ith_tiff = Tiff(fullfile(data_folder,tif_file_slice_list(i).name));
    ith_tiffImg = single(read(ith_tiff));
    ith_tiffImg = edensity(ith_tiffImg,metadata);
    % ith_tiffImg = density(ith_tiffImg, 60.08, 30); % assuming quartz
    %ith_tiffImg = ith_tiffImg.*mask;
    ith_tiffImg = imcrop(ith_tiffImg,[xc-r yc-r 2.*r 2.*r]);
    setTag(tiff_vol,tag_struct);
    write(tiff_vol,ith_tiffImg);
    writeDirectory(tiff_vol);
end
close(tiff_vol)
clear tiff_vol;

%%
tif_vol_file = 'S4_2_edensity.tif';
bim = blockedImage(tiffreadVolume(tif_vol_file), BlockSize=[100 100 100]);
mbim = makeMultiLevel3D(bim);
clear bim;

%% compute histogram from full tif volume
edges = 0:1e4;

[hbim, ~] = apply(mbim, ...
             @(bs)histcounts(bs.Data,edges));
hdata = gather(hbim);
hdata = sum(hdata,3);
hdata = sum(hdata,1);
hdata = reshape(hdata,length(edges)-1,[]);
hdata = sum(hdata,2)';

writematrix([edges./1e4; [hdata, NaN]],[tif_vol_file(1:end-4),'_histcounts.csv']);

%%
hdata = readmatrix('S4_1_edensity_histcounts.csv');
edges = hdata(1,:).*1e4;
hdata = hdata(2,1:end-1);

%%
stackedhistdata = [readmatrix('S4_1_edensity_histcounts.csv'); ...
                   readmatrix('S4_2_edensity_histcounts.csv'); ...
                   readmatrix('S8_1_edensity_histcounts.csv'); ...
                   readmatrix('S8_2_edensity_histcounts.csv')];
stackedhistdata = stackedhistdata([1 2:2:end], 2:end);
%%
hist_fig = figure;
for i = 2:size(stackedhistdata,1)
    hist_subplot(i-1) = subplot(size(stackedhistdata,1)-1, 1, i-1);
    stackedhist(i-1) = histogram('BinEdges',stackedhistdata(1, :),...
                            'BinCounts',stackedhistdata(i, 1:end-1));
    stackedhist(i-1).Normalization = "probability";
    stackedhist(i-1).DisplayStyle = "stairs";
    %ylim([0 7]);
    hist_subplot(i-1).XAxis.Visible = false;
    hist_subplot(i-1).YAxis.Visible = false;
end

hist_subplot(4).XAxis.Visible = true;
hist_fig.Position(3) = 400; % width
hist_fig.Position(4) = 800; % height

%% randomly sample densities and save
tifdata = tiffreadVolume(tif_vol_file);
tifdata = reshape(tifdata, [], 1);

tifsample = datasample(tifdata, 1e6, 'Replace',false); % randomly subsample tif volume to fit in memory
tifsmpsngl = single(tifsample)./1e4; % convert to floating precision electron densities
writematrix(tifsmpsngl,[tif_vol_file(1:end-4),'_sample.csv']);
clear tifdata tifsample tifsmpsngl;

%% fit sampled distribution of densities
tifsmpsngl = readmatrix([tif_vol_file(1:end-4),'_sample.csv']);

tifsmpsngl = tifsmpsngl(tifsmpsngl>0.2); % don't fit pores

clear gmPDF fitopts endmembers GMModel;

% try fitting without a starting guess
% might have to run a few times to get endmembers correctly
endmembers = 5;
fitopts = statset('Display','final','MaxIter',500,'TolFun',1e-6);
GMModel = fitgmdist(tifsmpsngl,...
    endmembers, ...
    'Options',fitopts, 'SharedCovariance', false, 'CovarianceType','diagonal');

% save([tif_vol_file(1:end-4),'_GMModel.mat'],"GMModel");
disp(GMModel.mu);
gmPDF = @(x) arrayfun(@(x0) pdf(GMModel,x0),x);

figure;
    h = histogram(tifsmpsngl);
    h.Normalization = "pdf";
    hold on;
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    for p = 1:endmembers
        plot(h.BinEdges,...
            normpdf(h.BinEdges,GMModel.mu(p),sqrt(GMModel.Sigma(:,:,p))).*GMModel.ComponentProportion(p),...
            'LineWidth',1);
    end
    %xlim([.01 6.55])
    ylim([0 10])
    hold off;
    xlabel('electron density');
    ylabel('probability (%)')

%%
figure;
    hS41 = histogram(readmatrix('S4_1_edensity_sample.csv'), 'DisplayName', 'S4-1');
    hS41.Normalization = 'pdf';
    hS41.DisplayStyle = 'stairs';
    hS41.morebins;
    hold on;
    hS42 = histogram(readmatrix('S4_2_edensity_sample.csv'), 'DisplayName', 'S4-2');
    hS42.Normalization = 'pdf';
    hS42.DisplayStyle = 'stairs';
    hS42.morebins;
    hS81 = histogram(readmatrix('S8_1_edensity_sample.csv'), 'DisplayName', 'S8-1');
    hS81.Normalization = 'pdf';
    hS81.DisplayStyle = 'stairs';
    hS81.morebins;
    hS82 = histogram(readmatrix('S8_2_edensity_sample.csv'), 'DisplayName', 'S8-2');
    hS82.Normalization = 'pdf';
    hS82.DisplayStyle = 'stairs';
    hS81.morebins;
    hold off;
    legend();
    xlim([0.01 3]);
    %ylim([0 6]);
    %set(gca,'YScale','log');
    xlabel('electron density');
    ylabel('probability (%)');
