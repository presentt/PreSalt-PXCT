% Visualize microfossil PXCT data
% Ted Present, Feb 2024

data_folder = "E:\PXCT\PXCT_data\test" ;
tif_file_slice_list = dir(fullfile(data_folder,"*.tif"));
firstTiff = Tiff(fullfile(data_folder,tif_file_slice_list(1).name));

% convert tif slices to multipage tif volume
tag_struct.ImageLength = firstTiff.getTag("ImageLength"); %todo: loop through Tiff.getTagNames
tag_struct.ImageWidth = firstTiff.getTag("ImageWidth");
tag_struct.SampleFormat = Tiff.SampleFormat.IEEEFP; %firstTiff.getTag("SampleFormat");
tag_struct.Photometric = firstTiff.getTag("Photometric");
tag_struct.BitsPerSample = 32; %firstTiff.getTag("BitsPerSample");
tag_struct.SamplesPerPixel = 1; % firstTiff.getTag("SamplesPerPixel");
tag_struct.Compression = Tiff.Compression.None; % firstTiff.getTag("Compression");
tag_struct.PlanarConfiguration = firstTiff.getTag("PlanarConfiguration");
tag_struct.ResolutionUnit = firstTiff.getTag("ResolutionUnit");
tag_struct.XResolution = firstTiff.getTag("XResolution");
tag_struct.YResolution = firstTiff.getTag("YResolution");
tag_struct.Orientation = firstTiff.getTag("Orientation");

%%
firstTiffImg = single(read(firstTiff));

%%
metadata = readMetadata("E:\PXCT\PXCT_data\metadata\S4_2\TIFF_delta_FBP_ram-lak_freqscl_1.00_cutoffs.txt");
metadata.photon_energy = 6.2; % keV       TODO: append these experimental info to the metadata structure
metadata.photon_flux = 8e5; % photons/s, for 5.7mm working distance
metadata.field_of_view = [80 83]; % um
metadata.step_size = 2.5; % um
metadata.exposure_time = 0.05; % s
metadata.projections = 640;
metadata.mu = NaN; % TODO: pull from Henke 1993 for quartz
metadata.rho = 2.65; % g/cc, bulk density of specimen assuming quartz; quartz GLASS is 2.2g/cc, and this is what the rendering picked up, too; also what Diaz uses for their capillary
metadata.N0 = 3.480e+06; % photons/um^2 % TODO: calculate from flux, FOV, step size, exposure time, projections
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

densityImg = (firstTiffImg.*(metadata.high_cutoff - metadata.low_cutoff)/(2^16-1) + metadata.low_cutoff).*metadata.factor_edensity.*3.08;

% s8_2 but not great crop...
% xc = 1060;
% yc = 1010;
% r = 850;

xc = 1025;
yc = 950;
r = 875;

figure;
    subplot(1,2,1)
        histogram(densityImg)
    subplot(1,2,2)
        imshow(densityImg)
        viscircles([xc yc],r,'EdgeColor','r', 'LineWidth', 1);

%%
[xDim,yDim] = size(densityImg);
[xx,yy] = meshgrid(1:yDim,1:xDim);
mask = false(xDim,yDim);
mask = mask | hypot(xx - xc, yy - yc) < r;

croppedImg = densityImg.*mask;
croppedImg = imcrop(croppedImg,[xc-r yc-r 2.*r 2.*r]);

figure;
    imshow(croppedImg)

%%
figure;
    h = histogram(croppedImg);
    h.Normalization = "pdf";
    hold on;
    
    guess.mu = [0; 0; 1.2; 1.5; 1.7; 2.3; 4];
    guess.Sigma(1,1,1) = 1e-4;
    guess.Sigma(1,1,2) = 0.5;
    guess.Sigma(1,1,3) = 0.5;
    guess.Sigma(1,1,4) = 0.5;
    guess.Sigma(1,1,5) = 0.5;
    guess.Sigma(1,1,6) = 0.5;
    guess.Sigma(1,1,7) = 0.5;
    guess.ComponentProportion = [0.05 0.04 0.3 0.1 0.1 0.4 0.01];
    endmembers = length(guess.mu);
    fitopts = statset('Display','final','MaxIter',500,'TolFun',1e-6);

    GMModel = fitgmdist(reshape(croppedImg,[],1),...
        endmembers,'RegularizationValue',0.00001, ...
        'Start',guess, 'Options',fitopts);
    disp(GMModel.mu);
    gmPDF = @(x) arrayfun(@(x0) pdf(GMModel,x0),x);
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    for p = 1:endmembers
        plot(h.BinEdges,...
            normpdf(h.BinEdges,GMModel.mu(p),sqrt(GMModel.Sigma(:,:,p))).*GMModel.ComponentProportion(p),...
            'LineWidth',1);
    end
    hold off;

%%
% overwrite tag_struct with file-specific metadata
tag_struct.ImageLength = size(croppedImg,1);
tag_struct.ImageWidth = size(croppedImg,2);
tag_struct.ResolutionUnit = 3; % cm
tag_struct.XResolution = metadata.pixel_size;
tag_struct.YResolution = metadata.pixel_size;

tiff_vol = Tiff('S4_2_density.tif', 'w8'); % prepare a BigTIFF for writing
for i = 1:length(tif_file_slice_list)
    ith_tiff = Tiff(fullfile(data_folder,tif_file_slice_list(i).name));
    ith_tiffImg = single(read(ith_tiff));
    ith_tiffImg = (ith_tiffImg.*(metadata.high_cutoff - metadata.low_cutoff)/(2^16-1) + metadata.low_cutoff).*metadata.factor_edensity.*3.08;
    ith_tiffImg = ith_tiffImg.*mask;
    ith_tiffImg = imcrop(ith_tiffImg,[xc-r yc-r 2.*r 2.*r]);
    setTag(tiff_vol,tag_struct);
    write(tiff_vol,ith_tiffImg);
    writeDirectory(tiff_vol);
end
close(tiff_vol)
clear tiff_vol;

%%
%reducedImg = imresize3(tiffreadVolume('S4_2_cropped.tif'), 0.25);
%%

bim = blockedImage(tiffreadVolume('testDensityVol.tif'), BlockSize=[300 300 300]);

%%
figure;
subplot(2,1,1)
    h = histogram(tiffreadVolume('testDensityVol.tif'));
    h.Normalization = "probability";
    range = [h.BinEdges(1) h.BinEdges(end)];
    xlim(range)

%%
voldisp = volshow(bim);

%%
alpha = [0 0 0 .02 .8 1 .8 0.02 .005 .005];
% color = [0 0 0;
%         200 140 75;
%         231 208 141;
%         231 208 141;
%         255 255 255;
%         255 255 255] ./ 255;
color = [255 255 255;
        255 255 255;
        255 255 255;
        200 140 75;
        200 140 75;
        200 140 75;
        200 140 75;
        200 140 75;
        255 255 255;
        255 255 255] ./ 255;
% intensity = [0 11300 11300.1 11800 11801 40000]; % S4_2
% intensity = [0 13500 13500.1 18000 18001 40000]; % S8_1
% intensity = [0 15500 15500.1 21500 21501 40000]; % S8_2
% intensity = [0 14000 14000.1 19000 19000.1 25000]; % test cropped
p = 5; % component to render
intensity = [range(1) 0 GMModel.mu(p)-3.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p) GMModel.mu(p)+sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+3.*sqrt(GMModel.Sigma(:,:,p)) range(end)]; % test density
%intensity = [0 1.1 1.11 1.9 1.91 h.BinEdges(end)]; % test cropped
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);

voldisp.Alphamap = alphamap;
voldisp.Colormap = colormap;

subplot(2,1,2)
    plot(intensity,alpha,'-o'); hold on;
    plot(intensity,color(:,3),'-x'); hold off;
    xlim(range)

%voldisp.RenderingStyle ="GradientOpacity";

%%
% use DataReadFinished event to record a video...

%% functions
function metadataStruct = readMetadata(filename)
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [1, 5];
    opts.Delimiter = ["#", "="];
    opts.VariableNames = ["Var1", "name", "value"];
    opts.SelectedVariableNames = ["name", "value"];
    opts.VariableTypes = ["string", "string", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Var1", "name"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "name"], "EmptyFieldRule", "auto");
    raw_metadata = readtable(filename, opts);
    
    raw_metadata.Properties.RowNames = raw_metadata.name;

    metadataStruct.low_cutoff = raw_metadata.value("low_cutoff");
    metadataStruct.high_cutoff = raw_metadata.value("high_cutoff");
    metadataStruct.factor = raw_metadata.value("factor");
    metadataStruct.pixel_size = raw_metadata.value("pixel size");
    metadataStruct.factor_edensity = raw_metadata.value("factor_edensity");
end