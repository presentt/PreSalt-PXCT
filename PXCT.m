% Visualize microfossil PXCT data
% Ted Present, Feb 2024

data_folder = "E:\PXCT\PXCT_data\test" ;
tif_file_slice_list = dir(fullfile(data_folder,"*.tif"));
firstTiff = Tiff(fullfile(data_folder,tif_file_slice_list(end).name));

% convert tif slices to multipage tif volume
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

%%
firstTiffImg = read(firstTiff);

% s8_2 but not great crop...
% xc = 1060;
% yc = 1010;
% r = 850;
xc = 1040;
yc = 1010;
r = 750;

figure;
subplot(1,2,1)
    imhist(firstTiffImg)
subplot(1,2,2)
    imshow(firstTiffImg)
    viscircles([xc yc],r,'EdgeColor','r', 'LineWidth', 1);

%%
[xDim,yDim] = size(firstTiffImg);
[xx,yy] = meshgrid(1:yDim,1:xDim);
mask = false(xDim,yDim);
mask = mask | hypot(xx - xc, yy - yc) <= r;
mask = uint16(mask);

croppedImg = firstTiffImg.*mask;
croppedImg = imcrop(croppedImg,[xc-r yc-r 2.*r 2.*r]);

imshow(croppedImg)

%%
% overwrite tag_struct with file-specific metadata
tag_struct.ImageLength = size(croppedImg,1);
tag_struct.ImageWidth = size(croppedImg,2);
tag_struct.ResolutionUnit = 3; % cm
tag_struct.XResolution = 3.484412e-06; % todo: read in metadata file
tag_struct.YResolution = tag_struct.XResolution;

tiff_vol = Tiff('testCroppedVol.tif', 'w8'); % prepare a BigTIFF for writing
for i = 1:length(tif_file_slice_list)
    ith_tiff = Tiff(fullfile(data_folder,tif_file_slice_list(i).name));
    ith_tiffImg = read(ith_tiff);
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

bim = blockedImage(tiffreadVolume('testCroppedVol.tif'), BlockSize=[300 300 300]);

imhist(tiffreadVolume('testCroppedVol.tif'));

%%
voldisp = volshow(bim);

%%
alpha = [0 0 1 1 0.002 0.002];
color = [0 0 0;
        200 140 75;
        231 208 141;
        231 208 141;
        255 255 255;
        255 255 255] ./ 255;
% intensity = [0 11300 11300.1 11800 11801 40000]; % S4_2
% intensity = [0 13500 13500.1 18000 18001 40000]; % S8_1
% intensity = [0 15500 15500.1 21500 21501 40000]; % S8_2
intensity = [0 14000 14000.1 19000 19000.1 25000]; % test
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);

voldisp.Alphamap = alphamap;
voldisp.Colormap = colormap;

%voldisp.RenderingStyle ="GradientOpacity";

%%
% use DataReadFinished event to record a video...

%% functions
test = readMetadata("E:\PXCT\PXCT_data\metadata\S4_2\TIFF_delta_FBP_ram-lak_freqscl_1.00_cutoffs.txt");

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