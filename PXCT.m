% Visualize microfossil PXCT data
% Ted Present, Feb 2024

data_folder = "./PXCT_data/test" ;
tif_file_slice_list = dir(fullfile(data_folder,"*.tif"));

% convert tif slices to multipage tif volume
firstTiff = Tiff(fullfile(data_folder,tif_file_slice_list(1).name));
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


tiff_vol = Tiff('testVol.tif', 'w8'); % prepare a BigTIFF for writing
for i = 1:length(tif_file_slice_list)
    ith_tiff = Tiff(fullfile(data_folder,tif_file_slice_list(i).name));
    setTag(tiff_vol,tag_struct);
    write(tiff_vol,read(ith_tiff));
    writeDirectory(tiff_vol);
end
close(tiff_vol)


%%
bim = blockedImage(tiffreadVolume('S4_2.tif'), BlockSize=[500 500 500]);

%%
imhist(read(firstTiff))

voldisp = volshow(bim);

%voldisp.DataLimits(1) = 14000; % crop out low-count edges and pores
%voldisp.DataLimits = [14000 17000]; % organics?
voldisp.DataLimits = [20000 25000]; % chert?
voldisp.Alphamap(end) = 0;  

%voldisp.RenderingStyle ="GradientOpacity";
%voldisp.Alphamap(1:100) = 0;


%% functions

