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
tif_vol_file = 'S4_1_edensity.tif';
bim = blockedImage(tiffreadVolume(tif_vol_file), BlockSize=[100 100 100]);
mbim = makeMultiLevel3D(bim);
clear bim;

%%
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
%%
% fit phases as gaussians
X = hdata(2000:end)./sum(hdata(2000:end));
xvals = edges(2000:end-1)';

% guess/define initial params
zero_amp = X(1);
void_amp = X(2);
void_sd = 200;
n_comps = 3; % number of non-void components to fit as gaussians
guess_gaus_amps = 1e-6 .* ones(1, n_comps);
guess_guas_mus = [4000 6000 8000];
guess_guas_sds = 200.*ones(1, n_comps);

guess_params = [zero_amp(:); void_amp(:); void_sd(:); ...
    guess_gaus_amps(:); guess_guas_mus(:); guess_guas_sds(:)];
mixGausMdl = @(b,x) [zero_amp; zeros(length(x)-1, 1)] + ...
                    void_amp .* exp(-(x(:, 1)).^2 ./ (2 .*b(3).^2)) + ...
                    sum( ...
                        b(4:4+n_comps-1)' .* exp(-(x(:, 1) - b(4+n_comps:4+2*n_comps-1)').^2 ./...
                        (2.*b(4+2*n_comps:4+3*n_comps-1)'.^2)), ...
                    2);

lower_bounds = zeros(size(guess_params));
high_bounds = [];

% Poisson log-likelihood; fmincon finds the minimum of this function. Thus, 
% will find the max of the poisson loglikelihood.
poissonLogLikelihood = @(r,lambda) r(:)'*log(lambda(:)) - sum(lambda);
% loss functions for fitting
nllFcn = @(params) -poissonLogLikelihood(X, mixGausMdl(params, xvals(:)));

test = mixGausMdl(guess_params, xvals(:));
test2 = poissonLogLikelihood(X, test);
test3 = nllFcn(guess_params);

opts = optimset('TolX', 1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4, 'Display', 'on');
[paramEst, ~, exitFlag, outStruct] = fmincon(nllFcn, guess_params, [], [], [], [], lower_bounds, high_bounds, [],opts);

histGMModel = gmdistribution(paramEst(4+n_comps:4+2*n_comps-1), ...
                            reshape(paramEst(4+2*n_comps:4+3*n_comps-1).^2, [1 1 n_comps]), ...
                            paramEst(4:4+n_comps-1));
gmPDF = @(x) arrayfun(@(x0) pdf(histGMModel,x0),x);

figure;
    h = histogram('BinEdges',edges,'BinCounts',hdata);
    h.Normalization = "pdf";
    hold on;
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    plot(h.BinEdges,mixGausMdl(paramEst, h.BinEdges), 'LineWidth', 3);
    for p = 1:histGMModel.NumComponents
        plot(h.BinEdges,...
            normpdf(h.BinEdges, histGMModel.mu(p), ...
                sqrt(histGMModel.Sigma(:,:,p))).*histGMModel.ComponentProportion(p),...
                'LineWidth',1);
    end
    plot(h.BinEdges, void_amp .* exp(-(h.BinEdges).^2 /(2*paramEst(3)^2)))
    xlim([2 1e4])
    ylim([0 1e-3])
    hold off;

%%
% fit phases as gaussians, modified from
% https://www.mathworks.com/matlabcentral/answers/720430-how-do-i-use-fitgmdist-with-a-set-of-data-that-is-not-in-a-histogram-format
X = hdata(2:end)./sum(hdata(2:end));
xvals = edges(2:end-1)';

% fit a half-gaussian for the pore space, and then several gaussians
mixGausMdl = @(b,x) b(1) .* exp(-(x(:, 1)).^2 /(2*b(2)^2)) + ...
                    b(3) .* exp(-(x(:, 1) - b(4)).^2 /(2*b(5)^2)) + ...
                    b(6) .* exp(-(x(:, 1) - b(7)).^2 /(2*b(8)^2)) + ...
                    b(9) .* exp(-(x(:, 1) -b(10)).^2/(2*b(11)^2)) + ...
                    b(12) .*exp(-(x(:, 1) -b(13)).^2/(2*b(14)^2));
%               amp,   mu,  sd 
guess_params = [500,       200, ...
                0.1,    1, 200, ...
                0.1, 4000, 200, ...
                0.1, 6000, 200, ...
                0.2, 7500, 500];
lower_bounds = zeros(size(guess_params));
high_bounds = [];

% Poisson log-likelihood; fmincon finds the minimum of this function. Thus, 
% will find the max of the poisson loglikelihood.
poissonLogLikelihood = @(r,lambda) r(:)'*log(lambda(:)) - sum(lambda);
% loss functions for fitting
nllFcn = @(params) -poissonLogLikelihood(X, mixGausMdl(params, xvals(:)));

test = mixGausMdl(guess_params, xvals(:));
test2 = poissonLogLikelihood(X, test);
test3 = nllFcn(guess_params);

opts = optimset('TolX', 1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4, 'Display', 'on');
[paramEst, ~, exitFlag, outStruct] = fmincon(nllFcn, guess_params, [], [], [], [], lower_bounds, high_bounds, [],opts);

histGMModel = gmdistribution(paramEst(4:3:end)', ...
                            reshape(paramEst(5:3:end).^2, [1 1 4]), ...
                            paramEst(3:3:end));
gmPDF = @(x) arrayfun(@(x0) pdf(histGMModel,x0),x);

figure;
    h = histogram('BinEdges',edges,'BinCounts',hdata);
    h.Normalization = "pdf";
    hold on;
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    for p = 1:histGMModel.NumComponents
        plot(h.BinEdges,...
            normpdf(h.BinEdges, histGMModel.mu(p), ...
                sqrt(histGMModel.Sigma(:,:,p))).*histGMModel.ComponentProportion(p),...
                'LineWidth',1);
    end
    plot(h.BinEdges, paramEst(1) .* exp(-(h.BinEdges).^2 /(2*paramEst(2)^2)))
    xlim([2 1e4])
    ylim([0 1e-3])
    hold off;

%% try least squares fit
paramEst = lsqcurvefit(@(params,x) mixGausMdl(params,x), guess_params, xvals, X');

histGMModel = gmdistribution(paramEst(4:3:end)', ...
                            reshape(paramEst(5:3:end).^2, [1 1 4]), ...
                            paramEst(3:3:end));
gmPDF = @(x) arrayfun(@(x0) pdf(histGMModel,x0),x);

figure;
    h = histogram('BinEdges',edges,'BinCounts',hdata);
    h.Normalization = "pdf";
    hold on;
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    for p = 1:histGMModel.NumComponents
        plot(h.BinEdges,...
            normpdf(h.BinEdges, histGMModel.mu(p), ...
                sqrt(histGMModel.Sigma(:,:,p))).*histGMModel.ComponentProportion(p),...
                'LineWidth',1);
    end
    plot(h.BinEdges, paramEst(1) .* exp(-(h.BinEdges).^2 /(2*paramEst(2)^2)))
    xlim([2 1e4])
    ylim([0 1e-3])
    hold off;

%% fit distribution of densities (not binned density counts)
tifdata = tiffreadVolume('S4_1_edensity.tif');
tifdata = reshape(tifdata, [], 1);
%%
tifsample = datasample(tifdata, 1e6, 'Replace',false); % randomly subsample tif volume to fit in memory
tifsmpsngl = single(tifsample)./1e4;
guess.mu = [0; .1; .4; .6; .75];
guess.Sigma(1,1,1) = .1;
guess.Sigma(1,1,2) = .3;
guess.Sigma(1,1,3) = .3;
guess.Sigma(1,1,4) = .3;
guess.Sigma(1,1,5) = .3;
guess.ComponentProportion = [0.5 .1 0.1 0.1 0.2];

endmembers = length(guess.mu);
fitopts = statset('Display','final','MaxIter',500,'TolFun',1e-9);
GMModel = fitgmdist(tifsmpsngl,...
    endmembers,'RegularizationValue',1e-9, ...
    'Start',guess, 'Options',fitopts);
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
    xlim([.01 1])
    ylim([0 10])
    hold off;
    xlabel('electron density');
    ylabel('probability (%)')


% Error using zeros
% Requested 6018274080x5 (112.1GB) array exceeds maximum array size preference (27.7GB). This might cause MATLAB
% to become unresponsive.
% 
% Error in wdensity (line 17)
%     log_lh = zeros(n,k,'like',X);
% 
% Error in gmcluster>gmcluster_learn (line 271)
%         log_lh=wdensity(X,S.mu, S.Sigma, S.PComponents, SharedCov, CovType);
% 
% Error in gmcluster (line 197)
%         [S0,ll0,  optimInfo0] = gmcluster_learn...
% 
% Error in gmdistribution.fit (line 102)
%         gmcluster(X,k,start,reps,CovType,SharedCov,RegV,options,probtol);
% 
% Error in fitgmdist (line 135)
% gm = gmdistribution.fit(X,k,varargin{:});
%%
voldisp = volshow(mbim);

%%
GMModel = histGMModel;
alpha = [0 0 .02 .8 1 .8 0.02 .005 .005];
colors = [0 0 0;
          200 140 75;
          231 208 141;
          255 255 255;
          128 0 0] ./ 255;
color = ones(9, 3);

p = 2; % component to render
color(3:7, :) = repmat(colors(p, :), 5, 1);
intensity = [0 GMModel.mu(p)-3.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p) GMModel.mu(p)+sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+3.*sqrt(GMModel.Sigma(:,:,p)) 1e4];
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);

voldisp.Alphamap = alphamap;
voldisp.Colormap = colormap;
% 
% p = 3; % component to render
% color(3:7, :) = repmat(colors(p, :), 5, 1);
% intensity = [0 GMModel.mu(p)-3.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p) GMModel.mu(p)+sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+3.*sqrt(GMModel.Sigma(:,:,p)) 1e4];
% queryPoints = linspace(min(intensity),max(intensity),256);
% alphamap = interp1(intensity,alpha,queryPoints)';
% colormap = interp1(intensity,color,queryPoints);
% 
% voldisp.Alphamap = (voldisp.Alphamap + alphamap) ./ max(voldisp.Alphamap + alphamap);
% voldisp.Colormap = (voldisp.Colormap + colormap) ./ 2;

subplot(2,1,1)
    h = histogram('BinEdges',edges,'BinCounts',hdata);
    hold on;
    h.Normalization = "pdf";
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    hold off;
    xlim([2 1e4])
    ylim([0 max(gmPDF(h.BinEdges))])
subplot(2,1,2)
    plot(queryPoints,voldisp.Alphamap,'k-o'); hold on;
    h = plot(queryPoints,voldisp.Colormap,'-x'); hold off;
    h(1).Color = 'red'; h(2).Color = 'green'; h(3).Color = 'blue';
    xlim([2 1e4])
    legend('alpha','red','green','blue');
    xlabel('e^- density x 10^4');

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

function Nr = edensity(int16data, varargin)
    % calculate electron density (electrons/cubic angstrom) 
    % from integer data and metadata structure or scaling factors

    floatingdata = single(int16data);

    switch length(varargin)
        case 3
            high_cutoff = varargin{1};
            low_cutoff = varargin{2};
            factor_edensity = varargin{3};
        case 1
            metadataStruct = varargin{1};
            high_cutoff = metadataStruct.high_cutoff;
            low_cutoff = metadataStruct.low_cutoff;
            factor_edensity = metadataStruct.factor_edensity;
    end

    Nr = (floatingdata.*(high_cutoff - low_cutoff)/(2^16-1) + low_cutoff).*factor_edensity;

    Nr = uint16(round(Nr.*1e4)); % to make integers meaningful?
    % WARNING: cannot handle electron densities higher than 6.5535
    % electrons/cubic angstrom (typically nonphysical)
end

function rho = density(n_e, A, Z)
    % calculates mass density from electron density (e- per Angstrom^3),
    % molar mass (g/mol), and electrons per mole

    N_A = 6.0221367e23; % mol^-1, Avogadro constant

    rho = (n_e .* A) ./ (N_A .* Z) .* 1e24; % g/cm^3
end