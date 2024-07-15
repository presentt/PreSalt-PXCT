%%
tif_vol_file = 'S8_2_edensity.tif';
bim = blockedImage(tiffreadVolume(tif_vol_file), BlockSize=[100 100 100]);
mbim = makeMultiLevel3D(bim);
clear bim;

%%
tifsmpsngl = readmatrix([tif_vol_file(1:end-4),'_sample.csv']);
tifsmpsngl = tifsmpsngl(tifsmpsngl>0.2); % don't fit pores

GMModel = load([tif_vol_file(1:end-4),'_GMModel.mat'], "GMModel");
GMModel = GMModel.GMModel;
gmPDF = @(x) arrayfun(@(x0) pdf(GMModel,x0),x);

%%
viewer = viewer3d();

info = imfinfo(tif_vol_file);
res = info.XResolution; % should be cm
px = 1/(res(1)*1e6); % pixel to nm ratio
T = [px 0 0 0;
     0 px 0 0;
     0 0 px 0;
     0 0 0 1];

voldisp = volshow(mbim, ...
                  'Transformation', affinetform3d(T), ...
                  parent=viewer);

viewer.ScaleBar = 'on';
viewer.ScaleBarUnits = 'nm';

%%
alpha = [0 0 .02 .8 1 .8 0.02 .005 .005];
colors = [200 140 75;
          200 140 75;
          231 208 141;
          255 255 255;
          0   128 0;
          128 64 0;] ./ 255;
color = ones(length(alpha), 3);

p = 3; % component to render
color(3:7, :) = repmat(colors(p, :), 5, 1);
intensity = [0 GMModel.mu(p)-3.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p) GMModel.mu(p)+sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+3.*sqrt(GMModel.Sigma(:,:,p)) max(tifsmpsngl)];
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);

voldisp.Alphamap = alphamap;
voldisp.Colormap = colormap;
% 
% p = 6; % component to render
% color(3:7, :) = repmat(colors(p, :), 5, 1);
% intensity = [0 GMModel.mu(p)-3.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p) GMModel.mu(p)+sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+3.*sqrt(GMModel.Sigma(:,:,p)) 6.55];
% alphamap = interp1(intensity,alpha,queryPoints)';
% colormap = interp1(intensity,color,queryPoints);
% 
% voldisp.Alphamap = (voldisp.Alphamap + alphamap) ./ max(voldisp.Alphamap + alphamap);
% voldisp.Colormap = (voldisp.Colormap + colormap) ./ 2;

subplot(2,1,1)
    h = histogram(tifsmpsngl); %histogram('BinEdges',edges,'BinCounts',hdata);
    hold on;
    h.Normalization = "pdf";
    plot(h.BinEdges,gmPDF(h.BinEdges),'LineWidth',3);
    for p = 1:GMModel.NumComponents
        plot(h.BinEdges,...
            normpdf(h.BinEdges,GMModel.mu(p),sqrt(GMModel.Sigma(:,:,p))).*GMModel.ComponentProportion(p),...
            'LineWidth',1);
    end
    hold off;
    xlim([0 max(tifsmpsngl)])
    ylim([0 10])
subplot(2,1,2)
    plot(queryPoints,voldisp.Alphamap,'k-o'); hold on;
    h = plot(queryPoints,voldisp.Colormap,'-x'); hold off;
    h(1).Color = 'red'; h(2).Color = 'green'; h(3).Color = 'blue';
    xlim([0 max(tifsmpsngl)])
    legend('alpha','red','green','blue');
    xlabel('e^- density');

%voldisp.RenderingStyle ="GradientOpacity";

%%
% use DataReadFinished event to record a video...