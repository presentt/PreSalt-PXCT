%%
tif_vol_file = 'S8_1_edensity.tif';
bim = blockedImage(tiffreadVolume(tif_vol_file), BlockSize=[100 100 100]);
mbim = makeMultiLevel3D(bim);
clear bim;

%%
tifsmpsngl = readmatrix([tif_vol_file(1:end-4),'_sample.csv']);
%tifsmpsngl = tifsmpsngl(tifsmpsngl>0.2); % don't fit pores
scale = numel(tifsmpsngl(tifsmpsngl>0.2))/numel(tifsmpsngl); % probability correction factor because we didn't fit the pores

GMModel = load([tif_vol_file(1:end-4),'_GMModel.mat'], "GMModel");
GMModel = GMModel.GMModel;
gmPDF = @(x) arrayfun(@(x0) pdf(GMModel,x0),x);

disp(GMModel);
disp(GMModel.mu');
disp(sqrt(reshape(GMModel.Sigma(1,1,:),1,GMModel.NumComponents)));
disp(GMModel.ComponentProportion.*100);

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

voldisp.RenderingStyle = "SlicePlanes";

viewer.ScaleBar = 'on';
viewer.ScaleBarUnits = 'nm';

%%
alpha = [0 0 0.005 1 1 1 1 1 .005 .005];
color = ones(length(alpha), 3);

colors = [231 208 141;
          150 100 50; % kerogen
          231 208 141;
          255 0 0;
          255 255 255;] ./ 255;

p = 1; % component to render
color(4:8, :) = repmat(colors(p, :), 5, 1);
intensity = [0 min(tifsmpsngl)+.001 GMModel.mu(p)-3.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)-sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p) GMModel.mu(p)+sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+2.*sqrt(GMModel.Sigma(:,:,p)) GMModel.mu(p)+3.*sqrt(GMModel.Sigma(:,:,p)) max(tifsmpsngl)];
queryPoints = linspace(min(intensity), max(intensity), 256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);

voldisp.Alphamap = alphamap;
voldisp.Colormap = colormap;

subplot(2,1,1)
    h = histogram(tifsmpsngl);
    hold on;
    h.Normalization = "pdf";
    plot(h.BinEdges,gmPDF(h.BinEdges).*scale,'LineWidth',3);
    for p = 1:GMModel.NumComponents
        plot(h.BinEdges,...
            scale.*normpdf(h.BinEdges,GMModel.mu(p),sqrt(GMModel.Sigma(:,:,p))).*GMModel.ComponentProportion(p),...
            'LineWidth',1);
    end
        plot(queryPoints,voldisp.Alphamap+1,'k-o');
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
    hold off;

voldisp.RenderingStyle = "VolumeRendering";

%%
% use DataReadFinished event to record a video...