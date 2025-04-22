% Video Generation Script

% Setup VideoWriter
v = VideoWriter('S4-2 slice');
v.FrameRate = 30;
v.Quality = 100;
open(v);

% Initialize variables for frame capture
totalFrames = 150; % Adjust based on how many frames you want

% Get volume dimensions in the original coordinate system
volumeSize = size(mbim, [1 2 3]);

% Ensure the viewer has at least one clipping plane
if isempty(viewer.ClippingPlanes)
    viewer.ClippingPlanes = [0 0 1 0];
end

% Get the initial clipping plane
initialClipPlane = viewer.ClippingPlanes(1,:);
normal = initialClipPlane(1:3);
normal = normal / norm(normal); % Ensure it's a unit vector

% Calculate the scale factor based on the initial D value
%scaleFactor = abs(initialClipPlane(4)) / max(volumeSize);

% Calculate total distance to move in the scaled coordinate system
% totalDistance = 2 * scaleFactor * dot(normal, volumeSize);

% Calculate start and end D values
% startD = initialClipPlane(4);
% endD = startD - totalDistance;
startD = 1.8447e3;
endD = 7.7050e4;

% Calculate increment per frame
increment = (endD - startD) / (totalFrames - 1);

% Display calculated values for verification
disp(['Start D: ', num2str(startD)]);
disp(['End D: ', num2str(endD)]);
disp(['Increment: ', num2str(increment)]);

% Main loop for video generation
tic;
for currentFrame = 1:totalFrames
    % Calculate new D value
    newD = startD + (currentFrame - 1) * increment;
    
    % Update only the fourth component of the clipping plane
    viewer.ClippingPlanes(1,4) = newD;
    
    % Force render update
    drawnow;
    pause(8); % Adjust this value if needed to ensure proper rendering
    
    try
        frame = getframe(viewer.Parent);
        writeVideo(v, frame);
    catch ME
        warning('Failed to capture or write frame %d. Error: %s', currentFrame, ME.message);
    end
    
    % Display progress
    disp(['Processed frame ' num2str(currentFrame) ' of ' num2str(totalFrames)]);
end
toc;

% Reset clipping plane to initial position
viewer.ClippingPlanes(1,4) = initialClipPlane(4);

% Clean up
close(v);
clear v