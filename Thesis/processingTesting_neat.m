source('mSphere.m');
% File: ProcessingTesting.m
% Author: Michael Thoreau
% Organisation: Commonwealth Scientific and Industrial Research Organisation - Autonommous Systems Lab
%
% note: code has been vectorised for performance > readability, performance increase by orders of magnitude
%

source('delay_calculation.m');

tic

%Constants
Fs = 48000; %Hz
F = 200;  %Hz
numPeriods = 1;
t=[0:1/Fs:(numPeriods/F)];  %s
numSamples = ceil((numPeriods/F)*Fs);
Nmicrophones = 8;
Nbeams = 1000;

results = []

% Calculate sample offsets for each beam - this is what would be provided to the beamformer on the fpga
[beams offsets] = delay_calculation(Nbeams);
buffer = max(max(offsets))+1
padding = buffer*3

% Generate testing sinusoid 
sinWave = sin(2*pi*F*t);

% Generate blackman window
w = blackman(length(sinWave));
wSine = zeros(numSamples);
wSine = w'.*sinWave;
wSine = [zeros(1, padding) wSine zeros(1, padding)];
t2 = [0:1/Fs:(numPeriods/F)+2*padding/Fs];


%The number of beams requested may not fit well on the sphere, so we may have slightly less
Nbeams = length(beams);

% pick a direction for the testing signal to originate from - for simulation this will be one of the beams
direction = 1;
receivedSig = zeros(length(wSine), Nmicrophones);

% Generate the signal at each microphone (with delays)
for j  = 1:Nmicrophones
    offset = offsets(direction, j);
    range = buffer:length(wSine)-buffer;
    receivedSig(range,j) = wSine(range-offset);
    receivedSig2(range,j) = receivedSig(range+offset, j);
end


% % Plot the signals at the microphone
figure(2);
plot(t2, receivedSig);
title('Signal as received by all microphones');
xlabel('time - s');
ylabel('Magnitude');



% declare empty array for energy vector
energy = zeros(Nbeams, 1);


% Compute the energy of the received signal in the direction of each beam
for direction = 1:Nbeams
    % Zero out the virtual signal
    range = 1:length(receivedSig)-buffer;


     % Elementwise sum and accumulate of the signals received at each microphone, with compensating delays added
    energy(direction) = sum(abs( receivedSig(range+offsets(direction, 1), 1)...
    + receivedSig(range+offsets(direction, 2), 2)...
    + receivedSig(range+offsets(direction, 3), 3)...
    + receivedSig(range+offsets(direction, 4), 4)...
    + receivedSig(range+offsets(direction, 5), 5)...
    + receivedSig(range+offsets(direction, 6), 6)...
    + receivedSig(range+offsets(direction, 7), 7)...
    + receivedSig(range+offsets(direction, 8), 8)));

end
time1 = toc



% convert to twos compliment magnitude for testing
energy = floor(energy*8388608);
% find Cartesian coordinates for the energy vectors
starts = zeros(Nbeams,3);
[x, y, z] = sph2cart(beams(:,1), beams(:,2), energy);

% scatter3(beams(:,1), beams(:,2), energy')


% Interpolate data over a grid of 360x360 points
resolution = 360;
[xx yy] = meshgrid(-pi:2*pi/resolution:pi, -pi/2:pi/resolution:pi/2);
ee = griddata(beams(:,1), beams(:,2), energy, xx, yy, 'linear')

maxIndex = find(energy==max(energy), 1)



figure(3);
surf(xx,yy,ee, 'EdgeColor','None');
hold on
scatter3(beams(maxIndex,1), beams(maxIndex,2), max(energy));
hold off
%view(2) %view as an image (optional)
title('Beamformer output image');
xlabel('Azimuth - radians');
ylabel('elevation - radians');
zlabel('Energy');
xlim([-pi pi]);
ylim([-pi/2 pi/2]);



% apply the optimal delay to the inputs
receivedSig2 = zeros(length(receivedSig), Nmicrophones);
for j  = 1:Nmicrophones
    offset = offsets(maxIndex, j);
    range = buffer:length(receivedSig)-buffer;
    receivedSig2(range,j) = receivedSig(range+offset,j);
end

% Plot the signals at the microphone
figure(4);
plot(receivedSig2);
title('Signal as received by all microphones - Estimated Position');
xlabel('time - s');
ylabel('Magnitude');
% xlim([0.0015 0.0027]);
ylim([-1 1]);





% Write samples to file - For use in VHDL testbench
fileID = fopen('input1.txt','w');
accumulator = 0;
a = 0
for row = receivedSig'
    tmp = 0;
	for j = 1:length(row)
		sample = int32(row(j)*8388608);
		fprintf(fileID, '%d\n', sample);
	end

end

fclose(fileID);


%Write offsets to file - For use in VHDL testbench (or in hardware)
fileID = fopen('offsets1.txt','w');
accumulator = 0;
a = 0
for row = offsets'
    tmp = 0;
    for j = 1:length(row)
        offset = row(j);
        fprintf(fileID, '%d\n', offset);
    end
end

fclose(fileID);