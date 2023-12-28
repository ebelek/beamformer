
% Beamformer for ultrasound Bmode imaging
close all; clc; clear all;
%% define data paths

rcvDataPath = 'C:\Users\path';
imgDataPath = 'C:\Users\path';
fileLabel = {'filename.mat'};

%% beamform

for nn = 1:length(fileLabel)

    load(fullfile(rcvDataPath,[fileLabel{nn},'RcvData.mat']))
    load(fullfile(imgDataPath,[fileLabel{nn},'ImgDataP.mat']))
    load(fullfile(imgDataPath,[fileLabel{nn},'workspace.mat']))

    disp('Done loading files!')

    x_ppw = 4;
    z_ppw = 4;

    % get relevant info from verasonics output
    c = 1540;                      % Speed of sound (m/s)
    fs = (Receive(1).ADCRate)*1e6; % [Hz]
    f = Trans.frequency*1e6;
    lambda = c/f;
    foc = P.focus;
    numElem = Trans.numelements;
    numRays = P.numRays;
    elemPos = Trans.ElementPos(:,1)*c/f;
    RcvData = RcvData{1,1};
    ImgData = ImgDataP{1,1};
    rf_raw = RcvData(:,Trans.Connector,1);
    dx = lambda/x_ppw;
    dz = lambda/z_ppw;

    % delay and sum

    %%%%%%%%%%% delay equation %%%%%%%%%%%%%%%%%%
    %%%% T(x_1,x,z)=(z+sqrt(z^2+(x−x_1)^2)/c %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic
    P.startDepth = 0e-3;
    t0 = 2*P.startDepth/c;

    % define imaging axes
    x = linspace(elemPos(1),elemPos(end),P.numRays); % lateral positions
    z = linspace(P.startDepth,P.endDepth,round((P.endDepth-P.startDepth)/dz)+1);    % axial positions
    x_interp = linspace(elemPos(1),elemPos(end),round((elemPos(end)-elemPos(1))/(dx))+1);  % interpolated lateral positions

    ds = zeros(length(x),length(z)); % delay and sum array
    for ii = 1:length(x) % for each Ray line
        rf = rf_raw(Receive(ii).startSample:Receive(ii).endSample,:);
        t0(ii) = 2*P.startDepth/c - max(TX(ii).Delay)*lambda./c;

        for jj = 1:length(z) 

            % calculate delays
            aperture = find(abs(elemPos-x(ii))<= z(jj)/P.fnum); % define aperture based on f-number

            % calculate 2-way delay
            delay2 = (z(jj) + sqrt((elemPos(aperture)-x(ii)).^2 + z(jj)^2))./c-t0(ii); 
            idt = round((delay2).*fs)+1;
            idt = idt(idt <= Receive(ii).endSample-Receive(ii).startSample);
            temp = [];
            for kk = 1:length(idt)
                temp(kk) = sum(rf(idt(kk),aperture(kk)));
            end
            % calculate sum
            %temp = abs(hilbert(temp));
            ds(ii,jj) = sum(temp); % delay and summed image
        end

    end
    toc

    env = abs(hilbert(ds)); % envelope detection
    if round(length(x)/length(x_interp)) > 1
        env = interp(env,round(length(x)/length(x_interp)));
    end
    logenv = mag2db(env/max(env(:)));   % log compression


    %% plots
    figure
    subplot(1,2,1)
    rawimg = squeeze(ImgData(2:end,:,:,1));
    logimg = mag2db(rawimg/max(rawimg(:)));
    interpfactor = 2;
    imagesc(x*1e3,z*1e3,logimg,[-30 0]);
    %interpimg = interp2(logimg,interpfactor);
    %imagesc(x*1e3,z*1e3,interpimg)
    colormap gray; colorbar
    title('ImgData from Verasonics (all pass filter)')

    subplot(1,2,2)
    imagesc(x*1e3,z*1e3,logenv',[-30 0])
    colormap gray; colorbar
    title('Beamformed Image (all pass filter)')
     
end

