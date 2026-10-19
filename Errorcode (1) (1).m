
port1 = readmatrix("subsonic_Port01_1_0.dat"); 
port2 = readmatrix("subsonic_Port02_1_0.dat");
port3 = readmatrix("subsonic_Port03_1_0.dat");
port4 = readmatrix("subsonic_Port04_1_0.dat");
port5 = readmatrix("subsonic_Port05_1_0.dat");
port6 = readmatrix("subsonic_Port06_1_0.dat");
port7 = readmatrix("subsonic_Port07_1_0.dat");

port11 = readmatrix("supersonic_Port01_1_0.dat"); 
port12 = readmatrix("supersonic_Port02_1_0.dat");
port13 = readmatrix("supersonic_Port03_1_0.dat");
port14 = readmatrix("supersonic_Port04_1_0.dat");
port15 = readmatrix("supersonic_Port05_1_0.dat");
port16 = readmatrix("supersonic_Port06_1_0.dat");
port17 = readmatrix("supersonic_Port07_1_0.dat");

gain = 1.7838e4; % Pa/V

startData = readmatrix('Group02/supersonic_Port07_offset_0.dat', ...
    'FileType', 'text', 'Delimiter', '\t', ...
    'NumHeaderLines', 18);

endData = readmatrix('Group02/offset_final_Port07_0.dat', ...
    'FileType', 'text', 'Delimiter', '\t', ...
    'NumHeaderLines', 18);

startVoltage = startData(:, 2:13);
endVoltage = endData(:, 2:13);

startOffset = mean(startVoltage, 1, 'omitnan');
endOffset = mean(endVoltage, 1, 'omitnan');

channelOffset = (startOffset + endOffset) / 2;
averageOffset = mean(channelOffset);
averagepress = mean(gain*channelOffset);

N_val = zeros(1,16);

% Store the available datasets for subsequent processing
dataSets = {port1, port2, port3, port4, port5, port6, port7, ...
    port11, port12, port13, port14, port15, port16, port17,...
    startData,endData};

for i = 1:16 %[output:group:1001cd00]

    time = dataSets{i}(:, 1);

    % Average the 12 channels at each time sample
    averageVoltage1 = mean(dataSets{i}(:, 2:13), 2, 'omitnan');

    correctedVoltage = averageVoltage1 - averageOffset;
    pressure = gain * correctedVoltage;

    % Mean pressure and pressure fluctuations
    meanPressure = mean(pressure, 'omitnan');
    pressureFluctuation = pressure - meanPressure;
%[output:9bd8f1df]
    % Time between samples
    dt = mean(diff(time));

    % Calculate normalized autocorrelation
    [Bxx, lags] = xcorr(pressureFluctuation, 'coeff');

    % Keep only zero and positive time lags
    positiveLags = lags >= 0;
    Bxx = Bxx(positiveLags);
    lags = lags(positiveLags);

    % Convert lags to seconds
    timeLag = lags * dt;

    % Find the first zero crossing
    zeroIndex = find(Bxx <= 0, 1);

    integralTime = 2*trapz(timeLag(1:zeroIndex), Bxx(1:zeroIndex));

    N = dataSets{i}(end,1)/integralTime;
    N = floor(N);

    N_val(1,i)=N;
end %[output:group:1001cd00]

% Plot autocorrelation
figure;
plot(timeLag, Bxx);
xlabel('Time lag (s)');
ylabel('Normalized autocorrelation');
title('Pressure Autocorrelation');
grid on;

% Integral time scale up to the first zero crossing


fprintf('First zero crossing = %.5f s\n', timeLag(zeroIndex)); %[output:5b43f644]
fprintf('Integral time scale = %.5f s\n', integralTime); %[output:3eb9f7fb]
fprintf('"N" independent data point = %.5f\n',N); %[output:8428c229]


bais_e = 0.01*10;

offset_e = sqrt((std(startOffset)/sqrt(N_val(1,15)))^2 + (std(endOffset)/sqrt(N_val(1,16)))^2 + bais_e^2);

stat_p_1_sub_e = gain * sqrt(offset_e^2 + (std(port1(:,2))./sqrt(N_val(1,1)).^2));
stat_p_2_sub_e = gain * sqrt(offset_e^2 + (std(port2(:,3))./sqrt(N_val(1,2))).^2);
stat_p_3_sub_e = gain * sqrt(offset_e^2 + (std(port3(:,4))./sqrt(N_val(1,3))).^2);
stat_p_4_sub_e = gain * sqrt(offset_e^2 + (std(port4(:,6))./sqrt(N_val(1,4))).^2);
stat_p_5_sub_e = gain * sqrt(offset_e^2 + (std(port5(:,5))./sqrt(N_val(1,5))).^2);
stat_p_6_sub_e = gain * sqrt(offset_e^2 + (std(port6(:,7))./sqrt(N_val(1,6))).^2);
stat_p_7_sub_e = gain * sqrt(offset_e^2 + (std(port7(:,8))./sqrt(N_val(1,7))).^2);
stat_p_8_sub_e = gain * sqrt(offset_e^2 + (std(port7(:,9))./sqrt(N_val(1,7))).^2);
stat_p_9_sub_e = gain * sqrt(offset_e^2 + (std(port7(:,10))./sqrt(N_val(1,7))).^2);
stat_p_10_sub_e = gain * sqrt(offset_e^2 + (std(port7(:,12))./sqrt(N_val(1,7))).^2);
stat_p_11_sub_e = gain * sqrt(offset_e^2 + (std(port7(:,11))./sqrt(N_val(1,7))).^2);

stat_p_sub_e   = [stat_p_1_sub_e, stat_p_2_sub_e, stat_p_3_sub_e, stat_p_4_sub_e, ...
    stat_p_5_sub_e, stat_p_6_sub_e, stat_p_7_sub_e, stat_p_8_sub_e, stat_p_9_sub_e, stat_p_10_sub_e, stat_p_11_sub_e];

stat_p_1_super_e = gain * sqrt(offset_e^2 + (std(port1(:,2))./sqrt(N_val(1,1))).^2);
stat_p_2_super_e = gain * sqrt(offset_e^2 + (std(port2(:,3))./sqrt(N_val(1,2))).^2);
stat_p_3_super_e = gain * sqrt(offset_e^2 + (std(port3(:,4))./sqrt(N_val(1,3))).^2);
stat_p_4_super_e = gain * sqrt(offset_e^2 + (std(port4(:,6))./sqrt(N_val(1,4))).^2);
stat_p_5_super_e = gain * sqrt(offset_e^2 + (std(port5(:,5))./sqrt(N_val(1,5))).^2);
stat_p_6_super_e = gain * sqrt(offset_e^2 + (std(port6(:,7))./sqrt(N_val(1,6))).^2);
stat_p_7_super_e = gain * sqrt(offset_e^2 + (std(port7(:,8))./sqrt(N_val(1,7))).^2);
stat_p_8_super_e = gain * sqrt(offset_e^2 + (std(port7(:,9))./sqrt(N_val(1,7))).^2);
stat_p_9_super_e = gain * sqrt(offset_e^2 + (std(port7(:,10))./sqrt(N_val(1,7))).^2);
stat_p_10_super_e = gain * sqrt(offset_e^2 + (std(port7(:,12))./sqrt(N_val(1,7))).^2);
stat_p_11_super_e = gain * sqrt(offset_e^2 + (std(port7(:,11))./sqrt(N_val(1,7))).^2);

stat_p_super_e = [stat_p_1_super_e, stat_p_2_super_e, stat_p_3_super_e, stat_p_4_super_e, ...
    stat_p_5_super_e, stat_p_6_super_e, stat_p_7_super_e, stat_p_8_super_e, stat_p_9_super_e, stat_p_10_super_e, stat_p_11_super_e];

tot_p_1_sub_e = gain*sqrt(((std(port1(:,13)./sqrt(7)).^2 + offset_e^2) +(std(port1(:,2)./sqrt(N_val(1,1))).^2)));
tot_p_2_sub_e = gain*sqrt(((std(port2(:,13)./sqrt(7)).^2 + offset_e^2) + (std(port2(:,3)./sqrt(N_val(1,2))).^2)));
tot_p_3_sub_e = gain*sqrt(((std(port3(:,13)./sqrt(7)).^2 + offset_e^2) + std((port3(:,4)./sqrt(N_val(1,3))).^2)));
tot_p_4_sub_e = gain*sqrt((std(port4(:,13))./sqrt(7)).^2 + offset_e^2 +(std(port4(:,6))./sqrt(N_val(1,4))).^2);
tot_p_5_sub_e = gain*sqrt((std(port5(:,13))./sqrt(7)).^2 + offset_e^2 + (std(port5(:,5))./sqrt(N_val(1,5))).^2);
tot_p_6_sub_e = gain*sqrt((std(port6(:,13))./sqrt(7)).^2 + offset_e^2 + (std(port6(:,7))./sqrt(N_val(1,6))).^2);
tot_p_7_sub_e = gain*sqrt(((std(port7(:,13))./sqrt(7)).^2 + offset_e^2) + (std(port6(:,8))./sqrt(N_val(1,6))).^2);
tot_p_1_super_e = gain*sqrt(((std(port11(:,13))./sqrt(7)).^2 + offset_e^2) + (std(port6(:,2))./sqrt(N_val(1,6))).^2);
tot_p_2_super_e = gain * sqrt((std(port12(:,13))/sqrt(7))^2 + offset_e^2 + (std(port6(:,3))/sqrt(N_val(1,6)))^2);
tot_p_3_super_e = gain * sqrt((std(port13(:,13))/sqrt(7))^2 + offset_e^2 + (std(port6(:,4))/sqrt(N_val(1,6)))^2);
tot_p_4_super_e = gain * sqrt((std(port14(:,13))/sqrt(7))^2 + offset_e^2 +(std(port6(:,6))/sqrt(N_val(1,6)))^2);
tot_p_5_super_e = gain * sqrt((std(port15(:,12))/sqrt(7))^2 + offset_e^2) + (std(port6(:,5))/sqrt(N_val(1,6)))^2;
tot_p_6_super_e = gain * sqrt((std(port16(:,12))/sqrt(7))^2 + offset_e^2 + (std(port6(:,7))/sqrt(N_val(1,6)))^2);
tot_p_7_super_e = gain * sqrt((std(port17(:,12))/sqrt(7))^2 + offset_e^2 + (std(port7(:,8))/sqrt(N_val(1,7)))^2);

tot_p_sub_e   = [tot_p_1_sub_e, tot_p_2_sub_e, tot_p_3_sub_e, tot_p_4_sub_e, ...
    tot_p_5_sub_e, tot_p_6_sub_e, tot_p_7_sub_e];
tot_p_super_e = [tot_p_1_super_e, tot_p_2_super_e, tot_p_3_super_e, tot_p_4_super_e, ...
    tot_p_5_super_e, tot_p_6_super_e, tot_p_7_super_e];

stat_p_sub_e_fM   = [stat_p_1_sub_e, stat_p_2_sub_e, stat_p_3_sub_e, stat_p_4_sub_e, ...
    stat_p_5_sub_e, stat_p_6_sub_e, stat_p_7_sub_e];

stat_p_super_e_fM = [stat_p_1_super_e, stat_p_2_super_e, stat_p_3_super_e, stat_p_4_super_e, ...
    stat_p_5_super_e, stat_p_6_super_e, stat_p_7_super_e];

pressure_ratio_sub_e = pressure_ratio_sub.*sqrt((stat_p_sub_e_fM ./ static_sub).^2 + (tot_p_sub_e ./ tot_p_sub).^2); %[output:4dadb30f]
pressure_ratio_super_e = pressure_ratio_super.*sqrt((stat_p_super_e_fM ./ static_super).^2 + (tot_p_super_e ./ tot_p_super).^2);

M_sub_e = M_subsonic .* sqrt((stat_p_sub_e_fM ./ static_sub).^2 + (tot_p_sub_e ./ tot_p_sub).^2) / 2;
M_super_e = M_supersonic .* sqrt((stat_p_super_e_fM ./ static_super).^2 + (tot_p_super_e ./ tot_p_super).^2) / 2;


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":9}
%---
%[output:9bd8f1df]
%   data: {"dataType":"text","outputData":{"text":"Mean pressure = -20146.38 Pa\nMean pressure = -20113.24 Pa\nMean pressure = -19887.78 Pa\nMean pressure = -19722.67 Pa\nMean pressure = -19592.98 Pa\nMean pressure = -19498.42 Pa\nMean pressure = -19674.15 Pa\nMean pressure = -50203.77 Pa\nMean pressure = -50476.37 Pa\nMean pressure = -50192.87 Pa\nMean pressure = -51436.57 Pa\nMean pressure = -51487.84 Pa\nMean pressure = -51262.21 Pa\nMean pressure = -51580.81 Pa\nMean pressure = 34.44 Pa\nMean pressure = -34.44 Pa\n","truncated":false}}
%---
%[output:5b43f644]
%   data: {"dataType":"text","outputData":{"text":"First zero crossing = 2.90000 s\n","truncated":false}}
%---
%[output:3eb9f7fb]
%   data: {"dataType":"text","outputData":{"text":"Integral time scale = 0.61288 s\n","truncated":false}}
%---
%[output:8428c229]
%   data: {"dataType":"text","outputData":{"text":"\"N\" independent data point = 14.00000\n","truncated":false}}
%---
%[output:4dadb30f]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'pressure_ratio_sub'."}}
%---
