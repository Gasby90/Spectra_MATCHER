%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SPECTRA MATCHER Ver 1.2 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Riccardo Gasbarrone, 2016 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script allows to identify an unknown spectrum comparing it with a 
%%% pure spectra library. When the euclidean distance calculated for
%%% each library spectrum returns a small value, then is possible that the
%%% spectrum is identified with one occurring in the library.
%% 
%%% Minimum euclidean distance spectrum find in the library is then 
%%% displayed and comparised with the original input spectrum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
%% 
%%% In database or spectra matrix must be inserted reflectance or
%%% absorbance values
%%% a= SPECTRA DATABASE = spectradatabase matrix (row = known spectra; 
%%% columns = reflectance or absorbance values;
%%% b= INPUT = unknow spectra matrix(column-reflectance or assorbance
%%% values);
%%% c= SPECTRA NAMES= must contain library spectra specimen names(column - 
%%% it must be a cell array);
%%% w= WAVELENGHT = must be the wavelenght of the spectra (column - matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%% 
%%% Read file
%   SPECTRA DATABASE
[stmfile, stmpath] = uigetfile('*.txt', 'Open reference database, a');
a = importdata(fullfile(stmpath, stmfile), '\t');
%   INPUT
[stmfile, stmpath] = uigetfile('*.txt', 'Open unknown spectrum, b');
b = importdata(fullfile(stmpath, stmfile), '\t');
%   WAVELENGHT
[stmfile, stmpath] = uigetfile('*.txt', 'Open name spectra database,w');
w = importdata(fullfile(stmpath, stmfile), '\t');
%   SPECTRA NAMES
[stmfile, stmpath] = uigetfile('*.txt', 'Open reference database, c');
c = importdata(fullfile(stmpath, stmfile), '\t');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%Calculating matrix and euclidean distance
%from distance.m by R. Bunschoten(2009)
%https://www.mathworks.com/matlabcentral/fileexchange/71-distance-m/content
%//distance.m
aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
%% 
%Finding the minimum distance
MIN = min(d);
%Display minimum distance
if MIN < 1
    display ('there is a match!');
    warndlg('There is a match!','!!Match found!!')
elseif MIN< 10
   warndlg('Several matches are possible!','Warning')
else   MIN = 100;
        h = msgbox('No match was found', 'ERROR','Error');
end
X=['Minimum distance is ',num2str(MIN)];
disp(X);
%% 
%Return the location of the minimum distance
LOCATION = find(d == MIN);
%Return the name of the identified spectra
SPECTRA = c(LOCATION,1);
X1=['Identified spectrum = ',SPECTRA];
disp(X1)
% Identified spectra
is=a(:,LOCATION);
% Display the original spectra & the identified spectra
% figure % new figure
%% For subplotting
% ax1 = subplot(2,1,1); % top subplot
% ax2 = subplot(2,1,2); % bottom subplot
% ax3 = subplot(2,1,3);
% 
% plot(w,b(:,1),'b')
% title(ax1,'Input spectrum')
% ylabel(ax1,'Reflectance')
% xlabel(ax1,'Wavelenght [um]')
% 
% plot(w,is, 'r')
% title(ax2,'Identified spectrum')
% ylabel(ax2,'Reflectance')
% xlabel(ax2,'Wavelenght [um]');
% Plot on one figure
%% 
figure
plot(w,b,'b')
hold on
plot(w,is,'r')
title('Input & identified spectra')
ylabel('Reflectance')
xlabel('Wavelenght [um]')
legend('Input spectrum','Identified spectrum');
%% 
%Euclidean distance score
figure
scatter(b,is)
title('Euclidean distance score')
ylabel('Identified spectrum')
xlabel('Input spectrum')
%% 
%Calculate: coefficient of determination (R-Squared)
%Ordinary — Ordinary (unadjusted) R-squared R^2 = SSR/SST
%SSE is the sum of squared error, SSR is the sum of squared regression, 
mdl = fitlm(b,is)
%find local max and min; From: http://www.billauer.co.il/peakdet.html
figure; plot(w,b);
[maxtab, mintab] = peakdet(b, 0.001, w);
title('Local Max & Min')
hold on; plot(mintab(:,1), mintab(:,2), 'g*');
plot(maxtab(:,1), maxtab(:,2), 'r*');
