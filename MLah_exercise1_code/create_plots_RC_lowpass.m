%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab-Script for Measurement Laboratory at Home
% RC Low-Pass Filter
% calculate response of RC circuit and plot bode diagram
% by Werner Hemmert, TUM, 18. Oct. 2016 - 16. Nov. 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          % ALWAYS start with clean workspace
close all;                  
clear all;
%%  Define figure size such that you can read the labels in a report/paper
figure                      % open figure 
fontSize=8;                 % define reasonable font size, not too small!
lw=1;                      % define linewidth, otherwise you get a hairline 
                       % which can get lost when you reduce the image size!
%% ======================= R/C Tiefpass ===================================

%---------------- Parameters: ALWAYS use SI units!! ----------------------
u1=1;                       % input voltage /V
R=10e6;                     % resistance /Ohm
C=16e-12;                   % capacity /F
f_c=1/(R*C*2*pi)            % R/C corner frequency /Hz

f=logspace(1,5,100);        % create frequency vector /Hz
u2=u1./(1+i*2*pi*f*R*C);    % calculate output voltage /V

%% --------------------------------- plot --------------------------------
a(1)=subplot(2,1,1);        % create plots for magnitude and phase 
loglog(f/1000,abs(u2),'LineWidth',lw);         % plot magnitude
%xlabel('Frequenz / kHz','FontSize',fontSize)  % omit x-label in upper plot
ylabel('|H(f)|','FontSize',fontSize)
axis([20/1000 20 3e-2 2])                      % define range x and y axes

%% ------------- some tricks for plotting --------------------------------
y_pos=[0.05 0.1 0.2 0.5 1];                   % position for y-axes labels
%  set(gca,'YtickLabel',y_pos,'FontSize',fontSize);
set(gca,'YTick',y_pos)
x_pos=[0.1 1 10];                             % position for x-axes labels
set(gca,'XtickLabel',[],'FontSize',fontSize);
set(gca,'XTick',x_pos)

H=line([0.01 f_c/1000],[1 1]);                % Bode limits: flat
set(H,'LineStyle','-.','Color','r')
H=line([f_c/1000 20],[1 0.05]);               % Bode limits: -6dB/oct
set(H,'LineStyle','-.','Color','r')
text(f_c/1000,1.3,'f_C','FontSize',fontSize,'Color','r'); % annotate
text(3,0.4,'-20 dB/Dek','FontSize',fontSize,'Color','r'); % annotate

%% --------------------------- plot phase --------------------------------
subplot(2,1,2);                               % plot phase below
semilogx(f/1000,angle(u2)*180/pi,'LineWidth',lw);
xlabel('Frequenz / kHz','FontSize',fontSize)
ylabel('Phase(H(f)) / Deg','FontSize',fontSize)
axis([20/1000 20 -100 10])                    % scale plot
y_pos=[-90 -45 0];                            % position for y-axes labels
set(gca,'YTick',y_pos)
%  set(gca,'YtickLabel',y_pos,'FontSize',fontSize);
x_pos=[0.1 1 10];                             % position for x-axes labels
set(gca,'XTick',x_pos)
x_pos=['0.1';' 1 ';' 10'];                    % define labels for x-axes
set(gca,'XtickLabel',x_pos,'FontSize',fontSize);
H=line([0.01 20],[-45 -45]);                  % Bode limits
set(H,'LineStyle','-.','Color','r')
text(f_c/1000,-35,'f_C','FontSize',fontSize,'Color','r'); % annotate

%%  Define figure size in original size 
%   then you can read the labels in the report/thesis/paper
set(gcf,'Units','Centimeters','Position',[0 0 8.4 9],'PaperPositionMode','auto')

print('rc_lp', '-depsc')                      % create scaleable figure
print('rc_lp', '-dtiff', '-r300')             % cretes pixel figure
print('rc_lp', '-dmeta')                      % windows only: emf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------- D O N E -----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
