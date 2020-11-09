clear; clf
addpath(genpath('C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\BCT\'))

%% APOE22
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
kmax=numel(myseq);
mydeg1=zeros(kmax, 332);
for k=1:kmax
    n=myseq(k);
    CIJ1= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE22_' , num2str(n), '.csv']);
    [id,od,deg1] = degrees_dir(CIJ1) ; 
    mydeg1(k,:)=deg1;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\','DEG1' ,'.csv'], mydeg1);

mypath='C:\Users\natal\OneDrive - Duke University\Research\CVNconnectivity\Independent Study Fall 2020-NW\'


    
%% APOE33
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
mydeg2=zeros(mmax, 332);
for m=1:mmax
    n=myseq(m);
    CIJ2= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE33_' , num2str(n), '.csv']);
    [id,od,deg2] = degrees_dir(CIJ2) ; 
    mydeg2(m,:)=deg2;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'DEG2' ,'.csv'], mydeg2);

[h, p, ci, stats]= ttest2(mydeg1', mydeg2')

%% APOE44
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
mydeg3=zeros(mmax, 332);
for x=1:mmax
    n=myseq(x);
    CIJ2= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE44_' , num2str(n), '.csv']);
    [id,od,deg3] = degrees_dir(CIJ2) ; 
    mydeg3(x,:)=deg3;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'DEG3' ,'.csv'], mydeg3);

%% Graphing Different Regions
clf
myindex=[42 43 51 59 62 65 81 91 92 119 120 121 122 124];
mylabels= {'Caudomedial Entorhinal Cortex', 'Dorsal Intermediate Entorhinal Cortex',...
    'Hippocampus', 'Hypothalamus', 'Septum', 'Amygdala', 'Superior Colliculus',...
    'Cerebellum', 'Dentate Cerebellum',    'Optic Tracts', 'Fimbria',...
    'Corpus Callosum', 'Fornix', 'Cingulum'}

for i=1:numel(myindex)
%     hcNOS=mydeg1(:,myindex(i))
%     hcAD=mydeg2(:,myindex(i))
%     plot(myseq', hcNOS,'ro')
%     hold on
%     plot(myseq', hcAD,'b*')
%     hold off

    figure(i)
    % Interpolation
    mydeg1new=mydeg1(:,i)
    mydeg2new=mydeg2(:,i)
    mydeg3new=mydeg3(:,i)
    plot(myseq, mydeg1new, 'r.-')
    hold on
    plot(myseq, mydeg2new, 'b-o')
    hold off
    hold on
    plot(myseq, mydeg3new, 'k.-')
    hold off


    % Title
    mytitle=char(mylabels(i))
    title(mytitle)
    xlabel('Network Densities')
    ylabel('Degrees of Connectivity')
    legend("apoe22", "apoe33", "apoe44")

end


%% Regional Correlations
apoe22_mat= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'APOE22_0.4.csv']);
myregions=[42 43 51 59 62 65 81 91 92 119 120 121 122 124];
for a = 1:length(myregions)
    reg_corr(a, a) = apoe22_mat(myregions(a), myregions(a))
    b = a
    while b~=length(myregions)
        b = b + 1
        reg_corr(a, b) = apoe22_mat(myregions(a), myregions(b))
        reg_corr(b, a) = apoe22_mat(myregions(a), myregions(b))
    end
end
figure(15)
imagesc(reg_corr)
colorbar
caxis([-1 1])
axis equal
mycolormap=jet(11)
mycolormap(6, :)=[1 1 1]
colormap(mycolormap)
title("Correlation of Selected APOE22 Regions at k = 0.40")


apoe33_mat= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'APOE33_0.4.csv']);
myregions=[42 43 51 59 62 65 81 91 92 119 120 121 122 124];
for a = 1:length(myregions)
    reg_corr2(a, a) = apoe33_mat(myregions(a), myregions(a))
    b = a
    while b~=length(myregions)
        b = b + 1
        reg_corr2(a, b) = apoe33_mat(myregions(a), myregions(b))
        reg_corr2(b, a) = apoe33_mat(myregions(a), myregions(b))
    end
end
figure(16)
imagesc(reg_corr2)
colorbar
caxis([-1 1])
axis equal
mycolormap=jet(11)
mycolormap(6, :)=[1 1 1]
colormap(mycolormap)
title("Correlation of Selected APOE33 Regions at k = 0.40")


apoe44_mat= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'APOE44_0.4.csv']);
myregions=[42 43 51 59 62 65 81 91 92 119 120 121 122 124];
for a = 1:length(myregions)
    reg_corr3(a, a) = apoe44_mat(myregions(a), myregions(a))
    b = a
    while b~=length(myregions)
        b = b + 1
        reg_corr3(a, b) = apoe44_mat(myregions(a), myregions(b))
        reg_corr3(b, a) = apoe44_mat(myregions(a), myregions(b))
    end
end
figure(17)
imagesc(reg_corr3)
colorbar
caxis([-1 1])
axis equal
mycolormap=jet(11)
mycolormap(6, :)=[1 1 1]
colormap(mycolormap)
title("Correlation of Selected APOE44 Regions at k = 0.40")


%% To Do
% Look in BCT tool box to find global network properties
% then focus specifically on four regions
