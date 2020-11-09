clear; clf
addpath(genpath('C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\BCT\'))

%% APOE22
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
kmax=numel(myseq);
MyClust_apoe22=zeros(kmax, 332);
for k=1:kmax
    n=myseq(k);
    CIJ1= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE22_' , num2str(n), '.csv']);
    Clust_apoe22=clustering_coef_bu(CIJ1);
    MyClust_apoe22(k,:)=Clust_apoe22;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\','MyClust_apoe22' ,'.csv'], MyClust_apoe22);

mypath='C:\Users\natal\OneDrive - Duke University\Research\CVNconnectivity\Independent Study Fall 2020-NW\'


    
%% APOE33
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
MyClust_apoe33=zeros(mmax, 332);
for m=1:mmax
    n=myseq(m);
    CIJ2= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE33_' , num2str(n), '.csv']);
    Clust_apoe33=clustering_coef_bu(CIJ2);
    MyClust_apoe33(m,:)=Clust_apoe33;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'MyClust_apoe33' ,'.csv'], MyClust_apoe33);

%[h, p, ci, stats]= ttest2(mydeg1', mydeg2')

%% APOE44
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
MyClust_apoe44=zeros(mmax, 332);
for x=1:mmax
    n=myseq(x);
    CIJ3= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE44_' , num2str(n), '.csv']);
    Clust_apoe44=clustering_coef_bu(CIJ3);
    MyClust_apoe44(x,:)=Clust_apoe44;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'MyClust_apoe44' ,'.csv'], MyClust_apoe44);

%% Graphing Different Regions
clf
myindex=[42 43 51 59 62 65 81 91 92 119 120 121 122 124];
mylabels= {'Caudomedial Entorhinal Cortex', 'Dorsal Intermediate Entorhinal Cortex',...
    'Hippocampus', 'Hypothalamus', 'Septum', 'Amygdala', 'Superior Colliculus',...
    'Cerebellum', 'Dentate Cerebellum',    'Optic Tracts', 'Fimbria',...
    'Corpus Callosum', 'Fornix', 'Cingulum'}

for i=1:numel(myindex)

    figure(i)
    % Interpolation
    MyClust_apoe22_new=MyClust_apoe22(:,i)
    MyClust_apoe33_new=MyClust_apoe33(:,i)
    MyClust_apoe44_new=MyClust_apoe44(:,i)
    plot(myseq, MyClust_apoe22_new, 'r.-')
    hold on
    plot(myseq, MyClust_apoe33_new, 'b-o')
    hold off
    hold on
    plot(myseq, MyClust_apoe44_new, 'k.-')
    hold off


    % Title
    mytitle=char(mylabels(i))
    title(mytitle)
    xlabel('Network Densities')
    ylabel('Clustering Coefficient')
    legend("apoe22", "apoe33", "apoe44")

end


%% Characteristic Path Length - CHECK IF CHARPATH INPUT IS THE CONNECTION MATRIX
%% APOE22
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
kmax=numel(myseq);
for k=1:kmax
    n=myseq(k);
    CIJ1= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE22_' , num2str(n), '.csv']);
    CharPath_apoe22=charpath(CIJ1);
    MyCharPath_apoe22(k,:)=CharPath_apoe22;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\','MyCharPath_apoe22' ,'.csv'], MyCharPath_apoe22);

mypath='C:\Users\natal\OneDrive - Duke University\Research\CVNconnectivity\Independent Study Fall 2020-NW\'


    
%% APOE33
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
for m=1:mmax
    n=myseq(m);
    CIJ2= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE33_' , num2str(n), '.csv']);
    CharPath_apoe33=charpath(CIJ2);
    MyCharPath_apoe33(m,:)=CharPath_apoe33;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'MyCharPath_apoe33' ,'.csv'], MyCharPath_apoe33);

%[h, p, ci, stats]= ttest2(mydeg1', mydeg2')

%% APOE44
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
for x=1:mmax
    n=myseq(x);
    CIJ3= csvread(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE44_' , num2str(n), '.csv']);
    CharPath_apoe44=charpath(CIJ3);
    MyCharPath_apoe44(x,:)=CharPath_apoe44;
end
csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
    'MyCharPath_apoe44' ,'.csv'], MyCharPath_apoe44);



%% Constructing 20 Random 332x332 Matrices
n_rand = 20;
tot_rand_mat = zeros(332, 332);
j = 0;
while j < n_rand
    j = j + 1;
    rand_mat = rand(332);
    for i = 1:length(rand_mat)
        rand_mat(i, i) = 1;
    end
    tot_rand_mat = tot_rand_mat + rand_mat;
end

avg_rand_mat = tot_rand_mat/n_rand;

% Rand Clustering Coefficient and Characteristic Path Length
myseq=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
mmax=numel(myseq);
MyClust_rand=zeros(mmax, 332);
for r=1:mmax
    Clust_rand=clustering_coef_bu(avg_rand_mat);
    MyClust_rand(r,:)=Clust_rand;
    CharPath_rand=charpath(avg_rand_mat)
    MyCharPath_rand(r,:)=CharPath_rand
end

%% Graphing Small World Index
myindex=[42 43 51 59 62 65 81 91 92 119 120 121 122 124];
mylabels= {'Caudomedial Entorhinal Cortex', 'Dorsal Intermediate Entorhinal Cortex',...
    'Hippocampus', 'Hypothalamus', 'Septum', 'Amygdala', 'Superior Colliculus',...
    'Cerebellum', 'Dentate Cerebellum',    'Optic Tracts', 'Fimbria',...
    'Corpus Callosum', 'Fornix', 'Cingulum'}

for i=1:numel(myindex)

    % Clustering Coeff
    MyClust_apoe22_new=MyClust_apoe22(:,i);
    MyClust_apoe33_new=MyClust_apoe33(:,i);
    MyClust_apoe44_new=MyClust_apoe44(:,i);
    % Characteristic Path Length
    MyCharPath_apoe22_new=MyCharPath_apoe22;
    MyCharPath_apoe33_new=MyCharPath_apoe33;
    MyCharPath_apoe44_new=MyCharPath_apoe44;
    % Clustering Coeff Rand
    MyClust_rand_new=MyClust_rand;
    % Characteristic Path Length Rand
    MyCharPath_rand_new=MyCharPath_rand;
    
    SWI_apoe22 = ((MyClust_apoe22_new)./(MyClust_rand_new))./((MyCharPath_apoe22_new)./(MyCharPath_rand_new));
    SWI_apoe33 = ((MyClust_apoe33_new)./(MyClust_rand_new))./((MyCharPath_apoe33_new)./(MyCharPath_rand_new));
    SWI_apoe44 = ((MyClust_apoe44_new)./(MyClust_rand_new))./((MyCharPath_apoe44_new)./(MyCharPath_rand_new));
    
    figure(i+14)
    a = plot(myseq, SWI_apoe22, 'r.-')
    hold on
    b = plot(myseq, SWI_apoe33, 'b-o')
    hold off
    hold on
    c = plot(myseq, SWI_apoe44, 'k.-')
    hold off
    % Title
    mytitle=char(mylabels(i))
    title(mytitle)
    xlabel('Network Densities')
    ylabel('Small World Index (SWI)')
    legend("apoe22", "apoe33", "apoe44", 'location', 'best')

end


