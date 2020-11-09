% Natalie Wilkinson, 9/3/2020

%% Retrieving Data
full_datatable = readtable("Volume_Data.csv");
genotype = table2array(full_datatable(1:29, 5));
tot_brain_vol = csvread("Volume_Data.csv", 1, 340, [1 340 29 340]);
mydata = csvread("Volume_Data.csv", 1, 8, [1 8 29 339]);

dimensions = size(mydata);
rows = dimensions(1);
columns = dimensions(2);

%% Normalizing Data
for i=1:columns
    normalized_data(:,i) = mydata(:,i)./tot_brain_vol(:,1);
end

%% Split Data According to Genotype
a = 0;
b = 0;
c = 0;
for j=1:rows
    if genotype(j, 1) == "APOE22"
        a = a + 1;
        index_apoe22(a) = j;
        apoe22_data(a, 1:332) = normalized_data(j, 1:end); % (g1)
    elseif genotype(j, 1) == "APOE33"
        b = b + 1;
        index_apoe33(b) = j;
        apoe33_data(b, 1:332) = normalized_data(j, 1:end); % (g2)
    elseif genotype(j, 1) == "APOE44"
        c = c + 1;
        index_apo44(c) = j;
        apoe44_data(c, 1:332) = normalized_data(j, 1:end);
    end
end

%% T-Test
myindex=[42 43 51 59 62 65 81 91 92 119 120 121 122 124]; % Regions correct?!
mylabels= {'Caudomedial Entorhinal Cortex', 'Dorsal Intermediate Entorhinal Cortex',...
    'Hippocampus', 'Hypothalamus', 'Septum', 'Amygdala', 'Superior Colliculus',...
    'Cerebellum', 'Dentate_Cerebellum',    'Optic Tracts', 'Fimbria',...
    'Corpus Callosum', 'Fornix', 'Cingulum'};


p1=0;
p2=0;
p1r=0;
p2r=0;
numperms=100 ;
NetThreshVals=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9]; %Usually Threshold between 13-45%
pparametrictotal=zeros(332,numel(NetThreshVals)); %?

degk1=zeros(1,numel(NetThreshVals));
degk2=zeros(1,numel(NetThreshVals));

for m=1:numel(myindex) %index of region
    for k=1:numel(NetThreshVals); %loop over net thresh
        NetThresh=NetThreshVals(k);
        [con1, con2, deg1, deg2]=PermConnectFunc(normalized_data', index_apoe22, index_apoe33, NetThresh);
        tcan=abs(deg1-deg2);
%number of non zero entries in teh correlation matrices - global values
        deg1_mat(k)=deg1;
        deg2_mat(k)=deg2;
    %global degrees
   
        reg_apoe22=con1(:,myindex(m));
        reg_apoe33=con2(:,myindex(m));


        [h pparametrictotal(myindex(m),k) ci stats]=ttest(reg_apoe22',reg_apoe33'); % finds p-values for the comparison of each regions connectivity and puts them in a 332x8 matrix with each region specified and each network density
    end
end
[h pparam ci stats]=ttest(deg1_mat,deg2_mat);  % global difference- completed at each network density



%% Finding P-Values more Manually... Review
allmyindices=[index_apoe22, index_apoe33]; % matrix of all the index values

p1_reg=0;
p2_reg=0;
numperms=250; % gets a more accurate dataset
NetThreshVals=[.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9];
degk1=zeros(1,numel(NetThreshVals));
degk2=zeros(1,numel(NetThreshVals));
for k=1:numel(NetThreshVals);
    NetThresh=NetThreshVals(k);
    [con1, con2, deg1, deg2]=PermConnectFunc(normalized_data', index_apoe22, index_apoe33, NetThresh);
    t_val=abs(deg1-deg2)
    deg1_matrix(k)=deg1;
    deg2_matrix(k)=deg2;
    for i=1:numperms %(eventually will go from 1:1000)
        g1=allmyindices(randperm(length(index_apoe22)+length(index_apoe33), 11)); % 11 integers randomly between 1 and 24
        g2=setdiff(allmyindices, g1); % returns the data in allmyindices that isn't in g1
        [con1, con2, deg1, deg2]=PermConnectFunc(normalized_data', g1, g2, NetThresh);
        t=abs(deg1-deg2);
        if t>t_val
            p1_reg=p1_reg+1; %-
        else 
            p2_reg=p2_reg+1; %-
        end
end
p1total(k)=p1_reg/numperms;
end
%Get rid of diagonal element in connectome 1. 
save('degk1')
save('degk2')