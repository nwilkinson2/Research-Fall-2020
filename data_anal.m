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
        index_apoe22(a) = j
        apoe22_data(a, 1:332) = normalized_data(j, 1:end);
    elseif genotype(j, 1) == "APOE33"
        b = b + 1;
        index_apoe22(b) = j
        apoe33_data(b, 1:332) = normalized_data(j, 1:end);
    elseif genotype(j, 1) == "APOE44"
        c = c + 1;
        index_apoe22(c) = j
        apoe44_data(c, 1:332) = normalized_data(j, 1:end);
    end
end

%% Generate Covariance Matrices - Cov
cov_apoe22 = cov(apoe22_data); %(g1)
cov_apoe33 = cov(apoe33_data); %(g2)
cov_apoe44 = cov(apoe44_data);


%% Plotting
% figure(1)
% imagesc(cov_apoe22)
% colormap(jet)
% axis equal
% caxis([0 30*10^-10])
% colorbar
% title("Covariance of APOE22")
% 
% figure(2)
% imagesc(cov_apoe33)
% colormap(jet)
% axis equal
% caxis([0 30*10^-10])
% colorbar
% title("Covariance of APOE33")
% 
% figure(3)
% imagesc(cov_apoe44)
% colormap(jet)
% axis equal
% caxis([0 30*10^-10])
% colorbar
% title("Covariance of APOE44")

%% Generate Covariance Matrices - Corr
corr_apoe22 = corr(apoe22_data);
corr_apoe33 = corr(apoe33_data);
corr_apoe44 = corr(apoe44_data);


%% Plotting
figure(4)
imagesc(corr_apoe22)
colorbar
caxis([-1 1])
axis equal
mycolormap=jet(11)
mycolormap(6, :)=[1 1 1]
colormap(mycolormap)
title("Correlation of APOE22")

figure(5)
imagesc(corr_apoe33)
colorbar
caxis([-1 1])
axis equal
mycolormap=jet(11)
mycolormap(6, :)=[1 1 1]
colormap(mycolormap)
title("Correlation of APOE33")

figure(6)
imagesc(corr_apoe44)
colorbar
caxis([-1 1])
axis equal
mycolormap=jet(11)
mycolormap(6, :)=[1 1 1]
colormap(mycolormap)
title("Correlation of APOE44")


% Remove bottom x(100)% of data.
Thresholds = [.15, .20, .25, .30, .35, .40, .45 .50 .55 .60 .65 .70 .75 .8 .85 .9]
count = 0
for i = 1:length(Thresholds)
    count = count + 3
    x = Thresholds(i)
    percentage = x*100
    MyThresh=ceil((x)*(332*332-332))
    inv_corr_apoe22=corr_apoe22(:); % creates column matrix
    [val,myindex]=sort(abs(inv_corr_apoe22), 'ascend'); % sorts matrix into ascending order and returns the index for each value
    MyNewThresh=inv_corr_apoe22(myindex(1:MyThresh)); % same thing as val?
    R=find(abs(inv_corr_apoe22)<=max(MyNewThresh));
    corr_apoe22_x=corr_apoe22; 
    corr_apoe22_x(R)=0; % sets x% of dat to be 0! Removea certain percentage of data
    corr_apoe22_x=reshape(corr_apoe22_x,332,332);
    csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE22_' , num2str(Thresholds(i)), '.csv'],corr_apoe22_x);
    figure(15 + count)
    imagesc(corr_apoe22_x)
    colorbar
    caxis([-1 1])
    axis equal
    mycolormap=jet(11)
    mycolormap(6, :)=[1 1 1]
    colormap(mycolormap)
    title("Correlation of " + percentage + "% of APOE22")

    MyThresh=ceil((x)*(332*332-332))
    inv_corr_apoe33=corr_apoe33(:);
    [val,myindex]=sort(abs(inv_corr_apoe33), 'ascend');
    MyNewThresh=inv_corr_apoe33(myindex(1:MyThresh));
    R=find(abs(inv_corr_apoe33)<=max(MyNewThresh));
    corr_apoe33_x=corr_apoe33; 
    corr_apoe33_x(R)=0;
    corr_apoe33_x=reshape(corr_apoe33_x,332,332);
    csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE33_' , num2str(Thresholds(i)), '.csv'],corr_apoe33_x);
    figure(16 + count)
    imagesc(corr_apoe33_x)
    colorbar
    caxis([-1 1])
    axis equal
    mycolormap=jet(11)
    mycolormap(6, :)=[1 1 1]
    colormap(mycolormap)
    title("Correlation of " + percentage + "% of APOE33")



    MyThresh=ceil((x)*(332*332-332))
    inv_corr_apoe44=corr_apoe44(:);
    [val,myindex]=sort(abs(inv_corr_apoe44), 'ascend');
    MyNewThresh=inv_corr_apoe44(myindex(1:MyThresh));
    R=find(abs(inv_corr_apoe44)<=max(MyNewThresh));
    corr_apoe44_x=corr_apoe44; 
    corr_apoe44_x(R)=0;
    corr_apoe44_x=reshape(corr_apoe44_x,332,332);
    csvwrite(['C:\Users\natal\OneDrive - Duke University\Research\Independent Study Fall 2020-NW\',...
        'APOE44_' , num2str(Thresholds(i)), '.csv'],corr_apoe44_x);
    figure(17 + count)
    imagesc(corr_apoe44_x)
    colorbar
    caxis([-1 1])
    axis equal
    mycolormap=jet(11)
    mycolormap(6, :)=[1 1 1]
    colormap(mycolormap)
    title("Correlation of " + percentage + "% of APOE44")
end


