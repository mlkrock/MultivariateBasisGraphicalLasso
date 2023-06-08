%% Introduction and Setup Data
% reading in the variables and standardizing them
% download CAM ensemble from https://gdex.ucar.edu/dataset/371_abaker.html (warning - it is 500 GB)
% run this script in the same folder as the 343 .nc files
% this will produce the data_array used in example_dataanalysis.m

clc
clear all

Files = dir(fullfile('*.nc'));
m = length(Files);

all_variable_names =["AEROD_v","ANRAIN","ANSNOW","AODDUST1","AODDUST3","AODVIS","AQRAIN","AQSNOW","AREI","AREL","AWNC","AWNI","BURDEN1","BURDEN2","BURDEN3","BURDENBC","BURDENDUST","BURDENPOM","BURDENSEASALT","BURDENSO4","BURDENSOA","CCN3","CDNUMC","CLDHGH","CLDICE","CLDLIQ","CLDLOW","CLDMED","CLDTOT","CLOUD","DCQ","DMS_SRF","DTCOND","DTV","DTWR_H2O2","DTWR_H2SO4","DTWR_SO2","EMISCLD","FICE","FLDS","FLNS","FLNSC","FLNT","FLNTC","FLUT","FLUTC","FREQI","FREQL","FREQR","FREQS","FSDS","FSDSC","FSNS","FSNSC","FSNT","FSNTC","FSNTOA","FSNTOAC","FSUTOA","H2O2_SRF","H2SO4_SRF","ICEFRAC","ICIMR","ICWMR","IWC","LANDFRAC","LHFLX","LWCF","NUMICE","NUMLIQ","OCNFRAC","OMEGA","OMEGAT","PBLH","PHIS","PRECC","PRECL","PRECSC","PRECSL","PS","PSL","Q","QFLX","QREFHT","QRL","QRS","RELHUM","SHFLX","SNOWHICE","SNOWHLND","SO2_SRF","SOAG_SRF","SOLIN","SWCF","T","TAUGWX","TAUGWY","TAUX","TAUY","TGCLDCWP","TGCLDIWP","TGCLDLWP","TMQ","TREFHT","TS","TSMN","TSMX","U","U10","UU","V","VD01","VQ","VT","VU","VV","WGUSTD","WSUB","Z3","bc_a1_SRF","dst_a1_SRF","dst_a3_SRF","ncl_a1_SRF","ncl_a2_SRF","ncl_a3_SRF","num_a1_SRF","num_a2_SRF","num_a3_SRF","pom_a1_SRF","so4_a1_SRF","so4_a2_SRF","so4_a3_SRF","soa_a1_SRF","soa_a2_SRF"];
variable_dimlist =[2,3,3,2,2,2,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,3,2,2,3,3,2,2,2,3,3,2,3,3,3,3,3,3,3,2,2,2,2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,2,2,2,3,3,2,3,3,2,2,2,2,2,2,2,2,3,2,2,3,3,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,3,2,3,3,3,3,3,3,3,2,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];

keep_these_vars = [6,13,14,15,16,18,19,20,21,23,24,28,29,40,41,42,43,44,51,52,53,54,56,57,67,68,74,80,84,88,94,98,99,100,101,102,103,104,109];
all_variable_names = all_variable_names(keep_these_vars);
variable_dimlist = variable_dimlist(keep_these_vars); %should be all 2's

ncvar_start_index= [1,2];
ncvar_count_index= [Inf,Inf];
ncvar_stride_index = [1,1];

% read in longitude and latitude
% may want to save as .csv file rather than loading in the .nc files like this
lon = ncread(Files(1).name,"lon");
lat = ncread(Files(1).name,"lat");
n = length(lon);

% PRECT variable done separately from the other 39 variables
prect_mean = zeros([n,1]);
prect_variance = zeros([n,1]);

for i = 1:m
    prect_mean = prect_mean + ncread(Files(i).name,"PRECC",ncvar_start_index,ncvar_count_index,ncvar_stride_index) + ncread(Files(i).name,"PRECL",ncvar_start_index,ncvar_count_index,ncvar_stride_index);
end
prect_mean = prect_mean/m;

variable_by_files = zeros([n,m]);
var_holder = zeros([n,1]);

for i = 1:m
    variable_by_files(:,i) = -prect_mean + ncread(Files(i).name,"PRECC",ncvar_start_index,ncvar_count_index,ncvar_stride_index) + ncread(Files(i).name,"PRECL",ncvar_start_index,ncvar_count_index,ncvar_stride_index);
end

for i =1:n
    prect_variance(i) = var(variable_by_files(i,:));
end

prect = zeros([n,m]);
for i = 1:m
    prect(:,i)=  (-prect_mean + ncread(Files(i).name,"PRECC",ncvar_start_index,ncvar_count_index,ncvar_stride_index) + ncread(Files(i).name,"PRECL",ncvar_start_index,ncvar_count_index,ncvar_stride_index))./(sqrt(prect_variance));
end

% Now for the rest of the variables

ncvar_start_index= cell([length(all_variable_names) 1]);
ncvar_count_index= cell([length(all_variable_names) 1]);
ncvar_stride_index= cell([length(all_variable_names) 1]);

for i = 1:length(all_variable_names)
    if variable_dimlist(i)==2
        ncvar_start_index{i} = [1,2];
        ncvar_count_index{i} = [Inf,Inf];
        ncvar_stride_index{i} = [1,1];
    else
        ncvar_start_index{i} = [1,30,2];
        ncvar_count_index{i} = [Inf,1,1];
        ncvar_stride_index{i} = [1,1,1];
    end
    
end

pixelwise_mean = zeros(n, length(all_variable_names));

for j = 1:length(all_variable_names)    
    for i = 1:m
        ncid = ncread(Files(i).name,all_variable_names(j),ncvar_start_index{j},ncvar_count_index{j},ncvar_stride_index{j});
        pixelwise_mean(:,j) = ncid + pixelwise_mean(:,j);
    end
end

pixelwise_mean = pixelwise_mean/m;

pixelwise_variance = zeros([n length(all_variable_names)]);
tempvarvec = zeros([n m]);
varholder = zeros([n 1]);

for j = 1:length(all_variable_names)
    vmn = pixelwise_mean(:,j);

    for i = 1:m
        ncid = ncread(Files(i).name,all_variable_names(j),ncvar_start_index{j},ncvar_count_index{j},ncvar_stride_index{j});
        tempvarvec(:,i) = (ncid - vmn);
    end
    
    for k = 1:n
        varholder(k) = var(tempvarvec(k,:));
    end
    
    pixelwise_variance(:,j) = varholder;
end

clear variable_by_files
clear tempvarvec

%set up data

all_variable_names(length(all_variable_names)+1)="PRECT";
p = length(all_variable_names);

data_array = zeros([p,n,m]);

for j = 1:(p-1)
    vmn = pixelwise_mean(:,j);
    vsd = sqrt(pixelwise_variance(:,j));
    
    for i = 1:m
        ncid = ncread(Files(i).name,all_variable_names(j),ncvar_start_index{j},ncvar_count_index{j},ncvar_stride_index{j});
        data_array(j,:,i) = (ncid - vmn)./vsd;
    end
end

data_array(p,:,:) = prect;

%% save as data_array object as .nc file which is about 6 GB

for i=1:p
    nccreate('data_array.nc',all_variable_names(i),...
            'Dimensions',{'r',size(data_array,2),'c',size(data_array,3)},...
            'Format','64bit')
        
    ncwrite(      'data_array.nc',all_variable_names(i),squeeze(data_array(i,:,:)))
end

% for future use, can read in data_array from the data_array.nc file directly:

% p=40;
% n=48602;
% m=343;
% data_array = zeros([p,n,m]);
% for j = 1:p
%     data_array(j,:,:)= ncread('data_array.nc',all_variable_names(j));
% end

% now can follow example_dataanalysis.m with this data_array