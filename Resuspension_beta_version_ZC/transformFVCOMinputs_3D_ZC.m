function transformFVCOMoutputs_3D(conf_name)
% TRANSFORM OGCM outputs (SARCCM FVCOM version) to TrackMPD format
% I.Jalon-Rojas  8 July 2019; based on LoadFVCOMFiles_3D


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (SARCCM FVCOM version)
% FVCOM output file
% numlat number of points in the new rectangular grid (latitude dimension)
% numlon number of points in the new rectangular gird (longitude dimension)
% name of time variable

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,w,E,depth,time,time_str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the model configuration and inputs files

conf=feval(conf_name);

file = conf.OGCM.FVCOMFile;
grid = conf.OGCM.FVCOMGrid;

% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)


%% Define model parameters (SARCCM FVCOM Version)

numlon=conf.OGCM.NumLonGrid; 
numlat=conf.OGCM.NumLatGrid;
time_name=conf.OGCM.NameTime;
w_name=conf.OGCM.NameW;

[NumGridPts,numlvl,NTimeStamps]=size(ncread(file,'u'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the Grid info from Casename_Geo.grd file ZC  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --load   *geo.grd file for grid information-----------
fid=fopen(grid);
[dims]=textscan(fid,'%d%d',1,'Headerlines',1);%Read header
[info_node]=textscan(fid,'%d%f%f%f',dims{1,2}(1)); %Read info of nodes
[info_cell]=textscan(fid,'%d%d%d%d%d',dims{1,1}(1)); %Read info of cells
cell_num=dims{1,1}(1);
lat_node=info_node{1,3};
lon_node=info_node{1,2};
cell_3nodeid(:,1)=info_cell{1,3};
cell_3nodeid(:,2)=info_cell{1,4};
cell_3nodeid(:,3)=info_cell{1,5};
%----calculate latitude and longtitude for each cell----
for idt_cell=1:cell_num
lat_cell(idt_cell)=lat_node(cell_3nodeid(idt_cell,1))+lat_node(cell_3nodeid(idt_cell,2))+lat_node(cell_3nodeid(idt_cell,3));
lat_cell(idt_cell)=lat_cell(idt_cell)/3;
lon_cell(idt_cell)=lon_node(cell_3nodeid(idt_cell,1))+lon_node(cell_3nodeid(idt_cell,2))+lon_node(cell_3nodeid(idt_cell,3));
lon_cell(idt_cell)=lon_cell(idt_cell)/3;
end

lat_v=lat_cell; % lat at cells
lon_v=lon_cell; %lon at cells
x=lon_node; %lon at nodes
y=lat_node; %lat at nodes

% For FVCOM using spherical coordinates

% lat_v=double(ncread(file,'latc')); % lat at cells
% lon_v=double(ncread(file,'lonc')); %lon at cells
% x=double(ncread(file,'lon')); %lon at nodes
% y=double(ncread(file,'lat')); %lat at nodes

nv=double(ncread(file,'nv')); %number of cells
TR = triangulation(nv,x,y); %triangular grid

% Tranformation to rectangular grid

Lat=linspace(min(lat_v),max(lat_v),numlat);
Lon=linspace(min(lon_v),max(lon_v),numlon);

[Lon_matrix,Lat_matrix]=meshgrid(Lon,Lat);

lon_nn=linspace(min(x),max(x),numlon); %ZC 03/02/20
lat_nn=linspace(min(y),max(y),numlat);
[Lonn,Latn]=meshgrid(lon_nn,lat_nn);

for i=1:numlat
    for j=1:numlon

        ti(i,j) = pointLocation(TR,[Lon_matrix(i,j),Lat_matrix(i,j)]); %ti=Nan-->Land point

    end
end

% (water/land) mask
mask_water=zeros(size(ti));
mask_water(~isnan(ti))=1;
mask_land = ~mask_water;

mask_land3D = mask_land;
for i=1:numlvl-1
    mask_land3D = cat(3,mask_land3D,mask_land);
end

%Bottom Depth
h=double(ncread(file,'h'));
BottomDepth=double(griddata(x,y,h,Lonn,Latn,'nearest')); % Interpolation of bottom depth ZC 03/02/20
BottomDepth(mask_water==0)=1; %Land=1

% Verification plot
figure;
mesh(Lon,Lat,BottomDepth);
hold on
plot3(x,y,h,'o');
title('Verification plot for grid transformation: Bottom depth in the new grid')

% save grid
save([conf.Data.BaseDir '\grid.mat'],'Lat','Lon','BottomDepth','mask_water');
fprintf('saving grid\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save Time information    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT=double(ncread(file,time_name)); %modified julian date
TT=datetime(TT(:),'convertfrom','modifiedjuliandate');
timestamps=datenum(TT);
save([conf.Data.BaseDir '\timestamps.mat'],'timestamps');
fprintf('saving timestamps\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the variables varing with time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read variables

UU=double(ncread(file,'u')); % units = 'm/s'
VV=double(ncread(file,'v')); % units = 'm/s'
WW=double(ncread(file,w_name)); % variable name of the vertical velocity from FVCOM

ele=double(ncread(file,'zeta')); %elevation (m)

siglay=double(ncread(file,'siglay')); 
h_con=double(griddata(x,y,h,Lonn,Latn,'nearest')); %ZC 03/02/20

% Resuspension ZC 03/02/20

tauc=double(ncread(file,'tauc')); 
kv_aux = double(ncread(file,'kh')); 
siglev=double(ncread(file,'siglev'));
%kh_aux = double(ncread(ncfile,'viscofh')); % ZC for horizontal diffusivity 03/02/20 

for i=1:length(siglay(1,:))
   siglayzc(:,:,i)=h_con.*siglay(1,i);
end

for i=1:length(siglev(1,:))
   dsiglev(:,:,i)=h_con.*siglev(1,i); % ZC 03/02/20
end

% siglayzc=double(ncread(file,'siglayzc')); %can be read from output ZC

% Loop for each time step

zeros_matrix=zeros(numlat,numlon,1);
depth=nan(numlat,numlon,numlvl);

for i=1:NTimeStamps
    
    u=nan(numlat,numlon,numlvl);
    v=nan(numlat,numlon,numlvl);
    w=nan(numlat,numlon,numlvl);
    E=nan(numlat,numlon);
    
    % 3D variables
    for j=1:numlvl
        Uaux=griddata(lon_v,lat_v,UU(:,j,i),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of U in the new grid
        Uaux(mask_water==0)=0; %Land point=0
        Uaux(isnan(Uaux))=0; 
        u(:,:,j)=Uaux;

        Vaux=griddata(lon_v,lat_v,VV(:,j,i),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of U in the new grid
        Vaux(mask_water==0)=0;
        Vaux(isnan(Vaux))=0; 
        v(:,:,j)=Vaux;
        
        Waux=griddata(lon_v,lat_v,WW(:,j,i),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of U in the new grid
        Waux(mask_water==0)=0; 
        Waux(isnan(Waux))=0; 
        w(:,:,j)=Waux;
        
        if i==1 %depth does not change with time in FVCOM (independent of elevation)
        Depth_aux=siglayzc(:,:,j);
        Depth_aux(mask_water==0)=0; 
        depth(:,:,j)=Depth_aux;
        end
        
        clear Uaux Vaux Waux
    end
    
    %2D variables: elevation
    ELEaux=double(griddata(x,y,ele(:,i),Lonn,Latn,'nearest')); % Interpolation of wlw in the new grid ZC 03/02/20
    ELEaux(mask_water==0)=0; %Land point=0
    ELEaux(isnan(ELEaux))=0; 
    E(:,:)=ELEaux;
   
    clear ELEaux
    
    %Resuspension variables: bottom stress ZC 03/02/20
    Tauc_aux=double(griddata(x,y,tauc(:,i),Lonn,Latn,'nearest')); % Interpolation of tauc in the new grid
    Tauc_aux(mask_water==0)=0; %Land point=0
    Tauc_aux(Tauc_aux==nan)=0;
    TAUC(:,:)=Tauc_aux;
    
    clear Tauc_aux
    
       
    %1D variables: time
    time=timestamps(i);
    time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
    fprintf('changing format for time %s\n',time_str);
    
    
    % From m/s to cm/s
    
    u=u*100;
    v=v*100;
    w=w*100;
    
    TAUC=TAUC.*1025; % FVCOM tauc=Cff2*Umag*Umag, no RHO(density) ZC 03/02/20
    
    % Add a layer for surface and for bottom (for interpolation purpose near the boundaries)
    
    if i==1
    depth=cat(3,zeros_matrix,depth,BottomDepth);
    end
    
    depth=depth-depth(:,:,1); %Change in the reference system (surface constant, varying bottom)
    
    %Vertical diffusivity
    for j=1:length(siglev(1,:))
        Kvaux=griddata(x,y,kv_aux(:,j,i),Lonn,Latn,'nearest'); 
        Kvaux(mask_water==0)=0; 
        Kvaux(Kvaux==nan)=0; 
        Kvv(:,:,j)=Kvaux;
    end
    
     for j=1:numlat
     for k=1:numlon

        Kv(j,k,:) = interp1(squeeze(dsiglev(j,k,:)),squeeze(Kvv(j,k,:)),squeeze(depth(j,k,:)),'linear'); %interpolate Kv at siglev into siglayzc ZC 03/02/20

     end
     end
    
   clear Kvaux Kvv
    
   %Horizontal diffusivity
%    for j=1:1:length(siglay(1,:))
%         Khaux=griddata(x,y,kh_aux(:,j),Lonn,Latn,'nearest'); 
%         Khaux(mask_water==0)=0; 
%         Khaux(Khaux==nan)=0; 
%         Khh(:,:,j)=Khaux;
%    end

%    for j=1:numlat
%    for k=1:numlon
% 
%         Kh(j,k,:) = interp1(squeeze(siglayzc(j,k,:)),squeeze(Khh(j,k,:)),squeeze(depth(j,k,:)),'linear'); %interpolate Kv at siglev into siglayzc ZC 03/02/20
% 
%    end
%    end

%    clear Khaux Khh
    
    u=cat(3,u(:,:,1),u,zeros_matrix);
    v=cat(3,v(:,:,1),v,zeros_matrix);
    w=cat(3,w(:,:,1),w,zeros_matrix);
    
    % save data for each time step 
    save([conf.Data.BaseDir '\TrackMPDInput' num2str(i) '.mat'],'u','v','w','E','time','time_str','depth','TAUC','Kv');
    
end

end
